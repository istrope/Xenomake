'''
Author: Ivy Strope
Creation: 6/16/23
Contact: benjamin.strope@bcm.edu
Last Edit: 8/18/23
'''
#############################################
#            SPECIFY PARAMETERS
##############################################
cell_barcode_flags="BASE_RANGE=" + str(config['project']['cell_barcode']) + " " + "BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1"
umi_flags="BASE_RANGE=" + str(config['project']['umi']) + " " + "BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=MI NUM_BASES_BELOW_QUALITY=1"
star_mapping_flags="--readFilesCommand samtools view -f 4 --readFilesType SAM SE --outSAMtype BAM Unsorted --genomeLoad NoSharedMemory --outSAMprimaryFlag AllBestScore --outSAMattributes All --outSAMunmapped None --outStd BAM_Unsorted --limitOutSJcollapsed 5000000"
#############################################
#        SPECIFY WILDCARD VARIABLE
#############################################
OUTDIR=config['project']['outdir']
sample=config['project']['sample']
threads=config['project']['threads']
dropseq=config['project']['dropseq_tools']
picard=config['project']['picard']
repo=config['repository']
############################################
#     TAG,TRIM,AND PREPROCESS READS
############################################
rule Generate_uBAM:
    input:
        read1=config['project']['r1'],
        read2=config['project']['r1']
    output:
        bam=temp('{OUTDIR}/{sample}/preprocess/unaligned.bam')
    log:
        stdout='{OUTDIR}/{sample}/logs/Generate_uBAM.log'
    params:
        Fastq = "PLATFORM=illumina SORT_ORDER=queryname SAMPLE_NAME={sample}"
    shell:
        """
        java -jar {picard} FastqToSam F1={input.read1} F2={input.read2} O={output.bam} {params.Fastq}
        """

rule Tag:
    input:
        bam='{OUTDIR}/{sample}/preprocess/unaligned.bam'
    output:
        processed_bam = temp('{OUTDIR}/{sample}/preprocess/tagged.bam')
    params:
        barcode=config['spatial']['cell_barcode'],
        umi=config['spatial']['umi']
    log:
        stdout='{OUTDIR}/{sample}/logs/Preprocess.log'
    threads: config['project']['threads']
    shell:
        """
        python {repo}/scripts/tag_reads.py --input_bam {input.bam} --output_bam {output.processed_bam} \
         --threads {threads} --barcode_range {params.barcode} --umi_range {params.umi}
        """

rule Generate_SE_Ubam:
    input:
        bam=temp('{OUTDIR}/{sample}/preprocess/tagged.bam'),
    output:
        bam=temp('{OUTDIR}/{sample}/preprocess/SE.bam'),
        header=temp('{OUTDIR}/{sample}/preprocess/fixed_header.sam')
    threads: config['project']['threads']
    shell:
        """
        samtools view -@ {threads} -b -f 128 -o {output.bam} {input.bam}
        """

rule Trim:
    input:
        read='{OUTDIR}/{sample}/preprocess/SE.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/preprocess/unaligned_bc_umi_tagged.bam')
    threads: config['project']['threads']
    shell:
        "python {repo}/scripts/trim_sequence.py {input.bam} {output.bam} --num_threads {threads}"

############################################
#      RUN STAR INDEX AND ALIGNMENT
############################################

rule Index_Genomes:
    input:
        mouse_ref = 'species/mouse/genome.fa',
        mouse_annotation='species/mouse/annotation.gtf',
        human_ref='species/human/genome.fa',
        human_annotation = 'species/human/annotation.gtf',
    output:
        human_index=directory('species/human/star_index/'),
        human_index_file='species/human/star_index/SAindex',
        mouse_index=directory('species/mouse/star_index/'),
        mouse_index_file='species/mouse/star_index/SAindex'
    threads: config['project']['threads']
    shell:
        """
            mkdir {output.human_index} &
	    STAR --runMode genomeGenerate \
                 --runThreadN {threads} \
                 --genomeDir {output.human_index} \
                 --genomeFastaFiles {input.human_ref} \
                 --sjdbGTFfile {input.human_annotation}

            mkdir {output.mouse_index} &
	    STAR --runMode genomeGenerate \
                 --runThreadN {threads} \
                 --genomeDir {output.mouse_index} \
                 --genomeFastaFiles {input.mouse_ref} \
                 --sjdbGTFfile {input.mouse_annotation}
        """
rule STAR_Human:
    input:
        ubam='{OUTDIR}/{sample}/preprocess/unaligned_bc_umi_tagged.bam',
        index='species/human/star_index/'
    output:
        aln=temp('{OUTDIR}/{sample}/mapping/human_Aligned.out.bam'),
    params:
        file_prefix = '{OUTDIR}/{sample}/mapping/human_',
        annotation='species/human/annotation.gtf',
        aln=star_mapping_flags,
        tmpdir = '{OUTDIR}/{sample}/mapping/human__STARgenome',
    log:
        log_progress='{OUTDIR}/{sample}/mapping/human_Log.progress.out',
        params='{OUTDIR}/{sample}/mapping/human_Log.out',
        sj='{OUTDIR}/{sample}/mapping/human_SJ.out.tab',
        log_final='{OUTDIR}/{sample}/mapping/human_Log.final.out',
        stdout='{OUTDIR}/{sample}/logs/human_TagGeneFunction.log'
    threads:
        config['project']['threads']
    shell:
        """
        STAR --genomeDir {input.index} \
            --readFilesIn {input.ubam} \
            --outFileNamePrefix {params.file_prefix} \
            --sjdbGTFfile {params.annotation} \
            {params.aln} \
            --runThreadN {threads} | \
            python {repo}/scripts/splice_bam_header.py --in-ubam {input.ubam} --in-bam /dev/stdin --out-bam /dev/stdout | \
            sambamba sort -t {threads} -n /dev/stdin -o /dev/stdout | \
            {dropseq}/TagReadWithGeneFunction I=/dev/stdin O={output.aln} ANNOTATIONS_FILE={params.annotation} &> {log.stdout}
        """

rule STAR_Mouse:
    input:
        ubam='{OUTDIR}/{sample}/preprocess/unaligned_bc_umi_tagged.bam',
        index='species/mouse/star_index/'
    output:
        aln='{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam',
    params:
        file_prefix = '{OUTDIR}/{sample}/mapping/mouse_',
        annotation='species/mouse/annotation.gtf',
        aln=star_mapping_flags,
        tmpdir = '{OUTDIR}/{sample}/mapping/mouse__STARgenome',
    log:
        log_progress='{OUTDIR}/{sample}/mapping/mouse_Log.progress.out',
        params='{OUTDIR}/{sample}/mapping/mouse_Log.out',
        sj='{OUTDIR}/{sample}/mapping/mouse_SJ.out.tab',
        log_final='{OUTDIR}/{sample}/mapping/mouse_Log.final.out',
        stdout='{OUTDIR}/{sample}/logs/mouse_TagGeneFunction.log'
    threads:
        config['project']['threads']
    shell:
        """
        STAR --genomeDir {input.index} \
            --readFilesIn {input.ubam} \
            --outFileNamePrefix {params.file_prefix} \
            --sjdbGTFfile {params.annotation} \
            {params.aln} \
            --runThreadN {threads} | \
            python {repo}/scripts/splice_bam_header.py --in-ubam {input.ubam} --in-bam /dev/stdin --out-bam /dev/stdout | 
            sambamba sort -t {threads} -n /dev/stdin -o /dev/stdout | \
            {dropseq}/TagReadWithGeneFunction I=/dev/stdin O={output.aln} ANNOTATIONS_FILE={params.annotation} &> {log.stdout}
        """




