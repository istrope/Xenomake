'''
Author: Ivy Strope
Creation: 6/16/23
Contact: benjamin.strope@bcm.edu
Last Edit: 6/20/23
'''
#############################################
#       SPECIFY PARAMETERS
##############################################
cell_barcode_flags="BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1"
umi_flags="BASE_RANGE=17-28 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=MI NUM_BASES_BELOW_QUALITY=1"
star_mapping_flags="--readFilesCommand samtools view -f 4 --readFilesType SAM SE --outSAMtype BAM Unsorted --genomeLoad NoSharedMemory --outSAMprimaryFlag AllBestScore --outSAMattributes All --outSAMunmapped None --outStd BAM_Unsorted --limitOutSJcollapsed 5000000"
#############################################
#       SPECIFY WILDCARD VARIABLE
#############################################
configfile: 'config.yaml'
OUTDIR=config['outdir']
sample=config['sample']
threads=config['threads']
dropseq=config['dropseq_tools']
picard=config['picard']

############################################
#     TAG,TRIM,AND PREPROCESS READS
############################################
rule FastqtoSam:
    input:
        read1=config['r1'],
        read2=config['r2']
    output:
        '{OUTDIR}/{sample}/preprocess/unaligned.bam'
    threads: config['threads']
    params:
        "PLATFORM=illumina SORT_ORDER=queryname SAMPLE_NAME={sample} "
    shell:
        """
        java -jar {picard} FastqToSam F1={input.read1} F2={input.read2} O={output} {params}
        """

rule TagCellBarcodes:
    input:
        '{OUTDIR}/{sample}/preprocess/unaligned.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/preprocess/unaligned_tagged_cell.bam'),
        summary='{OUTDIR}/{sample}/preprocess/tag_cell_barcodes.summary'
    params:
        cell_barcode_flags
    threads: config['threads']
    shell:
        """
        {dropseq}/TagBamWithReadSequenceExtended INPUT={input} OUTPUT={output.bam} \
        SUMMARY={output.summary} {params}
        """
rule TagUMI:
    input:
        '{OUTDIR}/{sample}/preprocess/unaligned_tagged_cell.bam',
    output:
        bam=temp('{OUTDIR}/{sample}/preprocess/unaligned_tagged_molecular.bam'),
        summary='{OUTDIR}/{sample}/preprocess/tag_umi.summary',
    params:
        umi_flags
    threads: config['threads']
    shell:
        """
        {dropseq}/TagBamWithReadSequenceExtended INPUT={input} OUTPUT={output.bam} SUMMARY={output.summary} {params}
        """
rule TrimAdapter:
    input:
        '{OUTDIR}/{sample}/preprocess/unaligned_tagged_molecular.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/preprocess/tagged_trimmed.bam'),
        summary='{OUTDIR}/{sample}/preprocess/trimmed_starting_sequence.summary'
    params:
        "SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"
    threads: config['threads']
    shell:
        """
        {dropseq}/TrimStartingSequence INPUT={input} OUTPUT={output.bam} OUTPUT_SUMMARY={output.summary} {params}
        """


rule PolyATrimmer:
    input:
        '{OUTDIR}/{sample}/preprocess/tagged_trimmed.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/preprocess/tagged_polyA_adapter_trimmed.bam'),
        summary='{OUTDIR}/{sample}/preprocess/polyATrimmer.summary'
    params:
        "MISMATCHES=0 NUM_BASES=6"
    threads: config['threads']
    shell:
        """
        {dropseq}/PolyATrimmer OUTPUT_SUMMARY={output.summary} INPUT={input} OUTPUT={output.bam} {params}
        """

rule Generate_SE_Ubam:
    input:
        bam='{OUTDIR}/{sample}/preprocess/tagged_polyA_adapter_trimmed.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/preprocess/unaligned_se.bam')
    threads: config['threads']
    shell:
        """
        samtools view {input.bam} | awk "NR%2==0" | samtools view -b - > {output.bam}
        """

rule Fix_Header:
    input:
        bam='{OUTDIR}/{sample}/preprocess/unaligned_se.bam',
        header='{OUTDIR}/{sample}/preprocess/unaligned.bam'
    output:
        bam='{OUTDIR}/{sample}/preprocess/unaligned_bc_umi_tagged.bam'
    threads: config['threads']
    shell:
        """
        samtools view -H {input.header} | \
        samtools reheader - {input.bam} > {output.bam}
        """
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
        human_index=directory('species/human/star_index'),
        human_index_file='species/human/star_index/SAindex',
        mouse_index=directory('species/mouse/star_index'),
        mouse_index_file='species/mouse/star_index/SAindex'
    threads: config['threads']
    shell:
        """
	    mkdir {output.human_index} &&
	    STAR --runMode genomeGenerate \
                 --runThreadN {threads} \
                 --genomeDir {output.human_index} \
                 --genomeFastaFiles {input.human_ref} \
                 --sjdbGTFfile {input.human_annotation}

	    mkdir {output.mouse_index} &&
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
        sj='{OUTDIR}/{sample}/mapping/human_SJ.out.tab',
        log_final='{OUTDIR}/{sample}/mapping/human_Log.final.out',
        log_progress='{OUTDIR}/{sample}/mapping/human_Log.progress.out',
        log='{OUTDIR}/{sample}/mapping/human_Log.out',
    params:
        file_prefix = '{OUTDIR}/{sample}/mapping/human_',
        annotation='species/human/annotation.gtf',
        aln=star_mapping_flags,
        tmpdir = '{OUTDIR}/{sample}mapping/mouse__STAR*'
    threads:
        config['threads']
    shell:
        """
        STAR --genomeDir {input.index} \
            --readFilesIn {input.ubam} \
            --outFileNamePrefix {params.file_prefix} \
            --sjdbGTFfile {params.annotation} \
            {params.aln} \
            --runThreadN {threads} | \
            python scripts/splice_bam_header.py --in-ubam {input.ubam} --in-bam /dev/stdin --out-bam /dev/stdout | \
            sambamba sort -n /dev/stdin -o /dev/stdout | \
            {dropseq}/TagReadWithGeneFunction I=/dev/stdin O={output.aln} ANNOTATIONS_FILE={params.annotation}
        """

rule STAR_Mouse:
    input:
        ubam='{OUTDIR}/{sample}/preprocess/unaligned_bc_umi_tagged.bam',
        index='species/mouse/star_index/'
    output:
        aln=temp('{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam'),
        sj='{OUTDIR}/{sample}/mapping/mouse_SJ.out.tab',
        log_final='{OUTDIR}/{sample}/mapping/mouse_Log.final.out',
        log_progress='{OUTDIR}/{sample}/mapping/mouse_Log.progress.out',
        log='{OUTDIR}/{sample}/mapping/mouse_Log.out'
    params:
        file_prefix = '{OUTDIR}/{sample}/mapping/mouse_',
        annotation='species/mouse/annotation.gtf',
        aln=star_mapping_flags,
        tmpdir = '{OUTDIR}/{sample}mapping/mouse__STAR*'
    threads:
        config['threads']
    shell:
        """
        STAR --genomeDir {input.index} \
            --readFilesIn {input.ubam} \
            --outFileNamePrefix {params.file_prefix} \
            --sjdbGTFfile {params.annotation} \
            {params.aln} \
            --runThreadN {threads} | \
            python scripts/splice_bam_header.py --in-ubam {input.ubam} --in-bam /dev/stdin --out-bam /dev/stdout | 
            sambamba sort -n /dev/stdin -o /dev/stdout | \
            {dropseq}/TagReadWithGeneFunction I=/dev/stdin O={output.aln} ANNOTATIONS_FILE={params.annotation}
        """




