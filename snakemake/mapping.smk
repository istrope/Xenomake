'''
Author: Ivy Strope
Creation: 6/16/23
Contact: benjamin.strope@bcm.edu
Last Edit: 6/20/23
'''

#############################################
#       SPECIFY WILDCARD VARIABLE
#############################################
#configfile: 'config.yaml'
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
        '{OUTDIR}/{sample}/mapping/unaligned.bam'
    threads: config['threads']
    params:
        "PLATFORM=illumina SORT_ORDER=queryname SAMPLE_NAME={sample}"
    shell:
        """
        java -jar {picard} FastqToSam F1={input.read1} F2={input.read2} O={output} {params}
        """

rule tagCellBarcodes:
    input:
        '{OUTDIR}/{sample}/mapping/unaligned.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/mapping/unaligned_tagged_cell.bam'),
        summary='{OUTDIR}/{sample}/mapping/tag_cell_barcodes.summary'
    params:
        "BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1"
    threads: config['threads']
    shell:
        """
        {dropseq}/TagBamWithReadSequenceExtended INPUT={input} OUTPUT={output.bam} \
        SUMMARY={output.summary} {params}
        """
rule tagUMI:
    input:
        '{OUTDIR}/{sample}/mapping/unaligned_tagged_cell.bam',
    output:
        bam=temp('{OUTDIR}/{sample}/mapping/unaligned_tagged_molecular.bam'),
        summary='{OUTDIR}/{sample}/mapping/tag_umi.summary',
    params:
        "BASE_RANGE=17-28 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=MI NUM_BASES_BELOW_QUALITY=1"
    threads: config['threads']
    shell:
        """
        {dropseq}/TagBamWithReadSequenceExtended INPUT={input} OUTPUT={output.bam} SUMMARY={output.summary} {params}
        """
rule TrimStartingSequence:
    input:
        '{OUTDIR}/{sample}/mapping/unaligned_tagged_molecular.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/mapping/tagged_trimmed.bam'),
        summary='{OUTDIR}/{sample}/mapping/trimmed_starting_sequence.summary'
    params:
        "SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"
    threads: config['threads']
    shell:
        """
        {dropseq}/TrimStartingSequence INPUT={input} OUTPUT={output.bam} OUTPUT_SUMMARY={output.summary} {params}
        """


rule PolyATrimmer:
    input:
        '{OUTDIR}/{sample}/mapping/tagged_trimmed.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/mapping/tagged_polyA_adapter_trimmed.bam'),
        summary='{OUTDIR}/{sample}/mapping/polyATrimmer.summary'
    params:
        "MISMATCHES=0 NUM_BASES=6"
    threads: config['threads']
    shell:
        """
        {dropseq}/PolyATrimmer OUTPUT_SUMMARY={output.summary} INPUT={input} OUTPUT={output.bam} {params}
        """


rule SamtoFastq:
    input:
        bam='{OUTDIR}/{sample}/mapping/tagged_polyA_adapter_trimmed.bam'
    output:
        read1=temp('{OUTDIR}/{sample}/mapping/{sample}_processed_R1.fastq'),
        read2=temp('{OUTDIR}/{sample}/mapping/{sample}_processed_R2.fastq')
    threads: config['threads']
    shell:
        """
        java -jar {picard} SamToFastq I={input} FASTQ={output.read1} SECOND_END_FASTQ={output.read2}
        """
############################################
#      RUN STAR INDEX AND ALIGNMENT
############################################

rule index_genome:
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
	i
	    mkdir {output.mouse_index} &&
	    STAR --runMode genomeGenerate \
                 --runThreadN {threads} \
                 --genomeDir {output.mouse_index} \
                 --genomeFastaFiles {input.mouse_ref} \
                 --sjdbGTFfile {input.mouse_annotation}
        """
rule gzip_fastq:
    input:
	    r1='{OUTDIR}/{sample}/mapping/{sample}_processed_R1.fastq',
	    r2='{OUTDIR}/{sample}/mapping/{sample}_processed_R2.fastq'
    output:
	    r1='{OUTDIR}/{sample}/mapping/{sample}_processed_R1.fastq.gz',
	    r2='{OUTDIR}/{sample}/mapping/{sample}_processed_R2.fastq.gz'
    threads: config['threads']
    shell:
	    """
	    gzip -nc {input.r1} > {output.r1}
	    gzip -nc {input.r2} > {output.r2}
	    """
rule STAR_aln_human:
    input:
        reads=expand('{OUTDIR}/{sample}/mapping/{sample}_processed_{read}.fastq.gz',OUTDIR=config['outdir'],sample=config['sample'],read=['R1','R2']),
        index='species/human/star_index/'
    output:
        aln=temp('{OUTDIR}/{sample}/mapping/human_Aligned.out.bam'),
        sj='{OUTDIR}/{sample}/mapping/human_SJ.out.tab',
        log_final='{OUTDIR}/{sample}/mapping/human_Log.final.out',
        log_progress='{OUTDIR}/{sample}/mapping/human_Log.progress.out',
        log='{OUTDIR}/{sample}/mapping/human_Log.out'
    params:
        file_prefix = '{OUTDIR}/{sample}/mapping/human_',
        annotation='species/human/annotation.gtf',
        aln="--readFilesCommand zcat --outSAMtype BAM Unsorted --genomeLoad NoSharedMemory --outSAMprimaryFlag AllBestScore --outSAMattributes All --outSAMunmapped Within --outStd BAM_Unsorted --limitOutSJcollapsed 5000000"
    threads:
        config['threads']
    shell:
        """
        STAR --genomeDir {input.index} \
            --readFilesIn {input.reads} \
            --outFileNamePrefix {params.file_prefix} \
            --sjdbGTFfile {params.annotation} \
            {params.aln} \
            --runThreadN {threads} > {output.aln}
        """

rule STAR_aln_mouse:
    input:
        reads=expand('{OUTDIR}/{sample}/mapping/{sample}_processed_{read}.fastq.gz',OUTDIR=config['outdir'],sample=config['sample'],read=['R1','R2']),
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
        aln="--readFilesCommand zcat --outSAMtype BAM Unsorted --genomeLoad NoSharedMemory --outSAMprimaryFlag AllBestScore --outSAMattributes All --outSAMunmapped Within --outStd BAM_Unsorted --limitOutSJcollapsed 5000000"
    threads:
        config['threads']
    shell:
        """
        STAR --genomeDir {input.index} \
            --readFilesIn {input.reads} \
            --outFileNamePrefix {params.file_prefix} \
            --sjdbGTFfile {params.annotation} \
            {params.aln} \
            --runThreadN {threads} > {output.aln}
        """
############################################
#  SORT AND MERGE TAGS INTO ALIGNED FILE
############################################
rule SortSam:
    input:
        human='{OUTDIR}/{sample}/mapping/human_Aligned.out.bam',
        mouse='{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_sorted_aligned.bam'),
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_sorted_aligned.bam')
    threads: config['threads']
    shell:
        """
        java -jar {picard} SortSam I={input.human} O={output.human} SORT_ORDER=queryname
        java -jar {picard} SortSam I={input.mouse} O={output.mouse} SORT_ORDER=queryname
        """

rule MergeAlignment:
    input:
        human_aln='{OUTDIR}/{sample}/mapping/human_sorted_aligned.bam',
        ubam='{OUTDIR}/{sample}/mapping/tagged_polyA_adapter_trimmed.bam',
        mouse_aln='{OUTDIR}/{sample}/mapping/mouse_sorted_aligned.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_merged.bam'),
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_merged.bam')
    threads: config['threads']
    params:
        human_ref='species/human/genome.fa',
        mouse_ref='species/mouse/genome.fa',
        other='INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false'
    shell:
        """
        java -jar {picard} MergeBamAlignment REFERENCE_SEQUENCE={params.human_ref}\
                                            UNMAPPED_BAM={input.ubam} \
                                            ALIGNED_BAM={input.human_aln} \
                                            OUTPUT={output.human} \
                                            {params.other}
        java -jar {picard} MergeBamAlignment REFERENCE_SEQUENCE={params.mouse_ref} \
                                            UNMAPPED_BAM={input.ubam} \
                                            ALIGNED_BAM={input.mouse_aln} \
                                            OUTPUT={output.mouse} \
                                            {params.other}
        """
############################################
# TAG READ with GENE and FILTER MULTIMAPPED
############################################
rule TagReadWithGeneFunction:
    input:
        human='{OUTDIR}/{sample}/mapping/human_merged.bam',
        mouse='{OUTDIR}/{sample}/mapping/mouse_merged.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_final.bam'),
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_final.bam')
    params:
        human_annotation='species/human/annotation.gtf',
        mouse_annotation='species/mouse/annotation.gtf'
    threads: config['threads']
    shell:
        """
        {dropseq}/TagReadWithGeneFunction I={input.human} O={output.human} ANNOTATIONS_FILE={params.human_annotation}
        {dropseq}/TagReadWithGeneFunction I={input.human} O={output.mouse} ANNOTATIONS_FILE={params.human_annotation}
        """

#rule filter_mm_reads:
 #   input:
  #      human='{OUTDIR}/{sample}/mapping/human_merged_tagged.bam',
   #     mouse='{OUTDIR}/{sample}/mapping/mouse_merged_tagged.bam'
   # output:
    #    human='{OUTDIR}/{sample}/mapping/human_final.bam',
     3   mouse='{OUTDIR}/{sample}/mapping/mouse_final.bam'
   # shell:
    #    """
     #   python scripts/filter_mm_reads.py --in-bam {input.human} --out-bam {output.human}
      #  python scripts/filter_mm_reads.py --in-bam {input.mouse} --out-bam {output.mouse}
      #  """

