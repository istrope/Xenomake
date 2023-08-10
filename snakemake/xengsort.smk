
'''
Author: Ivy Strope
Created: 5/31/23
Contact: benjamin.strope@bcm.edu
Last Edited: 6/20/23
'''
 ############################################
 # INCLUDE MAPPING MODULE AND POST PROCESSING
 ############################################
include: 'mapping.smk'
#configfile: 'config.yaml'
#############################################
#       SPECIFY WILDCARD VARIABLE
#############################################
OUTDIR=config['outdir']
sample=config['sample']
threads=config['threads']
dropseq=config['dropseq_tools']
picard=config['picard']


rule all:
    input:
        human=expand('{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz',sample=config['sample'],OUTDIR=config['outdir']),
        mouse=expand('{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz',sample=config['sample'],OUTDIR=config['outdir'])
    threads: config['threads']

#############################################
#           EXTRACT MAPPED READS
#############################################
rule mapped_reads:
    input:
        human='{OUTDIR}/{sample}/mapping/human_final.bam',
        mouse='{OUTDIR}/{sample}/mapping/mouse_final.bam'
    output:
        human=temp('{OUTDIR}/{sample}/human_mapped.bam'),
        mouse=temp('{OUTDIR}/{sample}/mouse_mapped.bam')
    threads: config['threads']
    shell:
        """
        samtools view -F 4 -b {input.human} > {output.human}
        samtools view -F 4 -b {input.mouse} > {output.mouse}
        """
rule mapped_readnames:
    input:
        human='{OUTDIR}/{sample}/human_mapped.bam',
        mouse='{OUTDIR}/{sample}/mouse_mapped.bam'
    output:
        human=temp('{OUTDIR}/{sample}/human_readnames.txt'),
        mouse=temp('{OUTDIR}/{sample}/mouse_readnames.txt')
    threads: config['threads']
    shell:
        """
        samtools view {input.human} | cut -f1 | sort > {output.human}
        samtools view {input.mouse} | cut -f1 | sort > {output.mouse}
        """
rule overlapped_reads:
    input:
        human='{OUTDIR}/{sample}/human_readnames.txt',
        mouse='{OUTDIR}/{sample}/mouse_readnames.txt'
    output:
        '{OUTDIR}/{sample}/xengsort/overlapped_reads.txt'
    threads: config['threads']
    shell:
        'comm -12 {input.mouse} {input.human} > {output}'

rule subset_reads:
    input:
        bam='{OUTDIR}/{sample}/mouse_mapped.bam',
        readlist='{OUTDIR}/{sample}/xengsort/overlapped_reads.txt'
    output:
        temp('{OUTDIR}/{sample}/xengsort/overlapped.bam')
    threads: config['threads']
    params:
        picard=config['picard']
    shell:
        'java -jar {params.picard} FilterSamReads I={input.bam} O={output} READ_LIST_FILE={input.readlist} FILTER=includeReadList'

rule bam_to_fastq:
    input:
        '{OUTDIR}/{sample}/xengsort/overlapped.bam'
    output:
        '{OUTDIR}/{sample}/xengsort/overlapped.fastq'
    threads: config['threads']
    shell:
        'samtools fastq {input} > {output}'

#############################################
#           RUN XENOGRAFT SORTING
#############################################
rule xengsort_index:
    params:
        human='species/human/genome.fa',
        mouse='species/mouse/genome.fa'
    output:
        directory('species/idx.zarr')
    threads: config['threads']
    shell:
        'xengsort index -H {params.mouse} -G {params.human} -n 4_500_000_000 -W {threads} --index {output}'


rule xengsort_clasify:
    input:
        fastq='{OUTDIR}/{sample}/xengsort/overlapped.fastq',
        index='species/idx.zarr'
    output:
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.fq',
        graft='{OUTDIR}/{sample}/xengsort/{sample}-graft.fq',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.fq',
        neither='{OUTDIR}/{sample}/xengsort/{sample}-neither.fq',
        ambigious='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.fq'
    threads: config['threads']
    log:
        "{OUTDIR}/{sample}/xengsort/{sample}_xengsort.log"
    params:
        outprefix=lambda wc: f"{wc.OUTDIR}/{wc.sample}/xengsort/{wc.sample}",
        depug='-DD',
        prefetch=1
    run:
        """
        xengsort {params.debug} classify \
        --index {input.index} \
        --fastq {input.fastq} \
        --out {params.outprefix} \
        --threads {threads} \
        --chunksize 16.0 \
        --prefetch {params.prefetch} &> {log}'
        """
rule filter_files:
    input:
        mouse='{OUTDIR}/{sample}/xengsort/{sample}-host.fq',
        human='{OUTDIR}/{sample}/xengsort/{sample}-graft.fq',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.fq',
        neither='{OUTDIR}/{sample}/xengsort/{sample}-neither.fq',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.fq'
    output:
        human=temp('{OUTDIR}/{sample}/xengsort/filter_human.fastq'),
        mouse=temp('{OUTDIR}/{sample}/xengsort/filter_mouse.fastq')
    threads: config['threads']
    shell:
        '''
        cat {input.ambiguous} {input.both} {input.mouse} {input.neither} > {output.human}
        cat {input.ambiguous} {input.both} {input.human} {input.neither} > {output.mouse}
        '''

rule extract_readnames:
    input:
        human='{OUTDIR}/{sample}/xengsort/filter_human.fastq',
        mouse='{OUTDIR}/{sample}/xengsort/filter_mouse.fastq'
    output:
        human=temp('{OUTDIR}/{sample}/xengsort/filter_human.txt'),
        mouse=temp('{OUTDIR}/{sample}/xengsort/filter_mouse.txt')
    threads: config['threads']
    shell:
        """
        awk 'NR%4==1 {{print substr($1,2)}}' {input.human} > {output.human}
        awk 'NR%4==1 {{print substr($1,2)}}' {input.mouse} > {output.mouse}
        """
#############################################
# REMOVE NON-ORGANISM READS FROM MAPPED FILES
#############################################
rule subset_xengsort:
    input:
        human_bam='{OUTDIR}/{sample}/human_mapped.bam',
        mouse_bam='{OUTDIR}/{sample}/mouse_mapped.bam',
        human_reads='{OUTDIR}/{sample}/xengsort/filter_human.txt',
        mouse_reads='{OUTDIR}/{sample}/xengsort/filter_mouse.txt'
    output:
        human='{OUTDIR}/{sample}/final/{sample}_human_final.bam',
        mouse='{OUTDIR}/{sample}/final/{sample}_mouse_final.bam'
    threads: config['threads']
    params:
        picard=config['picard']
    shell:
        '''
        java -jar {params.picard} FilterSamReads I={input.human_bam} O={output.human} READ_LIST_FILE={input.human_reads} FILTER=excludeReadList
        java -jar {params.picard} FilterSamReads I={input.mouse_bam} O={output.mouse} READ_LIST_FILE={input.mouse_reads} FILTER=excludeReadList
        '''

#############################################
#           CREATE COUNT MATRIX
#############################################
rule create_dge:
    # creates the dge. depending on if the dge has _cleaned in the end it will require the
    # topBarcodesClean.txt file or just the regular topBarcodes.txt
    input:
        barcodes=config['barcode'],
        human='{OUTDIR}/{sample}/final/{sample}_human_final.bam',
        mouse='{OUTDIR}/{sample}/final/{sample}_mouse_final.bam'
    output:
        human='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz',
        human_summary='{OUTDIR}/{sample}/final/human.dge_summary',
        mouse='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz',
        mouse_summary='{OUTDIR}/{sample}/final/mouse.dge_summary'
    threads: config['threads']
    params:
        dropseq_tools=config['dropseq_tools']
    shell:
        """
        {params.dropseq_tools}/DigitalExpression I= {input.human} O= {output.human} \
        SUMMARY= {output.human_summary} CELL_BC_FILE={input.barcodes} \
        CELL_BARCODE_TAG=CB MOLECULAR_BARCODE_TAG=MI

        {params.dropseq_tools}/DigitalExpression I= {input.mouse} O= {output.mouse} \
        SUMMARY= {output.mouse_summary} CELL_BC_FILE={input.barcodes} \
        CELL_BARCODE_TAG=CB MOLECULAR_BARCODE_TAG=MI
        """
