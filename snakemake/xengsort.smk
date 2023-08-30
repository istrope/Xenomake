
'''
Author: Ivy Strope
Created: 5/31/23
Contact: benjamin.strope@bcm.edu
Date Created: 6/20/23
Last Edited: 8/8/23
'''
 ############################################
 # INCLUDE MAPPING MODULE AND POST PROCESSING
 ############################################
configfile: 'config.yaml'
#############################################
#       SPECIFY WILDCARD VARIABLE
#############################################
OUTDIR=config['outdir']
sample=config['sample']
threads=config['threads']
dropseq=config['dropseq_tools']
picard=config['picard']

#############################################
#       SPECIFY PARAMETERS
#############################################
xengsort_outs = ['host','graft','both','ambiguous','neither']
#############################################
#           EXTRACT MAPPED READS
#############################################
rule Find_Overlap_Reads:
    input:
        mouse='{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam',
        human='{OUTDIR}/{sample}/mapping/human_Aligned.out.bam'
    output:
        temp('{OUTDIR}/{sample}/sorting/overlapped.txt')
    threads: config['threads']
    shell:
        "python scripts/overlapped_reads.py --human {input.human} --mouse {input.mouse} --out {output}"

rule Remove_Overlapped_Reads:
    input:
        mouse_bam='{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam',
        human_bam='{OUTDIR}/{sample}/mapping/human_Aligned.out.bam',
        overlapped_reads='{OUTDIR}/{sample}/sorting/overlapped.txt'
    output:
        human=temp('{OUTDIR}/{sample}/human_remaining.bam'),
        mouse=temp('{OUTDIR}/{sample}/mouse_remaining.bam')
    threads: config['threads']
    log:
        stdout='{OUTDIR}/{sample}/logs/subset_unique_reads.log'
    shell:
        """
        java -jar {picard} FilterSamReads I={input.mouse_bam} O={output.mouse} READ_LIST_FILE={input.overlapped_reads} FILTER=excludeReadList &>> {log.stdout}
        java -jar {picard} FilterSamReads I={input.human_bam} O={output.human} READ_LIST_FILE={input.overlapped_reads} FILTER=excludeReadList &>> {log.stdout}
        """

rule Subset_Overlapped:
    input:
        mouse_bam='{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam',
        human_bam='{OUTDIR}/{sample}/mapping/human_Aligned.out.bam',        
        readlist='{OUTDIR}/{sample}/sorting/overlapped.txt'
    output:
        mouse_bam=temp('{OUTDIR}/{sample}/sorting/mouse_overlapped.bam'),
        human_bam=temp('{OUTDIR}/{sample}/sorting/human_overlapped.bam')
    log:
        stdout='{OUTDIR}/{sample}/logs/subset_overlapped_reads.log'    
    shell:
        """
        java -jar {picard} FilterSamReads I={input.mouse_bam} O={output.mouse_bam} READ_LIST_FILE={input.readlist} FILTER=includeReadList &>> {log.stdout}
        java -jar {picard} FilterSamReads I={input.human_bam} O={output.human_bam} READ_LIST_FILE={input.readlist} FILTER=includeReadList &>> {log.stdout}
        """
#############################################
#           RUN XENOGRAFT SORTING
#############################################

rule Xengsort_Index:
    input:
        human='species/human/genome.fa',
        mouse='species/mouse/genome.fa'
    output:
        directory('species/idx.zarr')
    threads: config['threads']
    params:
        stdout='{OUTDIR}/{sample}/logs/xengsort_index.log'
    shell:
        """
        xengsort index -H {input.mouse} -G {input.human} \
        -n 4_500_000_000 -W {threads} --index {output} &> {params.stdout}
        """


rule Xengsort_Clasify:
    input:
        bam='{OUTDIR}/{sample}/sorting/mouse_overlapped.bam',
        index='species/idx.zarr'
    output: #only need to specify one file, if this succeeds then all files generated correctly
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.fq',
        graft='{OUTDIR}/{sample}/xengsort/{sample}-graft.fq',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.fq',
        neither='{OUTDIR}/{sample}/xengsort/{sample}-neither.fq',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.fq'
    threads: config['threads']
    params:
        stdout="{OUTDIR}/{sample}/logs/xengsort_classify.log",
        outprefix=lambda wc: f"{wc.OUTDIR}/{wc.sample}/xengsort/{wc.sample}",
        debug='-DD',
        prefetch=1
    shell:
        """
        java -jar {picard} SamToFastq I={input.bam} FASTQ=/dev/stdout | \
        xengsort {params.debug} classify \
        --index {input.index} \
        --fastq /dev/stdin \
        --out {params.outprefix} \
        --threads {threads} \
        --chunksize 16.0 \
        --prefetch {params.prefetch} &> {params.stdout}
        """

rule Compress_Xengsort_Fastqs:
    input:
        expand('{OUTDIR}/{sample}/xengsort/{sample}-{wildcard}.fq',OUTDIR=config['outdir'],sample=config['sample'],wildcard=xengsort_outs)
    output:
        expand('{OUTDIR}/{sample}/xengsort/{sample}-{wildcard}.fq.gz',OUTDIR=config['outdir'],sample=config['sample'],wildcard=xengsort_outs)
    shell:
        'gzip {input}'

rule Extract_Xenograft_Readnames:
    input:
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.fq.gz',
        graft='{OUTDIR}/{sample}/xengsort/{sample}-graft.fq.gz',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.fq.gz',
        neither='{OUTDIR}/{sample}/xengsort/{sample}-neither.fq.gz',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.fq.gz'

    output:
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.txt',
        graft='{OUTDIR}/{sample}/xengsort/{sample}-graft.txt',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.txt',
        neither='{OUTDIR}/{sample}/xengsort/{sample}-neither.txt',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.txt'
    shell:
        """
        zcat {input.host} |grep "@"|sed "s/@//g"| sort -u > {output.host}
        zcat {input.graft} |grep "@"|sed "s/@//g"| sort -u > {output.graft}
        zcat {input.ambiguous} |grep "@"|sed "s/@//g"| sort -u > {output.ambiguous}
        zcat {input.neither} |grep "@"|sed "s/@//g"| sort -u > {output.neither}
        zcat {input.both} |grep "@"|sed "s/@//g"| sort -u > {output.both}

        """

