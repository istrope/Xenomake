
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
#       SPECIFY tWILDCARD VARIABLE
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
rule Extract_Human_Readnames:
    input:
        bam='{OUTDIR}/{sample}/mapping/human_Aligned.out.bam',
    output:
        reads=temp('{OUTDIR}/{sample}/sorting/human_readnames.txt'),
    threads: config['threads']
    shell:
        "samtools view -F 4 {input.bam} | cut -f1 | sort > {output.reads}"
rule Extract_Mouse_Readnames:
    input:
        bam='{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam'
    output:
        reads=temp('{OUTDIR}/{sample}/sorting/mouse_readnames.txt')
    shell:
        "samtools view -F 4 {input.bam} | cut -f1 | sort > {output.reads}"

rule Find_Overlapped_Reads:
    input:
        human='{OUTDIR}/{sample}/sorting/human_readnames.txt',
        mouse='{OUTDIR}/{sample}/sorting/mouse_readnames.txt',
    output:
        temp('{OUTDIR}/{sample}/sorting/overlapped.txt')
    threads: config['threads']
    params:
        picard=config['picard']
    shell:
        "comm -12 {input.mouse} {input.human} > {output}"

rule Remove_Overlapped_Reads:
    input:
        mouse_bam='{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam',
        human_bam='{OUTDIR}/{sample}/mapping/human_Aligned.out.bam',
        overlapped_reads='{OUTDIR}/{sample}/sorting/overlapped.txt'
    output:
        human=temp('{OUTDIR}/{sample}/human_remaining.bam'),
        mouse=temp('{OUTDIR}/{sample}/mouse_remaining.bam')
    threads: config['threads']
    shell:
        """
        java -jar {picard} FilterSamReads I={input.mouse_bam} O={output.mouse} READ_LIST_FILE={input.overlapped_reads} FILTER=excludeReadList
        java -jar {picard} FilterSamReads I={input.human_bam} O={output.human} READ_LIST_FILE={input.overlapped_reads} FILTER=excludeReadList
        """

rule Mouse_Overlapped:
    input:
        bam='{OUTDIR}/{sample}/mapping/mouse_Aligned.out.bam',
        readlist='{OUTDIR}/{sample}/sorting/overlapped.txt'
    output:
        bam=temp('{OUTDIR}/{sample}/sorting/mouse_overlapped.bam')
    threads: config['threads']
    shell:
        "java -jar {picard} FilterSamReads I={input.bam} O={output.bam} READ_LIST_FILE={input.readlist} FILTER=includeReadList"

rule Human_Overlapped:
    input:
        bam='{OUTDIR}/{sample}/mapping/human_Aligned.out.bam',
        readlist='{OUTDIR}/{sample}/sorting/overlapped.txt'
    output:
        bam=temp('{OUTDIR}/{sample}/sorting/human_overlapped.bam')
    threads: config['threads']
    shell:
        "java -jar {picard} FilterSamReads I={input.bam} O={output.bam} READ_LIST_FILE={input.readlist} FILTER=includeReadList"
#############################################
#           RUN XENOGRAFT SORTING
#############################################

rule Xengsort_Index:
    params:
        human='species/human/genome.fa',
        mouse='species/mouse/genome.fa'
    output:
        directory('species/idx.zarr')
    threads: config['threads']
    shell:
        'xengsort index -H {params.mouse} -G {params.human} -n 4_500_000_000 -W {threads} --index {output}'


rule Xengsort_Clasify:
    input:
        bam='{OUTDIR}/{sample}/sorting/mouse_overlapped.bam',
        index='species/idx.zarr'
    output: 
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.fq',
        graft='{OUTDIR}/{sample}/xengsort/{sample}-graft.fq',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.fq',
        neither='{OUTDIR}/{sample}/xengsort/{sample}-neither.fq',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.fq'
    threads: config['threads']
    log:
        "{OUTDIR}/{sample}/xengsort/{sample}_xengsort.log"
    params:
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
        --prefetch {params.prefetch} &> {log}
        """

#############################################
#           FILTER_MULTIMAPPED
#############################################
rule Extract_Xenograft_Readnames:
    input:
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.fq',
        graft='{OUTDIR}/{sample}/xengsort/{sample}-graft.fq',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.fq',
        neither='{OUTDIR}/{sample}/xengsort/{sample}-neither.fq',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.fq'
    
    output:
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.txt',
        graft='{OUTDIR}/{sample}/xengsort/{sample}-graft.txt',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.txt',
        neither='{OUTDIR}/{sample}/xengsort/{sample}-neither.txt',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.txt'
    
    threads: config['threads']
    shell:
        """
        cat {input.host} |grep "@"|sed "s/@//g"| sort -u > {output.host}
        cat {input.graft} |grep "@"|sed "s/@//g"| sort -u > {output.graft}
        cat {input.ambiguous} |grep "@"|sed "s/@//g"| sort -u > {output.ambiguous}
        cat {input.neither} |grep "@"|sed "s/@//g"| sort -u > {output.neither}
        cat {input.both} |grep "@"|sed "s/@//g"| sort -u > {output.both}

        """

rule Mouse_Xenograft_BAM_Files:
    input:
        bam='{OUTDIR}/{sample}/sorting/mouse_overlapped.bam',
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.txt',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.txt',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.txt'

    output:
        host=temp('{OUTDIR}/{sample}/xengsort/mouse/{sample}-host.bam'),
        ambiguous=temp('{OUTDIR}/{sample}/xengsort/mouse/{sample}-ambiguous.bam'),
        both=temp('{OUTDIR}/{sample}/xengsort/mouse/{sample}-both.bam')
    threads:config['threads']
    shell:
        """
        java -jar {picard} FilterSamReads I={input.bam} O={output.host} READ_LIST_FILE={input.host} FILTER=includeReadList
        java -jar {picard} FilterSamReads I={input.bam} O={output.ambiguous} READ_LIST_FILE={input.ambiguous} FILTER=includeReadList
        java -jar {picard} FilterSamReads I={input.bam} O={output.both} READ_LIST_FILE={input.both} FILTER=includeReadList
        """

rule Human_Xenograft_BAM_Files:
    input:
        bam='{OUTDIR}/{sample}/sorting/human_overlapped.bam',
        graft='{OUTDIR}/{sample}/xengsort/{sample}-graft.txt',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.txt',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.txt'

    output:
        graft=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-graft.bam'),
        ambiguous=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-ambiguous.bam'),
        both=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-both.bam')
    threads:config['threads']
    shell:
        """
        java -jar {picard} FilterSamReads I={input.bam} O={output.graft} READ_LIST_FILE={input.graft} FILTER=includeReadList
        java -jar {picard} FilterSamReads I={input.bam} O={output.ambiguous} READ_LIST_FILE={input.ambiguous} FILTER=includeReadList
        java -jar {picard} FilterSamReads I={input.bam} O={output.both} READ_LIST_FILE={input.both} FILTER=includeReadList
        """

rule filter_mm_both:
    input:
        mouse='{OUTDIR}/{sample}/xengsort/mouse/{sample}-both.bam',
        human='{OUTDIR}/{sample}/xengsort/human/{sample}-both.bam'
    output:
        mouse=temp('{OUTDIR}/{sample}/xengsort/mouse/{sample}-both.filtered.bam'),
        human=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-both.filtered.bam')
    shell:
        "python scripts/filter_mm_reads_both.py --in1 {input.human} --in2 {input.mouse} --out1 {output.human} --out2 {output.mouse}"

rule filter_mm_ambiguous:
    input:
        mouse='{OUTDIR}/{sample}/xengsort/mouse/{sample}-ambiguous.bam',
        human='{OUTDIR}/{sample}/xengsort/human/{sample}-ambiguous.bam'
    output:
        mouse=temp('{OUTDIR}/{sample}/xengsort/mouse/{sample}-ambiguous.filtered.bam'),
        human=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-ambiguous.filtered.bam')
    shell:
        "python scripts/filter_mm_reads_both.py --in1 {input.human} --in2 {input.mouse} --out1 {output.human} --out2 {output.mouse}"


#############################################
        # MERGE BAM ALIGNMENTS
#############################################
rule Merge_Mouse:
    input:
        unique='{OUTDIR}/{sample}/mouse_remaining.bam',
        host='{OUTDIR}/{sample}/xengsort/mouse/{sample}-host.bam',
        both='{OUTDIR}/{sample}/xengsort/mouse/{sample}-both.filtered.bam',
        ambiguous='{OUTDIR}/{sample}/xengsort/mouse/{sample}-ambiguous.filtered.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/mouse_merged.bam')
    shell:
        "sambamba merge -t {threads} {output.bam} {input.unique} {input.host} {input.both} {input.ambiguous}"
rule Merge_Human:
    input:
        unique='{OUTDIR}/{sample}/human_remaining.bam',
        graft='{OUTDIR}/{sample}/xengsort/human/{sample}-graft.bam',
        both='{OUTDIR}/{sample}/xengsort/human/{sample}-both.filtered.bam',
        ambiguous='{OUTDIR}/{sample}/xengsort/human/{sample}-ambiguous.filtered.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/human_merged.bam')
    shell:
        "sambamba merge -t {threads} {output.bam} {input.unique} {input.graft} {input.both} {input.ambiguous}"

rule Filter_MultiMapped_Reads:
    input:
        human='{OUTDIR}/{sample}/human_merged.bam',
        mouse='{OUTDIR}/{sample}/mouse_merged.bam'
    output:
        human='{OUTDIR}/{sample}/human_final.bam',
        mouse='{OUTDIR}/{sample}/mouse_final.bam'
    shell:
        """
        python scripts/filter_mm_reads.py --in-bam {input.human} --out-bam {output.human}
        python scripts/filter_mm_reads.py --in-bam {input.mouse} --out-bam {output.mouse}
        """
