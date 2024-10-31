 '''
Author: Ivy Strope
Created: 8/18/23
Contact: u247529@bcm.edu
Last Edited: 8/18/23 
'''
 
 ############################################
 # INCLUDE MAPPING MODULE AND POST PROCESSING
 ############################################
sample=config['project']['sample']
OUTDIR=config['project']['outdir']
repo=config['repository']
#############################################
#           FILTER_MULTIMAPPED
#############################################

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
    shell:
        """
        java -jar {picard} FilterSamReads I={input.bam} O={output.host} READ_LIST_FILE={input.host} FILTER=includeReadList SORT_ORDER=queryname
        java -jar {picard} FilterSamReads I={input.bam} O={output.ambiguous} READ_LIST_FILE={input.ambiguous} FILTER=includeReadList SORT_ORDER=queryname
        java -jar {picard} FilterSamReads I={input.bam} O={output.both} READ_LIST_FILE={input.both} FILTER=includeReadList SORT_ORDER=queryname
        """

rule Human_Xenograft_BAM_Files:
    input:
        bam='{OUTDIR}/{sample}/sorting/human_overlapped.bam',
        host='{OUTDIR}/{sample}/xengsort/{sample}-host.txt',
        ambiguous='{OUTDIR}/{sample}/xengsort/{sample}-ambiguous.txt',
        both='{OUTDIR}/{sample}/xengsort/{sample}-both.txt'
    output:
        graft=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-graft.bam'),
        ambiguous=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-ambiguous.bam'),
        both=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-both.bam')
    shell:
        """
        java -jar {picard} FilterSamReads I={input.bam} O={output.graft} READ_LIST_FILE={input.host} FILTER=includeReadList SORT_ORDER=queryname
        java -jar {picard} FilterSamReads I={input.bam} O={output.ambiguous} READ_LIST_FILE={input.ambiguous} FILTER=includeReadList SORT_ORDER=queryname
        java -jar {picard} FilterSamReads I={input.bam} O={output.both} READ_LIST_FILE={input.both} FILTER=includeReadList SORT_ORDER=queryname
        """

rule filter_mm_both:
    input:
        mouse='{OUTDIR}/{sample}/xengsort/mouse/{sample}-both.bam',
        human='{OUTDIR}/{sample}/xengsort/human/{sample}-both.bam'
    output:
        mouse=temp('{OUTDIR}/{sample}/xengsort/mouse/{sample}-both.filtered.bam'),
        human=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-both.filtered.bam')
    shell:
        "python {repo}/scripts/filter_mm_reads_both.py --in1 {input.human} --in2 {input.mouse} --out1 {output.human} --out2 {output.mouse}"

rule filter_mm_ambiguous:
    input:
        mouse='{OUTDIR}/{sample}/xengsort/mouse/{sample}-ambiguous.bam',
        human='{OUTDIR}/{sample}/xengsort/human/{sample}-ambiguous.bam'
    output:
        mouse=temp('{OUTDIR}/{sample}/xengsort/mouse/{sample}-ambiguous.filtered.bam'),
        human=temp('{OUTDIR}/{sample}/xengsort/human/{sample}-ambiguous.filtered.bam')
    shell:
        "python {repo}/scripts/filter_mm_reads_both.py --in1 {input.human} --in2 {input.mouse} --out1 {output.human} --out2 {output.mouse}"


#############################################
        # MERGE BAM ALIGNMENTS
#############################################
input_mouse_a = ['{OUTDIR}/{sample}/mouse_remaining.bam',
                '{OUTDIR}/{sample}/xengsort/mouse/{sample}-host.bam',
                '{OUTDIR}/{sample}/xengsort/mouse/{sample}-both.filtered.bam',
                '{OUTDIR}/{sample}/xengsort/mouse/{sample}-ambiguous.filtered.bam']
input_mouse_b = ['{OUTDIR}/{sample}/mouse_remaining.bam',
                '{OUTDIR}/{sample}/xengsort/mouse/{sample}-host.bam']
rule Merge_Mouse:
    input:
        input_mouse_a if (config['run']['ambiguous'] in [True,'True','true']) else input_mouse_b
    output:
        bam=temp('{OUTDIR}/{sample}/mouse_merged.bam')
    threads: config['project']['threads']
    shell:
        "sambamba merge -t {threads} {output.bam} {input}"

input_human_a = ['{OUTDIR}/{sample}/human_remaining.bam',
                '{OUTDIR}/{sample}/xengsort/human/{sample}-graft.bam',
                '{OUTDIR}/{sample}/xengsort/human/{sample}-both.filtered.bam',
                '{OUTDIR}/{sample}/xengsort/human/{sample}-ambiguous.filtered.bam']
            
input_human_b = ['{OUTDIR}/{sample}/human_remaining.bam',
                '{OUTDIR}/{sample}/xengsort/human/{sample}-graft.bam']
rule Merge_Human:
    input:
        input_human_a if (config['run']['ambiguous'] in [True,'True','true']) else input_human_b
    output:
        bam=temp('{OUTDIR}/{sample}/human_merged.bam')
    shell:
        "sambamba merge -t {threads} {output.bam} {input}"

rule Filter_MultiMapped_Human:
    input:
        '{OUTDIR}/{sample}/human_merged.bam',
    output:
        '{OUTDIR}/{sample}/human_final.bam'
    log:
        stdout='{OUTDIR}/{sample}/logs/filter_mm.log'
    shell:
        "python {repo}/scripts/filter_mm_reads.py --in-bam {input} --out-bam {output} &>> {log.stdout}"

rule Filter_MultiMapped_Mouse:
    input:
        '{OUTDIR}/{sample}/mouse_merged.bam'
    output:
        '{OUTDIR}/{sample}/mouse_final.bam'
    log:
        stdout='{OUTDIR}/{sample}/logs/filter_mm.log'
    shell:
        "python {repo}/scripts/filter_mm_reads.py --in-bam {input} --out-bam {output} &>> {log.stdout}"
