'''
Author: Ivy Strope
Date: 8/16/23
Contact: u247529@bcm.edu
Last Edited: 8/17/23
'''
 #Determine how to change this based upon input
sample=config['sample']
OUTDIR=config['outdir']
repo=config['repository']
configfile: 'config.yaml'

rule Human_Downstream:
    input:
        '{OUTDIR}/{sample}/final/{sample}_human.h5ad'
    output:
        '{OUTDIR}/{sample}/final/{sample}_human_processed.h5ad'
    shell:
        'python {repo}/scripts/Downstream.py --input {input} --output {output}'

rule Mouse_Downstream:
    input:
        '{OUTDIR}/{sample}/final/{sample}_mouse.h5ad'
    output:
        '{OUTDIR}/{sample}/final/{sample}_mouse_processed.h5ad'
    shell:
        'python {repo}/scripts/Downstream.py --input {input} --output {output}'
