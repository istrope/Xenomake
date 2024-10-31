'''
Author: Ivy Strope
Date: 8/16/23
Contact: u247529@bcm.edu
Last Edited: 8/17/23
'''
 #Determine how to change this based upon input
sample=config['project']['sample']
OUTDIR=config['project']['outdir']
repo=config['repository']


#Create Anndata object for data storage and subsequent processing
rule H5ad_Human:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz',
        tissue_file=config['spatial']['coordinate_file']
    output:
        h5ad='{OUTDIR}/{sample}/final/{sample}_human.h5ad'
    threads: config['project']['threads']
    shell:
        """
        python {repo}/scripts/make_h5ad.py --counts {input.counts} \
        --output {output.h5ad} \
        --spatial_coordinates {input.tissue_file}
        """

rule H5ad_Mouse:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz',
        tissue_file=config['spatial']['coordinate_file']
    output:
        h5ad='{OUTDIR}/{sample}/final/{sample}_mouse.h5ad'
    threads: config['project']['threads']
    shell:
        """
        python {repo}/scripts/make_h5ad.py --counts {input.counts} \
        --output {output.h5ad} \
        --spatial_coordinates {input.tissue_file}
        """

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
