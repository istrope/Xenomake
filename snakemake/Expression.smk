'''
Author: Ivy Strope
Created: 5/31/23
Contact: benjamin.strope@bcm.edu
Date Created: 8/8/23
Last Edited: 
'''

############################################
 # INCLUDE MAPPING MODULE AND POST PROCESSING
 ############################################
include: 'mapping.smk'
include: 'xengsort.smk'
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
##############################################
out_prefix='{OUTDIR}/{sample}/final/{sample}_'
expression_outs = [
    out_prefix + 'human_counts.tsv.gz',
    out_prefix + 'mouse_counts.tsv.gz',
    out_prefix + 'human.hdf5',
    out_prefix + 'mouse.hdf5',
    out_prefix + 'human.h5ad',
    out_prefix + 'mouse.hdf5'
]

rule all:
    input:
        expand(expression_outs,sample=config['sample'],OUTDIR=config['outdir'])
    threads: config['threads']

#############################################
#           CREATE COUNT MATRIX
#############################################
rule DGE_Mouse:
    input:
        barcodes=config['barcode'],
        bam='{OUTDIR}/{sample}/mouse_final.bam'
    output:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz',
        summary='{OUTDIR}/{sample}/final/mouse.dge_summary'
    threads: config['threads']
    params:
        'CELL_BARCODE_TAG=CB MOLECULAR_BARCODE_TAG=MI READ_MQ=0'
    shell:
        "{dropseq}/DigitalExpression I= {input.bam} O={output.counts} SUMMARY= {output.summary} CELL_BC_FILE={input.barcodes} {params}"

rule DGE_Human:
    input:
        bam='{OUTDIR}/{sample}/human_final.bam',
        barcodes=config['barcode']
    output:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz',
        summary='{OUTDIR}/{sample}/final/human.dge_summary',
    params: 'CELL_BARCODE_TAG=CB MOLECULAR_BARCODE_TAG=MI READ_MQ=0'
    shell:
        "{dropseq}/DigitalExpression I={input.bam} O={output.counts} SUMMARY={output.summary} CELL_BC_FILE={input.barcodes} {params}"

#Create Visium Tissue Positions File
rule create_tissue_file:
    input:
        dge='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz'
    output:
        spatial='{OUTDIR}/{sample}/final/{sample}_tissue_positions_list.csv'
    params:
        coordinates='barcodes/visium_barcode_positions.csv'
    threads: config['threads']
    run:
        """
        python scripts/spatial_barcode.py --input {snakemake.input.dge} \
        --output {snakemake.output.spatial} \
        --spatial_coordinates {snakemake.params.coordinates}
        """
#Replicate hdf5 architecture used by 10x genomics
rule HDF5_Human:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz'
    output:
        hdf5='{OUTDIR}/{sample}/final/{sample}_human.hdf5'
    threads: config['threads']
    shell:
        'Rscript scripts/create_h5.R --input {input.counts} --out {output.hdf5} --species human --assembly_version hg38'
        
rule HDF5_Mouse:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz'
    output:
        hdf5='{OUTDIR}/{sample}/final/{sample}_mouse.hdf5'
    threads: config['threads']
    shell:
        'Rscript scripts/create_h5.R --input {input.counts} --out {output.hdf5} --species mouse --assembly_version mm10'

#Create Anndata object for data storage and subsequent processing
rule H5ad_Human:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz',
    output:
        h5ad='{OUTDIR}/{sample}/final/{sample}_human.h5ad'
    params:
        tissue_file='{OUTDIR}/{sample}/final/tissue_positions_list.csv'
    threads: config['threads']
    run:
        """
        python scripts/make_h5ad.py --input {snakemake.input.counts} \
        --output {snakemake.output.h5ad} \
        --spatial_coordinates {snakemake.params.tissue_file}
        """

rule H5ad_Mouse:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz',
    output:
        h5ad='{OUTDIR}/{sample}/final/{sample}_mouse.h5ad'
    params:
        tissue_file='{OUTDIR}/{sample}/final/tissue_positions_list.csv'
    threads: config['threads']
    run:
        """
        python scripts/make_h5ad.py --input {snakemake.input.counts} \
        --output {snakemake.output.h5ad} \
        --spatial_coordinates {snakemake.params.tissue_file}
        """
