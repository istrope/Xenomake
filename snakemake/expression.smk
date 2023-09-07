'''
Author: Ivy Strope
Created: 5/31/23
Contact: u247529@bcm.edu
Last Edited: 8/18/23 
'''
#############################################
#       SPECIFY WILDCARD VARIABLE
#############################################
OUTDIR=config['outdir']
sample=config['sample']
threads=config['threads']
dropseq=config['dropseq_tools']
picard=config['picard']
repo=config['repository']
#############################################
#           CREATE COUNT MATRIX
#############################################
locus_1 = "LOCUS_FUNCTION_LIST=CODING LOCUS_FUNCTION_LIST=UTR"
locus_2 = ['INTERGENIC', 'INTRONIC', 'UTR', 'CODING', 'RIBOSOMAL']
rule DGE_Mouse:
    input:
        barcodes=config['repository'] +'/' + config['barcode'],
        bam='{OUTDIR}/{sample}/mouse_final.bam'
    output:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz',
        summary='{OUTDIR}/{sample}/qc/mouse.dge_summary'
    threads: config['threads']
    log:
        stdout='{OUTDIR}/{sample}/logs/mouse_dge.log'
    params:
        default = 'CELL_BARCODE_TAG=CB MOLECULAR_BARCODE_TAG=MI READ_MQ=0 TMP_DIR=/tmp',
        locus = locus_1 if config['genic_only'] else locus_2
    shell:
        """
        {repo}/{dropseq}/DigitalExpression I= {input.bam} O={output.counts} \
        SUMMARY={output.summary} CELL_BC_FILE={input.barcodes} {params.default} \
        {params.locus} &> {log.stdout}
        """

rule DGE_Human:
    input:
        bam='{OUTDIR}/{sample}/human_final.bam',
        barcodes=config['repository'] + '/' + config['barcode']
    output:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz',
        summary='{OUTDIR}/{sample}/qc/human.dge_summary',
    params:
        default = 'CELL_BARCODE_TAG=CB MOLECULAR_BARCODE_TAG=MI READ_MQ=0 TMP_DIR=/tmp',
        locus = locus_1 if (config['genic_only'] == True) else locus_2
    log:
        stdout='{OUTDIR}/{sample}/logs/human_dge.log'
    shell:
        """
        {repo}/{dropseq}/DigitalExpression I={input.bam} O={output.counts} \
        SUMMARY={output.summary} CELL_BC_FILE={input.barcodes} {params.default} \
        {params.locus} &> {log.stdout}
        """

#Create Visium Tissue Positions File
rule create_tissue_file:
    input:
        dge='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz'
    output:
        spatial='{OUTDIR}/{sample}/final/tissue_positions_list.csv'
    params:
        coordinates='barcodes/visium_barcode_positions.csv'
    threads: config['threads']
    shell:
        """
        python {repo}/scripts/spatial_barcode.py --counts {input.dge} \
        --output {output.spatial} \
        --spatial_coordinates {params.coordinates}
        """
#Replicate hdf5 architecture used by 10x genomics
rule HDF5_Human:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz'
    output:
        hdf5='{OUTDIR}/{sample}/final/{sample}_human.hdf5'
    threads: config['threads']
    shell:
        """
        cat species/human/annotation.gtf | {repo}/scripts/gtf_parse.sh > species/human/ids.txt
        python {repo}/scripts/create_hdf5.py --counts {input.counts} --output {output.hdf5} --species human --assembly hg38
        """
        
rule HDF5_Mouse:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz'
    output:
        hdf5='{OUTDIR}/{sample}/final/{sample}_mouse.hdf5'
    threads: config['threads']
    shell:
        """
        cat species/mouse/annotation.gtf | {repo}/scripts/gtf_parse.sh > species/mouse/ids.txt
        python {repo}/scripts/create_hdf5.py --counts {input.counts} --output {output.hdf5} --species mouse --assembly mm10
        """

#Create Anndata object for data storage and subsequent processing
rule H5ad_Human:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz',
        tissue_file='{OUTDIR}/{sample}/final/tissue_positions_list.csv'
    output:
        h5ad='{OUTDIR}/{sample}/final/{sample}_human.h5ad'
    threads: config['threads']
    shell:
        """
        python {repo}/scripts/make_h5ad.py --counts {input.counts} \
        --output {output.h5ad} \
        --spatial_coordinates {input.tissue_file}
        """

rule H5ad_Mouse:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz',
        tissue_file='{OUTDIR}/{sample}/final/tissue_positions_list.csv'
    output:
        h5ad='{OUTDIR}/{sample}/final/{sample}_mouse.h5ad'
    threads: config['threads']
    shell:
        """
        python {repo}/scripts/make_h5ad.py --counts {input.counts} \
        --output {output.h5ad} \
        --spatial_coordinates {input.tissue_file}
        """
