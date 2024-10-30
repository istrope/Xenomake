'''
Author: Ivy Strope
Created: 5/31/23
Contact: u247529@bcm.edu
Last Edited: 11/13/23
'''
#############################################
#       SPECIFY WILDCARD VARIABLE
#############################################
OUTDIR=config['project']['outdir']
sample=config['project']['sample']
threads=config['project']['threads']
dropseq=config['project']['dropseq_tools']
picard=config['project']['picard']
repo=config['repository']

#############################################
#           CREATE COUNT MATRIX
#############################################
locus_1 = "LOCUS_FUNCTION_LIST=CODING LOCUS_FUNCTION_LIST=UTR"
locus_2 = "LOCUS_FUNCTION_LIST=INTERGENIC LOCUS_FUNCTION_LIST=INTRONIC LOCUS_FUNCTION_LIST=UTR LOCUS_FUNCTION_LIST=UTR LOCUS_FUNCTION_LIST=CODING LOCUS_FUNCTION_LIST=RIBOSOMAL

if config['project']['spatial_mode'] == 'seq-scope':
    rule CreateBarcodes:
        input:
            r1=config['project']['r1'],
            r2=config['project']['r2']
        output:
            coordinates='spatial_coordinates.txt',
            whitelist='whitelist.txt'
        params:
            hdmilength=20
        shell:
            "{repo}/{scripts}/extract_seqscope_coordinates.sh {input.r1} {input.r2} {params.hdmilength}"

if config['project']['spatial_mode'] == 'sc10x_v3':
    rule uncompress_whitelist:
        input:
            '{repo}/data/barcodes/3M-febrary-2018.txt.gz'    
        output:
            '{repo}/data/barcodes/3M-febrary-2018.txt'
        shell:
            "gunzip {input}"

rule DGE_Mouse:
    input:
        barcodes=config['repository'] +'/' + config['spatial']['barcode_file'],
        bam='{OUTDIR}/{sample}/mouse_final.bam'
    output:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz',
        summary='{OUTDIR}/{sample}/qc/mouse.dge_summary',
        long='{OUTDIR}/{sample}/final/{sample}_CellCell_counts.tsv.gz'
    threads: config['project']['threads']
    log:
        stdout='{OUTDIR}/{sample}/logs/mouse_dge.log'
    params:
        default = 'CELL_BARCODE_TAG=CB MOLECULAR_BARCODE_TAG=MI READ_MQ=0 TMP_DIR=/tmp',
        locus = locus_1 if (config['run']['genic_only'] in ['True','true',True]) else locus_2
    shell:
        """
        {dropseq}/DigitalExpression I={input.bam} O={output.counts} \
        SUMMARY={output.summary} CELL_BC_FILE={input.barcodes} {params.default} \
        {params.locus} &> {log.stdout}

        python {repo}/scripts/convert_long.py --counts {output.counts} --output {output.long}
        """

rule DGE_Human:
    input:
        bam='{OUTDIR}/{sample}/human_final.bam',
        barcodes=config['repository'] + '/' + config['spatial']['barcode_file']
    output:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz',
        summary='{OUTDIR}/{sample}/qc/human.dge_summary',
        long='{OUTDIR}/{sample}/final/{sample}_CellInteract_counts.tsv.gz'
    params:
        default = 'CELL_BARCODE_TAG=CB MOLECULAR_BARCODE_TAG=MI READ_MQ=0 TMP_DIR=/tmp',
        locus = locus_1 if (config['run']['genic_only'] in ['True','true',True]) else locus_2
    log:
        stdout='{OUTDIR}/{sample}/logs/human_dge.log'
    shell:
        """
        {dropseq}/DigitalExpression I={input.bam} O={output.counts} \
        SUMMARY={output.summary} CELL_BC_FILE={input.barcodes} {params.default} \
        {params.locus} &> {log.stdout}

        python {repo}/scripts/convert_long.py --counts {output.counts} --output {output.long}
        """


#Replicate hdf5 architecture used by 10x genomics
rule HDF5_Human:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_human_counts.tsv.gz'
    output:
        hdf5='{OUTDIR}/{sample}/final/{sample}_human.hdf5'
    threads: config['project']['threads']
    shell:
        """
        cat species/human/annotation.gtf | bash {repo}/scripts/gtf_parse.sh > species/human/ids.txt
        python {repo}/scripts/create_hdf5.py --counts {input.counts} --output {output.hdf5} --species human --assembly hg38
        """
        
rule HDF5_Mouse:
    input:
        counts='{OUTDIR}/{sample}/final/{sample}_mouse_counts.tsv.gz'
    output:
        hdf5='{OUTDIR}/{sample}/final/{sample}_mouse.hdf5'
    threads: config['project']['threads']
    shell:
        """
        cat species/mouse/annotation.gtf | bash {repo}/scripts/gtf_parse.sh > species/mouse/ids.txt
        python {repo}/scripts/create_hdf5.py --counts {input.counts} --output {output.hdf5} --species mouse --assembly mm10
        """
