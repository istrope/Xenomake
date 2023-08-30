'''
Author: 8/18/23
Contact: u247529@bcm.edu
Last Edited: 8/18/23 
'''

############################################
 # INCLUDE MAPPING MODULE AND POST PROCESSING
 ############################################
include: 'mapping.smk'
include: 'xengsort.smk'
include: 'expression.smk'
include: 'multimapped.smk'
include: 'annData.smk'
include: 'qc.smk'
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
default = [
    out_prefix + 'human_counts.tsv.gz',
    out_prefix + 'mouse_counts.tsv.gz',
    out_prefix + 'human.hdf5',
    out_prefix + 'mouse.hdf5',
    out_prefix + 'human.h5ad',
    out_prefix + 'mouse.hdf5',
    out_prefix + 'mouse.h5ad'
]

downstream = [
    out_prefix + 'human_processed.h5ad',
    out_prefix + 'mouse_processed.h5ad'
]

qc_prefix = '{OUTDIR}/{sample}/qc/'
qc = [qc_prefix + 'human_qc.png']
#############################################
#       SPECIFY FINAL OUTPUT FILES
##############################################
rule all:
    input:
        expand(default,sample=config['sample'],OUTDIR=config['outdir']),
        expand(downstream,sample=config['sample'],OUTDIR=config['outdir']) if config['downstream'] else next,
        expand(qc,sample=config['sample'],OUTDIR=config['outdir'])
    threads: config['threads']
