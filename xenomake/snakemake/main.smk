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
#       SPECIFY PARAMETERS
##############################################
out_prefix='{OUTDIR}/{sample}/final/{sample}_'
default = [
    out_prefix + 'human_counts.tsv.gz',
    out_prefix + 'mouse_counts.tsv.gz',
    out_prefix + 'human.hdf5',
    out_prefix + 'mouse.hdf5',
]

spatial = [
    out_prefix + 'mouse.h5ad',
    out_prefix + 'human.h5ad'
]
downstream = [
    out_prefix + 'human_processed.h5ad',
    out_prefix + 'mouse_processed.h5ad',
]

qc_prefix = '{OUTDIR}/{sample}/qc/'
qc = [qc_prefix + 'human_qc.png']


#############################################
#       SPECIFY FINAL OUTPUT FILES
##############################################
run_args = default + qc
if config['project']['spatial_mode'] not in ['sc10x_v2','sc10x_v3']:
    run_args += spatial

if config['run']['downstream'] in ('true',True,'True'):
    run_args += downstream
#############################################
#               RUN RULE ALL
#############################################
rule all:
    input:
        expand(run_args,sample=config['project']['sample'],OUTDIR=config['project']['outdir'])
