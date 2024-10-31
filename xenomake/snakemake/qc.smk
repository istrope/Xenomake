'''
Author: Ivy Strope
Date: 8/22/23
Contact: u247529@bcm.edu
Last Edited: 8/22/23
'''
#############################################
#       SPECIFY WILDCARD VARIABLES
#############################################
sample=config['project']['sample']
OUTDIR=config['project']['outdir']
repo=config['repository']
#############################################
#       SPECIFY PARAMETERS
##############################################
xengsort_outs = ['host','graft','both','ambiguous','neither']
in_prefix = '{OUTDIR}/{sample}/final/{sample}_'
out_prefix = '{OUTDIR}/{sample}/qc/'
#############################################
#       CREATE PIPELINE RULES
##############################################
rule Expression_QC:
    input:
        human = in_prefix + 'human.h5ad',
        mouse = in_prefix + 'mouse.h5ad'
    output: #only need to specify one file to make sure outputs correctly
        out_prefix + 'human_qc.png'
    params:#prefix for multiple qc plots
        human = out_prefix + 'human',
        mouse = out_prefix + 'mouse'
    shell:
        """
        python {repo}/scripts/qc.py --input {input.human} --out {params.human}
        python {repo}/scripts/qc.py --input {input.mouse} --out {params.mouse}
        """
