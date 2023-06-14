import os
import subprocess
import argparse
def get_args(required=True):
    parser=argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False
    )

    parser.add_argument(
        '--input_type',
        help='enter the string "fastq" or "bam" as the input type',
        required=True,
        type=str
    )
    parser.add_argument(
        '--r1',
        help='Path to the paired end r1 file generated from bcl2fastq',
        type=str,
        required=False,
        default=None
    )
    parser.add_argument(
        '--r2',
        help='Path to the paired end r2 file genereated from bcl2fastq',
        type=str,
        required=False,
        default=None
    )
    parser.add_argument(
        '--human_bam',
        help='Path to the trimmed,annotated,and mapped reads to the human genome',
        type=str,
        required=False,
        default=None
    )
    parser.add_argument(
        '--mouse_bam',
        help='Path to the trimmed,annotated,and mapped reads to the mouse genome',
        type=str,
        required=False,
        default=None
    )
    parser.add_argument(
        '--outdir',
        help='Desitnation directory',
        required=True,
        type=str
    )
    parser.add_argument(
        '--sample',
        help='sample name used to prepend filenames',
        required=True,
        type=str
    )
    parser.add_argument(
        '--barcode',
        help='visium barcode file provided by 10X genomics',
        default='barcodes/visium-v2.txt',
        required=False
    )
    parser.add_argument(
        '--threads',
        help='number of cores to use',
        required=False,
        default=8
    )
    args = parser.parse_args()
    return args
fn = get_args()
if fn['input_type'] == 'bcl'& fn['bcl'] == None:
    raise Exception('please provide the path for the bcl file with --bcl flag')

if fn['input_type'] == 'fastq' & ((fn['r1'] == None) or (fn['r2'] == None)):
    raise Exception('please provide the fastq file with --r1 and --r2 flags')

#Write parser object to config.yaml file to be used in snakemake commands
import yaml
file=open('workflow/config.yaml','w')
yaml.dump(fn,file)
file.close()
