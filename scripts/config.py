import os
import subprocess
import argparse
def get_args(required=True):
    parser=argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False
    )
    parser.add_argument(
        '--r1',
        help='Path to the paired end r1 file generated from bcl2fastq',
        type=str,
        required=True,
        default=None
    )
    parser.add_argument(
        '--r2',
        help='Path to the paired end r2 file genereated from bcl2fastq',
        type=str,
        required=True,
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
        '--ambiguous',
        help = 'handle multimap from "both" and "ambiguous" output reads (Must be boolean True/False)',
        required=False,
        default=True,
        type=bool
    )
    parser.add_argument(
        '--downstream',
        help = 'Performs Scanpy Processing of Data (Must be boolean True/False)',
        required = False,
        default = True,
        type = bool
    )
    parser.add_argument(
        '--barcode',
        help='visium barcode file provided by 10X genomics',
        default='barcodes/visium-v2.txt',
        required=False
    )
    parser.add_argument(
        '--picard',
        help='Path to picard jar executable',
        default='tools/picard.jar',
        required=False
    )
    parser.add_argument(
        '--dropseq_tools',
        help='Path to Dropseq tools directory',
        default='tools/Drop-seq_tools-2.5.3',
        required=False
    )
    parser.add_argument(
        '--threads',
        help='number of cores to use',
        required=False,
        default=8,
        type=int
    )
    parser.add_argument(
        '--reference',
        help='name of reference (default: genome)',
        type=str,
        default='genome'
    )
    parser.add_argument(
        '--mouse_assembly',
        help='name of mouse assembly version used (default=mm10)',
        default='mm10',
        type=str
    )
    parser.add_argument(
        '--human_assembly',
        help='human_assembly version (default=hg38)',
        type=str,
        default='hg38'
    )
    args = parser.parse_args()
    return args
fn = get_args()

#Write parser object to config.yaml file to be used in snakemake commands
import yaml
file=open('config.yaml','w')
yaml.dump(fn,file)
file.close()

with open('config.yaml', 'r') as fin:
    data = fin.read().splitlines(True)
with open('config.yaml', 'w') as fout:
    fout.writelines(data[1:])
