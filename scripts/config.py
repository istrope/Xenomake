import subprocess
import argparse
import os
def get_args(required=True):
    parser=argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False
    )
    parser.add_argument(
        '--repository',
        help='path to xenomake cloned directory',
        type=str,
        required=True,
        default=None
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
    return parser.parse_args()
fn = get_args()

if fn.repository.endswith('/'): #format directory correctly for snakemake scripts
    fn.repository = fn.repository[:-1]

#Write parser object to config.yaml file to be used in snakemake commands
import yaml
file=open('config.yaml','w')
yaml.dump(fn,file)
file.close()

with open('config.yaml', 'r') as fin:
    data = fin.read().splitlines(True)
with open('config.yaml', 'w') as fout:
    fout.writelines(data[1:])

#ENSURE THAT TOOLS ARE EXECUTABLE
if fn.dropseq_tools == 'tools/Drop-seq_tools-2.5.3tools/drop':
    dropseq_dir = fn.repository + fn.dropseq_tools
    os.system('chmod +x -R %s' % dropseq_dir)

if fn.picard == 'tools/picard.jar':
    picard_dir = fn.repository + fn.picard
    os.system('chmod +x %s' % picard_dir)

print('assembly versions %s' % [fn.human_assembly,fn.mouse_assembly])
print('downstream processing: %s' % fn.downstream)
print('ambiguous reads handling: %s' % fn.ambiguous)
print('cores: %s' % fn.threads)
print('project "%s" initialized, proceed to snakemake execution' % fn.sample)