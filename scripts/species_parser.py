import os
import subprocess
import argparse

def species_parser(required=True):
    parser=argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False
    )
    parser.add_argument(
        '--mouse_ref',
        help='Path to the .fa file for mouse reference to be used in alignment and organism assignment',
        type=str,
        required=required
    )
    parser.add_argument(
        '--human_ref',
        help='Path to the .fa file for human reference to be used in alignment and organism assignment',
        type=str,
        required=True
    )
    parser.add_argument(
        '--human_annotation',
        help='Path to the .gtf file matching the human reference build',
        type=str,
        required=True
    )
    parser.add_argument(
        '--mouse_annotation',
        help='Path to the .gtf file matching the mouse reference build',
        required=True,
        type=str
    )
    parser.add_argument(
        '--xengsort_index',
        help='xengsort host and graft index directory',
        required=False,
        default=None
    )
    parser.add_argument(
        '--human_index',
        help='path to star index under human  genome build',
        required=False,
        default=None
    )
    parser.add_argument(
        '--mouse_index',
        help='path to star index under mouse genome build',
        required=False,
        default=None
    )
    args = parser.parse_args()
    return args

args = species_parser()

human_dir='species/human'
mouse_dir='species/mouse'

if os.path.isdir('species'):
    raise Exception('Species Directory Already Exists')
if not os.path.isdir('species'):
    os.system('mkdir -p %s' % human_dir)
    os.system('mkdir -p %s' % mouse_dir)

#symlink reference data annotations to new directory
os.symlink(args.mouse_ref,'species/mouse/genome.fa')
os.symlink(args.mouse_annotation,'species/mouse/annotation.gtf')

os.symlink(args.human_ref,'species/human/genome.fa')
os.symlink(args.human_annotation,'species/human/annotation.gtf')

if args.xengsort_index != None:
    os.symlink(args.xengsort_index,'species/idx.zarr')

if args.human_index != None:
    os.symlink(args.human_index,'species/human/star_index')

if args.mouse_index != None:
    os.symlink(args.mouse_index,'species/mouse/star_index')
