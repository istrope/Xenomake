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
    args = parser.parse_arg()
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
os.symlink(args.mouse_ref,'species/mouse/ref.fa')
os.symlink(args.mouse_annotation,'species/mouse/ref.gtf')

os.symlink(args.human_ref,'species/human/ref.fa')
os.symlink(args.human_annotation,'species/human/ref.gtf')

if args.xengsort_index != None:
    os.symlink(args.xengsort_index,'species/idx.zarr')
