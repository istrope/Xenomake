import sys
import os
import snakemake
import argparse
import yaml
import pandas as pd
import numpy as np
import scanpy as sc


from xenomake.errors import XenomakeError, FileWrongExtensionError
from xenomake.utils import check_index
from xenomake.config import ConfigFile, xenomake_init, ProjectDF, project_df
config_path = 'config.yaml'

class Xenomake:
    """Xenomake
    load processed data
    """

    def __init__(self):
        self.root = os.getcwd()
        self.config = ConfigFile.from_yaml('f{root}/config.yaml')

    def load_processed_adata(self,project_id,sample_id,organism):
        adata = sc.read(f'{self.root}/{project_id}/{sample_id}/final/{sample_id}_{organism}_processed.h5ad')
        return adata
    def load_raw_adata(self,project_id,sample_id,organism):
        adata = sc.read(f'{self.root}/{project_id}/{sample_id}/final/{sample_id}_{organism}.h5ad')
        return adata
    
def get_run_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)

    parser.add_argument(
        "--cores", default=1, type=int, help="number of cores to be used in total"
    )
    parser.add_argument(
        "--dryrun",
        "-n",
        action="store_true",
        help="invokes a dry snakemake run, printing only commands",
    )
    parser.add_argument(
        "--rerun-incomplete",
        "--ri",
        action="store_true",
        help="forces snakemake to rerun incompletely generated files",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="if a job fails, keep executing independent jobs",
    )
    return parser

def setup_init_parser(parent_parser):
    """setup_init_parser.

    :param parent_parser_subparsers:
    """
    parser_init = parent_parser.add_parser(
        "init",
        help="initialise spacemake: create config files, download genomes and annotations",
    )
    parser_init.add_argument(
        "--root_dir",
        default="",
        help="where to output the results of the spacemake run. defaults to .",
    )
    parser_init.add_argument(
        "--temp_dir",
        default="/tmp",
        help="temporary directory used when running spacemake. defaults to /tmp",
    )
    parser_init.add_argument(
        "--download_species",
        default=False,
        help="if set, upon initialisation, spacemake will download the mouse and human genome and index",
        action="store_true",
    )
    parser_init.add_argument(
        "--dropseq_tools",
        help="absolute path to dropseq_tools directory",
        required=True,
    )
    parser_init.set_defaults(func=xenomake_init)

    return parser_init

def setup_run_parser(parent_parser):
    parser_run = parent_parser.add_parser(
        'run',help='run xenomake',parents=[get_run_parser()]
    )
    parser_run.set_

def xenomake_config(args):
    
def xenomake_run(args):
    if not os.path.isfile(config_path):
        msg = 'please initialize xenomake run firts.\n'
        msg += 'please run xenomake setup.'
        raise XenomakeError(msg)
    samples = []
    projects = []
    targets = ['run_analysis']
    config_vars = {
        'project_df': pdf.file_path,
        'samples': samples,
        'projects': projects,
        'pwd': os.getcwd()
    }

    pdf.assert_valid()

    snakefile = os.path.join(os.path.dirname(__file__),'snakemake/main.smk')

    analysis_finished = snakemake.snakemake(
        snakefile,
        configfiles=[config_path],
        cores=args['cores'],
        dryrun=args['dryrun'],
        config=config_vars
    )
    if analysis_finished is False:
        raise XenomakeError('error occured while running xenomake')
    
parser_main = argparse.ArgumentParser(allow_abbrev=False,
                                      description='Xenomake: Pipeline for Spatial Xenograft Processing')
parser_main.add_argument('--version', action='store_true')
parser_main_subparsers = parser_main.add_subparsers(help='sub command help',dest='subcommand')
parser_run=None
parser_projects=None
parser_config=None
parser_init=None
parser_spatial=None

#Xenomake Init
parser_init = setup_init_parser(parser_main_subparsers)
#Xenomake Config
from xenomake.config import setup_config_parser
if os.path.isfile(config_path):
    cf = ConfigFile.from_yaml(config_path)
    parser_config = setup_config_parser(cf,parser_main_subparsers)

#Xenomake Species Setup
from xenomake.config import setup_species_parser
if os.path.isfile(config_path):
    pdf = ProjectDF(project_df,cf)
    parser_species = setup_species_parser(pdf,parser_main_subparsers)
    #Run Xenomake
    parser_run = setup_run_parser(pdf,parser_main_subparsers)


def cmdline():
    import importlib.metadata
    args = parser_main.parse_args()
    if args.version and args.subcommand is None:
        print(importlib.metadata.version('xenomake'))
        return 0
    else:
        del args.version
    
    parser_dict = {
        'config': parser_config,
        'species': parser_species,
        'run': parser_run,
        'main': parser_main
    }

    if 'func' in args:
        func = args.func

    else:
        if args.subcommaind is not None:
            parser_dict[args.subcommand].print_help()
        else:
            parser_dict['main'].print_help()
        return 0
    
    args = {key: value for key, value in vars(args).items() if value is  not None}
    args.pop('func',None)
    args.pop('subcommand',None)
    func(args)

if __name__ == '__main__':
    cmdline()