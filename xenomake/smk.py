import os
import snakemake
import argparse
import yaml
import scanpy as sc

import xenomake
from xenomake.errors import *
from xenomake.utils import *
from xenomake.config import ConfigFile, ConfigMainVariable, RunMode,Spatial_Setup
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
    
def setup_init_parser(parent_parser):
    """setup_init_parser.

    :param parent_parser_subparsers:
    """
    parser_init = parent_parser.add_parser(
        "init",
        help="initialise spacemake: create config files, download genomes and annotations",
    )
    parser_init.add_argument(
        '--read1',
        '-r1',
        help='Path to the paired end r1 file',
        type=str,
        required=True,
        default=None
    )
    parser_init.add_argument(
        '--read2',
        '-r2',
        help='Path to the paired end r2 file',
        type=str,
        required=True,
        default=None
    )
    parser_init.add_argument(
        '--spatial_mode',
        help='use preset run modes to run samples based upon technology used [custom,visium,slide-seq,dbit-seq,hdst,sc10x,seq-scope]',
        required=True,
        default=None,
        type=str
    )
    parser_init.add_argument(
        '--run_mode',
        help='use preset run modes to determine how ambiguous reads, multimapped reads, polyA tails, and intergenic sequences get processed'+
        '\nuse: [prude,lenient,custom]',
        default='lenient',
        required=True,
        type=str
    )
    parser_init.add_argument(
        '--sample',
        '-s',
        help='sample name used to prepend filenames and directories',
        required=True,
        type=str
    )
    parser_init.add_argument(
        '--outdir',
        '-o',
        help='output project directory',
        required=True,
        type=str
    )
    parser_init.add_argument(
        "--root_dir",
        default="./",
        help="where to output the results of the xenomake run. defaults to ."
    )
    parser_init.add_argument(
        "--temp_dir",
        default="/tmp",
        help="temporary directory used when running xenomake. defaults to /tmp",
    )
    parser_init.add_argument(
        "--download_species",
        default=False,
        help="if set, upon initialization, xenomake will download the mouse and human genome and annotations",
    )
    parser_init.add_argument(
        '--download_index',
        default=False,
        help='if set, upon initialization, xenomake will download the xengsort index and star index for genomes',
    )
    parser_init.add_argument(
        "--dropseq_tools",
        help="absolute path to dropseq_tools directory",
        required=False,
        default = os.path.join(xenomake.__path__[0],'data/Drop-seq_tools-2.5.3')
    )
    parser_init.add_argument(
        '--picard',
        help='absolute path to picard jar file',
        required=False,
        default = os.path.join(xenomake.__path__[0],'data/picard.jar')
    )
    parser_init.set_defaults(func=xenomake_init)

    return parser_init
def setup_config_parser(parent_parser):
    parser_config = parent_parser.add_parser(
        'config',
        help = 'create custom method for spatial technology'
    )
    parser_config.add_argument(
        '--ambiguous',
        help = 'handle multimap from "both" and "ambiguous" output reads (Must be boolean True/False)',
        required=False,
        default=None,
    )
    parser_config.add_argument(
        '--mm_reads',
        help='assign multimapped reads to highest confidence position',
        required=False,
        default=None,
    )
    parser_config.add_argument(
        '--downstream',
        help = 'Performs Scanpy Processing of Data (Must be boolean True/False)',
        required = False,
        default = None,
    )
    parser_config.add_argument(
        '--genic_only',
        help='use only genic/exonic reads and excluse flags for intergenic,intronic',
        default=None,
        required = False
    )
    parser_config.add_argument(
        '--barcode_file',
        type=str,help='barcode file used to demultiplex samples in digitial gene expression. default is visium',
        default=None,
        required=False,
    )
    parser_config.add_argument(
        '--umi',
        type = str,
        help='umi structure defining the positions of umi in units of bases'
        +'Example: Visium Cell Barcode is contained in bases 17-28 and umi flag would be'
        + '--umi 17-28',
        required = False
        )
    parser_config.add_argument(
        '--cell_barcode',
        required=False,
        help = 'cell barcode structure defining the positions of umi in units of bases'
        +'Example: Visium Cell Barcode is contained in bases 1-16 and umi flag would be'
        + '--cell_barcode 1-16',
        default=None
    )
    parser_config.add_argument(
        '--polyA_trimming',
        required=False,
        help='trim poly A tails off of reads before mapping',
        default=False
    )
    parser_config.set_defaults(func=xenomake_config)
    return parser_config

def setup_species_parser(parent_parser):
    parser_species= parent_parser.add_parser(
        'species',
        help = 'setup genome directory for xenomake run'
    )
    parser_species.add_argument(
        '-mr',
        '--mouse_reference',
        help='Path to the .fa file for mouse reference to be used in alignment and organism assignment',
        type=str,
        required=True
    )
    parser_species.add_argument(
        '-hr',
        '--human_reference',
        help='Path to the .fa file for human reference to be used in alignment and organism assignment',
        type=str,
        required=True
    )
    parser_species.add_argument(
        '-ha',
        '--human_annotation',
        help='Path to the .gtf file matching the human reference build',
        type=str,
        required=True
    )
    parser_species.add_argument(
        '-ma',
        '--mouse_annotation',
        help='Path to the .gtf file matching the mouse reference build',
        required=True,
        type=str
    )
    parser_species.add_argument(
        '--xengsort_hash',
        help='path to xengsort hash table file required for classification',
        required=False,
        default=None,
        
    )
    parser_species.add_argument(
        '--xengsort_info',
        help='path to xengsort info file requred for classification',
        required=False,
        default=None,
        type=str
    )
    parser_species.add_argument(
        '--human_index',
        help='path to star index under human  genome build',
        required=False,
        default=None,
        type=str
    )
    parser_species.add_argument(
        '--mouse_index',
        help='path to star index under mouse genome build',
        required=False,
        default=None,
        type=str
    )
    return parser_species

def get_run_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)

    parser.add_argument(
        "--cores", default=1, type=int, help="number of cores to be used in total",required=True
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
    parser_species.set_defaults(func=xenomake_species)
    return parser


def setup_run_parser(parent_parser):
    parser_run = parent_parser.add_parser(
        'run',help='run xenomake',parents=[get_run_parser()]
    )
    parser_run.set_defaults(func=lambda args: xenomake_run(args))
    parser_run_subparsers = parser_run.add_subparsers()

    return parser_run

#########################
# UTILITY FUNCTIONS
########################
def populate_default_spatmode(args):
    spatmode = args['spatial_mode']
    spat_config = ConfigFile.from_yaml(os.path.join(xenomake.__path__[0],'data/defaults_spatmode.yaml'))
    spat_args = match_default(spatmode,spat_config)
    return spat_args

def populate_default_runmode(args):
    runmode = args['run_mode']
    defaults = ConfigFile.from_yaml(os.path.join(xenomake.__path__[0],'data/defaults_runmode.yaml'))
    run_args = match_default(runmode,defaults)    
    return run_args

##################################
# RUNNING XENOMAKE SPECIES
##################################
def process_species_args(args):
        #check files have correct extensions
        if not os.path.exists('species'):
            os.system('mkdir species')
            os.system('mkdir species/human')
            os.system('mkdir species/mouse')
        
        assert_file(args['mouse_reference'],default_value=None,extension=[".fa",".fna",".fa.gz",".fna.gz"])
        assert_file(args['human_reference'],default_value=None,extension=[".fa",".fna",".fa.gz",".fna.gz"])

        assert_file(args['mouse_annotation'],default_value=None,extension=[".gtf"])
        assert_file(args['human_annotation'],default_value=None,extension=[".gtf"])

        #create symlinked directory for standard access in pipeline execution
        os.system('ln --force --symbolic "$(readlink --canonicalize %s)" species/human/annotation.gtf' % args['human_annotation'])
        os.system('ln --force --symbolic "$(readlink --canonicalize %s)" species/human/genome.fa' % args['human_reference'])
        
        os.system('ln --force --symbolic "$(readlink --canonicalize %s)" species/mouse/annotation.gtf' % args['mouse_annotation'])
        os.system('ln --force --symbolic "$(readlink --canonicalize %s)" species/mouse/genome.fa' % args['mouse_reference'])
        if 'human_index' in args:
            check_index(args['human_index'])
            os.system('ln --force --symbolic "$(readlink --canonicalize %s)" species/human/star_index' % args['human_index'])
        if 'mouse_index' in args:
            check_index(args['mouse_index'])
            os.system('ln --force --symbolic "$(readlink --canonicalize %s)" species/mouse/star_index' % args['mouse_index'])
        if ('xengsort_hash' in args) and ('xengsort_info' in args):
            assert_file(args['xengsort_hash'],default_value=None,extension=['.hash'])
            assert_file(args['xengsort_info'],default_value=None,extension=['.info'])
            os.system('ln --force --symbolic "$(readlink --canonicalize %s)" species/idx.hash' % args['xengsort_hash'])
            os.system('ln --force --symbolic "$(readlink --canonicalize %s)" species/idx.info' % args['xengsort_info'])

        print('\n species directory successfully completed, keep project execution in current directory')

##################################
#  INITIALIZE CONFIGFILE
#################################
import sys
def xenomake_init(args):
    #check inputs
    if os.path.isfile(config_path):
        raise Exception('config file exists, project has been initialized in current directory')
    try:
        validate_init_inputs(args)
    except XenomakeError as e:
        raise e
        sys.exit(1)
    initial_config = os.path.join(xenomake.__path__[0],'data/config_structure.yaml')
    cf = ConfigFile.from_yaml(initial_config)
    cf.set_file_path(config_path)

    #add variables from argument parser
    cf.variables['repository'] = xenomake.__path__[0]
    cf.variables['project']['sample'] = args['sample']
    cf.variables['project']['outdir'] = args['outdir']
    cf.variables['project']['r1'] = args['read1']
    cf.variables['project']['r2'] = args['read2']
    cf.variables['project']['run_mode'] = args['run_mode']
    cf.variables['project']['spatial_mode'] = args['spatial_mode']
    cf.variables['project']['root_dir'] = args['root_dir']
    cf.variables['project']['temp_dir'] = args['temp_dir']
    cf.variables['project']['dropseq_tools'] = args['dropseq_tools']
    cf.variables['project']['picard'] = args['picard']
    cf.variables['project']['download_index'] = args['download_index']
    cf.variables['project']['download_species'] = args['download_species']

    if args['download_species']:
        species_urls = os.path.join(xenomake.__path__[0],'data/urls/species_urls.txt')
        os.system('wget -i %s' %species_urls)
        os.system('gunzip *.gtf.gz')

    if args['download_index']:
        #set paths to url files
        hindex_urls = os.path.join(xenomake.__path__[0],'data/urls/human_index_urls.txt')
        mindex_urls = os.path.join(xenomake.__path__[0],'data/urls/mouse_index_urls.txt')
        xengsort_urls = os.path.join(xenomake.__path__[0],'data/urls/xengsort_index_urls.txt')
        #download
        os.system('mkdir human_index/')
        os.system('wget -P human_index -i %s' %hindex_urls)
        os.system('mkdir mouse_index')
        os.system('wget -P mouse_index -i %s' %mindex_urls)
        os.system('wget -i %s' %xengsort_urls)

    if args['spatial_mode'] in ['visium','seq-scope','dbit-seq']:
        spat_config = populate_default_spatmode(args)
        cf.variables['spatial'] = spat_config

    if args['run_mode'] in ['lenient','prude']:
        run_config = populate_default_runmode(args)
        cf.variables['run'] = run_config

    print('\nsample name: %s' % args['sample'])
    print('spatial chemistry: %s' % args['spatial_mode'])
    print('run mode: %s' % args['run_mode'])
    print('project %s initialized, proceed to species setup and project execution' % args['outdir'])
    cf.dump()


def xenomake_species(args):
    process_species_args(args)

################################ 
# HANDLING OF CUSTOM INPUTS FOR RUN MODE OR SPATIAL MODE 
###############################

def check_config(args,required_keys):
    missing_keys = [key for key in required_keys if key not in args or args[key] in (None,'')]
    if missing_keys:
        missing_keys_str = ', '.join(missing_keys)
        raise EmptyConfigVariableError(missing_keys_str)
def process_run_mode(args):
    # THIS FUNCTION IN UTILIZED IF CUSTOM METHOD IS UTILIZED AND WILL POPULATE CONFIG YAML FILE FROM ARGS
    cf = ConfigFile.from_yaml(config_path)
    check_config(args,['mmreads','ambiguous','downstream','genic_only','polyA_trimming'])

    cf.variables['run']['mmreads'] = args['mm_reads']
    cf.variables['run']['ambiguous'] = args['ambiguous']
    cf.variables['run']['downstream'] = args['downstream']
    cf.variables['run']['genic_only'] = args['genic_only']
    cf.variables['run']['polyA_trimming'] = args['polyA_trimming']

    cf.dump() 

def process_spatial_mode(args):
    #THIS FUNCTION IS UTILIZED IF CUSTOM METHOD IS DESCRIBED AND WILL POPULATE CONFIG FILE FROM ARGS
    cf = ConfigFile.from_yaml(config_path)
    check_config(args,['barcode_file','umi','cell_barcode'])

    cf.variables['spatial']['barcode_file'] = args['barcode_file']
    cf.variables['spatial']['umi'] = args['umi']
    cf.variables['spatial']['cell_barcode'] = args['cell_barcode']

    cf.dump()

def xenomake_config(args):
    cf = ConfigFile.from_yaml(config_path)
    if cf.variables['project']['run_mode'] == 'custom':
        process_run_mode(args)
    if cf.variables['project']['spatial_mode'] == 'custom':
        process_spatial_mode(args)

############################## 
#    RUN SNAKEMAKE
###############################
def xenomake_run(args):
    if not os.path.isfile(config_path):
        msg = 'please initialize xenomake run first.\n'
        msg += 'please run xenomake setup.'
        raise XenomakeError(msg)
    
    #load in config file and assure that all variables have been set up properly
    cf = ConfigFile.from_yaml(config_path)
    cf.variables['project']['threads'] = args['cores']
    cf.dump()
    cf.check_project()


    snakefile = os.path.join(xenomake.__path__[0],'snakemake/main.smk')

    os.system('snakemake ')
    analysis_finished = snakemake.snakemake(
        snakefile,
        configfiles=[config_path],
        cores=args['cores'],
        dryrun=args['dryrun'],
        keepgoing=args['keep_going']
    )
    if analysis_finished is False:
        raise XenomakeError('error occured while running xenomake')


##########################
        # Parser
##########################    
parser_main = argparse.ArgumentParser(allow_abbrev=False,description='Xenomake: Pipeline for Spatial Xenograft Processing')
from importlib.metadata import version
parser_main.add_argument('-v','--version',action='version',version=version('xenomake'))
parser_main_subparsers = parser_main.add_subparsers(help='sub command help',dest='subcommand')
parser_run=None
parser_species=None
parser_init=None
parser_config=None

#Xenomake Init
parser_init = setup_init_parser(parser_main_subparsers)
parser_species = setup_species_parser(parser_main_subparsers)

#Xenomake Config
parser_config = setup_config_parser(parser_main_subparsers)
#Run Xenomake
parser_run = setup_run_parser(parser_main_subparsers)


def cmdline():
    args = parser_main.parse_args()
    parser_dict = {
        'init' : parser_init,
        'species': parser_species,
        'config':parser_config,
        'run': parser_run,
        'main': parser_main
    }

    if 'func' in args:
        func = args.func

    else:
        if args.subcommand is not None:
            parser_dict[args.subcommand].print_help()
        else:
            parser_dict['main'].print_help()
        return 0

    args = {key: value for key, value in vars(args).items() if value is not None}
    args.pop('func',None)
    args.pop('subcommand',None)
    func(args)

if __name__ == '__main__':
    cmdline()
