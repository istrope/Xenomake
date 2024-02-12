import argparse
import os
import yaml
from xenomake.utils import assert_file, check_index, str2bool
from xenomake.errors import *
def get_species_parser(required=True):
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
        '--xengsort_hash',
        help='path to xengsort hash table file required for classification',
        required=False,
        default=None
    )
    parser.add_argument(
        '--xengsort_info',
        help='path to xengsort info file requred for classification'
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
    return parser

def get_run_parser(required=True):
    parser=argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False
    )
    parser.add_argument(
        '--project',
        help='project directory',
        required=True,
        type=str
    )
    parser.add_argument(
        '--sample',
        help='sample name used to prepend filenames and directories',
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
        '--mm_reads',
        help='assign multimapped reads to highest confidence position',
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
        '--genic_only',
        help='use only genic/exonic reads and excluse flags for intergenic,intronic',
        type=bool,
        default=True
    )
    parser.add_argument(
        '--download_species',
        help='download hg38 and mm10 genomes, annotations, STAR indices, and xengsort index',
        required=False,
        default=False,
        type=bool
    )
    parser.add_argument(
        '--barcode',type=str,help='barcode file used to demultiplex samples in digitial gene expression. default is visium',
        default='barcodes/visium-v2.txt',
        required=False
    )
    return parser
def get_spatial_parser(requred=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='spatial chemistry parser',
        add_help=False
    )
    parser.add_argument(
        '--name',type=str,help='name of run mode used option',required=True
    )
    parser.add_argument(
        '--umi',type = str,help='umi structure using python list syntax'
        + "13-20 NT of Read1, use --umi r1[12:20]. It is also possible to use the first 8nt of "
        + "Read2 as UMI: --umi r2[0:8]",
        required = True)
    parser.add_argument(
        '--beads',type=int,help='expeted number of spots/beads for run'
    )
    parser.add_argument(
        '--spot_diameter_um',type =float,required=False,help='diameter of spot in spatial array: visium = 55um'
    )
    parser.add_argument(
        '--spot_diameter_um',type=float,required=False,help='distance between spots in um, visium = 100um'
    )
    parser.add_argument(
        '--cell_barcode',
        help = 'structure of CELL BARCODE, '
    )
    return parser
def get_file_parser(required=True):
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
    return parser

def setup_config_parser(config,parent_parser):
    parser_config = parent_parser.add_parser(
        'config',help='configure xenomake run'
    )
    parser_config_subparsers = parser_config.add_subparsers(
        help = 'config sub-command help'
    )
    return parser_config

class ConfigMainVariable:
    def __init__(self,name,**kwargs):
        self.name = name
        self.variables = {}

        for variable_key, variable in kwargs.items():
            if variable_key not in self.variable_types:
                raise UnrecognizedConfigVariable(
                    f"{variable_key}",list(self.variable_types.keys())
                )
            elif self.variable_types[variable_key] == 'int_list':
                self.variables[variable_key] = [int(x) for x in kwargs[variable_key]]
            else:
                new_type = self.variable_types[variable_key]
                self.variables[variable_key] = new_type(kwargs[variable_key])
        def __str__(self):
            class_name = self.__class__.__name__
            return f"{class_name}: {self.name}. variables: {self.variables}"
        
        def update(self,other):
            self.variables.update(other.variables)

class RunMode(ConfigMainVariable):
    variable_types = {
        'genic_only': bool,
        'downstream': bool,
        'ambiguous': bool,
        'mm_reads': bool
    }

class Spatial_Setup(ConfigMainVariable):
    variable_types = {'barcodes':str,'spot_diameter_um':float,'width_um': int}

    @property
    def has_barcodes(self):
        return (
            "barcodes" in self.variables
            and self.variables["barcodes"]
            and self.variables["barcodes"] != "None"
        )
    

class ConfigFile:
    initial_path = os.path.join(os.path.dirname(__file__),'data/config/config.yaml')

    main_variables_pl2sg = {
        'species': 'species',
        'run_modes': 'run_mode',
        'spatial_setups': 'spatial_setup'
    }

    main_variables_sg2pl = {value: key for key,value in main_variables_pl2sg.items()}

    main_variable_sg2type = {
        'species': str,
        'run_mode': 'str_list',
        'spatial_setup': str,
        'file': str
    }
    def __init__(self):
        self.variables = {
            'root_dir': "",
            'temp_dir': '/tmp',
            'dropseq_tools': 'tools/Drop-seq_tools-2.5.3',
            'picard': 'tools/picard.jar',
            'species': {},
            'run_modes': {},
            'spatial_setups': {},
        }
        self.file_path = 'config.yaml'

    @classmethod
    def from_yaml(cls,file_path):
        cf=cls()

        config_yaml_variables = None
        with open(file_path,'r') as f:
            config_yaml_variables = yaml.load(f,Loader=yaml.FullLoader)

        if config_yaml_variables is not None:
            cf.variable.update(config_yaml_variables)

        cf.file_path = file_path

        if file_path != cf.initial_config_path:
            initial_config = ConfigFile.from_yaml(cf.initial_path)
            cf.correct()

            for main_variable in cf.main_variables_pl2sg:
                if main_variable not in cf.variables:
                    cf.variables[main_variable] = initial_config.variables[main_variable]
            
            cf.vars_with_default = [
                key
                for key,value in initial_config.variables.items()
                if 'default' in value
            ]
            for var_with_default in cf.vars_with_default:
                default_val = initial_config.variables[var_with_default]['default']
                if 'default' not in cf.variables[var_with_default]:
                    cf.variables[var_with_default]['default'] = default_val
                else:
                    default_val.update(cf.variables[var_with_default]['default'])
                    cf.variables[var_with_default]['default'] = default_val
        return cf
    
    def dump(self):
        with open(self.file_path,'w') as f:
            f.write(yaml.dump(self.variables))
    
    def check_main_variable(self,variable):
        if variable not in self.main_variables_pl2sg.keys():
            raise UnrecognizedConfigVariable(
                variable,list(self.main_variables_pl2sg.keys())
            )
        
    def spaital_data(self):
        return self.variables['spatial_data']
    
    def check_spatial_data(self,cell_barcode = None,umi=None,name=None):
        if umi is not None:

        
    def get_variable(self,variable,name):
        if not self.variable_exists(variable,name):
            raise ConfigVariableNotFoundError(variable,name)
        else:
            return self.variables[variable][name]
    
    def set_file_path(self,file_path):
        self.file_path = file_path
    
    @property
    def spatial_data(self):
        return self.variables['spatial_data']
    
    def variable_exists(self,variable_name,variable_key):
        return variable_key in self.variables[variable_name]

    def assert_variable(self,variable_name,variable_key):
        self.assert_main_variable(variable_name)
        if not isinstance(variable_key,list):
            variable_key = [variable_key]

        for key in variable_key:
            if not self.variable_exiests(variable_name,key):
                variable_name_sg = self.main_variables_pl2sg[variable_name]
                raise ConfigVariableNotFoundError(variable_name_sg,key)

    def process_run_mode_args(self,**kwargs):
        default_run_mode = self.get_variable('run_modes','default')
        for key,value in default_run_mode.items():
            if isinstance(value,bool) and key in kwargs.keys():
                kwargs[key] = str2bool(kwargs[key])

    def process_spatial_args(self,**kwargs):
        
    def process_species_args(self,
                         mouse_ref=None,
                         human_ref=None,
                         mouse_annotation=None,
                         human_annotation=None,
                         human_index=None,
                         mouse_index=None,
                         xengsort_hash=None,
                         xengsort_info=None):
        assert_file(mouse_ref,default_value=None,extension=[".fa",".fna",".fa.gz",".fna.gz"])
        assert_file(human_ref,default_value=None,extension=[".fa",".fna",".fa.gz",".fna.gz"])

        assert_file(mouse_annotation,default_value=None,extension=[".gtf",".gtf.gz"])
        assert_file(human_annotation,default_value=None,extension=[".gtf",".gtf.gz"])

        if human_index:
            check_index(human_index)
        if mouse_index:
            check_index(mouse_index)
        if xengsort_hash:
            assert_file(xengsort_hash,default_value=None,extension=['.hash'])
        if xengsort_info:
            assert_file(xengsort_info,default_value=None,extension=['.info'])

        d = dict()
        return(d)

