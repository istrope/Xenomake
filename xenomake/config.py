import argparse
import os
import yaml
from xenomake.utils import str2bool
from xenomake.errors import *

def setup_spatial_parser(parent_parser):
    parser_spatial = parent_parser.add_parser(
        'spatial',
        help = 'create custom method for spatial technology'
    )
    parser_spatial.add_argument(
        '--ambiguous',
        help = 'handle multimap from "both" and "ambiguous" output reads (Must be boolean True/False)',
        required=False,
        default=True,
        type=bool
    )
    parser_spatial.add_argument(
        '--mm_reads',
        help='assign multimapped reads to highest confidence position',
        required=False,
        default=True,
        type=bool
    )
    parser_spatial.add_argument(
        '--downstream',
        help = 'Performs Scanpy Processing of Data (Must be boolean True/False)',
        required = False,
        default = True,
        type = bool
    )
    parser_spatial.add_argument(
        '--genic_only',
        help='use only genic/exonic reads and excluse flags for intergenic,intronic',
        type=bool,
        default=True,
        required = False
    )
    parser_spatial.add_argument(
        '--barcode',type=str,help='barcode file used to demultiplex samples in digitial gene expression. default is visium',
        default='barcodes/visium-v2.txt',
        required=False
    )
    parser_spatial.add_argument(
        '--umi',
        type = str,
        help='umi structure defining the positions of umi in units of bases'
        +'Example: Visium Cell Barcode is contained in bases 17-28 and umi flag would be'
        + '--umi 17-28',
        required = True
        )
    parser_spatial.add_argument(
        '--beads',type=int,help='expeted number of spots/beads for run'
    )
    parser_spatial.add_argument(
        '--spot_diameter_um',type =float,required=False,help='diameter of spot in spatial array: visium = 55um'
    )
    parser_spatial.add_argument(
        '--spot_diameter_um',type=float,required=False,help='distance between spots in um, visium = 100um'
    )
    parser_spatial.add_argument(
        '--cell_barcode',
        required=True,
        help = 'cell barcode structure defining the positions of umi in units of bases'
        +'Example: Visium Cell Barcode is contained in bases 1-16 and umi flag would be'
        + '--cell_barcode 1-16'
    )
    return parser_spatial

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
    mapping_types = {'cell_barcode':str,'umi':str}

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
            'run_mode': {},
            'spatial_setup': {},
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

        return kwargs
    

    def process_spatial_data(self,cell_barcode = None,umi=None,name=None):
        bam_tags = "CR:{cell},CB:{cell},MI:{UMI},RG:{assigned}"
        if umi is not None:
            self.variables['umi'] = umi
        if cell_barcode is not None:
            self.variables['cell'] = cell_barcode
        if name is not None:
            self.variables[name] = name

        barcodes = {'bam_tags': bam_tags}

        self.variables['alignment_tags'] = bam_tags
        return barcodes
    