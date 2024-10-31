import os
import yaml
import xenomake
from xenomake.utils import str2bool
from xenomake.errors import *

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
        'mm_reads': bool,
        'polyA_trimming':bool
    }

class Spatial_Setup(ConfigMainVariable):
    variable_types = {'barcode_file':str,
                      'cell_barcode':str,
                      'umi':str}
    @property
    def has_barcodes(self):
        return (
            "barcodes" in self.variables
            and self.variables["barcodes"]
            and self.variables["barcodes"] != "None"
        )

class ConfigFile:
    initial_path = os.path.join(xenomake.__path__[0],'data/config_structure.yaml')

    main_variables_pl2sg = {
        'project',
        'run',
        'spatial'
    }

    #main_variables_sg2pl = {value: key for key,value in main_variables_pl2sg.items()}

    main_variable_sg2type = {
        'project':str,
        'run': str,
        'spatial':str
    }
    def __init__(self):
        self.variables = {
            'project':{},
            'run': {},
            'spatial': {}
        }
        self.file_path = 'config.yaml'

    @classmethod
    def from_yaml(cls,file_path):
        cf=cls()

        config_yaml_variables = None
        with open(file_path,'r') as f:
            config_yaml_variables = yaml.load(f,Loader=yaml.FullLoader)

        if config_yaml_variables is not None:
            cf.variables.update(config_yaml_variables)

        cf.file_path = file_path

        if file_path != cf.initial_path:
            initial_config = ConfigFile.from_yaml(cf.initial_path)
            #cf.correct()

            for main_variable in cf.main_variables_pl2sg:
                if main_variable not in cf.variables:
                    cf.variables[main_variable] = initial_config.variables[main_variable]
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

    def check_project(self):
        defaults = {
            'run_mode': {'lenient','prude','custom'},
            'spatial_mode': {'visium','seq-scope','dbit-seq','hdst','sc10x','slide-seq','custom'}

        }
        project_vars = self.variables['project']
        
        for key,value in project_vars.items():
            if (value  == '') or (value == None):
                raise EmptyConfigVariableError(key)
            
        if project_vars['spatial_mode'] not in defaults['spatial_mode']:
            raise UnrecognizedConfigVariable('spatial_mode')
        
        if project_vars['run_mode'] not in defaults['run_mode']:
            raise UnrecognizedConfigVariable('run_mode')
        
