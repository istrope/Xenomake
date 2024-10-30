import os

from xenomake.errors import *
#FUNCTION TAKEN FROM SPACEMAKE REPOSITORY https://github.com/rajewsky-lab/spacemake AUTHOR: Rajewsky-Lab
def check_index(star_index):
    star_version = os.popen("STAR --version").read().strip()

    with open(os.path.join(star_index, "Log.out"), "r") as f:
        first_line = f.readline().strip()
        index_version = first_line.split("=")[-1].split("_")[-1]
        
        if index_version != star_version:
             raise Exception(
                  f"STAR index version ({index_version}) is" + f" incompatible with your STAR version ({star_version})"
                  )
        
        

def assert_file(file_path,default_value=None,extension=['all']):
    if file_path == default_value:
        # file doesn't exist but has the default value,
        # so we do not need to assert anything
        return False
    
    if isinstance(extension, str):
        extension = [extension]

    if not isinstance(file_path, list):
        file_path = [file_path]

    for fp in file_path:
        # check if file exists, raise error if not
        if not os.path.isfile(fp):
            raise FileNotFoundError(fp)

        # check for all extensions
        if extension != ["all"]:
            extension_check = [fp.endswith(ex) for ex in extension]
            if not any(extension_check):
                raise FileWrongExtensionError(fp, extension)

    # return true if file exists and every test was good
    return True

def str2bool(var):
    bool_in_str = ["True", "true", "False", "false"]
    if isinstance(var,bool):
        return var
    if var in ['True','true']:
        return True
    elif var in ['False','false']:
        return False
    else:
        raise ValueError(f'variable should be boolean or one of {bool_in_str}')

def match_default(key,cf):
    return cf.variables[key]


def load_yaml_config(yaml_file):
    """Loads the YAML configuration file."""
    with open(yaml_file, 'r') as file:
        config = yaml.safe_load(file)
    return config