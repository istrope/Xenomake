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
        
        

def validate_directory(directory_path):
    """Check if the directory exists."""
    if not os.path.isdir(directory_path):
        raise FileNotFoundError(directory_path)

def validate_mode(mode, valid_modes, mode_type):
    """Validate if the mode is within the predefined set of valid modes."""
    if mode not in valid_modes:
        raise UnrecognizedConfigVariable(mode_type, mode)

def validate_init_inputs(init_args):
    VALID_RUN_MODES = ['prude', 'lenient', 'custom']
    VALID_SPATIAL_MODES = ['custom', 'visium', 'slide-seq', 'dbit-seq', 'hdst', 'sc10x', 'seq-scope']
    """
    Validate all input files, directories, and modes in the `xenomake init` command.

    :param init_args: Dictionary with keys: 'read1', 'read2', 'spatial_mode', 'run_mode', 'sample', 'outdir', 'root_dir', 'temp_dir'.
    """
    # Validate input files (read1 and read2) for existence and correct extension
    assert_file(init_args['read1'],default_value=None,extension=['.fq','.fq.gz','.fastq','.fastq.gz'])
    assert_file(init_args['read2'],default_value = None,extension = ['.fq','.fq.gz','.fastq','.fastq.gz'])

    # Validate directories
    if 'root_dir' in init_args:
        validate_directory(init_args['root_dir'])
    if 'temp_dir' in init_args:
        validate_directory(init_args['temp_dir'])

    # Validate modes
    validate_mode(init_args['spatial_mode'], VALID_SPATIAL_MODES, 'spatial_mode')
    validate_mode(init_args['run_mode'], VALID_RUN_MODES, 'run_mode')

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