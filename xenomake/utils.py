import sys
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
    

def match_default(run_mode):
    defaults = ['visium','slide_seq','hdst','stereo_seq','pixel_seq','dbit_seq']
    if run_mode not in defaults:
             raise XenomakeError(f'current run mode is {run_mode},please choose spatial chemistry matching one of the following {defaults}')
    
    visium = {
        'barcodes': 'barcodes/visium_barcodes.txt',
        'spot_diameter_um': 55,
        'width_um': 6500,
        'umi':17-28,
        'cell': 1-16,
        'polyA_trimming':False
        }
    slide_seq = {
        'spot_diameter_um': 10,
        'width_um': 3000,
        'umi': 15-23,
        'cell': 1-14,
        'polyA_trimming':True
        }
    hdst = {}
    stereo_seq = {}
    pixel_seq = {}
    slide_seq = {}
    dbit_seq = {
        'barcodes': 'barcodes/dbit_seq_barcodes.txt',
        'spot_diameter_um': 10,
        'width_um':1000,
        'umi': 23-32,
        'cell': 33-40,
        'polyA_trimming':False
    }
    methods = {'visium': visium , 'slide_seq':slide_seq , 'hdst':hdst,
               'stereo_seq':stereo_seq , 'pixel_seq':pixel_seq , 'dbit_seq':dbit_seq}
    return methods[run_mode]
    


    
