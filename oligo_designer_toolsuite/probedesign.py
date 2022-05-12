############################################
# imports
############################################

import time
import argparse
import logging
from datetime import datetime

timestamp = datetime.now()
logging.basicConfig(filename='log_probe_design_{}-{}-{}-{}-{}.txt'.format(timestamp.year, timestamp.month, timestamp.day, timestamp.hour, timestamp.minute), level=logging.NOTSET)

from oligo_designer_toolsuite.utils import get_config, print_config, rm_intermediate_files
from oligo_designer_toolsuite.pipeline import download_and_filter, find_non_overlapping_sets, design_padlock_probes



############################################
# functions
############################################

def args():
    '''Argument parser for command line arguments.

    :return: Command line arguments with their values.
    :rtype: dict
    '''
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('-c','--config',help='path to config yaml file',type=str,required=True)
    args_parser.add_argument('-o','--output',help='path of output folder',type=str,required=True)
    args_parser.add_argument('-d','--download_only',help='dry run that only downloads the necessary files',action='store_true',required=False)
    return args_parser.parse_args()

  
def main():
    '''Main function of probe designer.
    '''
    # get comman line arguments
    parameters = args()

    dir_output = parameters.output
    config = get_config(parameters.config)
    download_only = parameters.download_only
    print_config(config, logging)

    logging.info('#########Start Pipeline#########')
    print('#########Start Pipeline#########')

    t_pipeline = time.time()
    
    download_and_filter(config, dir_output, logging, download_only=download_only)
    find_non_overlapping_sets(config, dir_output, logging)
    design_padlock_probes(config, dir_output, logging)
    
    t_pipeline = (time.time() - t_pipeline)/60

    logging.info('Time Pipeline: {} min'.format(t_pipeline))
    logging.info('#########End Pipeline#########')

    print('Time Pipeline: {} min \n'.format(t_pipeline))
    print('#########End Pipeline#########')

    #rm_intermediate_files(dir_output)
