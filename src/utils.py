############################################
# imports
############################################

import yaml


############################################
# helper functions
############################################

def get_config(config):
    """
    Loads config file in yaml format. 
    Parameters
    ----------
        config: string
            File path to yaml config file.
    Returns
    -------
        yaml: dict
            User-defined parameters, where keys are the parameter names and values are the paremeter values.
    """
    with open(config, 'r') as ymlfile:
        return yaml.safe_load(ymlfile)


############################################
    
def print_config(config, logging):
    """
    Logs formatted config parameters as <parameter_name>: <parameter_value>. 
    Parameters
    ----------
        config: string
            File path to yaml config file.
        logging: logger
            Logger object to store important information.
    Returns
    -------
        --- none ---
    """
    logging.info('\nParameter settings:')
    for item, value in config.items(): 
        logging.info("{}: {}".format(item, value))

    
############################################
    