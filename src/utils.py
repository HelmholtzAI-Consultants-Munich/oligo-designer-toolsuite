############################################
# imports
############################################

import yaml

import subprocess as sp


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
    
def write_output(dir_output, mapping_gene_to_probes_blastn):

    dir_output = '{}results/'.format(dir_output)
    cmd = 'mkdir {}'.format(dir_output)
    sp.run(cmd, shell=True)

    # write results gene-wise to csv file
    for gene_id, mapping_probe_to_exon in mapping_gene_to_probes_blastn.items():
        file_output = '{}probes_{}'.format(dir_output, gene_id)

        with open(file_output, 'w') as handle:
            handle.write('gene_id\texon_ids\tprobe_sequence\tGC_content\tmelting_temperature\n')
            
            for probe, probe_attributes in mapping_probe_to_exon.items():
                exons = ';'.join(probe_attributes['exon_id'])
                output = '{}\t{}\t{}\t{}\t{}\n'.format(gene_id, exons, probe, round(probe_attributes['GC'],2), round(probe_attributes['Tm'],2))
                handle.write(output)


############################################