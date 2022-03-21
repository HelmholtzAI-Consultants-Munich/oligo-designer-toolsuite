import os
import random
import itertools
from pathlib import Path
import yaml
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
import src.utils as utils

def design_padlocks(config,dir_in_probes,dir_in_probesets,dir_out):
    """Design final padlock probe sequences
    
    config: dict
        Configuration dict    
    
    Arguments
    ---------
    config: dict
        Configuration dictionary. Here only needed for detection oligo design
    dir_in_probes: str
        Directory where probes_<gene>.txt files are located
    dir_in_probesets: str
        Directory where ranked_probesets_<gene>.txt files are located
    dir_out: str
        Results directory
    
        
    Saves table "padlock_probes_<gene>.txt" with final padlock probe sequences, detection oligo sequences, and infos
    ("detection_oligo_<gene>.txt" <-> not needed, let's save detection_oligo in padlock_probes..txt)
    
    """

    Path(dir_out).mkdir(parents=True, exist_ok=True)
    
    #config = get_config(config)
        
    probes_files = [f for f in os.listdir(dir_in_probes) if f.startswith("probes_")]
    genes = [f.split("_")[1].split(".")[0] for f in probes_files]
    probeset_files = [f"ranked_probesets_{gene}.txt" for gene in genes]
    
    yaml_dict = {}
    
    for gene_idx, (gene, probe_f, probeset_f) in enumerate(zip(genes,probes_files,probeset_files)):
        yaml_dict[gene] = {}
        
        probes = pd.read_csv(os.path.join(dir_in_probes,probe_f),index_col=0,sep="\t")
        df_probeset = pd.read_csv(os.path.join(dir_in_probesets,probeset_f),index_col=0,sep="\t")
        probeset = [df_probeset[col].iloc[0] for col in df_probeset.columns if col.startswith("probe")]
        
        for probe_idx, probe in enumerate(probeset):
            complementary_seq = probes.loc[probe,"probe_sequence"]
            ligation_idx = probes.loc[probe,"ligation_site"]
            #print(gene, gene_idx, probe_idx, "### padlock probe")
            full_seq, sub_seqs = get_padlock_probe(gene_idx,complementary_seq,ligation_idx,barcode_seed=0,barcode_length=4)
            #print(gene, gene_idx, probe_idx, "### detection oligo", complementary_seq)
            det_oligo_seq, det_oligo_Tm = get_detection_oligo(complementary_seq,ligation_idx,config,max_no_U_length=9)
            #print(det_oligo_seq, type(det_oligo_seq), det_oligo_Tm)
            
            yaml_dict[gene][f"{gene}_probe{probe_idx+1}"] = {}
        
            # potentially nice to also save organism, full gene name, reference genome
            yaml_dict[gene][f"{gene}_probe{probe_idx+1}"]["probe_id"] = str(probe)
            for key in ["gene_id", "transcript_id", "exon_id", "chromosome", "start", "end", "strand"]:
                yaml_dict[gene][f"{gene}_probe{probe_idx+1}"][key] = str(probes.loc[probe,key])
            yaml_dict[gene][f"{gene}_probe{probe_idx+1}"].update({
                "padlock_probe_full_sequence" : str(full_seq),
                "detection_oligo_sequence"    : str(det_oligo_seq + "[fluorophore]"),
                "padlock_arm1_sequence"       : str(sub_seqs["arm1"]),
                "padlock_accessory1_sequence" : str(sub_seqs["accessory1"]),
                "padlock_ISS_anchor_sequence" : str(sub_seqs["ISS_anchor"]),
                "padlock_barcode_sequence"    : str(sub_seqs["barcode"]),
                "padlock_accessory2_sequence" : str(sub_seqs["accessory1"]),
                "padlock_arm2_sequence"       : str(sub_seqs["arm2"]),
                "complementary_sequence"      : str(complementary_seq),
                "target_mRNA_sequence"        : str(complementary_seq[::-1]),
            })
            for key in ["GC_content", "melting_temperature", "melt_temp_arm1", "melt_temp_arm2", "melt_temp_dif_arms"]:
                yaml_dict[gene][f"{gene}_probe{probe_idx+1}"][key] = float(probes.loc[probe,key])
            for key in ["length", "ligation_site"]:
                yaml_dict[gene][f"{gene}_probe{probe_idx+1}"][key] = int(probes.loc[probe,key])
            yaml_dict[gene][f"{gene}_probe{probe_idx+1}"]["melt_temp_detection_oligo"] = float(det_oligo_Tm)
            
    with open(os.path.join(dir_out,'padlock_probes.yml'), 'w') as outfile:
        yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)
        
    yaml_order = {}
    for gene in yaml_dict:
        yaml_order[gene] = {}
        for probe_id in yaml_dict[gene]:
            yaml_order[gene][probe_id] = {}
            yaml_order[gene][probe_id]["padlock_probe_full_sequence"] = yaml_dict[gene][probe_id]["padlock_probe_full_sequence"]
            yaml_order[gene][probe_id]["detection_oligo_sequence"] = yaml_dict[gene][probe_id]["detection_oligo_sequence"]
            
    with open(os.path.join(dir_out,'padlock_probes_order.yml'), 'w') as outfile:
        yaml.dump(yaml_order, outfile, default_flow_style=False, sort_keys=False)


def get_barcode(gene_idx,length=4,seed=0):
    """Get barcode sub sequence of padlock probe for in situ sequencing
    
    For SCRINSHOT padlock probes this could be constant, however it makes sense to have
    different barcodes so that the probe set could also be used for ISS experiments.
    
    Arguments
    ---------
    gene_idx: int
        Identifier for a given gene. The identifier makes sure to return the same bar code
        for the different padlock probes of a given gene.
    length: int
        Length of barcode sequence
    seed: int
        Defines the random assignment of barcodes to each gene_idx.
        
    Returns
    -------
    str: barcode sequence (5' to 3')
    
    """
    bases = ["A","C","T","G"]
    
    #barcodes = list(itertools.product(bases, repeat=length))
    barcodes = ["".join(nts) for nts in itertools.product(bases, repeat=length)]
    random.seed(seed)
    random.shuffle(barcodes)
    
    if gene_idx >= len(barcodes):
        raise ValueError("Barcode index exceeds number of possible combinations of barcodes. Increase barcode length?")
    
    return barcodes[gene_idx]
    
    
def SCRINSHOT_or_ISS_backbone_sequence(gene_idx,barcode_length=4,barcode_seed=0):
    """Get backbone sequence of padblock probes for SCRINSHOT or ISS
    
    Arguments
    ---------
    gene_idx: int
        Identifier for a given gene. The identifier makes sure to return the same bar code
        for the different padlock probes of a given gene.
    barcode_length: int
        Length of barcode sequence
    barcode_seed: int
        Defines the random assignment of barcodes to each gene_idx.
        
    Returns
    -------
    str: backbone sequence (5' to 3')
    dict of strs:
        Individual parts of the backbone sequence
        
    """
    accessory1 = "TCCTCTATGATTACTGAC" 
    ISS_anchor = "TGCGTCTATTTAGTGGAGCC"
    barcode = get_barcode(gene_idx,length=barcode_length,seed=barcode_seed)
    accessory2 = "CTATCTTCTTT"
    
    sub_seqs = {"accessory1":accessory1 ,"ISS_anchor":ISS_anchor ,"barcode":barcode ,"accessory2":accessory2}
    full_seq = accessory1 + ISS_anchor + barcode + accessory2
    
    return full_seq, sub_seqs

def convert_complementary_seq_to_5to3(complementary_seq, ligation_idx):
    """Convert the complementary sequence of padlock probes to two arms with 5' to 3' convention
    
    E.g.
    complementary_seq = "AAAATGCTTAAGC" ligation_idx = 7
    --> cut after 7th base: "AAAATGC|TTAAGC"
    --> final result: ["CGTAAAA","TTAAGC"]
    
    Arguments
    ---------
    complementary_seq: str
        Sequence that hybridises with the target RNA
    ligation_idx: int
        Site where complementary_seq is cut in two arms according padlock probe design
        
    Returns
    -------
    list of strs: first on second arm sequences (both 5' to 3')
    
    """
    
    arm1 = complementary_seq[:ligation_idx][::-1]
    arm2 = complementary_seq[ligation_idx:]
    
    return [arm1,arm2]
    
def get_padlock_probe(gene_idx,complementary_seq,ligation_idx,barcode_seed=0,barcode_length=4):
    """Get full padlock probe for a given gene and complementary sequence
    
    Arguments
    ---------
    gene_idx: int
        Identifier for a given gene. The identifier makes sure to return the same bar code
        for the different padlock probes of a given gene.    
    complementary_seq: str
        Sequence that hybridises with the target RNA
    ligation_idx: int
        Site where complementary_seq is cut in two arms according padlock probe design  
    barcode_length: int
        Length of barcode sequence
    barcode_seed: int
        Defines the random assignment of barcodes to each gene_idx.        
    
    Returns
    -------
    str: padlock probe sequence (5' to 3')
    dict of strs:
        Individual parts of the padlock sequence
        
    """
    arms = convert_complementary_seq_to_5to3(complementary_seq, ligation_idx)
    sub_seqs = {"arm1":arms[0]}
    
    backbone_seq, backbone_sub_seqs = SCRINSHOT_or_ISS_backbone_sequence(gene_idx,barcode_seed=barcode_seed,barcode_length=barcode_length)
    
    sub_seqs.update(backbone_sub_seqs)
    sub_seqs.update({"arm2":arms[1]})
    
    full_seq = arms[0] + backbone_seq + arms[1]
    
    return full_seq, sub_seqs



def Ts_are_close_enough(probe,max_no_T_length=9):
    """Check if sequence has maximally `max_no_T_length` nucleotides without T in a row
    
    E.g. (max_no_T_length=9):
        "AAAAAAAAAAAAAAATAAAA" -> Error
        "TAAAAAAAAAAAAAAAAAAT" -> Error
        "AAAATAAAAATAAAATAAAA" -> no Error
    
    """
    if not "T" in probe:
        return False    
    for idx in range(len(probe)):
        if probe[idx:].find("T") > max_no_T_length:
            return False
    for idx in range(len(probe)):
        if probe[idx::-1].find("T") > max_no_T_length:
            return False        
    return True
            

def exchange_T_with_U(probe,max_no_U_length=9):
    """Exchange minimal number of T(hymines) with U(racils) in probe
    
    Arguments
    ---------
    probe: str
        Sequence
    max_no_U_length: int
        Maximal number of nucleotides in a row without U
        
    Returns
    -------
    str:
        Sequence with the minimal number of exchanged Us
    
    """
    
    if not Ts_are_close_enough(probe,max_no_T_length=max_no_U_length):
        return "NOT-ENOUGH-THYMINES-FOR-DETECTION-OLIGO"
        raise ValueError(f"Sequence {probe} has subsequences with more than {max_no_U_length} nucleotides without T.")
    
    probe_length = len(probe)
    
    idx = 0
    idxs = []
    while idx < (probe_length - 1 - max_no_U_length):
        idxs.append(idx + probe[idx:idx+max_no_U_length+int(idx>0)+1].rfind("T"))
        idx = idxs[-1]
    probe_option1 = "".join(["U" if i in idxs else nt for i,nt in enumerate(probe)])
    
    probe_rev = probe[::-1]
    idx = 0
    idxs = []
    while idx < (probe_length - 1 - max_no_U_length):
        idxs.append(idx + probe_rev[idx:idx+max_no_U_length+int(idx>0)+1].rfind("T"))
        idx = idxs[-1]        
    probe_option2 = "".join(["U" if i in idxs else nt for i,nt in enumerate(probe_rev)])[::-1]
    
    #print(probe_option1, probe_option1.count('U'))
    #print(probe_option2, probe_option2.count('U'))
    
    if probe_option1.count('U') >= probe_option2.count('U'):
        return probe_option1
    else:
        return probe_option2


def _get_oligo_Tm(probe_sequence,Tm_parameters,Tm_correction_parameters):
    """Compute the melting temperature for the detection oligo sequence
    
    #TODO: this function is not placed very nicely. Also it uses utils.get_Tm_parameters everytime we search a new oligo
           is calculated. In datamodule.py we have the same function for the probe sequence instead of the oligo. Idk
           atm, it could be handled nicer. Not rly important atm since we only compute a handful of detection oligos.
    
    In the first step the melting temperature is calculated based on sequence and salt concentrations. In the
    second step the temperature is corrected by formamide (and DMSO) percentage.
    
    :param probe_sequence: Sequence of probe
    :type probe_sequence: string
    """
    
    Tm = mt.Tm_NN(probe_sequence, **Tm_parameters)
    Tm_corrected = round(mt.chem_correction(Tm, **Tm_correction_parameters),2)
    return Tm_corrected


def get_detection_oligo(probe_sequence,ligation_site,config,max_no_U_length=9):
    """Get detection oligo sequence for a given probe
    
    Detection oligos have the same sequence as the complementary sequence (i.e. `probe_sequence`) but shortend and 
    reversed. The ligation site is placed in the middle of the detection oligo. The detection oligo is shortend to get
    its melting temperature as close as possible to config["detect_oligo_Tm_opt"] but not shorter than 
    config["detect_oligo_length_min"].
    
    
    Arguments
    ---------
    probe_sequence: str
        The sequence of the complementary probe
    ligation_site: int
        Ligation site
    config: dict
        Configuration for oligo design. 
    max_no_U_length: int
        Maximal number of nucleotides in a row without U(racil) in the final detection oligo
    
    Returns
    -------
    str:
        Sequence of detection oligo
    float:
        Melting temperature of detection oligo
    
    """
    
    probe_length = len(probe_sequence)
    oligo_length_max = max(ligation_site,probe_length-ligation_site)
    oligo_length_min = config["detect_oligo_length_min"]
    
    Tm_parameters = utils.get_Tm_parameters(config['Tm_parameters'], sequence='detection_oligo')
    Tm_correction_parameters = utils.get_Tm_correction_parameters(
        config['Tm_correction_parameters'], 
        sequence='detection_oligo'
    )
    optimal_Tm = config["detect_oligo_Tm_opt"]
    get_Tm_dif = lambda seq: abs(_get_oligo_Tm(seq,Tm_parameters,Tm_correction_parameters) - optimal_Tm)
    
    if ligation_site > (probe_length-ligation_site):
        best_oligo = probe_sequence[probe_length-oligo_length_max:]
    else:
        best_oligo = probe_sequence[:oligo_length_max]
    best_Tm_dif = get_Tm_dif(best_oligo)
        
    # Find detection oligo with best melting temperature
    oligo = best_oligo
    oligo_length = len(oligo)
    Tm_dif = best_Tm_dif
    
    count = 0
    while oligo_length > oligo_length_min:
        
        if bool(count % 2):
            oligo = oligo[1:]
        else:
            oligo = oligo[:-1]
            
        Tm_dif = get_Tm_dif(oligo)
        
        #TODO: if there are no sequences with proper distances between Ts (very unlikely) than we run into an error
        #      this is really bad since the pipeline will be interrupted. 
        if (Tm_dif < best_Tm_dif) and Ts_are_close_enough(oligo,max_no_T_length=max_no_U_length):
            best_Tm_dif = Tm_dif
            best_oligo = oligo
            
        count += 1
        oligo_length = len(oligo)
    
    # exchange T's with U (for enzymatic degradation of oligos)
    oligo_seq = exchange_T_with_U(best_oligo,max_no_U_length=max_no_U_length)
    if oligo_seq != "NOT-ENOUGH-THYMINES-FOR-DETECTION-OLIGO":
        oligo_Tm = _get_oligo_Tm(best_oligo,Tm_parameters,Tm_correction_parameters)
    else:
        oligo_Tm = 0
    
    return oligo_seq, oligo_Tm
        
    
    