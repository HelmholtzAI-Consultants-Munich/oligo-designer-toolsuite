############################################
# imports
############################################

import os
import yaml
import random
import logging
import itertools

from pathlib import Path

import pandas as pd

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

import oligo_designer_toolsuite.utils as utils

############################################
# probe set generator class
############################################

class ProbeSequenceDesigner:
    '''This class is used to design the final padlock sequences.
    '''

    def __init__(self, config, dir_output, dir_probes = None, dir_probesets=None):
        """Constructor method
        """
        # set logger
        self.logging = logging.getLogger('padlock_probe_designer')
    
        # set directory
        if dir_probes == None:
            self.dir_probes = os.path.join(dir_output, 'probes')
            Path(self.dir_probes).mkdir(parents=True, exist_ok=True)
        else:
            self.dir_probes = dir_probes

        if dir_probesets == None:
            self.dir_probesets = os.path.join(dir_output, 'probesets')
            Path(self.dir_probesets).mkdir(parents=True, exist_ok=True)
        else:
            self.dir_probesets = dir_probesets

        self.dir_padlock_probes = os.path.join(dir_output, 'padlock_probes')
        Path(self.dir_padlock_probes).mkdir(parents=True, exist_ok=True)
        
        # set parameters
        self.detect_oligo_length_min = config["detect_oligo_length_min"]
        self.detect_oligo_length_max = config["detect_oligo_length_max"]
        self.detect_oligo_Tm_opt = config["detect_oligo_Tm_opt"]
        self.Tm_parameters = utils.get_Tm_parameters(config['Tm_parameters'], sequence='detection_oligo')
        self.Tm_correction_parameters = utils.get_Tm_correction_parameters(config['Tm_correction_parameters'], sequence='detection_oligo')


    def design_padlocks(self):
        """Design final padlock probe sequences
        
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
        
        Saves
        -----
        - table at dir_out+"padlock_probes.yml" with final padlock probe sequences, detection oligo sequences, and infos
        - table at dir_out+"padlock_probes_order.yml" with final padlock probe sequences and detection oligo sequences
        
        """ 
        probeset_files = [f for f in os.listdir(self.dir_probesets) if f.startswith("ranked_probesets_")]
        genes = [f.split("_")[-1].split(".")[0] for f in probeset_files]
        probes_files = [f"probes_{gene_id}.txt" for gene_id in genes]
        
        yaml_dict = {}
        
        for gene_idx, (gene, probe_f, probeset_f) in enumerate(zip(genes, probes_files, probeset_files)):
            yaml_dict[gene] = {}
            
            probes = pd.read_csv(os.path.join(self.dir_probes, probe_f),index_col=0,sep="\t")
            table_probesets = pd.read_csv(os.path.join(self.dir_probesets, probeset_f),index_col=0,sep="\t")
            probeset_idx = self.best_probeset_with_possible_detection_oligos(table_probesets, probes, minT=2)
            probeset = [table_probesets[col].iloc[probeset_idx] for col in table_probesets.columns if col.startswith("probe")]
            
            for probe_idx, probe in enumerate(probeset):
                # NOTE: so far what we called "probe" is actually the sequence on the target mRNA, it was straightforward
                #       to adjust everything from here, but it would have been cleaner to take the reverse complement
                #       from the beginning. The previous ligation_site for example needed to be adjusted which might be 
                #       confusing when looking into files generated in previous steps. Also the arm Tm's could have 
                #       changed, but the naming (arm1, arm2) still fits and Tm(seq)==Tm(rev_compl(seq)) (properly tested).
                target_mRNA = probes.loc[probe,"probe_sequence"]
                complementary_seq = str(Seq(target_mRNA).reverse_complement())
                ligation_idx = len(target_mRNA) - probes.loc[probe,"ligation_site"]
                full_seq, sub_seqs = self.get_padlock_probe(gene_idx,complementary_seq,ligation_idx,barcode_seed=0,barcode_length=4)
                det_oligo_seq, det_oligo_Tm = self.get_detection_oligo(complementary_seq, ligation_idx, minT=2)
                
                yaml_dict[gene][f"{gene}_probe{probe_idx+1}"] = {}
            
                #TODO: potentially nice to also save organism, full gene name, reference genome
                yaml_dict[gene][f"{gene}_probe{probe_idx+1}"]["probe_id"] = str(probe)
                for key in ["gene_id", "transcript_id", "exon_id", "chromosome", "start", "end", "strand"]:
                    yaml_dict[gene][f"{gene}_probe{probe_idx+1}"][key] = str(probes.loc[probe,key])
                yaml_dict[gene][f"{gene}_probe{probe_idx+1}"].update({
                    "padlock_probe_full_sequence" : str(full_seq),
                    "detection_oligo_sequence"    : str(det_oligo_seq),
                    "padlock_arm1_sequence"       : str(sub_seqs["arm1"]),
                    "padlock_accessory1_sequence" : str(sub_seqs["accessory1"]),
                    "padlock_ISS_anchor_sequence" : str(sub_seqs["ISS_anchor"]),
                    "padlock_barcode_sequence"    : str(sub_seqs["barcode"]),
                    "padlock_accessory2_sequence" : str(sub_seqs["accessory1"]),
                    "padlock_arm2_sequence"       : str(sub_seqs["arm2"]),
                    "complementary_sequence"      : str(complementary_seq),
                    "recognised_mRNA_sequence"    : str(target_mRNA),#str(Seq(complementary_seq).complement()),
                })
                for key in ["GC_content", "melting_temperature", "melt_temp_arm1", "melt_temp_arm2", "melt_temp_dif_arms"]:
                    yaml_dict[gene][f"{gene}_probe{probe_idx+1}"][key] = float(probes.loc[probe,key])
                for key in ["length"]:#, "ligation_site"]:
                    yaml_dict[gene][f"{gene}_probe{probe_idx+1}"][key] = int(probes.loc[probe,key])
                yaml_dict[gene][f"{gene}_probe{probe_idx+1}"]["ligation_site"] = int(ligation_idx)
                yaml_dict[gene][f"{gene}_probe{probe_idx+1}"]["melt_temp_detection_oligo"] = float(det_oligo_Tm)
                
        with open(os.path.join(self.dir_padlock_probes,'padlock_probes.yml'), 'w') as outfile:
            yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)
            
        yaml_order = {}
        for gene in yaml_dict:
            yaml_order[gene] = {}
            for probe_id in yaml_dict[gene]:
                yaml_order[gene][probe_id] = {}
                yaml_order[gene][probe_id]["padlock_probe_full_sequence"] = yaml_dict[gene][probe_id]["padlock_probe_full_sequence"]
                yaml_order[gene][probe_id]["detection_oligo_sequence"] = yaml_dict[gene][probe_id]["detection_oligo_sequence"]
                
        with open(os.path.join(self.dir_padlock_probes,'padlock_probes_order.yml'), 'w') as outfile:
            yaml.dump(yaml_order, outfile, default_flow_style=False, sort_keys=False)


    def best_probeset_with_possible_detection_oligos(self, df_probeset, probes, minT=2):
        """Get row index of best probeset for which all detection oligos can be designed
        
        Criterions for detection oligo design are 1. a length constraint and 2. to have at least 2 T(hymines).
        
        TODO: The detection oligo thing is a bit annoying. Currently I don't want to set this as a constrained for
            finding probes for a gene. Otherwise we could directly filter probes before searching for probesets 
            which would be great. Maybe something to think about in the future. For now we only check after having 
            the probe sets per gene.
        
        Arguments
        ---------
        df_probeset: pd.DataFrame
            Dataframe with ranked probesets. Output of get_nonoverlapping_sets in nonoverlapping_sets.py.
        probes: pd.DataFrame
            Table with probe infos.
        minT: int
            Minimal number of T(hymines) in detection oligo.
            
        Returns
        -------
        int: 
            Numerical row index of first probeset that only contains probes for which a detection oligo can be designed.
            If no probeset can be found the first probeset index i.e. 0 is returned.
        
        """
        for probeset_idx in range(len(df_probeset)):
            probeset = [df_probeset[col].iloc[probeset_idx] for col in df_probeset.columns if col.startswith("probe")]
            for probe_idx, probe in enumerate(probeset):                
                target_mRNA = probes.loc[probe,"probe_sequence"]
                complementary_seq = str(Seq(target_mRNA).reverse_complement())
                ligation_idx = len(target_mRNA) - probes.loc[probe,"ligation_site"]
                
                start_oligo, start_oligo_long_left, start_oligo_long_right = self._get_initial_oligos_for_search(complementary_seq, ligation_idx)
                if (start_oligo_long_left is not None) and (start_oligo_long_left.count("T") >= minT):
                    return probeset_idx
                elif (start_oligo_long_right is not None) and (start_oligo_long_right.count("T") >= minT):
                    return probeset_idx
                elif (start_oligo.count("T") >= minT):
                    return probeset_idx
            
        return 0
        
        
    def get_padlock_probe(self, gene_idx, complementary_seq, ligation_idx, barcode_seed=0, barcode_length=4):
        """Get full padlock probe for a given gene and a given complementary sequence
        
        Arguments
        ---------
        gene_idx: int
            Identifier for a given gene. The identifier makes sure to return the same bar code
            for the different padlock probes of a given gene.    
        complementary_seq: str
            Sequence that hybridises with the target RNA
        ligation_idx: int
            Site where complementary_seq is cut in two arms according padlock probe design  
        barcode_seed: int
            Defines the random assignment of barcodes to each gene_idx.        
        barcode_length: int
            Length of barcode sequence        
        
        Returns
        -------
        str: 
            padlock probe sequence (5' to 3')
        dict of strs:
            Individual parts of the padlock sequence
            
        """
        def _convert_complementary_seq_to_arms(complementary_seq, ligation_idx):
            """Convert the complementary sequence of padlock probes to two arms with 5' to 3' convention
            
            E.g.
            complementary_seq = "AAAATGCTTAAGC" ligation_idx = 7
            --> cut after 7th base: "AAAATGC|TTAAGC"
            --> final result: ["CGTAAAA","CGAATT"]
            
            Arguments
            ---------
            complementary_seq: str
                Sequence that hybridises with the target RNA
            ligation_idx: int
                Site where complementary_seq is cut in two arms according padlock probe design
                
            Returns
            -------
            list of strs: first and second arm sequences (both 5' to 3')
            
            """
            arm1 = complementary_seq[ligation_idx:]
            arm2 = complementary_seq[:ligation_idx]
            
            return [arm1,arm2]

        def _get_barcode(gene_idx,length=4,seed=0):
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
            
            barcodes = ["".join(nts) for nts in itertools.product(bases, repeat=length)]
            random.seed(seed)
            random.shuffle(barcodes)
            
            if gene_idx >= len(barcodes):
                raise ValueError("Barcode index exceeds number of possible combinations of barcodes. Increase barcode length?")
            
            return barcodes[gene_idx]

        def _SCRINSHOT_or_ISS_backbone_sequence(gene_idx, barcode_length=4, barcode_seed=0):
            """Get backbone sequence of padlock probes for SCRINSHOT or ISS
            
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
            str: 
                backbone sequence (5' to 3')
            dict of strs:
                Individual parts of the backbone sequence
                
            """
            accessory1 = "TCCTCTATGATTACTGAC" 
            ISS_anchor = "TGCGTCTATTTAGTGGAGCC"
            barcode = _get_barcode(gene_idx, length=barcode_length, seed=barcode_seed)
            accessory2 = "CTATCTTCTTT"
            
            sub_seqs = {"accessory1":accessory1 ,"ISS_anchor":ISS_anchor ,"barcode":barcode ,"accessory2":accessory2}
            full_seq = accessory1 + ISS_anchor + barcode + accessory2
            
            return full_seq, sub_seqs
        
        arms = _convert_complementary_seq_to_arms(complementary_seq, ligation_idx)
        sub_seqs = {"arm1":arms[0]}
        
        backbone_seq, backbone_sub_seqs = _SCRINSHOT_or_ISS_backbone_sequence(gene_idx,barcode_seed=barcode_seed,barcode_length=barcode_length)
        
        sub_seqs.update(backbone_sub_seqs)
        sub_seqs.update({"arm2":arms[1]})
        
        full_seq = arms[0] + backbone_seq + arms[1]
        
        return full_seq, sub_seqs


    def get_detection_oligo(self, probe_sequence, ligation_site, minT=2):
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
        minT: int
            Minimal number of T(hymines) in detection oligo
        
        Returns
        -------
        str:
            Sequence of detection oligo
        float:
            Melting temperature of detection oligo
        
        """
        def _get_oligo_Tm(probe_sequence,Tm_parameters,Tm_correction_parameters):
            """Compute the melting temperature for the detection oligo sequence
            
            #TODO: this function is not placed very nicely. Also it uses utils.get_Tm_parameters everytime we calculate a new 
                oligo. In datamodule.py we have the same function for the probe sequence instead of the oligo. Idk atm, 
                it could be handled nicer. Not rly important atm since we only compute a handful of detection oligos.
            
            In the first step the melting temperature is calculated based on sequence and salt concentrations. In the
            second step the temperature is corrected by formamide (and DMSO) percentage.
            
            :param probe_sequence: Sequence of probe
            :type probe_sequence: string
            :param Tm_parameters: Parameters for melting temperature calculation
            :type Tm_parameters: dict
            :param Tm_correction_parameters: Parameters for melting temperature formamide correction
            :type Tm_correction_parameters: dict    
            """
            Tm = mt.Tm_NN(probe_sequence, **Tm_parameters)
            Tm_corrected = round(mt.chem_correction(Tm, **Tm_correction_parameters),2)
            return Tm_corrected

        def _find_best_oligo(start_oligo, best_oligo, best_Tm_dif, minT, get_Tm_dif):
            """Find detection oligo with best melting temperature
            
            We shorten the start_oligo by 1 nucleotide each cycle and test if the shortened oligo is better.
            The procedure runs through two times: 1. starting with cutting from left, 2. from right. (We
            don't care too much about computation time ...this procedure tests each even length substring
            two times.)
            
            Arguments
            ---------
            start_oligo: str
                Longest possible oligo with even length
            best_oligo: str
                Initial best observed oligo sequence (can be start_oligo or a 1 nt longer oligo with odd length)
            best_Tm_dif: float
                Best observed melting temperature difference
            oligo_length_min: int
                Minimal length of detection oligo
            minT: int
                Minimal number of T(hymines) in detection oligo
            get_Tm_dif: fct
                Function that calculates the melting temperature difference based on the sequence
                
            Returns
            -------
            str:
                Best observed detection oligo sequence
            float:
                Corresponding best melting temperature difference
            
            """

            #  The while loop shortens the oligo 1 nt each step
            
            # Start cutting from right
            oligo = start_oligo
            oligo_length = len(oligo)
            Tm_dif = best_Tm_dif
            count = 0
            while oligo_length > self.detect_oligo_length_min:
                
                if bool(count % 2):
                    oligo = oligo[1:]
                else:
                    oligo = oligo[:-1]
                    
                Tm_dif = get_Tm_dif(oligo)
                
                # only oligos with at least minT thymines are allowed
                if (Tm_dif < best_Tm_dif) and (oligo.count("T") >= minT):
                    best_Tm_dif = Tm_dif
                    best_oligo = oligo
                    
                count += 1
                oligo_length = len(oligo)
                
            # Start cutting from left
            oligo = start_oligo
            oligo_length = len(oligo)
            Tm_dif = best_Tm_dif
            count = 0
            while oligo_length > self.detect_oligo_length_min:
                
                if bool(count % 2):
                    oligo = oligo[:-1]
                else:
                    oligo = oligo[1:]
                    
                Tm_dif = get_Tm_dif(oligo)
                
                if (Tm_dif < best_Tm_dif) and (oligo.count("T") >= minT):
                    best_Tm_dif = Tm_dif
                    best_oligo = oligo
                    
                count += 1
                oligo_length = len(oligo)
                
            return best_oligo, best_Tm_dif

        def exchange_T_with_U(probe, minT=2, U_distance=5):
            """Exchange 2 T(hymines) with U(racils) and find best side for fluorophore (closest U)
            
            Arguments
            ---------
            probe: str
                Sequence
            minT: int
                Minimal number of T(hymines) in probe
            U_distance: int
                Preferred minimal distance between U(racils)
                
            Returns
            -------
            str:
                Oligo sequence with exchanged U(racils) and fluorophore
            str:
                Info if fluorophore should be placed left or right
            
            """
            
            if probe.count("T") < minT:
                return "NOT-ENOUGH-THYMINES-FOR-DETECTION-OLIGO", None
            
            if probe.find("T") < probe[::-1].find("T"):
                fluorophor_pos = "left"
                p = probe
            else:
                fluorophor_pos = "right"
                p = probe[::-1]
            
            pos = 0
            new_pos = 1
            for i in range(minT):
                T_not_found = True
                while T_not_found:
                    shift = 0 if (pos == 0 and (new_pos != 0)) else U_distance
                    new_pos = p[min(pos+shift,len(p)):].find("T")
                    if new_pos == -1:
                        pos = p.rfind("T") - U_distance # 0
                    else:
                        p = "".join(["U" if (i == pos+shift+new_pos) else nt for i,nt in enumerate(p)])
                        pos = pos+shift+new_pos
                        T_not_found = False
                        
            if fluorophor_pos == "right":
                p = p[::-1]
            return p, fluorophor_pos

        get_Tm_dif = lambda seq: abs(_get_oligo_Tm(seq, self.Tm_parameters, self.Tm_correction_parameters) - self.detect_oligo_Tm_opt)
        
        # Search for best oligos
        start_oligo, start_oligo_long_left, start_oligo_long_right = self._get_initial_oligos_for_search(probe_sequence, ligation_site)

        # Check which of the three initial oligos is the best one
        best_oligo = start_oligo
        # The 10000 is for the case that start_oligo doesn't contain minT: The longer sequences could still contain minT
        # and Tm_dif of the longer sequences are definitely below 10000.
        best_Tm_dif = get_Tm_dif(start_oligo) if (start_oligo.count("T") >= minT) else 10000 
        for tmp_oligo in [start_oligo_long_left,start_oligo_long_right]:
            if tmp_oligo is not None:
                Tm_dif = get_Tm_dif(tmp_oligo)
                if (Tm_dif < best_Tm_dif) and (tmp_oligo.count("T") >= minT):
                    best_Tm_dif = Tm_dif
                    best_oligo = tmp_oligo

        # Iterative search through shorter oligos
        best_oligo, best_Tm_dif = _find_best_oligo(start_oligo, best_oligo, best_Tm_dif, minT, get_Tm_dif)
        
        # exchange T's with U (for enzymatic degradation of oligos)
        oligo_seq, fluorophor_pos = exchange_T_with_U(best_oligo, minT=minT, U_distance=5)

        if oligo_seq != "NOT-ENOUGH-THYMINES-FOR-DETECTION-OLIGO":
            oligo_Tm = _get_oligo_Tm(best_oligo, self.Tm_parameters, self.Tm_correction_parameters)
        else:
            oligo_Tm = 0
            
        # Add fluorophore
        if fluorophor_pos == "left":
            oligo_seq = "[fluorophore]"+oligo_seq
        elif fluorophor_pos == "right":
            oligo_seq = oligo_seq+"[fluorophore]"
        else:
            oligo_seq = oligo_seq
        
        return oligo_seq, oligo_Tm
    

    def _get_initial_oligos_for_search(self, probe_sequence, ligation_site):  
        """Get initial oligos for best oligo search
        
        We only allow a difference of 1 nt for the sequences left and right of the ligation site.
        The search will start with the even length oligo. However, if the parameters allow for odd lengths we might also
        find oligos with an additional nucleotide on the left or right or both sides.
        
        In a firt step we find the oligo sequence that is possible based on the location of the ligation site. 
        E.g. (the "|" is only for marking the ligation site):
            - probe_sequence = AAA|CTGCTG -> oligo = AAA|CTGC 
            - probe_sequence = AAA|CTG    -> oligo = AAA|CTG
            - probe_sequence = AAAAA|CTG  -> oligo = AAAA|CTG
        
        Then the following parameter scenarios can occur:
        1. The length constraint is smaller than the length of the oligo
            1.1 the maximal length is even:
                --> only an even length oligo
            1.2 the maximal length is odd:
                --> three different oligos: even, longer left, longer right
        2. The length of the oligo is smaller than the length constraint
            2.1 the length of the oligo is even (this only happens when the ligation site is exactly in the middle of an even length probe)
                --> only an even length oligo
            2.2 the length of the oligo is odd
                2.2.1 ligation site is closer to the left
                    --> two different oligos: even, long right
                2.2.2 ligation site is closter to the right
                    --> two different oligos: even, long left
        
        Arguments
        ---------
        probe_sequence: str
            Sequence of probe for which a detection oligo is designed
        ligation_site: int
            Position of ligation site. E.g. probe_sequence="AACTG", ligation_site = 2: AA|CTG 
        oligo_length_max_constraint: int
            Maximal length of oligo sequence.
            
        Returns
        -------
        str:
            oligo with even length (ligation site is in the center of the oligo)
        str or None:
            oligo with odd length (ligation site is in the center + 1 of the oligo)
        str or None
            oligo with odd length (ligation site is in the center - 1 of the oligo)
        
        """
        if self.detect_oligo_length_max ==  0:
            return None, None, None
        
        max_len_constraint_is_even = (self.detect_oligo_length_max % 2) == 0
        constraint_half_len = self.detect_oligo_length_max//2
        
        probe_length = len(probe_sequence)
        oligo_half_length_max = min(ligation_site,probe_length-ligation_site)
        
        if ligation_site == (probe_length-ligation_site):
            oligo = probe_sequence
        elif ligation_site > (probe_length-ligation_site):
            oligo = probe_sequence[probe_length-2*oligo_half_length_max-1:]
        else:
            oligo = probe_sequence[:2*oligo_half_length_max+1]
        
        # Different scenarios
        if self.detect_oligo_length_max < len(oligo):
            # 1.1
            if max_len_constraint_is_even:
                start_oligo = probe_sequence[ligation_site-constraint_half_len:ligation_site+constraint_half_len]
                start_oligo_long_left = None
                start_oligo_long_right = None    
            # 1.2
            else:
                start_oligo = probe_sequence[ligation_site-constraint_half_len:ligation_site+constraint_half_len]
                start_oligo_long_left = probe_sequence[ligation_site-constraint_half_len-1:ligation_site+constraint_half_len]
                start_oligo_long_right = probe_sequence[ligation_site-constraint_half_len:ligation_site+constraint_half_len+1]
        else:
            # 2.1
            if (len(oligo) % 2) == 0:
                start_oligo = oligo
                start_oligo_long_left = None
                start_oligo_long_right = None    
            else:
                # 2.2.1
                if (len(oligo) - ligation_site) > ligation_site:
                    start_oligo_long_right = oligo
                    start_oligo_long_left = None
                    start_oligo = oligo[:-1]
                # 2.2.2
                else:
                    start_oligo_long_left = oligo
                    start_oligo_long_right = None
                    start_oligo = oligo[1:]
        
        return start_oligo, start_oligo_long_left, start_oligo_long_right