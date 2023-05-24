############################################
# imports
############################################

import os
import yaml
import random
import itertools
import pandas as pd

from pathlib import Path
from Bio.SeqUtils import MeltingTemp as mt

from ..database import OligoDatabase
from ..utils._sequence_design import SCRINSHOT_or_ISS_backbone_sequence, get_barcode, convert_complementary_seq_to_arms

############################################
# oligo set generator class
############################################

class PadlockSequence:
    """
    This class is used to design the final padlock sequences.

    :param detect_oligo_length_min: min lenght of the detection oligo
    :type detect_oligo_length_min: int
    :param detect_oligo_length_max: max lenght of the detection oligo
    :type detect_oligo_length_max: int
    :param detect_oligo_Tm_opt: optimal melting temperature of the detection oligo
    :type detect_oligo_Tm_opt: float
    :param Tm_parameters: parameters to compute the melting temperature, for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_parameters: dict
    :param Tm_correction_parameters: parameters to correct the melting temperature, for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
    :type Tm_correction_parameters: dict
    """

    def __init__(
        self,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        detect_oligo_Tm_opt: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        dir_output: str = "output",
    ):
        """Constructor method"""

        # set parameters
        self.detect_oligo_length_min = detect_oligo_length_min
        self.detect_oligo_length_max = detect_oligo_length_max
        self.detect_oligo_Tm_opt = detect_oligo_Tm_opt
        self.Tm_parameters = Tm_parameters
        self.Tm_correction_parameters = Tm_chem_correction_parameters

        self.dir_output = os.path.join(dir_output, "padlock_sequences")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

    def design_final_padlock_sequence(self, oligo_database: OligoDatabase):
        """Design final padlock oligo sequences for the oligosets in the given database.

        Saves in a subfolder of ``dir_output`` called ``padlock_sequences`` the following files:
        - table at dir_out+"padlock_oligos.yml" with final padlock oligo sequences, detection oligo sequences, and infos
        - table at dir_out+"padlock_oligos_order.yml" with final padlock oligo sequences and detection oligo sequences


        :param oligo_database: database containging all teh infromation on the oligo sequences
        :type oligo_database: OligoDatabase
        """
        database = oligo_database.database
        oligosets = oligo_database.oligosets
        yaml_dict = {}
        regions = list(oligosets.keys())

        for region_idx, region in enumerate(regions):
            yaml_dict[region] = {}

            database_region = database[region]
            oligosets_region = oligosets[region]
            oligoset = self._best_oligoset_with_possible_detection_oligos(
                oligosets_region, database_region, minT=2
            )

            for oligo_idx, oligo_id in enumerate(oligoset):
                # NOTE: so far what we called "oligo" is actually the sequence on the target mRNA, it was straightforward
                #       to adjust everything from here, but it would have been cleaner to take the reverse complement
                #       from the beginning. The previous ligation_site for example needed to be adjusted which might be
                #       confusing when looking into files generated in previous steps. Also the arm Tm's could have
                #       changed, but the naming (arm1, arm2) still fits and Tm(seq)==Tm(rev_compl(seq)) (properly tested).
                target_mRNA = database_region[oligo_id]["sequence"]
                complementary_seq = str(target_mRNA.reverse_complement())
                ligation_idx = len(target_mRNA) - int(
                    database_region[oligo_id]["ligation_site"]
                )
                full_seq, sub_seqs = self._get_padlock_oligo(
                    region_idx,
                    complementary_seq,
                    ligation_idx,
                    barcode_seed=0,
                    barcode_length=4,
                )
                det_oligo_seq, det_oligo_Tm = self._get_detection_oligo(
                    complementary_seq, ligation_idx, minT=2
                )

                yaml_dict[region][f"{region}_oligo{oligo_idx+1}"] = {}

                # TODO: potentially nice to also save organism, full region name, reference genome
                yaml_dict[region][f"{region}_oligo{oligo_idx+1}"]["id"] = oligo_id
                yaml_dict[region][f"{region}_oligo{oligo_idx+1}"]["region"] = region
                for key in [
                    "additional_information_fasta",
                    "chromosome",
                    "start",
                    "end",
                    "strand",
                ]:
                    if len(database_region[oligo_id][key]) == 1:
                        yaml_dict[region][f"{region}_oligo{oligo_idx+1}"][key] = str(
                            database_region[oligo_id][key][0]
                        )
                    else:
                        yaml_dict[region][f"{region}_oligo{oligo_idx+1}"][
                            key
                        ] = ",".join(str(database_region[oligo_id][key]))
                yaml_dict[region][f"{region}_oligo{oligo_idx+1}"].update(
                    {
                        "padlock_probe_full_sequence": str(full_seq),
                        "detection_probe_sequence": str(det_oligo_seq),
                        "padlock_arm1_sequence": str(sub_seqs["arm1"]),
                        "padlock_accessory1_sequence": str(sub_seqs["accessory1"]),
                        "padlock_ISS_anchor_sequence": str(sub_seqs["ISS_anchor"]),
                        "padlock_barcode_sequence": str(sub_seqs["barcode"]),
                        "padlock_accessory2_sequence": str(sub_seqs["accessory1"]),
                        "padlock_arm2_sequence": str(sub_seqs["arm2"]),
                        "complementary_sequence": str(complementary_seq),
                        "recognised_mRNA_sequence": str(
                            target_mRNA
                        ),  # str(Seq(complementary_seq).complement()),
                    }
                )
                for key in [
                    "GC_content",
                    "melting_temperature",
                    "melt_temp_arm1",
                    "melt_temp_arm2",
                    "dif_melt_temp_arms",
                ]:
                    yaml_dict[region][f"{region}_oligo{oligo_idx+1}"][key] = float(
                        database_region[oligo_id][key]
                    )
                for key in ["length"]:  # , "ligation_site"]:
                    yaml_dict[region][f"{region}_oligo{oligo_idx+1}"][key] = int(
                        database_region[oligo_id][key]
                    )
                yaml_dict[region][f"{region}_oligo{oligo_idx+1}"][
                    "ligation_site"
                ] = int(ligation_idx)
                yaml_dict[region][f"{region}_oligo{oligo_idx+1}"][
                    "melt_temp_detection_probe"
                ] = float(det_oligo_Tm)

        with open(os.path.join(self.dir_output, "padlock_probes.yml"), "w") as outfile:
            yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)

        yaml_order = {}
        for region in yaml_dict:
            yaml_order[region] = {}
            for oligo_id in yaml_dict[region]:
                yaml_order[region][oligo_id] = {}
                yaml_order[region][oligo_id]["padlock_probe_full_sequence"] = yaml_dict[
                    region
                ][oligo_id]["padlock_probe_full_sequence"]
                yaml_order[region][oligo_id]["detection_probe_sequence"] = yaml_dict[
                    region
                ][oligo_id]["detection_probe_sequence"]

        with open(
            os.path.join(self.dir_output, "padlock_probes_order.yml"), "w"
        ) as outfile:
            yaml.dump(yaml_order, outfile, default_flow_style=False, sort_keys=False)

    def _best_oligoset_with_possible_detection_oligos(
        self, oligosets_region: pd.DataFrame, oligos_DB_gene: dict, minT: int = 2
    ):
        """Get row index of best oligoset for which all detection oligos can be designed

        Criterions for detection oligo design are 1. a length constraint and 2. to have at least 2 T(hymines).

        TODO: The detection oligo thing is a bit annoying. Currently I don't want to set this as a constrained for
            finding oligos for a region. Otherwise we could directly filter oligos before searching for oligosets
            which would be great. Maybe something to think about in the future. For now we only check after having
            the oligo sets per region.

        Arguments
        ---------
        oligosets_region: pd.DataFrame
            Dataframe with ranked oligosets. Output of get_nonoverlapping_sets in nonoverlapping_sets.py.
        oligos_DB_gene: dict
            Dictionary with oligo infos.
        minT: int
            Minimal number of T(hymines) in detection oligo.

        Returns
        -------
        int:
            Numerical row index of first oligoset that only contains oligos for which a detection oligo can be designed.
            If no oligoset can be found the first oligoset index i.e. 0 is returned.

        """
        for _, row in oligosets_region.iterrows():
            oligoset = [row[col] for col in row.index if col.startswith("oligo_")]

            for oligo_id in oligoset:
                target_mRNA = oligos_DB_gene[oligo_id]["sequence"]
                complementary_seq = str(target_mRNA.reverse_complement())
                ligation_idx = len(target_mRNA) - int(
                    oligos_DB_gene[oligo_id]["ligation_site"]
                )

                (
                    start_oligo,
                    start_oligo_long_left,
                    start_oligo_long_right,
                ) = self._get_initial_oligos_for_search(complementary_seq, ligation_idx)
                if (start_oligo_long_left is not None) and (
                    start_oligo_long_left.count("T") >= minT
                ):
                    return oligoset
                elif (start_oligo_long_right is not None) and (
                    start_oligo_long_right.count("T") >= minT
                ):
                    return oligoset
                elif start_oligo.count("T") >= minT:
                    return oligoset

        return [
            oligosets_region[col].iloc[0]
            for col in oligosets_region.columns
            if col.startswith("oligo_")
        ]  # return the first set

    def _get_padlock_oligo(
        self,
        region_idx,
        complementary_seq,
        ligation_idx,
        barcode_seed=0,
        barcode_length=4,
    ):
        """Get full padlock oligo for a given region and a given complementary sequence

        Arguments
        ---------
        region_idx: int
            Identifier for a given region. The identifier makes sure to return the same bar code
            for the different padlock oligos of a given region.
        complementary_seq: str
            Sequence that hybridises with the target RNA
        ligation_idx: int
            Site where complementary_seq is cut in two arms according padlock oligo design
        barcode_seed: int
            Defines the random assignment of barcodes to each region_idx.
        barcode_length: int
            Length of barcode sequence

        Returns
        -------
        str:
            padlock oligo sequence (5' to 3')
        dict of strs:
            Individual parts of the padlock sequence

        """

        
        arms = convert_complementary_seq_to_arms(complementary_seq, ligation_idx)
        sub_seqs = {"arm1": arms[0]}

        backbone_seq, backbone_sub_seqs = SCRINSHOT_or_ISS_backbone_sequence(
            region_idx, barcode_seed=barcode_seed, barcode_length=barcode_length
        )

        sub_seqs.update(backbone_sub_seqs)
        sub_seqs.update({"arm2": arms[1]})

        full_seq = arms[0] + backbone_seq + arms[1]

        return full_seq, sub_seqs

    def _get_detection_oligo(self, oligo_sequence, ligation_site, minT=2):
        """Get detection oligo sequence for a given oligo

        Detection oligos have the same sequence as the complementary sequence (i.e. `oligo_sequence`) but shortend and
        reversed. The ligation site is placed in the middle of the detection oligo. The detection oligo is shortend to get
        its melting temperature as close as possible to config["detect_oligo_Tm_opt"] but not shorter than
        config["detect_oligo_length_min"].


        Arguments
        ---------
        oligo_sequence: str
            The sequence of the complementary oligo
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

        def _get_oligo_Tm(oligo_sequence, Tm_parameters, Tm_chem_correction_parameters):
            """Compute the melting temperature for the detection oligo sequence

            #TODO: this function is not placed very nicely. Also it uses utils.get_Tm_parameters everytime we calculate a new
                oligo. In datamodule.py we have the same function for the oligo sequence instead of the oligo. Idk atm,
                it could be handled nicer. Not rly important atm since we only compute a handful of detection oligos.

            In the first step the melting temperature is calculated based on sequence and salt concentrations. In the
            second step the temperature is corrected by formamide (and DMSO) percentage.

            :param oligo_sequence: Sequence of oligo
            :type oligo_sequence: string
            :param Tm_parameters: Parameters for melting temperature calculation
            :type Tm_parameters: dict
            :param Tm_chem_correction_parameters: Parameters for melting temperature formamide correction
            :type Tm_chem_correction_parameters: dict
            """
            Tm = mt.Tm_NN(oligo_sequence, **Tm_parameters)
            Tm_corrected = round(
                mt.chem_correction(Tm, **Tm_chem_correction_parameters), 2
            )
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

        def _exchange_T_with_U(oligo, minT=2, U_distance=5):
            """Exchange 2 T(hymines) with U(racils) and find best side for fluorophore (closest U)

            Arguments
            ---------
            oligo: str
                Sequence
            minT: int
                Minimal number of T(hymines) in oligo
            U_distance: int
                Preferred minimal distance between U(racils)

            Returns
            -------
            str:
                Oligo sequence with exchanged U(racils) and fluorophore
            str:
                Info if fluorophore should be placed left or right

            """

            if oligo.count("T") < minT:
                return "NOT-ENOUGH-THYMINES-FOR-DETECTION-OLIGO", None

            if oligo.find("T") < oligo[::-1].find("T"):
                fluorophor_pos = "left"
                p = oligo
            else:
                fluorophor_pos = "right"
                p = oligo[::-1]

            pos = 0
            new_pos = 1
            for i in range(minT):
                T_not_found = True
                while T_not_found:
                    shift = 0 if (pos == 0 and (new_pos != 0)) else U_distance
                    new_pos = p[min(pos + shift, len(p)) :].find("T")
                    if new_pos == -1:
                        pos = p.rfind("T") - U_distance  # 0
                    else:
                        pos = pos + shift + new_pos
                        p = p[:pos] + "U" + p[pos + 1 :]
                        T_not_found = False

            if fluorophor_pos == "right":
                p = p[::-1]
            return p, fluorophor_pos

        get_Tm_dif = lambda seq: abs(
            _get_oligo_Tm(seq, self.Tm_parameters, self.Tm_correction_parameters)
            - self.detect_oligo_Tm_opt
        )

        # Search for best oligos
        (
            start_oligo,
            start_oligo_long_left,
            start_oligo_long_right,
        ) = self._get_initial_oligos_for_search(oligo_sequence, ligation_site)

        # Check which of the three initial oligos is the best one
        best_oligo = start_oligo
        # The 10000 is for the case that start_oligo doesn't contain minT: The longer sequences could still contain minT
        # and Tm_dif of the longer sequences are definitely below 10000.
        best_Tm_dif = (
            get_Tm_dif(start_oligo) if (start_oligo.count("T") >= minT) else 10000
        )
        for tmp_oligo in [start_oligo_long_left, start_oligo_long_right]:
            if tmp_oligo is not None:
                Tm_dif = get_Tm_dif(tmp_oligo)
                if (Tm_dif < best_Tm_dif) and (tmp_oligo.count("T") >= minT):
                    best_Tm_dif = Tm_dif
                    best_oligo = tmp_oligo

        # Iterative search through shorter oligos
        best_oligo, best_Tm_dif = _find_best_oligo(
            start_oligo, best_oligo, best_Tm_dif, minT, get_Tm_dif
        )

        # exchange T's with U (for enzymatic degradation of oligos)
        oligo_seq, fluorophor_pos = _exchange_T_with_U(
            best_oligo, minT=minT, U_distance=5
        )

        if oligo_seq != "NOT-ENOUGH-THYMINES-FOR-DETECTION-OLIGO":
            oligo_Tm = _get_oligo_Tm(
                best_oligo, self.Tm_parameters, self.Tm_correction_parameters
            )
        else:
            oligo_Tm = 0

        # Add fluorophore
        if fluorophor_pos == "left":
            oligo_seq = "[fluorophore]" + oligo_seq
        elif fluorophor_pos == "right":
            oligo_seq = oligo_seq + "[fluorophore]"
        else:
            oligo_seq = oligo_seq

        return oligo_seq, oligo_Tm

    def _get_initial_oligos_for_search(self, oligo_sequence, ligation_site):
        """Get initial oligos for best oligo search

        We only allow a difference of 1 nt for the sequences left and right of the ligation site.
        The search will start with the even length oligo. However, if the parameters allow for odd lengths we might also
        find oligos with an additional nucleotide on the left or right or both sides.

        In a firt step we find the oligo sequence that is possible based on the location of the ligation site.
        E.g. (the "|" is only for marking the ligation site):
            - oligo_sequence = AAA|CTGCTG -> oligo = AAA|CTGC
            - oligo_sequence = AAA|CTG    -> oligo = AAA|CTG
            - oligo_sequence = AAAAA|CTG  -> oligo = AAAA|CTG

        Then the following parameter scenarios can occur:
        1. The length constraint is smaller than the length of the oligo
            1.1 the maximal length is even:
                --> only an even length oligo
            1.2 the maximal length is odd:
                --> three different oligos: even, longer left, longer right
        2. The length of the oligo is smaller than the length constraint
            2.1 the length of the oligo is even (this only happens when the ligation site is exactly in the middle of an even length oligo)
                --> only an even length oligo
            2.2 the length of the oligo is odd
                2.2.1 ligation site is closer to the left
                    --> two different oligos: even, long right
                2.2.2 ligation site is closter to the right
                    --> two different oligos: even, long left

        Arguments
        ---------
        oligo_sequence: str
            Sequence of oligo for which a detection oligo is designed
        ligation_site: int
            Position of ligation site. E.g. oligo_sequence="AACTG", ligation_site = 2: AA|CTG
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
        if self.detect_oligo_length_max == 0:
            return None, None, None

        max_len_constraint_is_even = (self.detect_oligo_length_max % 2) == 0
        constraint_half_len = self.detect_oligo_length_max // 2

        oligo_length = len(oligo_sequence)
        oligo_half_length_max = min(ligation_site, oligo_length - ligation_site)

        if ligation_site == (oligo_length - ligation_site):
            oligo = oligo_sequence
        elif ligation_site > (oligo_length - ligation_site):
            oligo = oligo_sequence[oligo_length - 2 * oligo_half_length_max - 1 :]
        else:
            oligo = oligo_sequence[: 2 * oligo_half_length_max + 1]

        # Different scenarios
        if self.detect_oligo_length_max < len(oligo):
            # 1.1
            if max_len_constraint_is_even:
                start_oligo = oligo_sequence[
                    ligation_site
                    - constraint_half_len : ligation_site
                    + constraint_half_len
                ]
                start_oligo_long_left = None
                start_oligo_long_right = None
            # 1.2
            else:
                start_oligo = oligo_sequence[
                    ligation_site
                    - constraint_half_len : ligation_site
                    + constraint_half_len
                ]
                start_oligo_long_left = oligo_sequence[
                    ligation_site
                    - constraint_half_len
                    - 1 : ligation_site
                    + constraint_half_len
                ]
                start_oligo_long_right = oligo_sequence[
                    ligation_site
                    - constraint_half_len : ligation_site
                    + constraint_half_len
                    + 1
                ]
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
