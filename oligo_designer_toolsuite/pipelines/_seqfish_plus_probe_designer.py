from ..database import OligoDatabase, ReferenceDatabase
from ..oligo_property_filter import (
    PropertyFilter,
    GCContent,
    ConsecutiveRepeats,
    MaskedSequences,
)
from ..oligo_specificity_filter import (
    SpecificityFilter,
    ExactMatches,
    Blastn,
)
from ..oligo_efficiency_filter import (
    SeqFISHOligoScoring,
    AverageSetScoring,
)
from ..oligo_selection import (
    padlock_heuristic_selection,
    OligosetGenerator,
)
from ..utils import create_seqfish_plus_barcodes

from ._base_probe_designer import BaseProbeDesigner

import os
import logging
import shutil

from typing import List, Tuple

from Bio.Seq import Seq


class SeqfishPlusProbeDesigner(BaseProbeDesigner):
    """_summary_

    Args:
        BaseProbeDesigner (_type_): _description_
    """

    def filter_probes_by_property(
        self, oligo_database, GC_content_min, GC_content_max, number_consecutive, n_jobs
    ):
        masked_sequences = MaskedSequences()
        gc_content = GCContent(
            GC_content_min=GC_content_min, GC_content_max=GC_content_max
        )
        proh_seq = ConsecutiveRepeats(num_consecutive=number_consecutive)
        filters = [masked_sequences, proh_seq, gc_content]
        property_filter = PropertyFilter(filters=filters)
        oligo_database = property_filter.apply(
            oligo_database=oligo_database, n_jobs=n_jobs
        )
        if self.write_intermediate_steps:
            file_database = oligo_database.write_database(
                filename="oligo_database_property_filter.txt"
            )
        return oligo_database, file_database

    def filter_probes_by_specificity(
        self,
        oligo_database: OligoDatabase,
        probe_length_max: int,
        word_size: int,
        percent_identity: float,
        coverage: float,
        strand: str,
        n_jobs: int,
    ) -> Tuple[OligoDatabase, str]:
        """_summary_

        Args:
            oligo_database (OligoDatabase): _description_
            probe_length_max (int): _description_
            word_size (int): _description_
            percent_identity (float): _description_
            coverage (float): _description_
            strand (str): _description_
            n_jobs (int): _description_

        Returns:
            Tuple[OligoDatabase, str]: _description_
        """

        num_genes_before, num_probes_before = self._get_probe_database_info(
            oligo_database.database
        )
        dir_specificity = os.path.join(
            self.dir_output, "specificity_temporary"
        )  # folder where the temporary files will be written

        file_transcriptome = (
            self.region_generator.generate_transcript_reduced_representation(
                include_exon_junctions=True, exon_junction_size=2 * probe_length_max
            )
        )
        self.reference = ReferenceDatabase(
            file_fasta=file_transcriptome,
            files_source=self.region_generator.files_source,
            species=self.region_generator.species,
            annotation_release=self.region_generator.annotation_release,
            genome_assembly=self.region_generator.genome_assembly,
            dir_output=self.dir_output,
        )
        self.reference.load_fasta_into_database()
        logging.info(f"Reference DB created")
        exact_mathces = ExactMatches(dir_specificity=dir_specificity)
        blastn = Blastn(
            dir_specificity=dir_specificity,
            word_size=word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand=strand,
            # strand='plus',
        )
        filters = [exact_mathces, blastn]
        specificity_filter = SpecificityFilter(filters=filters)

        self.ref_db_cross_hybroligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            reference_database=self.reference,
            n_jobs=n_jobs,
        )
        if self.write_intermediate_steps:
            file_database = oligo_database.write_database(
                filename="oligo_database_specificity_filter.txt"
            )

        num_genes_after, num_probes_after = self._get_probe_database_info(
            oligo_database.database
        )

        logging.info(
            f"Step - Filter Probes by Specificity: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        shutil.rmtree(dir_specificity)
        return oligo_database, file_database

    def design_readout_probes(
        self,
        GC_content_min: float,
        GC_content_max: float,
        number_consecutive: int,
        length: int,
        blast_filter: Blastn,
        sequence_alphabet: List[str],
        num_probes: int,
        reference_DB: ReferenceDatabase,
    ) -> OligoDatabase:
        """_summary_

        Args:
            GC_content_min (float): _description_
            GC_content_max (float): _description_
            number_consecutive (int): _description_
            length (int): _description_
            blast_filter (Blastn): _description_
            sequence_alphabet (List[str]): _description_
            num_probes (int): _description_
            reference_DB (ReferenceDatabase): _description_

        Returns:
            OligoDatabase: _description_
        """

        property_filters = [
            GCContent(GC_content_min=GC_content_min, GC_content_max=GC_content_max),
            ConsecutiveRepeats(num_consecutive=number_consecutive),
        ]
        specificity_filters = [blast_filter]

        readout_database = OligoDatabase(dir_output=self.dir_output)

        # TODO make sure to generate enough random probes
        readout_database.create_random_database(
            length, num_probes * 20, sequence_alphabet=sequence_alphabet
        )

        property_filter = PropertyFilter(property_filters)
        readout_database = property_filter.apply(readout_database)

        specificity_filter = SpecificityFilter(specificity_filters)
        readout_database = specificity_filter.apply(readout_database, reference_DB)
        return readout_database

    def design_final_SeqFishPlus_probes(
        self,
        oligo_database,
        readout_database,
        n_pseudocolors,
        seed,
    ):
        barcodes = create_seqfish_plus_barcodes(
            n_pseudocolors=n_pseudocolors,
            seed=seed,
            num_genes=len(oligo_database.database.keys()),
        )
        readout_sequences = readout_database.to_sequence_list()
        for i in oligo_database.keys():
            barcode = barcodes[i]
            left = readout_sequences[barcode[0]] + readout_sequences[barcode[1]]
            right = readout_sequences[barcode[2]] + readout_sequences[barcode[3]]

            for j in oligo_database[i].keys():
                seq = str(oligo_database[i][j]["sequence"])
                seq = left + seq
                seq = seq + right
                oligo_database[i][j]["sequence"] = Seq(seq)
        return oligo_database

    def check_cross_hybridization(
        self,
        oligo_database,
        dir_specificity,
        word_size,
        percent_identity,
        coverage,
        n_jobs,
    ):
        # CROSS_HYBRIDIZATION CHECK
        # Check for cross-hybridization is performed here, ExactMatches and Blastn filters are being used
        # Blastn filter is used with strand "minus"
        oligo_database.write_fasta_from_database(filename="fasta_from_our_db")
        ref_db_cross_hybr = ReferenceDatabase(
            file_fasta=self.dir_output + "/oligo_database/fasta_from_our_db.fna"
        )
        exact_mathces = ExactMatches(dir_specificity=dir_specificity)
        blastn_cross_hybr = Blastn(
            dir_specificity=dir_specificity,
            word_size=word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            # THE MOST IMPORTANT PART HERE IS STRAND = "MINUS"
            strand="minus",
        )
        filters = [exact_mathces, blastn_cross_hybr]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            reference_database=ref_db_cross_hybr,
            n_jobs=n_jobs,
        )
        if self.write_intermediate_steps:
            file_database = oligo_database.write_database(
                filename="cross_hybridation.txt"
            )

        logging.info(f"Cross-hybridization check is performed")
        return oligo_database, file_database

    def get_oligosets(
        self,
        oligo_database,
        GC_content_min,
        GC_content_opt,
        GC_content_max,
        GC_weight,
        oligoset_size,
        min_oligoset_size,
        oligos_scoring,
        n_sets,
        n_jobs,
    ):
        oligos_scoring = SeqFISHOligoScoring(
            GC_content_min=GC_content_min,
            GC_content_opt=GC_content_opt,
            GC_content_max=GC_content_max,
            GC_weight=GC_weight,
        )
        set_scoring = AverageSetScoring()
        oligoset_generator = OligosetGenerator(
            oligoset_size=oligoset_size,
            min_oligoset_size=min_oligoset_size,
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            heurustic_selection=padlock_heuristic_selection,
        )
        oligo_database = oligoset_generator.apply(
            oligo_database=oligo_database, n_sets=n_sets, n_jobs=n_jobs
        )
        if self.write_intermediate_steps:
            dir_oligosets = oligo_database.write_oligosets(folder="oligosets")
            file_database = oligo_database.write_database(
                filename="oligo_database_oligosets.txt"
            )
        logging.info(f"Oligoset generation is finished. Done!")

        return oligo_database, file_database, dir_oligosets
