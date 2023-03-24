from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_property_filter import (
    PropertyFilter,
    GCContent,
    ConsecutiveRepeats,
    MaskedSequences,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilter,
    ExactMatches,
    Blastn,
)
from oligo_designer_toolsuite.utils import get_barcode
import os


class SeqFishPlusProbeDesigner(BaseProbeDesigner):
    """_summary_

    Args:
        BaseProbeDesigner (_type_): _description_
    """

    def __init__(self, config) -> None:
        super().__init__(config)

    def filter_probes_by_property(
        self,
        GC_content_min,
        GC_content_max,
        number_consecutive,
    ):
        masked_sequences = MaskedSequences()
        gc_content = GCContent(
            GC_content_min=GC_content_min, GC_content_max=GC_content_max
        )
        proh_seq = ConsecutiveRepeats(num_consecutive=number_consecutive)
        filters = [masked_sequences, proh_seq, gc_content]
        property_filter = PropertyFilter(
            filters=filters,
            write_regions_with_insufficient_oligos=write_removed_genes,
        )
        self.oligo_database = property_filter.apply(
            oligo_database=self.oligo_database, n_jobs=n_jobs
        )
        self.oligo_database.remove_regions_with_insufficient_oligos(
            pipeline_step="after applying property filters"
        )
        if self.write_intermediate_steps:
            file_database = self.oligo_database.write_database(
                filename="oligo_database_property_filter.txt"
            )

    def filter_probes_by_specificity(
        self, probe_length_max, word_size, percent_identity, coverage, strand
    ):
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
        specificity_filter = SpecificityFilter(
            filters=filters,
            write_regions_with_insufficient_oligos=write_removed_genes,
        )

        self.ref_db_cross_hybroligo_database = specificity_filter.apply(
            oligo_database=self.oligo_database,
            reference_database=self.reference,
            n_jobs=n_jobs,
        )
        if self.write_intermediate_steps:
            file_database = self.oligo_database.write_database(
                filename="oligo_database_specificity_filter.txt"
            )

        num_genes_after, num_probes_after = self._get_probe_database_info(
            oligo_database.database
        )

        logging.info(
            f"Step - Filter Probes by Specificity: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        shutil.rmtree(dir_specificity)
        return self.oligo_database, file_database

    def design_readout_probes(
        self,
        GC_content_min,
        GC_content_max,
        number_consecutive,
        length,
        num_probes,
        reference_DB,
    ):
        """_summary_"""
        property_filters = [
            GCContent(GC_content_min=GC_content_min, GC_content_max=GC_content_max),
            ConsecutiveRepeats(num_consecutive=number_consecutive),
        ]
        specificity_filters = [self.blast_filter]

        ref = os.path.basename(reference_DB.file_fasta)

        readout_database = OligoDatabase(file_fasta=None, dir_output=self.dir_output)

        readout_database.create_random_database(
            length * 100, num_probes, sequence_alphabet=self.sequence_alphabet
        )

        property_filter = PropertyFilter(property_filters)
        readout_database = property_filter.apply(readout_database)

        specificity_filter = SpecificityFilter(specificity_filters)
        readout_database = specificity_filter.apply(readout_database, ref)
        return readout_database

    def write_barcodes(self):
        output_file_barcodes = os.path.join(self.dir_output, "region_to_barcode.txt")
        with open(output_file_barcodes, "w") as f:
            for i, region_id in enumerate(self.oligo_database.keys()):
                f.write(region_id + " : " + get_barcode(i))

        # logging.info(
        #     f"Barcodes are assigned"
        # )

    def design_final_SeqFishPlus_probes(self, oligos_DB, readout_probes, barcodes):
        for i in oligos_DB.keys():
            barcode = barcodes[i]
            s1 = readout_probes[barcode[0]]
            s2 = readout_probes[barcode[1]]
            s3 = readout_probes[barcode[2]]
            s4 = readout_probes[barcode[3]]
            left = str(s1) + str(s2)
            right = str(s3) + str(s4)
            for j in oligos_DB[i].keys():
                seq = str(oligos_DB[i][j]["sequence"])
                seq = left + seq
                seq = seq + right
                oligos_DB[i][j]["sequence"] = Seq(seq)
        return oligos_DB

    def check_cross_hybridation(
        self, dir_specificity, word_size, percent_identity, coverage, n_jobs
    ):
        # CROSS_HYBRIDIZATION CHECK
        # Check for cross-hybridization is performed here, ExactMatches and Blastn filters are being used
        # Blastn filter is used with strand "minus"
        self.oligo_database.write_fasta_from_database(filename="fasta_from_our_db")
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
        specificity_filter = SpecificityFilter(
            filters=filters,
            write_regions_with_insufficient_oligos=self.write_removed_genes,
        )
        self.oligo_database = specificity_filter.apply(
            oligo_database=self.oligo_database,
            reference_database=ref_db_cross_hybr,
            n_jobs=n_jobs,
        )

        logging.info(f"Cross-hybridization check is performed")
