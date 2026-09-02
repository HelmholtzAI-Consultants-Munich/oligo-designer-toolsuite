"""
seqFISH+ probe designer pipeline.

seqFISH+, or sequential fluorescence in situ hybridization plus, is an
image-based RNA detection method that reads combinatorial pseudocolor barcodes
over multiple hybridization rounds. This allows large gene panels to be imaged
with a limited number of fluorescence channels.

See :class:`SeqFishPlusProbeDesigner` for the full pipeline description and
probe structure. See :func:`seqfish_plus_probe_designer` for the config-driven
workflow.
"""

############################################
# imports
############################################

import os
import shutil
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.SeqUtils import Seq

from oligo_designer_toolsuite._exceptions import (
    ConfigurationError,
    FeatureNotImplementedError,
    FileFormatError,
)
from oligo_designer_toolsuite.config import SeqfishPlusProbeDesignerConfig
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    DeviationFromOptimalGCContentScorer,
    LowestSetScoring,
    OligoScoring,
    OverlapUTRScorer,
)
from oligo_designer_toolsuite.oligo_property_calculator import (
    IsoformConsensusProperty,
    PropertyCalculator,
    ReverseComplementSequenceProperty,
)
from oligo_designer_toolsuite.oligo_property_calculator._property_functions import calc_tm_nn
from oligo_designer_toolsuite.oligo_property_filter import (
    BasePropertyFilter,
    ComplementFilter,
    GCClampFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SelfComplementFilter,
    SoftMaskedSequenceFilter,
)
from oligo_designer_toolsuite.oligo_selection import IndependentSetsOligoSelection
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveAllFilterPolicy,
    RemoveByDegreeFilterPolicy,
    RemoveByLargerRegionFilterPolicy,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._utils import (
    apply_required_parameters,
    base_log_parameters,
    base_parser,
    check_content_oligo_database,
    format_sequence,
    load_config,
    pipeline_step_basic,
    preprocess_tm_parameters,
    validate_bit_mapping_table,
    validate_codebook,
    validate_primer_sequence,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import append_nucleotide_to_sequences, configure_root_logger, logger

############################################
# SeqFISH Plus Probe Designer
############################################


class SeqFishPlusProbeDesigner:
    """
    Design probe libraries for seqFISH+ experiments.

    This class runs the design workflow for seqFISH+ encoding probes, readout
    probe assignments, DNA template probes, and output files. The designed
    library can be used for multiplexed RNA detection with sequential imaging
    and combinatorial barcoding.

    Overview
    --------
    seqFISH+, or sequential fluorescence in situ hybridization plus, is an
    imaging-based method for measuring many RNA species in single cells while
    preserving their spatial positions in tissue.

    The method assigns each RNA species a barcode that is read over multiple
    rounds of hybridization and imaging. In each round, fluorescent readout
    probes bind to selected readout-binding sequences on the encoding probes.
    The pattern of fluorescence signals across rounds identifies the RNA
    molecule.

    seqFISH+ expands the number of possible barcodes by using pseudocolors.
    A pseudocolor is defined by a combination of hybridization round and
    fluorescence channel. This allows large gene panels to be imaged with a
    limited number of microscope channels.

    Probe Structure
    ---------------
    **Hybridization (Encoding) Probes**

    Encoding probes are single-stranded DNA oligos that bind directly to target
    mRNAs. Each target gene is represented by several encoding probes. Together,
    these probes provide enough signal to detect individual RNA molecules.

    Each encoding probe contains:

    - a target-binding sequence that is complementary to the RNA,
    - one readout-binding overhang per barcode round (``n_barcode_rounds``),
    - short spacer bases between the target-binding and readout-binding parts.

    A common seqFISH+ encoding probe layout (four barcode rounds) is::

        [Overhang I] + [Overhang II] + [spacer] + [target-binding sequence]
        + [spacer] + [Overhang III] + [Overhang IV]

    With a configurable number of rounds, the first half of the overhangs are
    placed 5' of the target and the second half 3' of the target. In the
    original seqFISH+ design, the target-binding sequence was typically
    28 nucleotides long, and each readout-binding overhang was 15 nucleotides
    long. The exact lengths depend on the design settings used in the pipeline.

    **Readout Probes**

    Readout probes are short fluorescent DNA oligos that bind to the overhangs
    on the encoding probes. Each readout probe reports one barcode position in a
    defined hybridization round and fluorescence channel.

    During the experiment, readout probes are added, imaged, and removed before
    the next imaging round. After all rounds are complete, the observed
    pseudocolor pattern is decoded to identify the RNA molecule.

    seqFISH+ codebooks can include an error-checking round. This helps identify
    barcode patterns that are inconsistent with the expected code and reduces
    false assignments.

    **Codebook**

    The codebook assigns each target gene to a barcode. Rows correspond to
    genes, and columns correspond to readout positions or pseudocolors,
    depending on the design.

    In seqFISH+, each barcode is read through several rounds. The use of
    pseudocolors increases the number of possible barcodes without requiring the
    same number of fluorophores.

    **DNA Template Probes**

    DNA template probes are synthesis-ready oligos used to prepare the final
    encoding probe library. They contain the encoding probe sequence flanked by
    primer binding sites.

    A simplified template layout is::

        [Forward primer] + [encoding probe] + [Reverse primer]

    The primer binding sites allow amplification of the oligo pool. Depending on
    the protocol, the amplified pool can be transcribed, reverse-transcribed,
    processed to remove primer regions, and purified to produce single-stranded
    DNA encoding probes.

    Probe Library Preparation
    -------------------------
    A seqFISH+ encoding probe library is usually prepared from a pooled DNA oligo
    library. The pool contains many template probes, each with primer binding
    sites and an encoding probe body.

    The pool is amplified by limited-cycle PCR. It can then be transcribed and
    reverse-transcribed to generate single-stranded DNA probes. If a
    uracil-containing primer is used, USER enzyme can be used to remove primer
    sequences before the final probe pool is purified.

    In the experiment, encoding probes are hybridized to RNA in fixed cells or
    tissue. The barcode is then read over multiple imaging rounds using
    fluorescent readout probes.

    Pipeline Overview
    -----------------
    The pipeline performs the main steps needed to design a seqFISH+ probe
    library:

    1. **Target probe design**

       Design gene-specific target-binding sequences that hybridize to RNA
       transcripts.

    2. **Readout probe design**

       Load or generate readout-binding sequences and create a codebook based on
       pseudocolors and imaging rounds.

    3. **Encoding probe assembly**

       Combine target-binding sequences with assigned readout-binding overhangs
       according to the codebook.

    4. **Primer handling**

       Load or generate the forward and reverse primer sequences used for
       template amplification.

    5. **DNA template assembly**

       Add primer binding sites to the encoding probe body to create
       synthesis-ready DNA template probes.

    6. **Output generation**

       Write the designed probes, codebook, readout information, and related
       probe properties to output files.

    References
    ----------
    Eng, C. H. L., Lawson, M., Zhu, Q., Dries, R., Koulena, N., Takei, Y.,
    Yun, J., Cronin, C., Karp, C., Yuan, G. C., & Cai, L. (2019).
    Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH+.
    Nature, 568, 235-239. https://doi.org/10.1038/s41586-019-1049-y

    :param write_intermediate_steps: If ``True``, save intermediate probe
        databases after pipeline steps. This can help with checking a design run
        or finding where probes were removed.
    :type write_intermediate_steps: bool
    :param dir_output: Directory where output files and intermediate results are
        saved. The directory is created if it does not exist.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(
        self,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the SeqFishPlusProbeDesigner class."""

        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        target_probes_parameters: dict,
    ) -> OligoDatabase:
        """
        Design the RNA-binding parts of the seqFISH+ hybridization probes.

        This step designs target-binding sequences for each target gene. The
        candidate probes are filtered for sequence quality and specificity, and a
        final probe set is selected for each target. Readout-binding sequences are
        added later when the hybridization probes are assembled.

        :param target_probes_parameters: Settings for target probe design from the
            ``target_probes`` section of the pipeline config. This includes probe
            generation, sequence filters, specificity filters, and probe set
            selection settings.
        :type target_probes_parameters: dict
        :return: Database containing the selected target-binding probes for each
            target region.
        :rtype: OligoDatabase
        """
        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)
        target_probes_database = target_probe_designer.generate_target_probes(
            target_probes_parameters=target_probes_parameters,
            write_intermediate_steps=self.write_intermediate_steps,
        )
        return target_probes_database

    def design_readout_probes(
        self,
        region_ids: list[str],
        readout_probe_parameters: dict,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Load or create the seqFISH+ readout probes and codebook.

        The codebook assigns each target region to a combinatorial barcode. Active
        bits in that barcode identify the readout probes used during imaging. These
        readout probes later bind to barcode sequences on the hybridization probes
        and report which RNA molecule is present.

        The codebook and readout probe table can either be loaded from files or
        generated from the config. Both tables are checked before they are returned,
        so missing targets or missing readout sequences are caught early.

        :param region_ids: Target regions that must be represented in the codebook,
            usually gene names or gene IDs.
        :type region_ids: list[str]
        :param readout_probe_parameters: Settings from the ``readout_probes`` section
            of the pipeline config. This includes the codebook settings and the
            readout probe table settings.
        :type readout_probe_parameters: dict
        :return: Codebook and readout probe table used to assign readout probes to
            target regions.
        :rtype: tuple[pd.DataFrame, pd.DataFrame]
        """
        readout_probe_designer = ReadoutProbeDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        if readout_probe_parameters["codebook"]["source"] == "load":
            codebook = readout_probe_designer.load_codebook(
                file_codebook=readout_probe_parameters["codebook"]["file"]
            )
            codebook_source = readout_probe_parameters["codebook"]["file"]
        else:
            codebook = readout_probe_designer.generate_codebook(
                region_ids=region_ids,
                codebook_parameters=readout_probe_parameters["codebook"],
            )
            codebook_source = readout_probe_parameters["codebook"]["source"]

        if readout_probe_parameters["readout_probe_table"]["source"] == "load":
            readout_probe_table = readout_probe_designer.load_readout_probe_table(
                file_readout_probe_table=readout_probe_parameters["readout_probe_table"]["file"]
            )
            readout_probe_table_source = readout_probe_parameters["readout_probe_table"]["file"]
        else:
            readout_probe_table = readout_probe_designer.generate_readout_probe_table(
                readout_probe_parameters=readout_probe_parameters["readout_probe_table"],
                codebook_parameters=readout_probe_parameters["codebook"],
                write_intermediate_steps=self.write_intermediate_steps,
            )
            readout_probe_table_source = readout_probe_parameters["readout_probe_table"]["source"]

        readout_probe_designer.validate(
            readout_probe_table=readout_probe_table,
            readout_probe_table_source=readout_probe_table_source,
            codebook=codebook,
            region_ids=region_ids,
            codebook_source=codebook_source,
            hamming_weight=readout_probe_parameters["codebook"]["n_barcode_rounds"],
        )

        return codebook, readout_probe_table

    def assemble_hybridization_probes(
        self,
        oligo_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
        readout_probe_parameters: dict,
    ) -> OligoDatabase:
        """
        Build the seqFISH+ hybridization probes.

        Combine each RNA-binding target sequence with the readout-binding sequences
        assigned by the codebook. Each gene gets one overhang for every barcode
        round (``n_barcode_rounds``). The target-binding part keeps the probe on
        the RNA. The overhangs are read out during imaging.

        Overhangs are split evenly around the target: the first half sit 5' of the
        target and the second half sit 3' of it. A simplified layout is::

            [rc(readout overhangs, 5' half)] + [T]
            + [target-binding sequence] + [T]
            + [rc(readout overhangs, 3' half)]

        The assembled sequences are added to the existing probe database.

        :param oligo_database: Database returned by
            :py:meth:`design_target_probes`. This database is updated with the
            assembled hybridization probe sequences.
        :type oligo_database: OligoDatabase
        :param codebook: Table assigning target regions to readout barcode bits.
            Rows are target regions and columns are barcode bits.
        :type codebook: pd.DataFrame
        :param readout_probe_table: Table linking each barcode bit to its readout
            probe sequence and related readout information.
        :type readout_probe_table: pd.DataFrame
        :param readout_probe_parameters: Settings from the ``readout_probes`` section
            of the pipeline config. This includes the codebook settings and the
            readout probe table settings.
        :type readout_probe_parameters: dict
        :return: Database with seqFISH+ hybridization probe sequences added.
        :rtype: OligoDatabase
        :raises ConfigurationError: If a region has a number of active bits that
            does not equal ``n_barcode_rounds``.
        """
        region_ids = list(oligo_database.database.keys())
        n_barcode_rounds = readout_probe_parameters["codebook"]["n_barcode_rounds"]
        readout_sequence_types = [f"sequence_readout_probe_{i}" for i in range(1, n_barcode_rounds + 1)]

        oligo_database.set_database_sequence_types(
            [
                "sequence_target",
                *readout_sequence_types,
                "sequence_hybridization_probe",
            ]
        )

        n_left = n_barcode_rounds // 2

        for region_id in region_ids:
            barcode = codebook.loc[region_id]
            bits = barcode[barcode == 1].index
            readout_probe_sequences = list(readout_probe_table.loc[bits, "readout_probe_sequence"])
            if len(readout_probe_sequences) != n_barcode_rounds:
                raise ConfigurationError(
                    f"Region '{region_id}' has {len(readout_probe_sequences)} active barcode bit(s), "
                    f"but n_barcode_rounds={n_barcode_rounds}."
                )

            left_readout_sequences = readout_probe_sequences[:n_left]
            right_readout_sequences = readout_probe_sequences[n_left:]

            probe_ids = list(oligo_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:

                new_properties[probe_id]["sequence_target"] = format_sequence(
                    database=oligo_database,
                    property="target",
                    region_id=region_id,
                    oligo_id=probe_id,
                )

                for i, sequence in enumerate(readout_probe_sequences, start=1):
                    new_properties[probe_id][f"sequence_readout_probe_{i}"] = sequence

                # RC so fluorescent readout oligos can anneal to these overhangs on the probe.
                new_properties[probe_id]["sequence_hybridization_probe"] = (
                    "".join(str(Seq(sequence).reverse_complement()) for sequence in left_readout_sequences)
                    + "T"
                    + format_sequence(
                        database=oligo_database,
                        property="oligo",
                        region_id=region_id,
                        oligo_id=probe_id,
                    )
                    + "T"
                    + "".join(str(Seq(sequence).reverse_complement()) for sequence in right_readout_sequences)
                )

            oligo_database.update_oligo_properties(new_properties)

        return oligo_database

    def design_primers(
        self,
        oligo_database: OligoDatabase,
        primer_parameters: dict,
    ) -> tuple[str, str]:
        """
        Load or design the PCR primers for DNA template probes.

        The forward and reverse primers are used to amplify the DNA template probe
        library before the final hybridization probes are prepared.

        The primers can either be loaded from the config or generated by the primer
        design workflow. Both primer sequences are checked before they are returned.
        When a forward primer is generated, the assembled hybridization probes are
        also used to avoid primers that bind the probe body.

        :param oligo_database: Database returned by
            :py:meth:`assemble_hybridization_probes`. The hybridization probe
            sequences are used when generating a new forward primer.
        :type oligo_database: OligoDatabase
        :param primer_parameters: Settings from the ``primers`` section of the
            pipeline config. This includes the forward and reverse primer entries
            and their sequences or design settings.
        :type primer_parameters: dict
        :return: Reverse primer sequence and forward primer sequence.
        :rtype: tuple[str, str]
        """
        # Hybridization probes are written as a FASTA reference so generated primers
        # that anneal to the probe body are rejected.
        file_fasta_hybridization_probes_database = oligo_database.write_database_to_fasta(
            filename="db_reference_hybridization_probes",
            save_description=False,
            region_ids=None,
            sequence_type="sequence_hybridization_probe",
        )

        primer_designer = PrimerDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        if primer_parameters["reverse_primer"]["source"] == "load":
            reverse_primer_sequence = primer_designer.load_reverse_primer(
                sequence=primer_parameters["reverse_primer"]["sequence"]
            )
        else:
            reverse_primer_sequence = primer_designer.generate_reverse_primer()

        if primer_parameters["forward_primer"]["source"] == "load":
            forward_primer_sequence = primer_designer.load_forward_primer(
                sequence=primer_parameters["forward_primer"]["sequence"]
            )
        else:
            forward_primer_sequence = primer_designer.generate_forward_primer(
                parameters=primer_parameters["forward_primer"],
                reverse_primer_sequence=reverse_primer_sequence,
                file_fasta_hybridization_probes_database=file_fasta_hybridization_probes_database,
                write_intermediate_steps=self.write_intermediate_steps,
            )

        os.remove(file_fasta_hybridization_probes_database)

        primer_designer.validate(
            forward_primer=forward_primer_sequence,
            reverse_primer=reverse_primer_sequence,
        )

        return reverse_primer_sequence, forward_primer_sequence

    def assemble_dna_template_probes(
        self,
        oligo_database: OligoDatabase,
        forward_primer_sequence: str,
        reverse_primer_sequence: str,
    ) -> OligoDatabase:
        """
        Build the synthesis-ready DNA template probes.

        This step adds the forward and reverse primer sequences to each
        hybridization probe. The resulting DNA template sequences are the sequences
        used for pooled oligo synthesis and library amplification.

        A simplified layout is::

            [forward primer] + [hybridization probe] + [reverse primer]

        :param oligo_database: Database returned by
            :py:meth:`assemble_hybridization_probes`. This database is updated with
            the DNA template sequences.
        :type oligo_database: OligoDatabase
        :param forward_primer_sequence: Forward primer sequence used to amplify the
            DNA template probe library.
        :type forward_primer_sequence: str
        :param reverse_primer_sequence: Reverse primer sequence used to amplify the
            DNA template probe library.
        :type reverse_primer_sequence: str
        :return: Database with DNA template probe sequences added.
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())
        oligo_database.set_database_sequence_types(
            [
                "sequence_forward_primer",
                "sequence_reverse_primer",
                "sequence_dna_template_probe",
            ]
        )

        for region_id in region_ids:
            probe_ids = list(oligo_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:
                new_properties[probe_id]["sequence_reverse_primer"] = reverse_primer_sequence
                new_properties[probe_id]["sequence_forward_primer"] = forward_primer_sequence

                new_properties[probe_id]["sequence_dna_template_probe"] = (
                    forward_primer_sequence
                    + format_sequence(
                        database=oligo_database,
                        property="sequence_hybridization_probe",
                        region_id=region_id,
                        oligo_id=probe_id,
                    )
                    + reverse_primer_sequence
                )

            oligo_database.update_oligo_properties(new_properties)

        return oligo_database

    def generate_output(
        self,
        oligo_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
        readout_probe_parameters: dict,
        output_properties: list[str] | None = None,
    ) -> None:
        """
        Write the completed seqFISH+ probe design to files.

        Save the final probe database, the codebook, and the readout probe table.
        Also write an order-ready file with the DNA template sequences and the
        readout sequences needed for synthesis.

        If no output properties are provided, a default set of annotations and
        sequence fields is written, including one readout sequence field for each
        barcode round.

        :param oligo_database: Database returned by
            :py:meth:`assemble_dna_template_probes`.
        :type oligo_database: OligoDatabase
        :param codebook: Table assigning target regions to readout barcode bits.
            Rows are target regions and columns are barcode bits.
        :type codebook: pd.DataFrame
        :param readout_probe_table: Table linking each barcode bit to its readout
            probe sequence and related readout information.
        :type readout_probe_table: pd.DataFrame
        :param readout_probe_parameters: Settings from the ``readout_probes`` section
            of the pipeline config. This includes the codebook settings and the
            readout probe table settings.
        :type readout_probe_parameters: dict
        :param output_properties: Probe properties to include in the detailed output
            files. If ``None``, a default set of annotations and sequences is used.
        :type output_properties: list[str] | None
        :return: None
        :rtype: None
        """
        n_barcode_rounds = readout_probe_parameters["codebook"]["n_barcode_rounds"]
        readout_sequence_properties = [f"sequence_readout_probe_{i}" for i in range(1, n_barcode_rounds + 1)]
        if output_properties is None:
            output_properties = [
                "source",
                "species",
                "annotation_release",
                "genome_assembly",
                "gene_id",
                "chromosome",
                "start",
                "end",
                "strand",
                "regiontype",
                "transcript_id",
                "exon_number",
                "sequence_target",
                *readout_sequence_properties,
                "sequence_hybridization_probe",
                "sequence_forward_primer",
                "sequence_reverse_primer",
                "sequence_dna_template_probe",
                "isoform_consensus",
            ]

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="gene_name")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        oligo_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="seqfish_plus_probes",
        )

        oligo_database.write_oligosets_to_table(
            properties=output_properties,
            ascending=True,
            filename="seqfish_plus_probes",
        )

        oligo_database.write_ready_to_order_yaml(
            properties=[
                "sequence_dna_template_probe",
                *readout_sequence_properties,
            ],
            ascending=True,
            filename="seqfish_plus_probes_order",
        )


############################################
# SeqFish Plus Target Probe Designer
############################################


class TargetProbeDesigner:
    """
    Design the RNA-binding part of seqFISH+ hybridization probes.

    This class designs the gene-specific target sequences that bind directly to
    RNA. Downstream steps later add the readout barcode sequences that encode the
    gene identity. Several target probes are usually selected for each gene so the
    combined signal is bright enough for reliable detection.

    The workflow has four main steps:

    1. **Candidate generation**

       Build candidate target probes from transcript or target-region FASTA
       files.

    2. **Sequence filtering**

       Remove candidates with unsuitable sequence properties, such as masked
       sequence, long single-base runs, unsuitable GC content, or strong
       secondary structure.

    3. **Specificity filtering**

       Remove candidates that may bind to unintended targets or cross-hybridize
       with other probes in the panel.

    4. **Probe set selection**

       Select final probe sets for each target region, using criteria such as GC
       content, UTR overlap, and probe spacing.

    :param dir_output: Directory where output files and intermediate results are
        saved.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the TargetProbeDesigner class."""

        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    def generate_target_probes(
        self,
        target_probes_parameters: dict,
        write_intermediate_steps: bool = False,
    ) -> OligoDatabase:
        """
        Run the full seqFISH+ target-probe design workflow.

        This method designs the RNA-binding target sequences used in seqFISH+
        hybridization probes. It starts from transcript or target-region sequences,
        creates candidate probes, filters them, checks their specificity, and
        selects final probe sets for each target region.

        :param target_probes_parameters: Settings from the ``target_probes`` section
            of the pipeline config. This includes candidate generation, sequence
            filters, specificity filters, and probe set selection.
        :type target_probes_parameters: dict
        :param write_intermediate_steps: If ``True``, save intermediate probe
            databases after each main step. This can help when checking where probes
            were removed.
        :type write_intermediate_steps: bool
        :return: Database containing selected seqFISH+ target probes for each target
            region.
        :rtype: OligoDatabase
        """
        oligo_generation_parameters = target_probes_parameters["oligo_generation"]
        property_filters_parameters = target_probes_parameters["property_filters"]
        specificity_filters_parameters = target_probes_parameters["specificity_filters"]
        probe_set_selection_parameters = target_probes_parameters["probe_set_selection"]

        oligo_database: OligoDatabase = self._create_oligo_database(
            region_ids=oligo_generation_parameters["region_ids"],
            oligo_length_min=oligo_generation_parameters["probe_length_min"],
            oligo_length_max=oligo_generation_parameters["probe_length_max"],
            files_fasta_oligo_database=oligo_generation_parameters["files_fasta_probe_database"],
            min_oligos_per_gene=probe_set_selection_parameters["independent_set_selection"]["set_size_min"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_target_probes_initial")
            logger.info(
                f"Saved target probe database for step 1 (Create Database) in directory {dir_database}"
            )

        oligo_database = self._filter_by_property(
            oligo_database=oligo_database,
            isoform_consensus_filter=property_filters_parameters["isoform_consensus_filter"],
            hard_masked_sequences_filter=property_filters_parameters["hard_masked_sequences_filter"],
            soft_masked_sequences_filter=property_filters_parameters["soft_masked_sequences_filter"],
            homopolymeric_runs_filter=property_filters_parameters["homopolymeric_runs_filter"],
            GC_content_filter=property_filters_parameters["GC_content_filter"],
            secondary_structure_filter=property_filters_parameters["secondary_structure_filter"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_target_probes_property_filter")
            logger.info(
                f"Saved target probe database for step 2 (Property Filters) in directory {dir_database}"
            )

        oligo_database = self._filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=specificity_filters_parameters["specificity_blastn_filter"],
            cross_hybridization_blastn_filter=specificity_filters_parameters[
                "cross_hybridization_blastn_filter"
            ],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_target_probes_specificity_filter")
            logger.info(
                f"Saved target probe database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        oligo_database = self._create_oligo_sets(
            oligo_database=oligo_database,
            independent_set_selection=probe_set_selection_parameters["independent_set_selection"],
            GC_content_score=probe_set_selection_parameters["GC_content_score"],
            UTR_score=probe_set_selection_parameters["UTR_score"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_target_probes_sets")
            logger.info(
                f"Saved target probe database for step 4 (Set Selection) in directory {dir_database}."
            )

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Create Database")
    def _create_oligo_database(
        self,
        region_ids: list | None,
        oligo_length_min: int,
        oligo_length_max: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
    ) -> OligoDatabase:
        """
        Create the first database of candidate seqFISH+ target probes.

        Candidate probes are generated by sliding windows across the input
        sequences. All probe lengths between the minimum and maximum length are
        considered. For each candidate, the transcript-facing sequence is stored,
        and the reverse complement is stored as the DNA probe sequence that will
        bind to the RNA.

        Regions with too few candidate probes are removed at this stage.

        :param region_ids: Target regions to design probes for, usually gene names
            or gene IDs. If ``None``, all regions in the input FASTA files are used.
        :type region_ids: list[str] | None
        :param oligo_length_min: Minimum candidate probe length in bases.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum candidate probe length in bases.
        :type oligo_length_max: int
        :param files_fasta_oligo_database: FASTA files containing the transcript or
            target-region sequences used for probe design.
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of candidate probes a region must
            have to remain in the database.
        :type min_oligos_per_gene: int
        :return: Database containing candidate probes with the target sequence and
            the DNA probe sequence.
        :rtype: OligoDatabase
        """
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            region_ids=region_ids,
            n_jobs=self.n_jobs,
        )

        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_gene,
            write_regions_with_insufficient_oligos=True,
            max_entries_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            database_overwrite=True,
            sequence_type="target",
            region_ids=region_ids,
        )
        oligo_database.set_database_sequence_types(["target", "oligo"])

        # Probe strand is the reverse complement of the transcript ("target") window.
        reverse_complement_sequence_property = ReverseComplementSequenceProperty(
            sequence_type_reverse_complement="oligo"
        )
        calculator = PropertyCalculator(properties=[reverse_complement_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Database Creation")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Property Filters")
    def _filter_by_property(
        self,
        oligo_database: OligoDatabase,
        isoform_consensus_filter: dict,
        hard_masked_sequences_filter: dict,
        soft_masked_sequences_filter: dict,
        homopolymeric_runs_filter: dict,
        GC_content_filter: dict,
        secondary_structure_filter: dict,
    ) -> OligoDatabase:
        """
        Remove candidate probes with unsuitable sequence properties.

        This step checks whether each candidate probe is likely to behave well in
        the experiment. It can remove probes that overlap masked sequence, contain
        long single-base runs, have unsuitable GC content, or are predicted to fold
        strongly onto themselves.

        If isoform consensus filtering is enabled, probes are also checked for how
        well they represent the annotated isoforms of the target gene.

        :param oligo_database: Candidate probe database returned by
            :py:meth:`_create_oligo_database`. This database is updated by the
            filtering step.
        :type oligo_database: OligoDatabase
        :param isoform_consensus_filter: Settings for keeping probes that target a
            sufficient fraction of annotated isoforms.
        :type isoform_consensus_filter: dict
        :param hard_masked_sequences_filter: Settings for removing probes that
            overlap hard-masked bases, such as ``N`` bases.
        :type hard_masked_sequences_filter: dict
        :param soft_masked_sequences_filter: Settings for removing probes that
            overlap soft-masked sequence, often used for repetitive or low-complexity
            regions.
        :type soft_masked_sequences_filter: dict
        :param homopolymeric_runs_filter: Settings for removing probes with long
            runs of the same base.
        :type homopolymeric_runs_filter: dict
        :param GC_content_filter: Settings for the allowed GC-content range.
        :type GC_content_filter: dict
        :param secondary_structure_filter: Settings for removing probes predicted to
            form stable self-structures.
        :type secondary_structure_filter: dict
        :return: Filtered database containing probes that passed the enabled
            sequence-property checks.
        :rtype: OligoDatabase
        """
        # Cheap property lookup first; drop weak isoform coverage before sequence work.
        if isoform_consensus_filter["enabled"]:
            isoform_consensus_property = IsoformConsensusProperty()
            calculator = PropertyCalculator(properties=[isoform_consensus_property])
            oligo_database = calculator.apply(
                oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
            )
            oligo_database.filter_database_by_property_threshold(
                property_name="isoform_consensus",
                property_thr=isoform_consensus_filter["isoform_consensus"],
                remove_if_smaller_threshold=True,
            )

        filters: list[BasePropertyFilter] = []
        if hard_masked_sequences_filter["enabled"]:
            hard_masked_sequences = HardMaskedSequenceFilter()
            filters.append(hard_masked_sequences)

        if soft_masked_sequences_filter["enabled"]:
            soft_masked_sequences = SoftMaskedSequenceFilter()
            filters.append(soft_masked_sequences)

        if homopolymeric_runs_filter["enabled"]:
            homopolymeric_runs = HomopolymericRunsFilter(
                base_n=homopolymeric_runs_filter["homopolymeric_base_n"],
            )
            filters.append(homopolymeric_runs)

        if GC_content_filter["enabled"]:
            gc_content = GCContentFilter(
                GC_content_min=GC_content_filter["GC_content_min"],
                GC_content_max=GC_content_filter["GC_content_max"],
            )
            filters.append(gc_content)

        if secondary_structure_filter["enabled"]:
            secondary_sctructure = SecondaryStructureFilter(
                T=secondary_structure_filter["T"],
                thr_DG=secondary_structure_filter["thr_DG"],
            )
            filters.append(secondary_sctructure)

        # Filters were queued cheapest-first so failing probes exit before thermodynamics.
        property_filter = PropertyFilter(filters=filters)
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Specificity Filters")
    def _filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        specificity_blastn_filter: dict,
        cross_hybridization_blastn_filter: dict,
    ) -> OligoDatabase:
        """
        Remove probes that may bind to the wrong place.

        This step checks whether candidate probes are specific to their intended RNA
        target. It removes exact duplicate probe sequences and, when enabled, uses
        BLASTN to find probes that may also bind to unintended reference sequences.

        The method can also remove probes that may cross-hybridize with other probes
        in the same panel.

        :param oligo_database: Probe database returned by
            :py:meth:`_filter_by_property`. This database is updated by the
            specificity filters.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Settings for checking probe specificity
            against reference sequences. This includes the reference FASTA files and
            BLASTN search settings.
        :type specificity_blastn_filter: dict
        :param cross_hybridization_blastn_filter: Settings for checking whether
            probes in the same panel may bind to each other or to unintended probe
            targets.
        :type cross_hybridization_blastn_filter: dict
        :return: Filtered database containing probes that passed the enabled
            specificity checks.
        :rtype: OligoDatabase
        """
        exact_matches = ExactMatchFilter(policy=RemoveAllFilterPolicy(), filter_name="exact_match")
        filters: list = [exact_matches]
        directories: list[str] = []

        if specificity_blastn_filter["enabled"]:
            reference_database = ReferenceDatabase(
                database_name=f"{self.subdir_db_reference}_sequences", dir_output=self.dir_output
            )
            reference_database.load_database_from_file(
                files=specificity_blastn_filter["files_fasta_reference_database"],
                file_type="fasta",
                database_overwrite=True,
            )
            specificity = BlastNFilter(
                remove_hits=True,
                search_parameters=specificity_blastn_filter["search_parameters"],
                hit_parameters=specificity_blastn_filter["hit_parameters"],
                filter_name="specificity_blastn_filter",
                dir_output=self.dir_output,
            )
            specificity.set_reference_database(reference_database=reference_database)
            filters.append(specificity)
            directories.append(specificity.dir_output)

        if cross_hybridization_blastn_filter["enabled"]:
            cross_hybridization_aligner = BlastNFilter(
                remove_hits=True,
                search_parameters=cross_hybridization_blastn_filter["search_parameters"],
                hit_parameters=cross_hybridization_blastn_filter["hit_parameters"],
                filter_name="cross_hybridization_blastn_filter",
                dir_output=self.dir_output,
            )
            cross_hybridization = CrossHybridizationFilter(
                policy=RemoveByLargerRegionFilterPolicy(),
                alignment_method=cross_hybridization_aligner,
                filter_name="cross_hybridization_blastn_filter",
                dir_output=self.dir_output,
            )
            filters.append(cross_hybridization)
            directories.append(cross_hybridization_aligner.dir_output)
            directories.append(cross_hybridization.dir_output)

        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        for directory in directories:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Set Selection")
    def _create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        independent_set_selection: dict,
        GC_content_score: dict,
        UTR_score: dict,
    ) -> OligoDatabase:
        """
        Select final seqFISH+ target-probe sets.

        This step chooses groups of probes from the filtered candidates. The selected
        probes should be well spaced along the target region and should meet the
        requested number of probes per target.

        Probe sets are scored using GC content and UTR overlap. The method can keep
        more than one possible probe set per target region, which gives users
        alternatives when several good designs are available. Regions without enough
        suitable probes are removed.

        :param oligo_database: Filtered probe database returned by
            :py:meth:`_filter_by_specificity`. This database is updated with the
            selected probe sets.
        :type oligo_database: OligoDatabase
        :param independent_set_selection: Settings that control how many probe sets
            are selected, how many probes each set should contain, and how far apart
            selected probes should be placed.
        :type independent_set_selection: dict
        :param GC_content_score: Settings for scoring probes by how close their GC
            content is to the desired value.
        :type GC_content_score: dict
        :param UTR_score: Settings for scoring probes by how well they overlap the
            target region annotation.
        :type UTR_score: dict
        :return: Database with selected probe sets attached to each remaining target
            region.
        :rtype: OligoDatabase
        """
        utr_scorer = OverlapUTRScorer(score_weight=UTR_score["weight"], property_name="regiontype")
        GC_scorer = DeviationFromOptimalGCContentScorer(
            GC_content_opt=GC_content_score["GC_content_opt"],
            score_weight=GC_content_score["weight"],
        )
        oligos_scoring = OligoScoring(scorers=[utr_scorer, GC_scorer])
        set_scoring = LowestSetScoring(ascending=True)

        base_log_parameters({"Set Selection": "Independent Sets"})
        oligoset_generator = IndependentSetsOligoSelection(
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            set_size_opt=independent_set_selection["set_size_opt"],
            set_size_min=independent_set_selection["set_size_min"],
            distance_between_oligos=independent_set_selection["distance_between_probes"],
            n_attempts_graph=independent_set_selection["n_attempts_graph"],
            n_attempts_clique_enum=independent_set_selection["n_attempts_clique_enum"],
            diversification_fraction=independent_set_selection["diversification_fraction"],
            jaccard_opt=independent_set_selection["jaccard_opt"],
            jaccard_step=independent_set_selection["jaccard_step"],
        )
        oligo_database = oligoset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=independent_set_selection["n_sets"],
            n_jobs=self.n_jobs,
        )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Oligo Selection")
        check_content_oligo_database(oligo_database)

        return oligo_database


############################################
# SeqFish Plus Readout Probe Designer
############################################


class ReadoutProbeDesigner:
    """
    Manage the readout probes and codebook used for seqFISH+ decoding.

    In seqFISH+, each target gene is assigned a combinatorial barcode. The
    barcode is read over several hybridization rounds. In each round, fluorescent
    readout probes bind to selected barcode sequences on the hybridization
    probes. The pattern of signals across rounds identifies the RNA molecule.

    This class handles the two tables needed for that assignment:

    - the **codebook**, which assigns each target region to a barcode built from
      pseudocolors across imaging rounds and fluorescence channels,
    - the **readout probe table**, which links each barcode bit to a readout probe
      sequence, barcode round, pseudocolor, and channel.

    Both tables can be loaded from files or generated from the config. When the
    readout probe table is generated, candidate sequences are drawn from a random
    DNA pool and filtered for sequence quality and specificity. Before the tables
    are used, they are checked together to make sure that all requested target
    regions and barcode bits are covered.

    :param dir_output: Directory where output files and intermediate results are
        saved.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the ReadoutProbeDesigner class."""

        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_readout_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    def load_codebook(self, file_codebook: str) -> pd.DataFrame:
        """
        Load a seqFISH+ codebook from a table file.

        The codebook file should contain one row per target region and one column
        per barcode bit. The target region column must be named ``gene_name``. Bit
        columns are expected to be named ``bit_1``, ``bit_2``, and so on.

        The table is loaded here, but the full content check is done later together
        with the readout probe table. This makes sure that the codebook bits and the
        available readout probes match.

        :param file_codebook: Path to the codebook file. The file must contain a
            ``gene_name`` column and one or more ``bit_*`` columns.
        :type file_codebook: str
        :return: Codebook table indexed by ``gene_name``.
        :rtype: pd.DataFrame
        """
        return pd.read_csv(file_codebook, sep=None, engine="python", index_col="gene_name")

    def generate_codebook(self, region_ids: list[str], codebook_parameters: dict) -> pd.DataFrame:
        """
        Generate a seqFISH+ codebook from the codebook settings.

        The generated codebook assigns each target region to a combinatorial
        barcode. Each row is a target region. Each column is a barcode bit. A value
        of ``1`` means that the target carries the readout-binding sequence for that
        bit.

        seqFISH+ barcodes are built from pseudocolors across several barcode rounds.
        In each round, one pseudocolor is selected and that choice is tied to one
        fluorescence channel. The first ``n_barcode_rounds - 1`` rounds are chosen
        freely. The last round is derived from the earlier rounds as a parity check.
        That final round helps catch barcode patterns that do not match the expected
        relation during decoding.

        Bit positions are arranged by barcode round, then pseudocolor, then channel.
        Unused bit columns are kept as all-zero columns so the codebook stays aligned
        with the full readout probe set.

        :param region_ids: Target regions that need barcode assignments, usually
            gene names or gene IDs.
        :type region_ids: list[str]
        :param codebook_parameters: Settings from the ``readout_probes.codebook``
            section of the pipeline config. This includes ``n_barcode_rounds``,
            ``n_pseudocolors``, and ``channels_ids``.
        :type codebook_parameters: dict
        :return: Codebook table with target regions as rows and barcode bits as
            columns. Each row has one active bit per barcode round.
        :rtype: pd.DataFrame
        :raises ConfigurationError: If there are not enough valid barcodes for all
            requested target regions.
        """

        def _generate_barcode(
            pseudocolors: list, channel: int, n_pseudocolors: int, n_channels: int
        ) -> np.ndarray:
            """
            Build one seqFISH+ barcode from free-round pseudocolors and a channel.

            The last round is filled in as a parity check over the free rounds. That
            extra round helps catch decoding errors later. The result has one active
            bit per barcode round.

            :param pseudocolors: Pseudocolor indices for the free barcode rounds.
            :type pseudocolors: list
            :param channel: Fluorescence channel index used for the barcode.
            :type channel: int
            :param n_pseudocolors: Total number of pseudocolors available.
            :type n_pseudocolors: int
            :param n_channels: Total number of fluorescence channels available.
            :type n_channels: int
            :return: Barcode vector with one active bit per barcode round.
            :rtype: np.ndarray
            :raises ConfigurationError: If a pseudocolor or channel index falls
                outside the configured ranges.
            """
            # Final-round pseudocolor is a parity checksum over the free rounds.
            pseudocolors = pseudocolors + [sum(pseudocolors) % n_pseudocolors]
            if max(pseudocolors) >= n_pseudocolors:
                raise ConfigurationError(
                    f"Barcode references pseudocolor index {max(pseudocolors)} but only "
                    f"{n_pseudocolors} pseudocolors are configured."
                )
            if channel >= n_channels:
                raise ConfigurationError(
                    f"Barcode references channel index {channel} but only "
                    f"{n_channels} channels are configured."
                )
            n_barcode_rounds = len(pseudocolors)
            barcode = np.zeros(n_channels * n_pseudocolors * n_barcode_rounds, dtype=np.int8)
            for i, pseudocolor in enumerate(pseudocolors):
                # Flat bit index: barcode round, then pseudocolor, then channel.
                barcode[i * n_pseudocolors * n_channels + n_channels * pseudocolor + channel] = 1
            return barcode

        n_regions = len(region_ids)
        n_pseudocolors = codebook_parameters["n_pseudocolors"]
        n_barcode_rounds = codebook_parameters["n_barcode_rounds"]
        n_channels = len(codebook_parameters["channels_ids"])

        codebook_list: list[np.ndarray] = []
        codebook_size = n_channels * (n_pseudocolors ** (n_barcode_rounds - 1))
        barcode_size = n_pseudocolors * n_barcode_rounds * n_channels

        if codebook_size < n_regions:
            raise ConfigurationError(
                f"The number of valid barcodes ({codebook_size}) is lower than the number of regions ({n_regions}). "
                f"Consider increasing the number of pseudocolors or barcoding rounds, or reducing the number of regions."
            )
        for pseudocolors in product(range(n_pseudocolors), repeat=n_barcode_rounds - 1):
            for channel in range(n_channels):
                barcode = _generate_barcode(
                    pseudocolors=list(pseudocolors),
                    channel=channel,
                    n_pseudocolors=n_pseudocolors,
                    n_channels=n_channels,
                )
                codebook_list.append(barcode)

        codebook: pd.DataFrame = pd.DataFrame(
            codebook_list[0:n_regions], index=region_ids, columns=[f"bit_{i+1}" for i in range(barcode_size)]
        )
        codebook.index.name = "gene_name"

        return codebook

    def load_readout_probe_table(self, file_readout_probe_table: str) -> pd.DataFrame:
        """
        Load a seqFISH+ readout probe table from a table file.

        The readout probe table links each barcode bit to a readout probe sequence.
        It must contain the bit name, barcode round, pseudocolor, fluorescence
        channel, and readout probe sequence.

        Column and sequence checks are done later by :py:meth:`validate`. When a
        codebook is also provided there, every codebook bit must appear in this
        table.

        :param file_readout_probe_table: Path to the readout probe table file. The
            file must contain ``bit``, ``barcode_round``, ``pseudocolor``,
            ``channel``, and ``readout_probe_sequence`` columns.
        :type file_readout_probe_table: str
        :return: Readout probe table indexed by ``bit``.
        :rtype: pd.DataFrame
        :raises FileFormatError: If the ``bit`` column is missing.
        """
        readout_probe_table = pd.read_csv(file_readout_probe_table, sep=None, engine="python")
        if "bit" not in readout_probe_table.columns:
            raise FileFormatError(
                f"Readout probe table '{file_readout_probe_table}' must contain a 'bit' column."
            )
        return readout_probe_table.set_index("bit")

    def generate_readout_probe_table(
        self,
        readout_probe_parameters: dict,
        codebook_parameters: dict,
        write_intermediate_steps: bool = False,
    ) -> pd.DataFrame:
        """
        Generate a seqFISH+ readout probe table.

        This method designs short DNA sequences that will serve as readout probes.
        Candidate sequences are drawn from a random DNA pool and filtered for
        sequence quality and specificity. The surviving probes are written into a
        bit-indexed table with barcode round, pseudocolor, and channel information.

        :param readout_probe_parameters: Settings from the
            ``readout_probes.readout_probe_table`` section of the pipeline config.
            This includes candidate generation, sequence filters, and specificity
            filters.
        :type readout_probe_parameters: dict
        :param codebook_parameters: Settings from the ``readout_probes.codebook``
            section of the pipeline config. The barcode-round, pseudocolor, and
            channel settings define how many readout probes are written into the
            final table.
        :type codebook_parameters: dict
        :param write_intermediate_steps: If ``True``, save intermediate probe
            databases after each main step. This can help when checking where probes
            were removed.
        :type write_intermediate_steps: bool
        :return: Readout probe table indexed by ``bit``, with barcode round,
            pseudocolor, channel, and sequence columns.
        :rtype: pd.DataFrame
        """
        oligo_database: OligoDatabase = self._create_oligo_database(
            oligo_length=readout_probe_parameters["oligo_generation"]["probe_length"],
            oligo_base_probabilities=readout_probe_parameters["oligo_generation"]["base_probabilities"],
            initial_num_sequences=readout_probe_parameters["oligo_generation"]["initial_num_sequences"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_readout_probes_initial")
            logger.info(
                f"Saved readout probe database for step 1 (Create Database) in directory {dir_database}"
            )

        oligo_database = self._filter_by_property(
            oligo_database=oligo_database,
            GC_content_filter=readout_probe_parameters["property_filters"]["GC_content_filter"],
            homopolymeric_runs_filter=readout_probe_parameters["property_filters"][
                "homopolymeric_runs_filter"
            ],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_readout_probes_property_filter")
            logger.info(
                f"Saved readout probe database for step 2 (Property Filters) in directory {dir_database}"
            )

        oligo_database = self._filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=readout_probe_parameters["specificity_filters"][
                "specificity_blastn_filter"
            ],
            cross_hybridization_blastn_filter=readout_probe_parameters["specificity_filters"][
                "cross_hybridization_blastn_filter"
            ],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_readout_probes_specificty_filter")
            logger.info(
                f"Saved readout probe database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        return self._format_readout_probe_table(
            readout_probe_database=oligo_database,
            channels_ids=codebook_parameters["channels_ids"],
            n_barcode_rounds=codebook_parameters["n_barcode_rounds"],
            n_pseudocolors=codebook_parameters["n_pseudocolors"],
        )

    def validate(
        self,
        readout_probe_table: pd.DataFrame | None = None,
        *,
        readout_probe_table_source: str | None = None,
        codebook: pd.DataFrame | None = None,
        region_ids: list[str] | None = None,
        codebook_source: str | None = None,
        hamming_weight: int | None = None,
    ) -> None:
        """
        Check the codebook and/or readout probe table.

        This method can check either table on its own, or both together. When a
        codebook is provided, it must contain all requested target regions, use
        ``gene_name`` as the row index, and contain only ``0`` and ``1`` values in
        its bit columns. Each target region must have exactly one active bit per
        barcode round. Pass that round count as ``hamming_weight`` (the same value
        as ``n_barcode_rounds`` in the codebook config).

        When a readout probe table is provided, it must include a valid DNA
        sequence, barcode round, pseudocolor, and channel for each bit. If a
        codebook is also provided, every bit used by the codebook must appear in
        the table.

        Running these checks before probe assembly helps catch mismatched files
        early, such as a codebook that refers to a bit missing from the readout
        probe table.

        :param readout_probe_table: Optional table linking each barcode bit to a
            readout probe sequence and its readout information. If ``None``, the
            table is not checked.
        :type readout_probe_table: pd.DataFrame | None
        :param readout_probe_table_source: File path or source label for the
            readout probe table. Used in error messages when the table is checked.
        :type readout_probe_table_source: str | None
        :param codebook: Optional codebook table assigning target regions to
            barcode bits. Rows are target regions and columns are ``bit_*``
            entries. If ``None``, the codebook is not checked.
        :type codebook: pd.DataFrame | None
        :param region_ids: Target regions that must be present in the codebook.
            Required when a codebook is checked.
        :type region_ids: list[str] | None
        :param codebook_source: File path or source label for the codebook. Used
            in error messages when the codebook is checked.
        :type codebook_source: str | None
        :param hamming_weight: Expected number of active bits in each codebook
            row. Required when a codebook is checked. For seqFISH+, pass
            ``n_barcode_rounds`` from the codebook config.
        :type hamming_weight: int | None
        :return: None
        :rtype: None
        :raises ValueError: If a table is checked without its required companion
            arguments (for example ``region_ids`` / ``codebook_source`` /
            ``hamming_weight`` with a codebook, or ``readout_probe_table_source``
            with a readout probe table).
        :raises FileFormatError: If the codebook or readout probe table is missing
            required information or contains invalid values.
        """
        if codebook is not None:
            if region_ids is None:
                raise ValueError("region_ids must be provided when validating a codebook.")
            if codebook_source is None:
                raise ValueError("codebook_source must be provided when validating a codebook.")
            if hamming_weight is None:
                raise ValueError("hamming_weight must be provided when validating a codebook.")
            validate_codebook(
                codebook=codebook,
                region_ids=region_ids,
                source=codebook_source,
                expected_hamming_weight=hamming_weight,
                index_name="gene_name",
            )
        if readout_probe_table is not None:
            if readout_probe_table_source is None:
                raise ValueError(
                    "readout_probe_table_source must be provided when validating a readout probe table."
                )
            validate_bit_mapping_table(
                table=readout_probe_table,
                codebook=codebook,
                source=readout_probe_table_source,
                required_columns=["barcode_round", "pseudocolor", "channel", "readout_probe_sequence"],
                sequence_columns=["readout_probe_sequence"],
            )

    @pipeline_step_basic(step_name="Readout Probe Generation - Create Oligo Database")
    def _create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
    ) -> OligoDatabase:
        """
        Create the first database of candidate readout probes.

        Candidate sequences are drawn from a random DNA pool. Readout probes do not
        bind endogenous RNA directly. They bind barcode sequences carried by the
        hybridization probes, so no transcript scaffold is needed here.

        :param oligo_length: Length of each candidate readout probe in bases.
        :type oligo_length: int
        :param oligo_base_probabilities: Probabilities used to draw random DNA
            sequences. Keys are ``A``, ``T``, ``C``, and ``G``.
        :type oligo_base_probabilities: dict
        :param initial_num_sequences: Number of random candidate sequences to create
            before filtering.
        :type initial_num_sequences: int
        :return: Database containing the random candidate readout probes.
        :rtype: OligoDatabase
        """
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_random(
            filename_out="readout_probe_sequences",
            length_sequences=oligo_length,
            num_sequences=initial_num_sequences,
            name_sequences="readout_probe",
            base_alphabet_with_probability=oligo_base_probabilities,
        )

        oligo_database = OligoDatabase(
            min_oligos_per_region=0,
            write_regions_with_insufficient_oligos=False,
            max_entries_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=None,
        )
        oligo_database.set_database_sequence_types(["oligo"])

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Property Filters")
    def _filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_filter: dict,
        homopolymeric_runs_filter: dict,
    ) -> OligoDatabase:
        """
        Remove candidate readout probes with unsuitable sequence properties.

        This step checks whether each candidate is likely to behave well during
        imaging. It can remove probes with unsuitable GC content or long runs of the
        same base. Keeping the remaining probes in a similar property range helps
        them hybridize under the same imaging conditions.

        :param oligo_database: Candidate probe database returned by
            :py:meth:`_create_oligo_database`. This database is updated by the
            filtering step.
        :type oligo_database: OligoDatabase
        :param GC_content_filter: Settings for the allowed GC-content range.
        :type GC_content_filter: dict
        :param homopolymeric_runs_filter: Settings for removing probes with long
            runs of the same base.
        :type homopolymeric_runs_filter: dict
        :return: Filtered database containing probes that passed the enabled
            sequence-property checks.
        :rtype: OligoDatabase
        """
        filters: list[BasePropertyFilter] = []

        if homopolymeric_runs_filter["enabled"]:
            homopolymeric_runs = HomopolymericRunsFilter(
                base_n=homopolymeric_runs_filter["homopolymeric_base_n"],
            )
            filters.append(homopolymeric_runs)

        if GC_content_filter["enabled"]:
            gc_content = GCContentFilter(
                GC_content_min=GC_content_filter["GC_content_min"],
                GC_content_max=GC_content_filter["GC_content_max"],
            )
            filters.append(gc_content)

        property_filter = PropertyFilter(filters=filters)
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Specificity Filters")
    def _filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        specificity_blastn_filter: dict,
        cross_hybridization_blastn_filter: dict,
    ) -> OligoDatabase:
        """
        Remove readout probes that may bind to the wrong place.

        This step checks whether candidate readout probes are specific enough for
        imaging. It removes exact duplicate sequences and, when enabled, uses BLASTN
        to find probes that may also bind to unintended reference sequences.

        This matters for seqFISH+ because a readout probe that matches endogenous
        RNA can create false signal during an imaging round. The method can also
        remove probes that may bind to each other instead of to the barcode
        sequences on the hybridization probes.

        :param oligo_database: Probe database returned by
            :py:meth:`_filter_by_property`. This database is updated by the
            specificity filters.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Settings for checking probe specificity
            against reference sequences. This includes the reference FASTA files and
            BLASTN search settings.
        :type specificity_blastn_filter: dict
        :param cross_hybridization_blastn_filter: Settings for checking whether
            readout probes in the same panel may bind to each other.
        :type cross_hybridization_blastn_filter: dict
        :return: Filtered database containing probes that passed the enabled
            specificity checks.
        :rtype: OligoDatabase
        """
        exact_matches = ExactMatchFilter(
            policy=RemoveAllFilterPolicy(), filter_name="readout_probes_exact_match"
        )
        filters: list = [exact_matches]
        directories: list[str] = []

        if specificity_blastn_filter["enabled"]:
            reference_database = ReferenceDatabase(
                database_name=self.subdir_db_reference, dir_output=self.dir_output
            )
            reference_database.load_database_from_file(
                files=specificity_blastn_filter["files_fasta_reference_database"],
                file_type="fasta",
                database_overwrite=True,
            )
            specificity = BlastNFilter(
                remove_hits=True,
                search_parameters=specificity_blastn_filter["search_parameters"],
                hit_parameters=specificity_blastn_filter["hit_parameters"],
                filter_name="readout_probes_blastn_specificity",
                dir_output=self.dir_output,
            )
            specificity.set_reference_database(reference_database=reference_database)
            filters.append(specificity)
            directories.append(specificity.dir_output)

        if cross_hybridization_blastn_filter["enabled"]:
            cross_hybridization_aligner = BlastNFilter(
                remove_hits=True,
                search_parameters=cross_hybridization_blastn_filter["search_parameters"],
                hit_parameters=cross_hybridization_blastn_filter["hit_parameters"],
                filter_name="readout_probes_blastn_crosshybridization",
                dir_output=self.dir_output,
            )
            cross_hybridization = CrossHybridizationFilter(
                policy=RemoveByDegreeFilterPolicy(),
                alignment_method=cross_hybridization_aligner,
                filter_name="readout_probes_blastn_crosshybridization",
                dir_output=self.dir_output,
            )
            filters.append(cross_hybridization)
            directories.append(cross_hybridization_aligner.dir_output)
            directories.append(cross_hybridization.dir_output)

        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        for directory in directories:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        check_content_oligo_database(oligo_database)
        return oligo_database

    def _format_readout_probe_table(
        self,
        readout_probe_database: OligoDatabase,
        channels_ids: list[str],
        n_barcode_rounds: int,
        n_pseudocolors: int,
    ) -> pd.DataFrame:
        """
        Build the bit-indexed seqFISH+ readout probe table.

        This step assigns each selected readout probe to one barcode bit. Bits are
        filled in a fixed order: barcode round first, then pseudocolor, then
        channel. That order keeps the table aligned with the codebook.

        :param readout_probe_database: Selected readout probe database returned by
            :py:meth:`_filter_by_specificity`.
        :type readout_probe_database: OligoDatabase
        :param channels_ids: Fluorescence channel names to assign across barcode
            bits.
        :type channels_ids: list[str]
        :param n_barcode_rounds: Number of barcode rounds in the experiment.
        :type n_barcode_rounds: int
        :param n_pseudocolors: Number of pseudocolors used in each barcode round.
        :type n_pseudocolors: int
        :return: Readout probe table indexed by ``bit``, with barcode round,
            pseudocolor, channel, and sequence columns.
        :rtype: pd.DataFrame
        :raises ConfigurationError: If fewer readout probes are available than are
            needed to fill all barcode bits.
        """
        n_channels = len(channels_ids)
        n_bits = n_barcode_rounds * n_pseudocolors * n_channels
        readout_probes = readout_probe_database.get_oligoid_sequence_mapping(
            sequence_type="oligo", sequence_to_upper=False
        )
        if len(readout_probes) < n_bits:
            raise ConfigurationError(
                f"Only {len(readout_probes)} readout probes are available but "
                f"{n_bits} bits are required. Increase ``initial_num_sequences`` or "
                f"relax the readout-probe filters to yield more candidates."
            )
        readout_probe_table = pd.DataFrame(
            columns=[
                "bit",
                "barcode_round",
                "pseudocolor",
                "channel",
                "readout_probe_sequence",
            ],
            index=list(range(n_bits)),
        )
        barcode_round = 0
        pseudocolor = 0
        channel = 0
        for i, readout_probe_sequence in enumerate(readout_probes.values()):
            readout_probe_table.iloc[i] = pd.Series(
                [
                    f"bit_{i + 1}",
                    barcode_round + 1,
                    pseudocolor + 1,
                    channels_ids[channel],
                    readout_probe_sequence,
                ],
                index=[
                    "bit",
                    "barcode_round",
                    "pseudocolor",
                    "channel",
                    "readout_probe_sequence",
                ],
            )
            channel = (channel + 1) % n_channels
            if channel == 0:
                pseudocolor = (pseudocolor + 1) % n_pseudocolors
                if pseudocolor == 0:
                    barcode_round = (barcode_round + 1) % n_barcode_rounds

            # Extra random probes may exist; stop once every barcode bit has a sequence.
            if i >= n_bits - 1:
                break
        readout_probe_table.set_index("bit", inplace=True)
        return readout_probe_table


############################################
# SeqFish Plus Primer Designer
############################################


class PrimerDesigner:
    """
    Load or design the PCR primers used for seqFISH+ DNA templates.

    seqFISH+ hybridization probes are prepared from DNA template sequences. These
    templates contain the probe body flanked by a forward primer and a reverse
    primer. The primer pair is used to amplify the pooled template library before
    the final single-stranded probes are prepared.

    The reverse primer is usually provided in the pipeline config. The forward
    primer can also be loaded from the config, or it can be designed automatically.
    Automatic forward-primer design draws random candidates, filters them for
    sequence quality and specificity, and then selects the candidate whose melting
    temperature best matches the reverse primer.

    :param dir_output: Directory where output files and intermediate results are
        saved.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the PrimerDesigner class."""

        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_primer"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    def load_reverse_primer(self, sequence: str) -> str:
        """
        Load the reverse primer sequence from the config.

        The reverse primer is used together with the forward primer to amplify the
        DNA template probe pool. This method trims surrounding whitespace. The
        sequence itself is checked later by :py:meth:`validate`.

        :param sequence: Reverse primer sequence used for PCR amplification of the
            DNA template probes.
        :type sequence: str
        :return: Cleaned reverse primer sequence.
        :rtype: str
        """
        return str(sequence).strip()

    def generate_reverse_primer(self) -> str:
        """
        Generate a reverse primer sequence.

        Automatic reverse primer design is not available yet. Provide the reverse
        primer sequence in the pipeline config instead.

        :return: Reverse primer sequence.
        :rtype: str
        :raises FeatureNotImplementedError: Always raised until reverse primer
            generation is implemented.
        """
        raise FeatureNotImplementedError(
            "Generation of the reverse primer is not yet implemented. "
            "Please provide a reverse_primer.sequence parameter and set reverse_primer.source to 'load'."
        )

    def load_forward_primer(self, sequence: str) -> str:
        """
        Load the forward primer sequence from the config.

        The forward primer is used together with the reverse primer to amplify the
        DNA template probe pool. This method trims surrounding whitespace. The
        sequence itself is checked later by :py:meth:`validate`.

        :param sequence: Forward primer sequence used for PCR amplification of the
            DNA template probes.
        :type sequence: str
        :return: Cleaned forward primer sequence.
        :rtype: str
        """
        return str(sequence).strip()

    def generate_forward_primer(
        self,
        parameters: dict,
        reverse_primer_sequence: str,
        file_fasta_hybridization_probes_database: str,
        write_intermediate_steps: bool = False,
    ) -> str:
        """
        Generate a forward primer matched to the reverse primer.

        This method designs a forward primer for amplifying the DNA template probe
        pool. Candidate sequences are drawn from a random DNA pool, filtered for
        sequence quality, and checked for unwanted binding to reference sequences
        or to the assembled hybridization probes. The surviving candidate whose
        melting temperature best matches the reverse primer is returned.

        :param parameters: Settings from the ``primers.forward_primer`` section of
            the pipeline config. This includes candidate generation, sequence
            filters, and specificity filters.
        :type parameters: dict
        :param reverse_primer_sequence: Reverse primer sequence used as the melting
            temperature reference for selecting the forward primer.
        :type reverse_primer_sequence: str
        :param file_fasta_hybridization_probes_database: FASTA file containing the
            assembled hybridization probe sequences. Used to avoid primers that bind
            the probe body.
        :type file_fasta_hybridization_probes_database: str
        :param write_intermediate_steps: If ``True``, save intermediate primer
            databases after each main step. This can help when checking where
            candidates were removed.
        :type write_intermediate_steps: bool
        :return: Selected forward primer sequence.
        :rtype: str
        """
        oligo_database: OligoDatabase = self._create_oligo_database(
            oligo_length=parameters["oligo_generation"]["probe_length"],
            oligo_base_probabilities=parameters["oligo_generation"]["base_probabilities"],
            initial_num_sequences=parameters["oligo_generation"]["initial_num_sequences"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_primers_initial")
            logger.info(f"Saved primer database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = self._filter_by_property(
            oligo_database=oligo_database,
            GC_content_filter=parameters["property_filters"]["GC_content_filter"],
            GC_clamp_filter=parameters["property_filters"]["GC_clamp_filter"],
            homopolymeric_runs_filter=parameters["property_filters"]["homopolymeric_runs_filter"],
            self_complementarity_filter=parameters["property_filters"]["self_complementarity_filter"],
            complement_reverse_primer_filter=parameters["property_filters"][
                "complement_reverse_primer_filter"
            ],
            Tm_filter=parameters["property_filters"]["Tm_filter"],
            secondary_structure_filter=parameters["property_filters"]["secondary_structure_filter"],
            reverse_primer_sequence=reverse_primer_sequence,
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_primer_property_filter")
            logger.info(f"Saved primer database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = self._filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=parameters["specificity_filters"]["specificity_blastn_filter"],
            hybridization_probes_blastn_filter=parameters["specificity_filters"][
                "hybridization_probes_blastn_filter"
            ],
            file_fasta_hybridization_probes_database=file_fasta_hybridization_probes_database,
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_primer_specificty_filter")
            logger.info(f"Saved primer database for step 3 (Specificity Filters) in directory {dir_database}")

        forward_primer_sequence = self._get_best_forward_primer(
            oligo_database=oligo_database,
            reverse_primer_sequence=reverse_primer_sequence,
            Tm_parameters=parameters["property_filters"]["Tm_filter"]["Tm_parameters"],
            Tm_chem_correction_parameters=parameters["property_filters"]["Tm_filter"][
                "Tm_chem_correction_parameters"
            ],
            Tm_salt_correction_parameters=parameters["property_filters"]["Tm_filter"][
                "Tm_salt_correction_parameters"
            ],
        )

        return forward_primer_sequence

    def validate(self, forward_primer: str, reverse_primer: str) -> None:
        """
        Check that both primer sequences are valid DNA sequences.

        Each primer must be a non-empty sequence containing only ``A``, ``C``,
        ``G``, and ``T``. This check catches missing primers and accidental
        characters before the template probes are assembled.

        :param forward_primer: Forward primer sequence to check.
        :type forward_primer: str
        :param reverse_primer: Reverse primer sequence to check.
        :type reverse_primer: str
        :return: None
        :rtype: None
        :raises FileFormatError: If either primer is empty or contains characters
            other than ``A``, ``C``, ``G``, and ``T``.
        """
        validate_primer_sequence(sequence=forward_primer, source="forward_primer")
        validate_primer_sequence(sequence=reverse_primer, source="reverse_primer")

    @pipeline_step_basic(step_name="Primer Generation - Create Oligo Database")
    def _create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
    ) -> OligoDatabase:
        """
        Create the first database of candidate forward primers.

        Candidate sequences are drawn from a random DNA pool. Each sequence is
        generated one base shorter than the requested primer length, and a terminal
        ``T`` is then added. This keeps the 3' end consistent across the candidate
        pool.

        :param oligo_length: Final length of each candidate primer in bases.
        :type oligo_length: int
        :param oligo_base_probabilities: Probabilities used to draw random DNA
            sequences. Keys are ``A``, ``T``, ``C``, and ``G``.
        :type oligo_base_probabilities: dict
        :param initial_num_sequences: Number of random candidate sequences to create
            before filtering.
        :type initial_num_sequences: int
        :return: Database containing the random candidate forward primers.
        :rtype: OligoDatabase
        """
        forward_primer_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        # Draw at length-1 because a fixed 3'-T is appended below.
        forward_primer_fasta_file = forward_primer_sequences.create_sequences_random(
            filename_out="forward_primer_sequences",
            length_sequences=oligo_length - 1,
            num_sequences=initial_num_sequences,
            name_sequences="forward_primer",
            base_alphabet_with_probability=oligo_base_probabilities,
        )
        forward_primer_fasta_file = append_nucleotide_to_sequences(forward_primer_fasta_file, nucleotide="T")

        oligo_database = OligoDatabase(
            min_oligos_per_region=0,
            write_regions_with_insufficient_oligos=False,
            max_entries_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=forward_primer_fasta_file,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=None,
        )
        oligo_database.set_database_sequence_types(["oligo"])

        dir = forward_primer_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

    @pipeline_step_basic(step_name="Primer Generation - Property Filters")
    def _filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_filter: dict,
        GC_clamp_filter: dict,
        homopolymeric_runs_filter: dict,
        self_complementarity_filter: dict,
        complement_reverse_primer_filter: dict,
        Tm_filter: dict,
        secondary_structure_filter: dict,
        reverse_primer_sequence: str,
    ) -> OligoDatabase:
        """
        Remove candidate primers with unsuitable sequence properties.

        This step checks whether each candidate is likely to work well during PCR.
        It can remove primers with unsuitable GC content, missing GC clamps, long
        single-base runs, strong self-complementarity, complementarity to the
        reverse primer, unsuitable melting temperature, or strong secondary
        structure.

        :param oligo_database: Candidate primer database returned by
            :py:meth:`_create_oligo_database`. This database is updated by the
            filtering step.
        :type oligo_database: OligoDatabase
        :param GC_content_filter: Settings for the allowed GC-content range.
        :type GC_content_filter: dict
        :param GC_clamp_filter: Settings for requiring G or C bases near the 3' end
            of the primer.
        :type GC_clamp_filter: dict
        :param homopolymeric_runs_filter: Settings for removing primers with long
            runs of the same base.
        :type homopolymeric_runs_filter: dict
        :param self_complementarity_filter: Settings for removing primers that can
            fold back on themselves.
        :type self_complementarity_filter: dict
        :param complement_reverse_primer_filter: Settings for removing primers that
            can bind strongly to the reverse primer.
        :type complement_reverse_primer_filter: dict
        :param Tm_filter: Settings for the allowed melting-temperature range and the
            conditions used for the calculation.
        :type Tm_filter: dict
        :param secondary_structure_filter: Settings for removing primers predicted to
            form stable self-structures.
        :type secondary_structure_filter: dict
        :param reverse_primer_sequence: Reverse primer sequence used when checking
            primer-primer complementarity.
        :type reverse_primer_sequence: str
        :return: Filtered database containing primers that passed the enabled
            sequence-property checks.
        :rtype: OligoDatabase
        """
        filters: list[BasePropertyFilter] = []

        if homopolymeric_runs_filter["enabled"]:
            homopolymeric_runs = HomopolymericRunsFilter(
                base_n=homopolymeric_runs_filter["homopolymeric_base_n"],
            )
            filters.append(homopolymeric_runs)

        if GC_content_filter["enabled"]:
            gc_content = GCContentFilter(
                GC_content_min=GC_content_filter["GC_content_min"],
                GC_content_max=GC_content_filter["GC_content_max"],
            )
            filters.append(gc_content)

        if GC_clamp_filter["enabled"]:
            gc_clamp = GCClampFilter(
                n_bases=GC_clamp_filter["number_three_prime_base_GCclamp"],
                n_GC=GC_clamp_filter["number_GC_GCclamp"],
            )
            filters.append(gc_clamp)

        if self_complementarity_filter["enabled"]:
            self_complement = SelfComplementFilter(
                max_len_selfcomplement=self_complementarity_filter["max_len_selfcomplement"],
            )
            filters.append(self_complement)

        if complement_reverse_primer_filter["enabled"]:
            complement = ComplementFilter(
                comparison_sequence=reverse_primer_sequence,
                max_len_complement=complement_reverse_primer_filter["max_len_complement"],
            )
            filters.append(complement)

        if Tm_filter["enabled"]:
            melting_temperature = MeltingTemperatureNNFilter(
                Tm_min=Tm_filter["Tm_min"],
                Tm_max=Tm_filter["Tm_max"],
                Tm_parameters=Tm_filter["Tm_parameters"],
                Tm_chem_correction_parameters=Tm_filter["Tm_chem_correction_parameters"],
                Tm_salt_correction_parameters=Tm_filter["Tm_salt_correction_parameters"],
            )
            filters.append(melting_temperature)

        if secondary_structure_filter["enabled"]:
            secondary_sctructure = SecondaryStructureFilter(
                T=secondary_structure_filter["T"],
                thr_DG=secondary_structure_filter["thr_DG"],
            )
            filters.append(secondary_sctructure)

        property_filter = PropertyFilter(filters=filters)
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Primer Generation - Specificity Filters")
    def _filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        specificity_blastn_filter: dict,
        hybridization_probes_blastn_filter: dict,
        file_fasta_hybridization_probes_database: str,
    ) -> OligoDatabase:
        """
        Remove primers that may bind to the wrong place.

        This step checks whether candidate primers are specific enough for library
        amplification. When enabled, it uses BLASTN to find primers that may also
        bind to unintended reference sequences.

        It can also check whether primers bind to the assembled hybridization
        probes. Primers that anneal inside the probe body can interfere with
        amplification of the intended template library.

        :param oligo_database: Primer database returned by
            :py:meth:`_filter_by_property`. This database is updated by the
            specificity filters.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Settings for checking primer specificity
            against reference sequences. This includes the reference FASTA files and
            BLASTN search settings.
        :type specificity_blastn_filter: dict
        :param hybridization_probes_blastn_filter: Settings for checking whether
            primers may bind to the assembled hybridization probes.
        :type hybridization_probes_blastn_filter: dict
        :param file_fasta_hybridization_probes_database: FASTA file containing the
            assembled hybridization probe sequences used for the probe-body check.
        :type file_fasta_hybridization_probes_database: str
        :return: Filtered database containing primers that passed the enabled
            specificity checks.
        :rtype: OligoDatabase
        """
        filters: list = []
        directories: list[str] = []

        if specificity_blastn_filter["enabled"]:
            reference_database = ReferenceDatabase(
                database_name=self.subdir_db_reference, dir_output=self.dir_output
            )
            reference_database.load_database_from_file(
                files=specificity_blastn_filter["files_fasta_reference_database"],
                file_type="fasta",
                database_overwrite=True,
            )
            specificity_refrence = BlastNFilter(
                search_parameters=specificity_blastn_filter["search_parameters"],
                hit_parameters=specificity_blastn_filter["hit_parameters"],
                filter_name="primer_blastn_specificity_reference",
                dir_output=self.dir_output,
            )
            specificity_refrence.set_reference_database(reference_database=reference_database)
            filters.append(specificity_refrence)
            directories.append(specificity_refrence.dir_output)

        if hybridization_probes_blastn_filter["enabled"]:
            hybridization_probes_database = ReferenceDatabase(
                database_name=self.subdir_db_reference, dir_output=self.dir_output
            )
            hybridization_probes_database.load_database_from_file(
                files=file_fasta_hybridization_probes_database,
                file_type="fasta",
                database_overwrite=True,
            )
            specificity_hybridization_probes = BlastNFilter(
                search_parameters=hybridization_probes_blastn_filter["search_parameters"],
                hit_parameters=hybridization_probes_blastn_filter["hit_parameters"],
                filter_name="primer_blastn_specificity_hybridization_probes",
                dir_output=self.dir_output,
            )
            specificity_hybridization_probes.set_reference_database(
                reference_database=hybridization_probes_database
            )
            filters.append(specificity_hybridization_probes)
            directories.append(specificity_hybridization_probes.dir_output)

        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        for directory in directories:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        check_content_oligo_database(oligo_database)

        return oligo_database

    def _get_best_forward_primer(
        self,
        oligo_database: OligoDatabase,
        reverse_primer_sequence: str,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
    ) -> str:
        """
        Choose the forward primer whose melting temperature best matches the reverse primer.

        After filtering, several forward-primer candidates may remain. This method
        compares their melting temperatures with the reverse primer and returns the
        closest match. Matching melting temperatures helps both primers anneal under
        the same PCR conditions.

        :param oligo_database: Filtered primer database returned by
            :py:meth:`_filter_by_specificity`.
        :type oligo_database: OligoDatabase
        :param reverse_primer_sequence: Reverse primer sequence used as the melting
            temperature reference.
        :type reverse_primer_sequence: str
        :param Tm_parameters: Settings used to calculate melting temperatures.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Optional chemical correction settings
            for the melting-temperature calculation.
        :type Tm_chem_correction_parameters: dict
        :param Tm_salt_correction_parameters: Optional salt correction settings for
            the melting-temperature calculation.
        :type Tm_salt_correction_parameters: dict
        :return: Forward primer sequence with the closest melting temperature match.
        :rtype: str
        """
        # Reverse-primer Tm is fixed; compute once and compare each forward candidate to it.
        Tm_reverse_primer = calc_tm_nn(
            sequence=reverse_primer_sequence,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        )

        min_dif_Tm = float("inf")
        forward_primer_sequence = ""
        for database_region in oligo_database.database.values():
            for primer_properties in database_region.values():
                Tm_forward_primer = calc_tm_nn(
                    sequence=primer_properties["oligo"],
                    Tm_parameters=Tm_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                )
                dif_Tm = abs(Tm_forward_primer - Tm_reverse_primer)
                if dif_Tm < min_dif_Tm:
                    min_dif_Tm = dif_Tm
                    forward_primer_sequence = primer_properties["oligo"]

        return forward_primer_sequence


############################################
# SeqFish Plus Probe Designer Pipeline
############################################


def _preprocess_config(config_validated: SeqfishPlusProbeDesignerConfig) -> dict[str, Any]:
    """
    Prepare the seqFISH+ config before the pipeline runs.

    This step converts the validated pydantic config to a plain dict and updates it
    so later design stages can read ready-to-use settings. Target-probe design in this
    pipeline does not use melting-temperature
    settings. When the forward primer is generated rather than loaded, this step also
    resolves melting-temperature tables, turns off unused temperature corrections,
    and copies the shared temperature settings into the forward-primer filter that
    needs them.

    It also expands an optional gene-list file into a concrete list of target
    regions. If no gene list is provided, all regions in the input FASTA files are
    used.

    Lastly, it inserts the parameters from required_parameters into the correct sections.

    :param config_validated: Validated pipeline configuration (pydantic model).
    :type config_validated: SeqfishPlusProbeDesignerConfig
    :return: The configuration converted to a dict, updated with the prepared settings.
    :rtype: dict
    """

    config = config_validated.model_dump()

    apply_required_parameters(config)

    # Readout probes and the forward primer only carry specificity filters when they
    # are generated rather than loaded.
    readout_probe_table_cfg = config["readout_probes"]["readout_probe_table"]
    if readout_probe_table_cfg["source"] == "generate":
        readout_probe_table_cfg["specificity_filters"]["specificity_blastn_filter"][
            "files_fasta_reference_database"
        ] = config["required_parameters"]["reference_genome"]

    file_region_ids = config["target_probes"]["oligo_generation"]["file_region_ids"]
    if file_region_ids is None:
        logger.warning(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. "
            "This choice can use a lot of resources."
        )
        config["target_probes"]["oligo_generation"]["region_ids"] = None
    else:
        with open(file_region_ids) as f:
            config["target_probes"]["oligo_generation"]["region_ids"] = sorted({line.rstrip() for line in f})

    forward_primer_cfg = config["primers"]["forward_primer"]
    if forward_primer_cfg["source"] == "generate":
        forward_primer_cfg["specificity_filters"]["specificity_blastn_filter"][
            "files_fasta_reference_database"
        ] = config["required_parameters"]["reference_genome"]

        # Resolve Tm table names and blank disabled chem/salt corrections to None so
        # downstream filters treat None as "no correction" without checking the flag.
        global_parameters = forward_primer_cfg["global_parameters"]
        global_parameters["Tm_parameters"] = preprocess_tm_parameters(global_parameters["Tm_parameters"])
        for correction in ["Tm_chem_correction_parameters", "Tm_salt_correction_parameters"]:
            correction_cfg = global_parameters[correction]
            if not correction_cfg["enabled"]:
                correction_cfg["parameters"] = None

        # Inline shared Tm settings into the forward-primer Tm filter.
        forward_primer_cfg["property_filters"]["Tm_filter"]["Tm_parameters"] = global_parameters[
            "Tm_parameters"
        ]
        forward_primer_cfg["property_filters"]["Tm_filter"]["Tm_chem_correction_parameters"] = (
            global_parameters["Tm_chem_correction_parameters"]["parameters"]
        )
        forward_primer_cfg["property_filters"]["Tm_filter"]["Tm_salt_correction_parameters"] = (
            global_parameters["Tm_salt_correction_parameters"]["parameters"]
        )

    return config


def seqfish_plus_probe_designer(config: SeqfishPlusProbeDesignerConfig) -> None:
    """
    Run the seqFISH+ probe design pipeline from a validated configuration (pydantic model).

    This function prepares the config with :func:`_preprocess_config`, then runs
    :class:`SeqFishPlusProbeDesigner` end to end. It designs target probes, loads or
    creates the readout probes and codebook, assembles hybridization probes and DNA
    templates with primers, and writes the final files under
    ``config['general']['dir_output']``. The caller should configure the library
    logger before calling this function (see :func:`main`).

    The config should follow ``data/configs/seqfish_plus_probe_designer.yaml``.

    Top-level config sections:

    - ``general``: output directory, intermediate-step writing, and worker count.
    - ``target_probes``: candidate generation, sequence filters, specificity filters,
      and probe set selection.
    - ``readout_probes``: codebook and readout probe table settings.
    - ``primers``: forward and reverse primer settings.

    Files written under ``dir_output``:

    - ``codebook.tsv``: barcode assignments for each target gene.
    - ``readout_probes.tsv``: readout probe sequences and related bit information.
    - ``seqfish_plus_probes.yml``: full probe records.
    - ``seqfish_plus_probes_order.yml``: sequences ready for synthesis.
    - ``seqfish_plus_probes.tsv`` / ``seqfish_plus_probes.xlsx``: probe sets as tables.

    Intermediate probe databases are also written when
    ``general.write_intermediate_steps`` is ``True``.

    See :class:`SeqFishPlusProbeDesigner` for the pipeline description and probe
    structure.

    :param config: Validated pipeline configuration. It is converted and prepared by
        :func:`_preprocess_config` before the pipeline runs.
    :type config: SeqfishPlusProbeDesignerConfig
    :return: None
    :rtype: None
    """

    config_dict = _preprocess_config(config)

    pipeline = SeqFishPlusProbeDesigner(
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        dir_output=config_dict["general"]["dir_output"],
        n_jobs=config_dict["general"]["n_jobs"],
    )

    target_probe_database = pipeline.design_target_probes(
        target_probes_parameters=config_dict["target_probes"],
    )

    codebook, readout_probe_table = pipeline.design_readout_probes(
        region_ids=list(target_probe_database.database.keys()),
        readout_probe_parameters=config_dict["readout_probes"],
    )

    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        oligo_database=target_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
        readout_probe_parameters=config_dict["readout_probes"],
    )

    reverse_primer_sequence, forward_primer_sequence = pipeline.design_primers(
        oligo_database=hybridization_probe_database,
        primer_parameters=config_dict["primers"],
    )

    dna_template_probe_database = pipeline.assemble_dna_template_probes(
        oligo_database=hybridization_probe_database,
        reverse_primer_sequence=reverse_primer_sequence,
        forward_primer_sequence=forward_primer_sequence,
    )

    pipeline.generate_output(
        oligo_database=dna_template_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
        readout_probe_parameters=config_dict["readout_probes"],
    )


def main() -> None:
    """
    Run the seqFISH+ probe design pipeline from the command line.

    Parses the required ``-c``/``--config`` argument, loads the YAML configuration
    file, and configures the library logger to write under the configured output
    directory. It then calls :func:`seqfish_plus_probe_designer`.

    :return: None
    :rtype: None
    """
    print("--------------START PIPELINE--------------")

    args = base_parser(
        prog="seqFISH+ Probe Designer",
        usage="seqfish_plus_probe_designer [options]",
        description=__doc__,
    )

    config_validated = load_config(args["config"], SeqfishPlusProbeDesignerConfig)

    # Configure logging only after dir_output is known so the log file lands there.
    configure_root_logger(
        dir_output=config_validated.general.dir_output,
        pipeline_name="seqfishplus_probe_designer",
    )

    seqfish_plus_probe_designer(config_validated)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
