import os
import sys

cwd=os.getcwd()
sys.path.append(cwd)

from oligo_designer_toolsuite.database import (
    CustomGenomicRegionGenerator,
    ReferenceDatabase,
    OligoDatabase,
)

annotation_file_ncbi = "tests/data/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
sequence_file_ncbi = "tests/data/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
genes = ["AARS1", "DECR2", "FAM234A", "RHBDF1", "WASIR2"]

region_generator_ncbi = CustomGenomicRegionGenerator(
    annotation_file_ncbi,
    sequence_file_ncbi,
    files_source="NCBI",
    species="Homo_sapiens",
    annotation_release="110",
    genome_assembly="GRCh38",
    dir_output='output',
)
file_ncbi_transcriptome = (
    region_generator_ncbi.generate_transcript_reduced_representation()
)

# file_ncbi_transcriptome = '/home/francesco/Desktop/Work/oligo-designer-toolsuite/tests/notebooks/output/annotation/transcriptome_source_NCBI_species_Homo_sapiens_annotation_release_110_genome_assemly_GRCh38_incl_exonjunctions_of_size_100.fna'
oligos = OligoDatabase(
    file_fasta=file_ncbi_transcriptome,
    min_oligos_per_region=0,
    files_source="NCBI",
    species="Homo_sapiens",
    annotation_release="110",
    genome_assembly="GRCh38",
    n_jobs=2,
)

oligos.create_database(oligo_length_min=90, oligo_length_max=90,region_ids=genes)
prob = 0
for key in oligos.database.keys():
    prob += len(list(oligos.database[key].values()))
    print(len(list(oligos.database[key].values())))
print(prob)