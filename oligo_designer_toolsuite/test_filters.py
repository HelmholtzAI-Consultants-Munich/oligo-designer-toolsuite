import os
import time

from Bio.SeqUtils import MeltingTemp as mt

os.chdir("../")
cwd = os.getcwd()
print(cwd)

from IO._database import CustomDB
from oligo_pre_filter._filter_base import GCContent, MaskedSequences, MeltingTemperature
from oligo_pre_filter._filter_padlock_probes import PadlockArms

dir_output = os.path.join(cwd, "output/annotations")

Tm_parameters = {
    "check": True,
    "strict": True,
    "c_seq": None,
    "shift": 0,
    "nn_table": getattr(mt, "DNA_NN3"),
    "tmm_table": getattr(mt, "DNA_TMM1"),
    "imm_table": getattr(mt, "DNA_IMM1"),
    "de_table": getattr(mt, "DNA_DE1"),
    "dnac1": 50,  # [nM]
    "dnac2": 0,
    "selfcomp": False,
    "dNTPs": 0,
    "saltcorr": 7,
    "Na": 1.25,  # [mM]
    "K": 75,  # [mM]
    "Tris": 20,  # [mM]
    "Mg": 10,  # [mM]
}
Tm_correction_parameters = {
    "DMSO": 0,
    "DMSOfactor": 0.75,
    "fmdfactor": 0.65,
    "fmdmethod": 1,
    "GC": None,
    "fmd": 20,
}
annotation = dir_output + "/GCF_000001405.40_GRCh38.p14_genomic.gtf"
sequence = dir_output + "/GCF_000001405.40_GRCh38.p14_genomic.fna"
genes = [
    "WASH7P",
    "DDX11L1",
    "TRNT",
    "NOC2L",
    "PLEKHN1",
    "AGRN",
    "UBE2J2",
    "DVL1",
    "MIB2",
    "LOC112268402_1",
]

masked_sequences = MaskedSequences()
GC_content = GCContent(GC_content_min=40, GC_content_max=60)
melting_temperature = MeltingTemperature(
    Tm_min=52,
    Tm_max=67,
    Tm_parameters=Tm_parameters,
    Tm_correction_parameters=Tm_correction_parameters,
)
arms_tm = PadlockArms(
    min_arm_length=10,
    max_Tm_dif=2,
    Tm_min=40,
    Tm_max=43,
    Tm_parameters=Tm_parameters,
    Tm_correction_parameters=Tm_correction_parameters,
)

filters = [masked_sequences, GC_content, melting_temperature, arms_tm]
start_time = time.time()
custom = CustomDB(
    probe_length_min=30,
    probe_length_max=40,
    file_annotation=annotation,
    file_sequence=sequence,
    filters=filters,
)
custom.create_reference_DB(dir_output=dir_output)
custom.create_oligos_DB(genes=genes, number_batchs=2, dir_output=dir_output, write=True)
t = (time.time() - start_time) / 60
print(" Generation the DB: --- %s min ---" % (t))
