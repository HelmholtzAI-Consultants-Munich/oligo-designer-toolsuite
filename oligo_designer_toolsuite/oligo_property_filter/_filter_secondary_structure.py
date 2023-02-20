from Bio.SeqUtils import Seq
from seqfold import dg

from . import PropertyFilterBase


class Secondary_struct(PropertyFilterBase):
    """ Filter the sequences by presence of a hairpin secondary structure
    :param T: Temperature at which the secondary structures are predicted (°C)
    :type T: float
    :param DG: Delta G threshold for hairpins in kcal/mol (Haipin DG should be weaker (more positive) than this threshold)
    :type DG: float
    
    """
    
    def __init__(self, T: float, DG: float) -> None:
        """Constructor"""
        super().__init__()
        self.T = T
        self.DG = DG

    def apply(self, sequence: Seq):
        """Applies the filter and returns True if there isn't any hairpin secondary structure.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrined is fulfilled
        :rtype: bool
        """
        
        if (dg(Seq, temp = self.T)< self.DG):
                return False,{}
        return True,{}