from Bio.SeqUtils import Seq
from seqfold import dg

from . import PropertyFilterBase


class Secondary_struct(PropertyFilterBase):
    """ Filter sequences by the minimum free energy of the folded sequence, i.e. secondary structure containing stacks, bulges, hairpins or interior loops.
    :param T: The temperature to fold at, i.e. temperature at which the sequence folding is predicted (°C)
    :type T: float
    :param DG: Delta G (minimum free energy) threshold for the folded sequence in kcal/mol (DG of folded sequence should be weaker (more positive) than this threshold)
    :type DG: float
    
    """
    
    def __init__(self, T: float, DG: float) -> None:
        """Constructor"""
        super().__init__()
        self.T = T
        self.DG = DG

    def apply(self, sequence: Seq):
        """Applies the filter and returns True if the minimum free energy of the folded sequence is below the given threshold.

        :param sequence: sequence to be filtered
        :type sequence: str
        :return: True if the constrain is fulfilled
        :rtype: bool
        """
        if (dg(sequence, temp = self.T)< self.DG):
                return False,{}
        return True,{}
