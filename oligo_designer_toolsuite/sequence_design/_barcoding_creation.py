import numpy as np 


class BarcodingCreator:
    """
    This class is used to create 4-digit barcode for each gene.
    The barcodes can manage the loss of 1 digit, as they are designed in the following way: (i, j, k, mod(i+j+k, num_pseudocolors))
    :param num_pseudocolors: number of pseudocolors, that are used in the experiment
    :type num_pseudocolors: int
    :param list_of_genes: list of names of genes, for which baracodes should be created
    :type list_of_genes: list
    """
    def __init__(
        self,
        num_pseudocolors: int,
        list_of_genes:list, 
        seed: int = 0
    ):
        self.pseudocolors = num_pseudocolors
        self.genes = list_of_genes
        self.num_genes = len(list_of_genes)
        self.seed = seed


    def create_barcodes(self):
        """
        Function to create barcodes for each gene, that is specified in self.genes
        :return: dictionary of barcodes for each gene; { gene : barcode }
        :rtype: dict (str : list of 4 int)
        """ 
        output = dict()
        all_possible_barcodes = [[i, (i+j+k)%self.pseudocolors,j, k] for i in range(self.pseudocolors) for j in range(self.pseudocolors) for k in range(self.pseudocolors)]
        l = len(all_possible_barcodes)
        arr = np.arange(0, l, 1)
        np.random.seed(self.seed)
        barcodes = np.random.choice(arr, size = self.num_genes ,replace = False)
        for i in range(self.num_genes):
            output[self.genes[i]] = all_possible_barcodes[barcodes[i]]
        return output

    