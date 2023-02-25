import numpy as np 


class BarcodingCreator:

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
        output = dict()
        all_possible_barcodes = [[i, (i+j+k)%self.pseudocolors,j, k] for i in range(self.pseudocolors) for j in range(self.pseudocolors) for k in range(self.pseudocolors)]
        l = len(all_possible_barcodes)
        arr = np.arange(0, l, 1)
        np.random.seed(self.seed)
        barcodes = np.random.choice(arr, size = self.num_genes ,replace = False)
        for i in range(self.num_genes):
            output[self.genes[i]] = all_possible_barcodes[barcodes[i]]
        return output

    