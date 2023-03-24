import numpy as np
import random
import itertools

def create_barcodes(pseudocolors=4, seed=33, num_genes=10, genes=['GENE_' + str(i) for i in range(10)]):
        """
        Function to create barcodes for each gene, that is specified in genes
        :return: dictionary of barcodes for each gene; { gene : barcode }
        :rtype: dict (str : list of 4 int)
        """ 
        output = dict()
        all_possible_barcodes = [[i, (i+j+k)%pseudocolors,j, k] for i in range(pseudocolors) for j in range(pseudocolors) for k in range(pseudocolors)]
        l = len(all_possible_barcodes)
        arr = np.arange(0, l, 1)
        np.random.seed(seed)
        barcodes = np.random.choice(arr, size = num_genes ,replace = False)
        for i in range(num_genes):
            output[genes[i]] = all_possible_barcodes[barcodes[i]]
        return output

def get_barcode(region_idx, length=2, seed=0, choices=["A", "C", "T", "G"]):

            barcodes = ["".join(nts) for nts in itertools.product(choices, repeat=length)]
            random.seed(seed)
            random.shuffle(barcodes)

            if region_idx >= len(barcodes):
                raise ValueError(
                    "Barcode index exceeds number of possible combinations of barcodes. Increase barcode length?"
                )

            return barcodes[region_idx]

print('v_students', create_barcodes(seed=0))
print('v_louiss', [get_barcode(i, choices=[str(i) for i in range(4)]) for i in range(10)])