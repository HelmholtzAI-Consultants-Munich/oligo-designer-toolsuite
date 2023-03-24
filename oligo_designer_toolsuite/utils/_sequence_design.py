import numpy as np 
import random
import itertools

from Bio.Seq import Seq


def generate_random_sequence(oligo_length, sequence_alphabet=['A', 'T', 'G', 'C']):
    seq_letters = np.random.choice(sequence_alphabet, oligo_length, replace = True)
    seq = ''.join(seq_letters)
    return Seq(seq)


def generate_binary_sequences(n_seq, n_bit=16, n_one=1):
    ''' 
    This generator randomly generates the codebook from user input, and returns a list of encoding sequences used for merfish.
    
    Params:
    n_bit: Number of bits in an encoding code
    n_one: Number of ones in an encoding code
    n_seq: Number of encoding codes in a list 
    '''

    if n_seq > 2 ** n_bit:
        raise Exception("Error: the number of sequences is tool large!")
    if n_one > n_bit:
        raise Exception("Error: the number of ones is tool large!")

    else:
        sequences = set()
        while len(sequences) < n_seq:
            candidate = format(random.getrandbits(n_bit), f"0{n_bit}b")
            if bin(int(candidate, 2)).count("1") == n_one and candidate not in sequences:
                sequences.add(candidate)
        
        return list(sequences)




def generate_codebook(num_seq: int, encoding_scheme: str):
    '''
    This function generates bit-codes

    return:
    A list of num_seq numbers of MHD2/MHD4 encoding codes

    Params:
    num_seq: Numer of Codes
    encoding_scheme: MHD2 or MHD4
    '''


    # initialize an empty list to store the codes
    codes = []
    pass_sequences = []
    if encoding_scheme == "MHD4":
        # loop through all possible combinations of 16 bits
        for i in range(2**16):
            # convert the integer i to a binary string of length 16
            bin_str = format(i, '016b')
            # count the number of ones in the binary string
            num_ones = bin_str.count('1')
            # if the number of ones is exactly 4
            if num_ones == 4:
                # append the binary string to the codes list
                codes.append(bin_str)
        # loop through all pairs of codes in the list

        for i in range(len(codes)):
            for j in range(i+1, len(codes)):

                if len(pass_sequences) == 140 or len(pass_sequences) == num_seq:
                    break

                # compute the Hamming distance between two codes
                ham_dist = sum(a != b for a, b in zip(codes[i], codes[j]))
                # if the Hamming distance is less than 4
                if ham_dist >= 4:
                    if len(pass_sequences) == 0:
                        pass_sequences.append(codes[i])
                    # remove one of the codes from the list
                    else:
                        ham_list = []
                        for pass_seq in pass_sequences:
                            ham_dist = sum(a != b for a, b in zip(pass_seq, codes[j]))
                            ham_list.append(ham_dist)

                        if all(4 <= ham for ham in ham_list):
                            pass_sequences.append(codes[j])
                            print(len(pass_sequences))
                            break
            if len(pass_sequences) == 140 or len(pass_sequences) == num_seq:
                break
                    # break out of the inner loop

    if encoding_scheme == "MHD2":
                # loop through all possible combinations of 14 bits
        for i in range(2**14):
            # convert the integer i to a binary string of length 16
            bin_str = format(i, '016b')
            # count the number of ones in the binary string
            num_ones = bin_str.count('1')
            # if the number of ones is exactly 4
            if num_ones == 4:
                # append the binary string to the codes list
                codes.append(bin_str)
        # loop through all pairs of codes in the list

        for i in range(len(codes)):
            for j in range(i+1, len(codes)):

                if len(pass_sequences) == 1001 or len(pass_sequences) == num_seq:
                    break
                # compute the Hamming distance between two codes
                ham_dist = sum(a != b for a, b in zip(codes[i], codes[j]))
                # if the Hamming distance is less than 2
                if ham_dist >= 2:
                    if len(pass_sequences) == 0:
                        pass_sequences.append(codes[i])
                    # remove one of the codes from the list
                    else:
                        ham_list = []
                        for pass_seq in pass_sequences:
                            ham_dist = sum(a != b for a, b in zip(pass_seq, codes[j]))
                            ham_list.append(ham_dist)
                        # check passed sequences
                        if all(2 <= ham for ham in ham_list):
                            pass_sequences.append(codes[j])
                            print(len(pass_sequences))
                            break
                    # break out of the inner loop
            if len(pass_sequences) == 1001 or len(pass_sequences) == num_seq:
                break
    return pass_sequences

def SCRINSHOT_or_ISS_backbone_sequence(
            region_idx, barcode_length=4, barcode_seed=0
        ):
            """Get backbone sequence of padlock oligos for SCRINSHOT or ISS

            Arguments
            ---------
            region_idx: int
                Identifier for a given region. The identifier makes sure to return the same bar code
                for the different padlock oligos of a given region.
            barcode_length: int
                Length of barcode sequence
            barcode_seed: int
                Defines the random assignment of barcodes to each region_idx.

            Returns
            -------
            str:
                backbone sequence (5' to 3')
            dict of strs:
                Individual parts of the backbone sequence

            """
            accessory1 = "TCCTCTATGATTACTGAC"
            ISS_anchor = "TGCGTCTATTTAGTGGAGCC"
            barcode = get_barcode(region_idx, length=barcode_length, seed=barcode_seed)
            accessory2 = "CTATCTTCTTT"

            sub_seqs = {
                "accessory1": accessory1,
                "ISS_anchor": ISS_anchor,
                "barcode": barcode,
                "accessory2": accessory2,
            }
            full_seq = accessory1 + ISS_anchor + barcode + accessory2

            return full_seq, sub_seqs

def get_barcode(region_idx, length=4, seed=0, choices=["A", "C", "T", "G"]):
            """Get barcode sub sequence of padlock oligo for in situ sequencing

            For SCRINSHOT padlock oligos this could be constant, however it makes sense to have
            different barcodes so that the oligo set could also be used for ISS experiments.

            Arguments
            ---------
            region_idx: int
                Identifier for a given region. The identifier makes sure to return the same bar code
                for the different padlock oligos of a given region.
            length: int
                Length of barcode sequence
            seed: int
                Defines the random assignment of barcodes to each region_idx.

            Returns
            -------
            str: barcode sequence (5' to 3')

            """

            barcodes = ["".join(nts) for nts in itertools.product(choices, repeat=length)]
            random.seed(seed)
            random.shuffle(barcodes)

            if region_idx >= len(barcodes):
                raise ValueError(
                    "Barcode index exceeds number of possible combinations of barcodes. Increase barcode length?"
                )

            return barcodes[region_idx]

def convert_complementary_seq_to_arms(complementary_seq, ligation_idx):
            """Convert the complementary sequence of padlock oligos to two arms with 5' to 3' convention

            E.g.
            complementary_seq = "AAAATGCTTAAGC" ligation_idx = 7
            --> cut after 7th base: "AAAATGC|TTAAGC"
            --> final result: ["CGTAAAA","CGAATT"]

            Arguments
            ---------
            complementary_seq: str
                Sequence that hybridises with the target RNA
            ligation_idx: int
                Site where complementary_seq is cut in two arms according padlock oligo design

            Returns
            -------
            list of strs: first and second arm sequences (both 5' to 3')

            """
            arm1 = complementary_seq[ligation_idx:]
            arm2 = complementary_seq[:ligation_idx]

            return [arm1, arm2]

