import unittest
from numpy import random
from scipy.spatial.distance import hamming


def generate_binary_sequences(n_bit, n_one, n_seq=None):
    ''' 
    This generator randomly generates the codebook from user input, and returns a list of encoding sequences used for merfish.
    
    Params:
    n_bit: Number of bits in an encoding code
    n_one: Number of ones in an encoding code
    n_seq: Number of encoding codes in a list 
    '''

    if n_seq > 2 ** n_bit:
        raise Exception("Sorry, the number of sequences were larger than possible")
    if n_one > n_bit:
        raise Exception("Sorry, the number of ones were were larger than possible")

    else:
        sequences = set()
        count = 0
        while count < n_seq:
            candidate = format(random.getrandbits(n_bit), f"0{n_bit}b")
            if bin(int(candidate, 2)).count("1") == n_one and candidate not in sequences:
                sequences.add(candidate)
                yield candidate
                count += 1


def get_binary_sequences(n_bit=16, n_one=4, n_seq=None):
    ''' 
    This function is using generator generate_binary_sequences() to return a list of encoding codes

    Params:
    n_bit: Number of bits in a code
    n_one: Number of ones in a code
    n_seq: Number of codes in a list 
    '''

    return list(generate_binary_sequences(n_bit, n_one, n_seq))


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

if __name__ == "__main__":
    unittest.main()
