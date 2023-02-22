import unittest
import random
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


def generate_codebook(num_seq, num_bits, num_ones, max_distance):
    # Initialize the list of sequences
    sequences = []

    # Generate sequences until the desired number of sequences is reached
    while len(sequences) < num_seq:
        # Generate a new sequence with the desired number of 1s
        sequence = [1] * num_ones + [0] * (num_bits - num_ones)
        random.shuffle(sequence)

        # check duplicates
        if sequence in sequences:
            continue

        # Check if the new sequence meets the minimum distance requirement
        if all(max_distance >= hamming(sequence, seq) * num_bits for seq in sequences):
            # Add the new sequence to the list of sequences
            sequences.append(sequence)

    # Convert the list of sequences to a list of strings and return it
    return ["".join(str(bit) for bit in sequence) for sequence in sequences]


if __name__ == "__main__":

    unittest.main()