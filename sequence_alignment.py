from Bio import SeqIO
import sys

orig_stdout = sys.stdout
f = open('bio_results.txt', 'w')
sys.stdout = f

def print_matrix(matrix):
    matrix_size=len(matrix)
    print(matrix_size)
    for i in range(0, matrix_size):
        sequence_name='seq_{}'.format(i)
        print(sequence_name,end=' ')
        for j in range(0, matrix_size):
            print('%.4f' % matrix[i][j], end=' ')
        print()
    print()

# check if the given (di)codon starts with 'start' codon
# return bool


def is_start_codon(codon):
    if(codon[0] == 'A' and codon[1] == 'T' and codon[2] == 'G'):
        return True
    return False

# check if the given (di)codon ends with 'stop' codon
# return bool


def is_stop_codon(codon):
    c=len(codon)-1
    b=c-1
    a=c-2
    if(codon[a] != 'T'):
        return False
    if(codon[b] == 'A'):
        if(codon[c] == 'A' or codon[c] == 'G'):
            return True
    if(codon[a] == 'G'):
        if(codon[b] == 'A'):
            return True
    return False

# find protein sequences in plasmide sequence
# return all found open reading frames


def find_protein_sequences(sequence, codon_length):
    orfs = []
    for index in range(0, 3):
        codons = []
        codon = []
        started = False
        counter = 0
        for i in sequence[index:]:
            codon.append(i)
            counter += 1
            if counter < codon_length:
                continue
            counter = 0
            if is_start_codon(codon):
                if started == False:
                    started = True
                codons.append(codon)
                codon=[]
                continue
            if started == False:
                pass
            if is_stop_codon(codon):
                codons.append(codon)
                orfs.append(codons)
                codons = []
                started = False
            else:
                codons.append(codon)
            codon = []
        # longest_codon = longest_protein(orf)
    return orfs


def print_codons(codons):
    for codon in codons:
        for i in codon:
            print(i)
        print('\n')

# drop sequences which are shorter than 100 base pairs.
# return list of sequences


def filter(sequence, limit):
    result = []
    for codon in sequence:
        if len(codon) > limit:
            result.append(codon)
    return result

# calculate frequencies of each codon in the given orf
# return dictionary: {codon: frequence}


def calculate_frequencies(orf):
    result = {}
    length = len(orf)

    for codon in orf:
        counter = 0
        for comparable_codon in orf:
            if comparable_codon == codon:
                counter += 1
        freq = counter/length
        result[str(codon)] = freq

    return result

# distance function


def function(arg1, arg2):
    return abs(arg1-arg2)**1.5  

# find the distance between two sequences calculated for each pair
# in two sequences having exact key and for distinct keys.
# return distance


def compare(seq1, seq2):
    found = []
    distance = 0
    for k1, v1 in seq1.items():
        is_found = False
        for k2, v2 in seq2.items():
            if k1 == k2:
                distance += function(v1, v2)
                is_found = True
        if is_found == False:
            distance += function(v1, 0)
        found.append(k1)

    for k2, v2 in seq2.items():
        if k2 not in found:
            distance += function(v2, 0)
            found.append(k2)

    return distance

# create a distance matrix from given sequences frequencies
# return an (n x n) matrix where x is the number of sequences.


def create_distance_matrix(frequencies):
    length = len(frequencies)
    matrix = [[0 for x in range(length)] for y in range(length)]  # init zeros

    for i in range(0, length):
        for j in range(0, length):
            distance = compare(frequencies[i], frequencies[j])
            matrix[i][j] = round(distance, 3)
    return matrix


for seq_record in SeqIO.parse('plazmide.fasta', 'fasta'):
    # dikodonų atvejis
    # codon_length=6
    # count_of_codons_lower_limit=16

    # kodonų atvejis
    codon_length=3
    count_of_codons_lower_limit=33
    
    codons = find_protein_sequences(seq_record,codon_length)
    longer = filter(codons,count_of_codons_lower_limit)

    index = 0
    frequencies = {}
    for i in longer:
        freq = calculate_frequencies(i)
        frequencies[index] = freq
        index += 1

    python = create_distance_matrix(frequencies)
    print_matrix(python)
    print()
    print('---------------------------------------')

sys.stdout = orig_stdout
f.close()
