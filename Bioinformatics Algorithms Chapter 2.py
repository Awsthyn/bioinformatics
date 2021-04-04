def distance_between_pattern_and_strings(pattern, dna):
    median_sum = 0
    k = len(pattern)
    for string in dna:
        distance_sum = float('inf')
        for kmer in range(len(string) - k):
            seq = string[kmer: kmer+k]
            d = hamming_distance(pattern,seq)
            if(distance_sum > d):
                
                distance_sum = d
        median_sum += distance_sum
    return median_sum

#BA2B 
def median_string(dna,k):
    distance = float('inf')
    median = ''
    nucleotides = ["A", "C", "G", "T"]
    for i in nucleotides:
        for j in nucleotides:
            for m in nucleotides:
                for n in nucleotides:
                    for o in nucleotides:
                        for p in nucleotides:
                            pattern = i+j+m+n+o+p
                            d = distance_between_pattern_and_strings(pattern, dna)
                            if (distance > d):
                                distance = d
                                median = pattern
    print( median)

from Bio import SeqIO

def parse_fasta(file):
    seqs = []
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs.append(str(record.seq))
    return seqs

def generate_nucleotide_matrix(seqs):
    nucleotides = {"A": [], "C": [], "G": [], "T": []}

    for j in range(len(seqs[0])):
        for n in nucleotides:
            nucleotides[n].append(0)
   
    for seq in seqs:
        for i in range(len(seq)):
            letter = seq[i]
            nucleotides[letter][i] += 1
    return nucleotides

def consensus(matrix):
    string = ""
    nucleotides = ["A","C","G", "T"]
    for p in range(len(matrix["A"])):
        a = matrix["A"][p]
        c = matrix["C"][p]
        g = matrix["G"][p]
        t = matrix["T"][p]
        count = [a,c,g,t]
        maximum = max(count)
        for n in range(len(count)):
            if(count[n] == maximum):
                string += nucleotides[n]
                break
    print(string)
    for i in list(matrix.keys()):
        print(i+": " + " ".join(str(item) for item in matrix[i]))

def find_most_probable_string(matrix):
    string = ""
    nucleotides = ["A","C","G", "T"]
    for p in range(len(matrix["A"])):
        a = matrix["A"][p]
        c = matrix["C"][p]
        g = matrix["G"][p]
        t = matrix["T"][p]
        count = [a,c,g,t]
        maximum = max(count)
        for n in range(len(count)):
            if(count[n] == maximum):
                string += nucleotides[n]
                break
    return string

def parse_matrix(matrix):
    parsed_matrix = []
    nucleotides = {"A": [], "C": [], "G": [], "T": []}

    for i in matrix:
        parsed_matrix.append(list(map(float ,i.split(' '))))
    nucleotides["A"] = parsed_matrix[0]
    nucleotides["C"] = parsed_matrix[1]
    nucleotides["G"] = parsed_matrix[2]
    nucleotides["T"] = parsed_matrix[3]
    return nucleotides

def find_nearest(text, pattern, k):
    distance = float('inf')
    string = ''
    for i in range(len(text) - len(pattern)):
        segment = text[i: i+k] 
        d = hamming_distance(segment, pattern)
        if (distance > d):
            distance = d
            string = segment
    return string

#BA2C
def profile_most_probable_k_mer(text, k , matrix):
    parsed_matrix = parse_matrix(matrix)
    most_probable_string = find_most_probable_string(parsed_matrix)
    nearest = find_nearest(text, most_probable_string, k)
    print(nearest)