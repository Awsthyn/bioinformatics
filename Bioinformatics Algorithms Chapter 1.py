#BA1A
def pattern_count(text, pattern):
    count = 0
    for i in range(len(text) - len(pattern)):
        if(pattern == text[i: i+len(pattern)]):
            count += 1   
    return count

#BA1B
def frequent_words(text, k):
    hash_table = {}
    for i in range(len(text) - k):
        try:
            hash_table[text[i: i+k]] += 1
        except KeyError:
            hash_table[text[i: i+k]] = 1
    maximum = max(hash_table.values())
    result = dict(filter(lambda elem: elem[1] == maximum, hash_table.items()))
    return result.keys()

#BA1C
def complement(data):
  reverse = data[::-1]
  reverse = reverse.replace("A", "t")
  reverse = reverse.replace("T", "a")
  reverse = reverse.replace("C", "g")
  reverse = reverse.replace("G", "c")
  return reverse.upper()

#BA1D
def pattern_matching(pattern, genome):
    list = []
    for i in range(len(genome) - len(pattern)):
        if(pattern == genome[i: i+len(pattern)]):
            list.append(i)   
    return " ".join(str(item) for item in list)

#BA1E
def frequent_words_clump_finding(text, k, t):
    hash_table = {}
    for i in range(len(text) - k):
        try:
            hash_table[text[i: i+k]] += 1
        except KeyError:
            hash_table[text[i: i+k]] = 1
    
    result = dict(filter(lambda elem: elem[1] == t-1, hash_table.items()))
    return list(result.keys())

def clump_finding(genome, k, L, t):
    clumps = set(())
    for i in  range(len(genome) - L):
        aux = frequent_words_clump_finding(genome[i: i+L], k, t)
        for i in aux:
            clumps.add(i)
    print(' '.join(clumps))

#BA1F
def minimum_skew(genome):
  count_c = 0
  count_g = 0
  counter = []
  oris = []
  for i in genome:
      counter.append(count_g - count_c)
      if(i == "C"): count_c += 1
      if(i == "G"): count_g += 1
  minimum = min(counter)
  for index, count in enumerate(counter):
      if(count == minimum): oris.append(index)
  return " ".join(str(item) for item in oris)

#BA1G
def hamming_distance(string1, string2):
    difference = 0
    for i, value in enumerate(string1):
        if(string1[i] != string2[i]): difference += 1
    return difference

#BA1H
def approximate_pattern_matching(pattern, text, d):
    positions = []
    for i in range(len(text) - len(pattern)):
        compare = text[i: i+len(pattern)]
        distance = hamming_distance(compare, pattern)
        if(distance <= d): positions.append(i)
    print(" ".join(str(item) for item in positions))
    
#BA1I
def frequent_words_with_mismatches(text, k, d):
    dictionary = {}
    for i in range(len(text) - k):  
        s = text[i: i+k]
        for index, l in enumerate(s):
            nucleotides = ["A", "T", "C", "G"]
            for nucleotid in nucleotides:
                mutation = s[:index] + nucleotid + s[index + 1:]
                matches = approximate_pattern_matching(mutation, text, d)
                dictionary[mutation] = matches
    maximum = max(dictionary.values())
    result = dict(filter(lambda elem: elem[1] == maximum, dictionary.items()))
    print(" ".join(list(result.keys())))
    
#BA1K
def frequency_array(text, y):
    hash_table = {}
    nucleotides = ["A", "C", "G", "T"]
    frequency_list = []
    for i in nucleotides:
        for j in nucleotides:
            for k in nucleotides:
                for l in nucleotides:
                    for m in nucleotides:
                        for n in nucleotides:
                            hash_table[i+j+k+l+m+n] = 0
    for x in range(len(text)-y+1):
        try:
            hash_table[text[x: x+y]] += 1
        except KeyError:
            None
    for key in sorted(hash_table):
       frequency_list.append(hash_table[key])
    print(" ".join(str(item) for item in frequency_list))

#BA1L
def pattern_to_number(pattern):
    if(pattern == ""):
        return 0
    symbol = pattern[len(pattern)-1:]
    prefix = pattern[:len(pattern)-1]
    symbol_to_number = {"A": 0, "C": 1, "G": 2, "T":3}
    return 4 * pattern_to_number(prefix) + symbol_to_number[symbol]


from math import floor
#BA1M
def number_to_pattern(index,k):
    number_to_symbol = ["A", "C", "G", "T"]
    if (k == 1):
        return number_to_symbol[index]
    prefix_index = floor(index / 4)
    r = index % 4
    symbol = number_to_symbol[r]
    return  number_to_pattern(prefix_index, k-1) + symbol
    
#BA1N
def neighbors(pattern, d):
    if(d == 0):
        return pattern;
    if(len(pattern) == 1):
        return {"A", "C", "G", "T"}
    neighborhood = set()
    suffix = pattern[1: len(pattern)]
    first_symbol = pattern[:1]
    suffix_neighbors = neighbors(suffix, d)
    for string in suffix_neighbors:
        if(hamming_distance(suffix, string) < d):
            neighborhood.add("A" + string)
            neighborhood.add("C" + string)
            neighborhood.add("G" + string)
            neighborhood.add("T" + string)
        else:
            neighborhood.add(first_symbol + string)
    return neighborhood