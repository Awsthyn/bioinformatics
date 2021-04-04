def counting_nucleotides(data):
  count_a = 0
  count_t = 0
  count_c = 0
  count_g = 0
  for letter in data:
    if (letter == "A"): count_a += 1
    if (letter == "T"): count_t += 1
    if (letter == "C"): count_c += 1
    if (letter == "G"): count_g += 1
  return str(count_a) + " " + str(count_c) + " " + str(count_g) + " " + str(count_t)


def transcription(data):
  return data.replace("T", "U")

def compute_gc_content(strings):
  split= strings.split(">")
  sample = ""
  max = 0
  for str in split:
    if(len(str) > 0):
      count_nucleotides = str.count("G") + str.count("C") + str.count("A") + str.count("T")
      count = ((str.count("G") + str.count("C")) / count_nucleotides * 100)
      if(count > max): 
        max = count
        sample = str[0:13]
  print(sample)
  print(max)

def genbank_search(term, start_date, end_date):
    from Bio import Entrez
    Entrez.email = "aguswagner008@gmail.com"
    handle = Entrez.esearch(db="nucleotide", term='"'+ term + '"[Organism] AND ("' +start_date +'"[Publication Date] : "' + end_date + '"[Publication Date])')
    record = Entrez.read(handle)
    print(record["Count"])


def convert_fastq_to_fasta(file):
    from Bio import SeqIO
    seqs = []
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seqs.append(record)
    return  SeqIO.write(seqs, "fastq_to_fasta.txt", "fasta")

def protein_translation(protein, coding_dna):
    from Bio.Seq import translate
    for i in range(1,16):
        if(i > 6 and i < 9): continue
        if(translate(coding_dna, table=i, stop_symbol="") == protein ): print(i)