def read_fasta_multi(filename): # Reads a FASTA with multiple sequences. Returns a dict: {header: sequence}.
    sequences = {}
    header = None
    seq_parts = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip() # removes empty spaces
            if line.startswith(">"):
                if header:  # vorige opslaan
                    sequences[header] = "".join(seq_parts)
                header = line[1:]  # safes header without '>'
                seq_parts = []
            else:
                seq_parts.append(line)

        # saves last sequence
        if header:
            sequences[header] = "".join(seq_parts)

    return sequences


def split_in_codons(dna_sequence):
    codons = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            codons.append(codon)
    return codons

# Codon dictionary
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'TTA':'L', 'TTG':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TGC':'C', 'TGT':'C', 'TGG':'W',
    'TAC':'Y', 'TAT':'Y',
    'TAA':'_', 'TAG':'_', 'TGA':'_'
}


def amino_seq(codons, codon_table):
    sequence = []
    for codon in codons:
        aa = codon_table.get(codon, 'X')  # X for unknown codons (not present in codon_table)
        if aa == '_':  # stopcodon
            break
        sequence.append(aa)
    return ''.join(sequence) # sequence changes from list to string


# Use
fasta_file = r"C:\file"
rna_dict = read_fasta_multi(fasta_file)

for header, seq in rna_dict.items():
    codons = split_in_codons(seq)
    protein = amino_seq(codons, codon_table)
    print(">", header)
    print(protein[:100], "...")   # first 100 amino acids 
    print("Length protein:", len(protein))
    print()
