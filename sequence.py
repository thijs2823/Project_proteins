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
    # Alanine
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    # Arginine
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    # Asparagine
    'AAT':'N', 'AAC':'N',
    # Aspartic acid
    'GAT':'D', 'GAC':'D',
    # Cysteine
    'TGT':'C', 'TGC':'C',
    # Glutamine
    'CAA':'Q', 'CAG':'Q',
    # Glutamic acid
    'GAA':'E', 'GAG':'E',
    # Glycine
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    # Histidine
    'CAT':'H', 'CAC':'H',
    # Isoleucine
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    # Leucine
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    # Lysine
    'AAA':'K', 'AAG':'K',
    # Methionine (start)
    'ATG':'M',
    # Phenylalanine
    'TTT':'F', 'TTC':'F',
    # Proline
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    # Serine
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    # Threonine
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    # Tryptophan
    'TGG':'W',
    # Tyrosine
    'TAT':'Y', 'TAC':'Y',
    # Valine
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    # Stop codons
    'TAA':'_', 'TAG':'_', 'TGA':'_'
}


def translate_from_first_ATG(dna_sequence, codon_table):
    """
    Translates a DNA sequence into a protein sequence starting from the first ATG codon (start codon).
    This ensures translation matches the actual protein (ignores UTRs before CDS).
    """
    dna_sequence = dna_sequence.upper().replace("\n","").replace(" ","")
    
    for i in range(0, len(dna_sequence)-2):
        codon = dna_sequence[i:i+3]
        if dna_sequence[i:i+3] == "ATG":  # find the first start codon
            protein = []
            for j in range(i, len(dna_sequence), 3):
                codon = dna_sequence[j:j+3]
                if len(codon) < 3:
                    break
                aa = codon_table.get(codon, 'X')  # unknown codons = 'X'
                if aa == '_':  # stop codon
                    return protein
                protein.append(aa)
            return protein
    return []  # if no ATG found


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
