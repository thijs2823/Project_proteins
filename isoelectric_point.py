def simple_pI(sequence): #Simple estimate of isoelectric point. Calculates approximate isoelectric point for the given protein sequence.
    """
    pKa_values = {
        'K': 10.8, 'R': 12.5, 'H': 6.5,
        'D': 3.9, 'E': 4.1, 'C': 8.5, 'Y': 10.1
    }

    ionisable = [pKa_values[aa] for aa in sequence if aa in pKa_values]

    if ionisable:
        return sum(ionisable) / len(ionisable)
    else:
        print('Warning: No ionisable residues found.')
        return float('nan')

def get_sequence_from_txt(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    seq_parts = []
    for line in lines:
        line = line.strip()
        if not line.startswith(">"):  # ignore headers
            seq_parts.append(line)

    return ''.join(seq_parts)

# Use
filepath = r'C:\file'
hlf_seq = get_sequence_from_txt(filepath)

print(hlf_seq)
print(f"Estimated pI (hlf): {simple_pI(hlf_seq):.2f}")
