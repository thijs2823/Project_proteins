# Proteins project
Various protein properties

---


## RNA to Protein Translator
This Python script reads one or multiple FASTA-formatted RNA/DNA sequences, splits them into codons, translates them into protein sequences using a codon table, and prints the first 100 amino acids along with the protein length. Unknown codons are marked as `X`.


DNA (template strand) --> Transcription (RNA polymerase) --> mRNA (messenger RNA) --> Translation (ribosome) --> Protein (amino acid sequence)


### Features

* Handles multiple sequences in a single FASTA file.
* Skips headers and concatenates sequence lines.
* Translates DNA/RNA sequences into amino acids.
* Starts translation from the first ATG codon (start codon); ignores UTRs before CDS.
* Stops translation at the first stop codon.
* Displays first 100 amino acids and protein length.

### Requirements

* Python 3.x


### Example Output (TP53)

FASTA file [1] contains:

```
>NM_000546.6 TP53 [organism=Homo sapiens] [GeneID=7157] [transcript=1]
CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGC...
```

Script output:

```
> NM_000546.6 TP53 [organism=Homo sapiens] [GeneID=7157] [transcript=1]
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD ...
Length protein: 393
```

* `>` indicates the FASTA header.
* The protein sequence shows amino acids, with `X` for unknown codons.
* Length is the total number of amino acids before the first stop codon.

1. National Center for Biotechnology Information. TP53 gene (Homo sapiens) [Internet]. Bethesda (MD): National Library of Medicine (US). Available from: https://www.ncbi.nlm.nih.gov/datasets/gene/7157/

## BiLSTM model for protein structure prediction
This project contains a Bi-directional LSTM (BiLSTM) model for predicting protein secondary structures (Q8 classification) using the NEW364 dataset.
The Q8 classification refers to a scheme for predicting protien secondary structures. Eeah amino adcid in a protien sequences is assigned to one of eight structured states.
- H: Alpha helix
- G: 3₁₀ helix
- I: Pi helix
- E: Beta strand
- B: Beta bridge
- T: Turn
- S: Bend
- C: Coil (loop or irregular structure) 

The BiLSTM model is used due to its ability to leverage both forward and backward sequence information.

### Dataset
- The model uses the NEW364 CSV dataset, formatted with two columns:
  - input: Amino acid sequence
  - dssp8: Corresponding secondary structure labels (8 classes
- Sequence length is padded or truncated to 700 residues (BiLSTM expects the same length for each input)
- Split: 80% training, 20% validation.

### Training Parameters
- Batch size: 16
- Learning rate: 0.001
- Epochs: 100
- Hidden dimensions: 128
- LSTM layers: 2

### Results

<img width="600" height="734" alt="image" src="https://github.com/user-attachments/assets/79a59975-8a9d-427a-a6cd-a3b3157e9b39" />


Figure 1. Training and validation loss, and validation accuracy

### Conclusion
After 50 epochs, the trainings loss and validation loss grow apart from earch other, detriment of validation loss; this indicates the beginning of overheating of the dataset. The model reaches validation accuracy of around 77%.


