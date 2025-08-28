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
