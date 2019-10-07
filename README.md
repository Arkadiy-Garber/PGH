# PseudoHunter
## Pseudogene Identification Program 
Identification of pseudogenes based on comparison between one or more genomes. To use this program, please provide contigs or gene-calls in FASTA format. You must also provide a reference dataset, which must consist of either contigs or gene-calls in FASTA format.

### easy-installation with conda
    git clone https://github.com/Arkadiy-Garber/PGH.git
    cd PGH
    ./setup.sh
    source activate pseudo
    PseudoHunter2.py -h

### quickstart with raw contigs
    PseudoHunter4.py -q contigs.fna -r referenceContigs.fna -out PseudoOutput/

### quickstart with annotated genes
    PseudoHunter4.py -n genesNucleicAcids.ffn -a genesAminoAcids.faa -rn referenceNucleicAcids.ffn -ra referenceAminoAcids.faa -gff genes.gff -out PseudoOutput/

