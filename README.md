# PGH
## Pseudogene Hunter 
Identification of pseudogenes based on pairwise comparison between one or more genomes. To use this program, please provide contigs or gene-calls in FASTA format. You must also provide a reference dataset, which must consist of either contigs or gene-calls in FASTA format.

### easy-installation with conda
    git clone https://github.com/Arkadiy-Garber/PGH.git
    cd PGH
    ./setup.sh
    source activate pseudo
    PseudoHunter2.py -h

### quickstart with raw contig
    PseudoHunter2.py -q contigs.fna -r referenceContigs.fna -out PseudoOutput/

### quickstart with annotated genes (required files: genes in nucleic acid and amino acid FASTA formats, as well as an affiliated gff file)
    PseudoHunter2.py -n genesNucleicAcids.ffn -a genesAminoAcids.faa -rn referenceNucleicAcids.ffn -ra referenceAminoAcids.faa -gff genes.gff -out PseudoOutput/

