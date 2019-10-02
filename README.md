# PGH
Identification of pseudogenes based on pairwise comparison between genomes

### easy-installation with conda
    git clone https://github.com/Arkadiy-Garber/PGH.git
    cd PGH
    ./setup.sh
    source activate pseudo
    PseudoHunter2.py -h

### quickstart with annotated genes (required files: genes in nucleic acid and amino acid FASTA formats, as well as a gff file)
    PseudoHunter2.py -n genesNucleicAcids.ffn -a genesAminoAcids.faa -rn referenceNucleicAcids.ffn -ra referenceAminoAcids.faa -gff genes.gff
    
### quickstart with raw contig (required files: genes in nucleic acid and amino acid formats, as well as a gff file)
