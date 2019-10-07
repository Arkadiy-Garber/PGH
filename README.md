# PseudoHunter
## Pseudogene Identification Program 
Identification of pseudogenes based on comparison between one or more genomes

### Inputs
To use this program, please provide contigs or gene-calls in FASTA format.If you would like to predict pseudogenes from previously-predicted gene calls, then you must supply a GFF file. Otherwise, PseudoHunter will take contigs in FASTA format and perform its own gene predictions using Prodigal.

### Reference dataset
To use PseudoHunter must provide a reference dataset, which can consist of either contigs or gene-calls in FASTA format. This reference dataset can be a single genome, or a collection of genomes. PseudoHunter will use these genomes as a benchmark to predict which genes are pseudogenized in your dataset-of-interest; thus, you must be sure that whatever pseudogenization has occurred in your dataset is not also present in your reference genomes.

### Intergenic regions
PseudoHunter can also looks in intergenic regions for pseudogenes. In this case, please be sure to provide contigs in FASTA format (regardless of whether you are providing gene calls + GFF inputs)

### Output
Inside the output directory, you will find a CSV file names "summary.csv". This file will contain information on all predicted ORFs in your dataset, including which are predicted to be pseudogenes, as well as dN/dS ratios, proportion of exected gene length, fragmentation due to stop mutations, etc.


### easy-installation with conda
    git clone https://github.com/Arkadiy-Garber/PseudoHunter.git
    cd PGH
    ./setup.sh
    source activate pseudo
    PseudoHunter4.py -h

### quickstart with raw contigs
    PseudoHunter4.py -q contigs.fna -r referenceContigs.fna

### quickstart with annotated genes
    PseudoHunter4.py -n genesNucleicAcids.ffn -a genesAminoAcids.faa -rn referenceNucleicAcids.ffn -ra referenceAminoAcids.faa -gff genes.gff
    
### quickstart with annotated genes and contigs
    PseudoHunter4.py -n genesNucleicAcids.ffn -a genesAminoAcids.faa -rn referenceNucleicAcids.ffn -ra referenceAminoAcids.faa -gff genes.gff -q contigs.fna -out PseudoHunter_output

### altering paramaters and skipping time-consuming steps by using previously-created output
    PseudoHunter4.py -n genesNucleicAcids.ffn -a genesAminoAcids.faa -rn referenceNucleicAcids.ffn -ra referenceAminoAcids.faa -gff genes.gff -q contigs.fna -out PseudoHunter_output --skip

