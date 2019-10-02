#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import numpy as np
import sys


def firstNonspace(ls):
    for i in ls:
        if i != "":
            break
    return i


def gc(seq):
    gc = 0
    for bp in seq:
        if bp == "C" or bp == "G":
            gc += 1
    return gc / len(seq)


def Dictparser(Dictionary):
    lowest = float(1000)
    for i in Dictionary:
        if float(Dictionary[i]) < float(lowest):
            lowest = Dictionary[i]
            key = i
    return [i, lowest]


def reverseComplement(seq):
    out = []
    for i in range(len(seq) - 1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def ribosome(seq):
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 4:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def ave(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count / len(ls)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count / len(seq)


def reject_outliers(data):
    m = 2
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
    return filtered


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls) - 1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


parser = argparse.ArgumentParser(
    prog="PseudoHunter.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
              @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
             @@@@@@@.........#@@@@@@@@@@@@@(.........@@@@@@@
            @@@@@@@@..,,.........@@@@@@@.........*,..@@@@@@@@
            @@@@@(...@@@@@@@@@@@.........@@@@@@@@@@@...,@@@@@
            @@@@@...@@@.......@@@&..%..@@@@@......@@@...@@@@@
            @@@@@@..@@@ PSEUDO @@.,&%%..@@@ HUNTER @@..@@@@@@
            @@@@@@@.@@@.......@@*.&&%%%.%@@@......@@@.@@@@@@@
            @@@@@@@.&@@@@@@@@@@..&&&%%%%..@@@@@@@@@@/.@@@@@@@
            @@@@@@@@..%@@@@@/..@&&&&%%%%%@..(@@@@@(.,@@@@@@@@
            @@@@@@@@@UGA@@@@@@&&&&&&%%%%%%%@@@@@@UAA@@@@@@@@@
            @@@@@@UAG@@@@@@@@@@...&&%%%...@@@@@@@@@@UGA@@@@@@
            @@@@UAA@@@@@@@@@.................@@@@@@@@UAG@@@@@
            @@@@@..@@@@@@@.....................@@@@@@%..@@@@@
            @@@@@..@@@@@............@............@@@@@..@@@@@
            @@@@@.................@@@@@.................@@@@@
             @@@@@@@..........&@@@@@@@@@@@#..........@@@@@@@
              @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    Developed by Arkadiy Garber and John McCutcheon
    University of Montana, Biological Sciences
    Please send comments and inquiries to arkadiy.garber@mso.umt.edu
    ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/ 

    You can provide this program with contigs via the \'-q and -r\'
    arguments, OR you can provide this program with pre-existing ORF-calls
    from Prodigal via the \'-n, -a, -rn, and -ra\' arguments. The ORF-calls
    must be from Prodigal, or formatted in such a way where the header name
    ends in an underscore and number (e.g. CP003881_1, CP003881_2, CP003881_3,...)

    The codeml control file is provided within the codeml program: paml4.8/codeml.ctl
    ************************************************************************
    '''))

# parser.add_argument('-blast', type=str, help="clu.tsv output from mmseqs", default="NA")

parser.add_argument('-q', type=str, help="input query genome (contigs in FASTA format)", default="NA")
parser.add_argument('-r', type=str,
                    help="input reference dataset (could be contigs from one genome or a set of closely-related genomes)",
                    default="NA")
parser.add_argument('-ref', type=str,
                    help="is the reference dataset one genome or multiple? (one/multiple) Default=one", default="one")

parser.add_argument('-n', type=str, help="genes in nucleic acid format", default="NA")
parser.add_argument('-a', type=str, help="genes in amino acid format", default="NA")
parser.add_argument('-rn', type=str, help="Reference genes in nucleic acid format", default="NA")
parser.add_argument('-ra', type=str, help="Reference genes in amino acid format", default="NA")

parser.add_argument('-gff', type=str, help="gff file for the query genome. Provide only if you are providing this program with"
                                           " genes, rather than contigs", default="NA")

parser.add_argument('-ctl', type=str, help="template control file for codeml", default="NA")

parser.add_argument('-out', type=str, help="name output directory", default="PseudoHunter_output")

parser.add_argument('-l', type=float,
                    help="minimum proportion of of target gene length that must be covered by alignment with query gene "
                         "for the query gene to be calssified as \'intact\' (default = 0.7)",
                    default=0.7)

parser.add_argument('-d', type=float, help="maximum dN/dS value for gene too be considered \'intact\' (default = 0.3)",
                    default=0.3)
parser.add_argument('-M', type=float, help="maximum dS value for dN/dS calculation (default = 3)", default=3)
parser.add_argument('-m', type=float, help="minimum dS value for dN/dS calculation (default = 0.001)", default=0.001)
parser.add_argument('-t', type=int, help="number of threads to use for BLAST", default=1)
parser.add_argument('-s', type=str, help="search engine to use (blast/diamond). Default = blast", default="blast")
parser.add_argument('--skip', type=str, help="By choosing this option, and providing pseudoHunter with the previously-created output directory, "
                                            "you are choosing to skip time-consuming steps of this pipeline "
                                            "(e.g. BLAST, Muscle, codeml), and would like to use the output files "
                                            "created from a previous run to re-do the anlaysis, "
                                            "perhaps with different parameters. All other arguments "
                                            "(e.g. \'-a\', \'-n\', \'-q\', \'-r\', and particulary \'-out\') still need to be provided as before", const=True, nargs="?")

args = parser.parse_args()


cwd = os.getcwd()

os.system("echo ${ctl} > ctl.txt")
file = open("ctl.txt")
for i in file:
    ctl = (i.rstrip())
os.system("rm ctl.txt")

# MAKING SURE ALL ARGUMENTS ARE THERE
# try:
#     test = os.listdir(cwd + "/dnds-analysis")
#     if args.o == "NA":
#         print("dnds-analysis directory already exists...please remove or rename prior to re-starting this pipeline")
#         raise SystemExit
# except FileNotFoundError:
#     os.system('mkdir ' + cwd + "/dnds-analysis")


# if you don't provide contigs, and use your own gene calls, but not provide a gff file:
if args.r == "NA" and args.q == "NA" and args.gff == "NA":
    if args.n != "NA" and args.rn != "NA" and args.a != "NA" and args.ra != "NA":
        print("Looks like you have provided your own gene calls to the program, "
              "but have not provided an accompanying gff file. PseudoHunter needs this gff file to get essential "
              "information about your gene calls. If you do not have the gff file, please provide raw contigs via the "
              "\'-q\' and \'-r\' arguments, so that this program can run Prodigal, and calculate this essential "
              "information itself. Alternatively, you can annotate your contigs using Prokka "
              "(https://github.com/tseemann/prokka), then come back here with "
              "all the essential files (gene calls and gff files)")
        print("")
        raise SystemExit
    else:
        print("Please check your command...it looks like one or more files needed to run PseudoHunter is missing.")
        print("")
        raise SystemExit

if args.r != "NA" and args.q != "NA" and args.gff != "NA":
    print("Looks like you have provided your genomes in raw contigs, and also a gff file."
          " If you would like to incorporate annotation data from the gff file, please provide associated gene calls "
          " via the \'-ra\', \'-rn\', \'-a\', and \'-r\' arguments. Otherwise, you can submit contigs as you have without the gff file.")
    print("")
    raise SystemExit

if not args.skip:
    print("Starting pipeline...")

    os.system("mkdir " + args.out)
    os.system('mkdir ' + args.out + "/dnds-analysis")

    if args.r != "NA" and args.q != "NA":
        query = args.q
        ref = args.r
        print("Predicting open-reading frames...")
        os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna" % (query, query, query))

        if args.ref == "one":
            os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna" % (ref, ref, ref))
        else:
            os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna -p meta" % (ref, ref, ref))

        # os.system("rm prodigal.out")
        # print("finished with prodigal")

        faaRef = open(ref + "-proteins.faa")
        faaRef = fasta2(faaRef)
        faa = open(query + "-proteins.faa")
        faa = fasta2(faa)

        fnaRef = open(ref + "-proteins.fna")
        fnaRef = fasta2(fnaRef)
        fna = open(query + "-proteins.fna")
        fna = fasta2(fna)

        if args.s == "blast":
            print("Running BLAST")
            os.system(
                "makeblastdb -dbtype prot -in %s-proteins.faa -out %s-proteins.faa" % (args.r, args.r))
            # os.system("rm makeblastdb.out")
            os.system("blastp -query %s-proteins.faa -db %s-proteins.faa "
                      "-outfmt 6 -out %s/pseudogene.blast -evalue 1E-6 -num_threads %s -max_target_seqs 1" % (
                          args.q, args.r, args.out, args.t))

            os.system("rm %s-proteins.faa.psq" % args.r)
            os.system("rm %s-proteins.faa.phr" % args.r)
            os.system("rm %s-proteins.faa.pin" % args.r)

        elif args.s == "diamond":
            print("Running DIAMOND")
            os.system(
                "diamond makedb --in %s-proteins.faa -d %s-proteins.faa &> makedb.out" % (args.r, args.r))
            os.system("rm makedb.out")
            os.system("diamond blastp --db %s-proteins.faa.dmnd --query %s-proteins.faa --outfmt 6 --out %s/pseudogene.blast "
                      "--max-target-seqs 1 --evalue 1E-6 --threads %d" % (args.r, args.q, args.out, args.t))

            os.system("rm %s-proteins.faa.dmnd" % args.r)

    else:
        faaRef = open(args.ra)
        faaRef = fasta2(faaRef)

        faa = open(args.a)
        faa = fasta2(faa)

        fnaRef = open(args.rn)
        fnaRef = fasta2(fnaRef)

        fna = open(args.n)
        fna = fasta2(fna)

        if args.s == "blast":
            print("Running BLAST")
            os.system("makeblastdb -dbtype prot -in %s -out %s &> makeblastdb.out" % (args.ra, args.ra))
            os.system("rm makeblastdb.out")
            os.system("blastp -query %s -db %s "
                      "-outfmt 6 -out %s/pseudogene.blast -evalue 1E-6 -num_threads %s -max_target_seqs 1" % (
                          args.a, args.ra, args.out, args.t))
        elif args.s == "diamond":
            print("Running DIAMOND")
            os.system(
                "diamond makedb --in %s -d %s &> makedb.out" % (args.ra, args.ra))
            os.system("rm makedb.out")
            os.system("diamond blastp --db %s.dmnd --query %s --outfmt 6 --out %s/pseudogene.blast "
                      "--max-target-seqs 1 --evalue 1E-6 --threads %d" % (args.ra, args.a, args.out, args.t))

            # os.system("rm %s.dmnd" % args.ra)

####################################################################################################################
    prescreened = []
    alnLengthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open("%s/pseudogene.blast" % args.out)
    for i in blast:
        ls = i.rstrip().split("\t")
        if ls[0] not in prescreened:

            alnLengthDict[ls[0]] = ls[3]

            outNUC = open(args.out + "/dnds-analysis/%s.faa.fna" % ls[0], "w")
            outNUC.write(">" + ls[1] + "\n")
            outNUC.write(fnaRef[ls[1]] + "\n")
            outNUC.write(">" + ls[0] + "\n")
            outNUC.write(fna[ls[0]] + "\n")
            outNUC.close()

            outAA = open(args.out + "/dnds-analysis/%s.faa" % ls[0], "w")
            outAA.write(">" + ls[1] + "\n")
            outAA.write(faaRef[ls[1]] + "\n")
            outAA.write(">" + ls[0] + "\n")
            outAA.write(faa[ls[0]] + "\n")
            outAA.close()

            if float(ls[3])/len(faaRef) > args.l:
                prescreened.append(ls[0])

    # ALIGNING PROTEIN SEQUENCES AND CREATING A CODON ALIGNMENT
    print("aligning files...")
    DIR = args.out + "/dnds-analysis"
    os.system("for i in %s/*faa; do"
              " muscle -in $i -out $i.aligned.fa &> muscle.out;"
              " rm muscle.out;"
              " pal2nal.pl $i.aligned.fa $i.fna -output fasta > $i.codonalign.fa;"
              " done" % DIR)

    # BUILDING CONTROL FILES
    print("preparing for codeml analysis")
    DIR = args.out + "/dnds-analysis"
    codealign = os.listdir(DIR)
    count = 0
    for file in codealign:
        if re.findall(r'codonalign', file):
            count += 1
            clu = file.split(".faa")[0]
            setup = open(ctl)
            out = open("%s/%s.ctl" % (DIR, str(clu)), "w")

            for i in setup:
                if re.findall('seqfile', i):
                    out.write(
                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' + args.out + '/dnds-analysis/' + file + ' ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                elif re.findall(r'outfile', i):
                    out.write('' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                              + args.out + '/dnds-analysis/mlcTree_' + str(
                        clu) + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' '
                              + '' + ' ' + '' + ' ' + '' + ' ' + '*' + ' ' + 'main' + ' ' + 'result' + ' ' + 'file' + ' ' + 'name\n')

                else:
                    out.write(i)
            out.close()
    print("")

    # RUNNING CODEML FOR DN/DS CALCULATION
    total = 0
    for i in codealign:
        if re.findall(r'codonalign', i):
            total += 1

    count = 0
    codealign = os.listdir(DIR)
    for file in codealign:
        if lastItem(file.split(".")) == "ctl":
            count += 1
            perc = (count / total) * 100
            sys.stdout.write("running codeml: %d%%   \r" % (perc))
            sys.stdout.flush()
            os.system("codeml %s/dnds-analysis/%s" % (args.out, file))
            # os.system("rm codeml.out")
    print("")


    # PARSING CODEML OUTPUT

    cwd = os.getcwd()
    DIR = args.out + "/dnds-analysis"


    if args.gff != "NA":
        gff = open(args.gff)
        gffDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'NA')))
        for i in gff:
            ls = i.rstrip().split("\t")
            if re.match(r'##FASTA', i):
                break
            else:
                if not re.match(r'#', i):
                    contig = ls[0]
                    orf = ls[8]
                    orf = orf.split(";")[0]
                    orf = orf.split("=")[1]

                    product = lastItem(ls[8].split(";")).split("=")[1]
                    product = replace(product, [","], ";")

                    gffDict[orf]["product"] = product
                    gffDict[orf]["contig"] = contig
                    gffDict[orf]["start"] = ls[3]
                    gffDict[orf]["end"] = ls[4]

    else:
        gffDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'NA')))
        faa = open(query + "-proteins.faa")
        faa = fasta(faa)

        for i in faa.keys():
            ls = i.split(" # ")
            start = ls[1]
            end = ls[2]
            contig = allButTheLast(ls[0], "_")
            gffDict[ls[0]]["contig"] = contig
            gffDict[ls[0]]["start"] = start
            gffDict[ls[0]]["end"] = end
            gffDict[ls[0]]["product"] = "NA"


    if args.r != "NA" and args.q != "NA":
        faaRef = open(ref + "-proteins.faa")
        faaRef = fasta2(faaRef)
        faa = open(query + "-proteins.faa")
        faa = fasta2(faa)

        fnaRef = open(ref + "-proteins.fna")
        fnaRef = fasta2(fnaRef)
        fna = open(query + "-proteins.fna")
        fna = fasta2(fna)
    else:
        faaRef = open(args.ra)
        faaRef = fasta2(faaRef)
        faa = open(args.a)
        faa = fasta2(faa)

        fnaRef = open(args.rn)
        fnaRef = fasta2(fnaRef)
        fna = open(args.n)
        fna = fasta2(fna)

    alnLengthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    alnIdDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open("%s/pseudogene.blast" % args.out)
    for i in blast:
        ls = i.rstrip().split("\t")
        alnLengthDict[ls[0]][ls[1]] = ls[3]
        alnIdDict[ls[0]][ls[1]] = ls[2]

    print("summarizing codeml output")
    codealign = os.listdir(DIR)
    dndsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in codealign:
        if re.findall(r'mlc', i):
            print(i)
            file = open(DIR + "/%s" % i, "r")
            for j in file:
                if re.search(r'#1', j):
                    ls = (j.rstrip().split(" "))
                    orf = ls[1]
                if re.search(r'#2', j):
                    ls = (j.rstrip().split(" "))
                    NODE = ls[1]
                line = j.rstrip()
            ls = line.split("  ")
            dS = remove(lastItem(ls), [" ", "=", "d", "S"])
            dN = remove(lastItem(ls[0:len(ls) - 1]), [" ", "=", "d", "N"])
            dndsDict[NODE]["orf"] = orf
            dndsDict[NODE]["dn"] = dN
            dndsDict[NODE]["ds"] = dS

    count = 0
    dndsList = []
    dndsDict2 = defaultdict(list)
    for i in sorted(dndsDict.keys()):
        count += 1
        if float(dndsDict[i]["dn"]) <= args.M and float(dndsDict[i]["ds"]) <= args.M and float(
                dndsDict[i]["ds"]) >= args.m and float(dndsDict[i]["dn"]) >= args.m:
            orf = dndsDict[i]["orf"]
            dn = dndsDict[i]["dn"]
            ds = dndsDict[i]["ds"]
            dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])
            dndsDict2[orf].append(i)
            dndsList.append(dnds)

    print("preparing final output file: summary.csv")
    dsList = []
    dnList = []
    dndsList2 = []
    out = open(args.out + "/summary.csv", "w")
    out.write(
        "ORF_calls" + "," + "Ortholog" + "," + "Pseudogene" + "," + "contig" + "," + "start" + "," + "end" + "," + "geneLength" + "," + "AlignmentLength" + "," +
        "OrthologLength" + "," + "Identity" + "," + "Annotation" + "," + "Pseudogene_confidence" + "," +
        "NumberOfGeneFrags" + "," + "AlignmentLength/OrthologLength" + "," + "geneLength/OrthologLength" + "," + "dN" + "," + "dS" + "," + "dN/dS" + "," +
        "Translation" + "," + "Sequences" + "\n")

    for i in dndsDict2.keys():
        if len(dndsDict2[i]) > 1:
            dndsDict3 = defaultdict(list)
            for k in dndsDict2[i]:
                dndsDict3[allButTheLast(k, "_")].append(int(lastItem(k.split("_"))))
            for l in dndsDict3.keys():
                listOfLists = (cluster(dndsDict3[l], 2))
                for m in listOfLists:
                    try:
                        if len(m) > 1:
                            seq = ''
                            seq2 = ''
                            dnLS = []
                            dsLS = []
                            annotations = ''
                            ORFs = ''
                            TotalAlnLength = 0

                            idLS = []

                            for n in m:
                                originalN = stabilityCounter(n)
                                # originalN = n
                                ORF = (l + "_" + str(originalN))

                                ORFs += ORF + "|"

                                annotation = gffDict[ORF]["product"]
                                annotations += annotation + "|"
                                seq += faa[ORF]
                                seq2 += fna[ORF]
                                dnLS.append(float(dndsDict[ORF]["dn"]))
                                dsLS.append(float(dndsDict[ORF]["ds"]))

                                alnLength = int(alnLengthDict[ORF][i])
                                TotalAlnLength += alnLength

                                identity = float(alnIdDict[ORF][i])
                                idLS.append(identity)

                            identity = ave(idLS)

                            dn = ave(dnLS)
                            ds = ave(dsLS)
                            dnds = dn / ds
                            fragments = len(m)
                            ORF = ORFs[0:len(ORFs) - 1]
                            annotation = annotations

                            contig = gffDict[ORF.split("|")[0]]["contig"]
                            start = gffDict[ORF.split("|")[0]]["start"]
                            end = gffDict[lastItem(ORF.split("|"))]["end"]

                        else:

                            fragments = 1
                            originalM0 = stabilityCounter(m[0])
                            # originalM0 = m[0]
                            ORF = (l + "_" + str(originalM0))

                            contig = gffDict[ORF]["contig"]
                            start = gffDict[ORF]["start"]
                            end = gffDict[ORF]["end"]

                            annotation = gffDict[ORF]["product"]
                            seq = faa[ORF]
                            seq2 = fna[ORF]

                            TotalAlnLength = int(alnLengthDict[ORF][i])

                            identity = float(alnIdDict[ORF][i])

                            dn = dndsDict[ORF]["dn"]
                            ds = dndsDict[ORF]["ds"]
                            dnds = float(dn) / float(ds)

                    except (TypeError, ValueError):
                        if len(m) > 1:
                            seq = ''
                            seq2 = ''
                            dnLS = []
                            dsLS = []
                            ORFs = ''
                            TotalAlnLength = 0

                            idLS = []

                            for n in m:
                                # originalN = stabilityCounter(n)
                                originalN = n
                                ORF = (l + "_" + str(originalN))


                                ORFs += ORF + "|"

                                annotation = gffDict[ORF]["product"]
                                annotations += annotation + "|"
                                seq += faa[ORF]
                                seq2 += fna[ORF]
                                dnLS.append(float(dndsDict[ORF]["dn"]))
                                dsLS.append(float(dndsDict[ORF]["ds"]))

                                alnLength = int(alnLengthDict[ORF][i])
                                TotalAlnLength += alnLength

                                identity = float(alnIdDict[ORF][i])
                                idLS.append(identity)

                            identity = ave(idLS)

                            dn = ave(dnLS)
                            ds = ave(dsLS)
                            dnds = dn / ds
                            fragments = len(m)
                            ORF = ORFs[0:len(ORFs)-1]
                            annotation = annotations

                            contig = gffDict[ORF.split("|")[0]]["contig"]
                            start = gffDict[ORF.split("|")[0]]["start"]
                            end = gffDict[lastItem(ORF.split("|"))]["end"]

                        else:

                            fragments = 1
                            # originalM0 = stabilityCounter(m[0])
                            originalM0 = m[0]
                            ORF = (l + "_" + str(originalM0))

                            contig = gffDict[ORF]["contig"]
                            start = gffDict[ORF]["start"]
                            end = gffDict[ORF]["end"]

                            annotation = gffDict[ORF]["product"]
                            seq = faa[ORF]
                            seq2 = fna[ORF]

                            TotalAlnLength = int(alnLengthDict[ORF][i])

                            identity = float(alnIdDict[ORF][i])

                            dn = dndsDict[ORF]["dn"]
                            ds = dndsDict[ORF]["ds"]
                            dnds = float(dn) / float(ds)

        else:
            fragments = 1
            ORF = dndsDict2[i][0]

            start = gffDict[ORF]["start"]
            end = gffDict[ORF]["end"]

            annotation = gffDict[ORF]["product"]
            seq = faa[ORF]
            seq2 = fna[ORF]
            dn = dndsDict[ORF]["dn"]
            ds = dndsDict[ORF]["ds"]
            dnds = float(dn) / float(ds)

            TotalAlnLength = int(alnLengthDict[ORF][i])

            identity = float(alnIdDict[ORF][i])

        ratio = TotalAlnLength/len(faaRef[i])

        out.write(ORF + "," + i + ',')

        if dnds > args.d or ratio < args.l or fragments > 1 or ratio / (len(seq)/len(faaRef[i])) < args.l:
            out.write("Y" + ",")
            prob = 1

            # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
            inflation = float(dnds) / ave(dndsList)

            prob = prob * inflation
            prob = prob * float(fragments)
            if ratio < 1:
                prob = prob / ratio
            else:
                prob = prob * ratio

        else:
            out.write("N" + ",")

            # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
            prob = 1
            inflation = float(dnds) / ave(dndsList)
            prob = prob * inflation
            prob = prob * float(fragments)
            if ratio < 1:
                prob = prob / ratio
            else:
                prob = prob * ratio

        # WRITING TO FILE
        out.write(str(contig) + "," + str(start) + "," + str(end) + "," + str(len(seq)) + "," + str(TotalAlnLength) + "," +
                  str(len(faaRef[i])) + "," + str(identity) + "," + str(annotation) + "," + str(prob) + "," +
                  str(fragments) + "," + str(ratio) + "," + str(len(seq)/len(faaRef[i])) + "," + str(dn) + "," +
                  str(ds) + "," + str(dnds) + "," + seq + "," + seq2 + "\n")

        dnList.append(dn)
        dsList.append(ds)
        dndsList2.append(dnds)

    out.close()

    if len(dsList) > 0:
        print("")
        print("Identified " + str(count) + " orthologs in reference dataset")
        print("Average dN among orthologs: " + str(ave(dnList)))
        print("Average dS among orthologs: " + str(ave(dsList)))
        print("Average dN/dS among orthologs: " + str(ave(dndsList2)))
        print("")
        os.system("rm 2NG.t 2NG.dN 2NG.dS rst1 rst 2ML.t 2ML.dN 2ML.dS 4fold.nuc rub")
        print("Pipeline finished without any crashes. Thanks for using pseudoHunter!")
    else:
        print("")
        print("Identified " + str(count) + " orthologs in reference dataset")
        print("It looks like no orthologs were identified below the specified dS threshold of: " + str(args.M))
        print("Please try running again, with a higher value for -M")
        print("Dont forget to add the \'-o %s\' flag, so that you don't need to wait for codeml to run again." % args.out)


else:
    cwd = os.getcwd()
    DIR = args.out + "/dnds-analysis"

    if args.gff != "NA":
        gff = open(args.gff)
        gffDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'NA')))
        for i in gff:
            ls = i.rstrip().split("\t")
            if re.match(r'##FASTA', i):
                break
            else:
                if not re.match(r'#', i):
                    contig = ls[0]
                    orf = ls[8]
                    orf = orf.split(";")[0]
                    orf = orf.split("=")[1]

                    product = lastItem(ls[8].split(";")).split("=")[1]
                    product = replace(product, [","], ";")

                    gffDict[orf]["product"] = product
                    gffDict[orf]["contig"] = contig
                    gffDict[orf]["start"] = ls[3]
                    gffDict[orf]["end"] = ls[4]

    else:
        gffDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'NA')))
        query = args.q
        ref = args.r

        faa = open(query + "-proteins.faa")
        faa = fasta(faa)

        for i in faa.keys():
            ls = i.split(" # ")
            start = ls[1]
            end = ls[2]
            contig = allButTheLast(ls[0], "_")
            gffDict[ls[0]]["contig"] = contig
            gffDict[ls[0]]["start"] = start
            gffDict[ls[0]]["end"] = end
            gffDict[ls[0]]["product"] = "NA"

    if args.r != "NA" and args.q != "NA":
        query = args.q
        ref = args.r
        faaRef = open(ref + "-proteins.faa")
        faaRef = fasta2(faaRef)
        faa = open(query + "-proteins.faa")
        faa = fasta2(faa)

        fnaRef = open(ref + "-proteins.fna")
        fnaRef = fasta2(fnaRef)
        fna = open(query + "-proteins.fna")
        fna = fasta2(fna)
    else:
        faaRef = open(args.ra)
        faaRef = fasta2(faaRef)
        faa = open(args.a)
        faa = fasta2(faa)

        fnaRef = open(args.rn)
        fnaRef = fasta2(fnaRef)
        fna = open(args.n)
        fna = fasta2(fna)

    alnLengthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    alnIdDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open(args.out + "/pseudogene.blast")
    for i in blast:
        ls = i.rstrip().split("\t")
        alnLengthDict[ls[0]][ls[1]] = ls[3]
        alnIdDict[ls[0]][ls[1]] = ls[2]

    print("summarizing codeml output")
    codealign = os.listdir(DIR)
    dndsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in codealign:
        if re.findall(r'mlc', i):
            file = open(DIR + "/%s" % i, "r")
            for j in file:
                if re.search(r'#1', j):
                    ls = (j.rstrip().split(" "))
                    orf = ls[1]
                if re.search(r'#2', j):
                    ls = (j.rstrip().split(" "))
                    NODE = ls[1]
                line = j.rstrip()
            ls = line.split("  ")
            dS = remove(lastItem(ls), [" ", "=", "d", "S"])
            dN = remove(lastItem(ls[0:len(ls) - 1]), [" ", "=", "d", "N"])
            dndsDict[NODE]["orf"] = orf
            dndsDict[NODE]["dn"] = dN
            dndsDict[NODE]["ds"] = dS

    count = 0
    dndsList = []
    dndsDict2 = defaultdict(list)
    for i in sorted(dndsDict.keys()):
        count += 1
        if float(dndsDict[i]["dn"]) <= args.M and float(dndsDict[i]["ds"]) <= args.M and float(
                dndsDict[i]["ds"]) >= args.m and float(dndsDict[i]["dn"]) >= args.m:
            orf = dndsDict[i]["orf"]
            dn = dndsDict[i]["dn"]
            ds = dndsDict[i]["ds"]
            dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])
            dndsDict2[orf].append(i)
            dndsList.append(dnds)

    print("preparing final output file: summary.csv")
    dsList = []
    dnList = []
    dndsList2 = []
    out = open(args.out + "/summary.csv", "w")
    out.write(
        "ORF_calls" + "," + "Ortholog" + "," + "Pseudogene" + "," + "contig" + "," + "start" + "," + "end" + "," + "geneLength" + "," + "AlignmentLength" + "," +
        "OrthologLength" + "," + "Identity" + "," + "Annotation" + "," + "Pseudogene_confidence" + "," +
        "NumberOfGeneFrags" + "," + "AlignmentLength/OrthologLength" + "," + "geneLength/OrthologLength" + "," + "dN" + "," + "dS" + "," + "dN/dS" + "," +
        "Translation" + "," + "Sequences" + "\n")

    for i in dndsDict2.keys():
        if len(dndsDict2[i]) > 1:
            dndsDict3 = defaultdict(list)
            for k in dndsDict2[i]:
                dndsDict3[allButTheLast(k, "_")].append(int(lastItem(k.split("_"))))
            for l in dndsDict3.keys():
                listOfLists = (cluster(dndsDict3[l], 2))
                for m in listOfLists:
                    try:
                        if len(m) > 1:
                            seq = ''
                            seq2 = ''
                            dnLS = []
                            dsLS = []
                            annotations = ''
                            ORFs = ''
                            TotalAlnLength = 0

                            idLS = []

                            for n in m:
                                originalN = stabilityCounter(n)
                                # originalN = n
                                ORF = (l + "_" + str(originalN))

                                ORFs += ORF + "|"

                                annotation = gffDict[ORF]["product"]
                                annotations += annotation + "|"
                                seq += faa[ORF]
                                seq2 += fna[ORF]
                                dnLS.append(float(dndsDict[ORF]["dn"]))
                                dsLS.append(float(dndsDict[ORF]["ds"]))

                                alnLength = int(alnLengthDict[ORF][i])
                                TotalAlnLength += alnLength

                                identity = float(alnIdDict[ORF][i])
                                idLS.append(identity)

                            identity = ave(idLS)

                            dn = ave(dnLS)
                            ds = ave(dsLS)
                            dnds = dn / ds
                            fragments = len(m)
                            ORF = ORFs[0:len(ORFs) - 1]
                            annotation = annotations

                            contig = gffDict[ORF.split("|")[0]]["contig"]
                            start = gffDict[ORF.split("|")[0]]["start"]
                            end = gffDict[lastItem(ORF.split("|"))]["end"]

                        else:

                            fragments = 1
                            originalM0 = stabilityCounter(m[0])
                            # originalM0 = m[0]
                            ORF = (l + "_" + str(originalM0))

                            contig = gffDict[ORF]["contig"]
                            start = gffDict[ORF]["start"]
                            end = gffDict[ORF]["end"]

                            annotation = gffDict[ORF]["product"]
                            seq = faa[ORF]
                            seq2 = fna[ORF]

                            TotalAlnLength = int(alnLengthDict[ORF][i])

                            identity = float(alnIdDict[ORF][i])

                            dn = dndsDict[ORF]["dn"]
                            ds = dndsDict[ORF]["ds"]
                            dnds = float(dn) / float(ds)

                    except (TypeError, ValueError):
                        if len(m) > 1:
                            seq = ''
                            seq2 = ''
                            dnLS = []
                            dsLS = []
                            ORFs = ''
                            TotalAlnLength = 0

                            idLS = []

                            for n in m:
                                # originalN = stabilityCounter(n)
                                originalN = n
                                ORF = (l + "_" + str(originalN))

                                ORFs += ORF + "|"

                                annotation = gffDict[ORF]["product"]
                                annotations += annotation + "|"
                                seq += faa[ORF]
                                seq2 += fna[ORF]
                                dnLS.append(float(dndsDict[ORF]["dn"]))
                                dsLS.append(float(dndsDict[ORF]["ds"]))

                                alnLength = int(alnLengthDict[ORF][i])
                                TotalAlnLength += alnLength

                                identity = float(alnIdDict[ORF][i])
                                idLS.append(identity)

                            identity = ave(idLS)

                            dn = ave(dnLS)
                            ds = ave(dsLS)
                            dnds = dn / ds
                            fragments = len(m)
                            ORF = ORFs[0:len(ORFs) - 1]
                            annotation = annotations

                            contig = gffDict[ORF.split("|")[0]]["contig"]
                            start = gffDict[ORF.split("|")[0]]["start"]
                            end = gffDict[lastItem(ORF.split("|"))]["end"]

                        else:

                            fragments = 1
                            # originalM0 = stabilityCounter(m[0])
                            originalM0 = m[0]
                            ORF = (l + "_" + str(originalM0))

                            contig = gffDict[ORF]["contig"]
                            start = gffDict[ORF]["start"]
                            end = gffDict[ORF]["end"]

                            annotation = gffDict[ORF]["product"]
                            seq = faa[ORF]
                            seq2 = fna[ORF]

                            TotalAlnLength = int(alnLengthDict[ORF][i])

                            identity = float(alnIdDict[ORF][i])

                            dn = dndsDict[ORF]["dn"]
                            ds = dndsDict[ORF]["ds"]
                            dnds = float(dn) / float(ds)

        else:
            fragments = 1
            ORF = dndsDict2[i][0]

            start = gffDict[ORF]["start"]
            end = gffDict[ORF]["end"]

            annotation = gffDict[ORF]["product"]
            seq = faa[ORF]
            seq2 = fna[ORF]
            dn = dndsDict[ORF]["dn"]
            ds = dndsDict[ORF]["ds"]
            dnds = float(dn) / float(ds)

            TotalAlnLength = int(alnLengthDict[ORF][i])

            identity = float(alnIdDict[ORF][i])

        ratio = TotalAlnLength / len(faaRef[i])

        out.write(ORF + "," + i + ",")

        if dnds > args.d or ratio < args.l or fragments > 1 or ratio / (len(seq) / len(faaRef[i])) < args.l:
            out.write("Y" + ",")
            prob = 1

            # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
            inflation = float(dnds) / ave(dndsList)

            prob = prob * inflation
            prob = prob * float(fragments)
            if ratio < 1:
                prob = prob / ratio
            else:
                prob = prob * ratio

        else:
            out.write("N" + ",")

            # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
            prob = 1
            inflation = float(dnds) / ave(dndsList)
            prob = prob * inflation
            prob = prob * float(fragments)
            if ratio < 1:
                prob = prob / ratio
            else:
                prob = prob * ratio

        # WRITING TO FILE
        out.write(
            str(contig) + "," + str(start) + "," + str(end) + "," + str(len(seq)) + "," + str(TotalAlnLength) + "," +
            str(len(faaRef[i])) + "," + str(identity) + "," + str(annotation) + "," + str(prob) + "," +
            str(fragments) + "," + str(ratio) + "," + str(len(seq) / len(faaRef[i])) + "," + str(dn) + "," +
            str(ds) + "," + str(dnds) + "," + seq + "," + seq2 + "\n")

        dnList.append(dn)
        dsList.append(ds)
        dndsList2.append(dnds)

    out.close()
    if len(dsList) > 0:
        print("")
        print("Identified " + str(count) + " orthologs in reference dataset")
        print("Average dN among orthologs: " + str(ave(dnList)))
        print("Average dS among orthologs: " + str(ave(dsList)))
        print("Average dN/dS among orthologs: " + str(ave(dndsList2)))
        print("")

        print("Pipeline finished without any crashes. Thanks for using pseudoHunter!")
    else:
        print("")
        print("Identified " + str(count) + " orthologs in reference dataset")
        print("It looks like no orthologs were identified below the specified dS threshold of: " + str(args.M))
        print("Please try running again, with a higher value for -M")
        print("Dont forget to add the \'-o %s\' flag, so that you don't need to wait for codeml to run again." % args.out)