import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--NoGaps", dest="NG",
                  help="logical parameter", action="store_true")
parser.add_option("--MinCov", dest="MC",
                  help="numerical parameter", default=4)
parser.add_option("--MaxPropGaps", dest="MG",
                  help="numerical parameter", default=1)
parser.add_option("--exclude", dest="EX", help="Input file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


GTs = d(lambda: d(list))
Del = d(lambda: d(list))
C = 1
if options.EX != None:
    EX = options.EX.split(",")
else:
    EX = []

Positions = d(str)
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    a = l.rstrip().split()

    # get names from headers
    if l.startswith("#"):
        header = [x.split("/")[-1].split(".bam")[0] for x in a[9:]]
        continue

    # obtain alleles
    REF = a[3]
    ALT = a[4]
    ALLELE = [REF, ALT]

    # ignore tri- and tetra-allelic SNPs
    if len(ALT) > 1 or len(REF) > 1:
        continue

    pops = a[9:]

    TestGT = []
    for i in range(len(header)):
        if header[i] in EX:
            continue
        GT, PLi, DP, AD = pops[i].split(":")

        if int(DP) < int(options.MC):
            GTs[header[i]][C] = "N"

        if GT == "0/0":
            GTs[header[i]][C] = REF
        elif GT == "0/1" or GT == "1/0":
            import random
            GTs[header[i]][C] = random.choice([ALT, REF])
        elif GT == "1/1":
            GTs[header[i]][C] = ALT
        elif GT == "./.":
            GTs[header[i]][C] = "N"

    C += 1

# test if length of missing poistions larger than MG. if yes, ignore sample
for k, v in list(GTs.items()):
    if list(v.values()).count("N")/len(v.values()) > float(options.MG):
        del (GTs[k])


# Test if header already printeds
TEST = ""
for i in header:

    # ignore samples that were excluded
    if i not in GTs:
        continue

    PL = []
    for j in range(1, C):
        PL.append(GTs[i][j])
    if TEST == "":
        print(len(GTs.keys()), len(PL))
        TEST = 1
    print(i+"\t"+"".join(PL))
