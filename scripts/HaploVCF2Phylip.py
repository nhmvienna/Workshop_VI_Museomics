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
parser.add_option("--MinAlt", dest="MA",
                  help="numerical parameter", default=4)
parser.add_option("--MinCov", dest="MC",
                  help="numerical parameter", default=4)
parser.add_option("--MaxPropGaps", dest="MG",
                  help="numerical parameter", default=1)
parser.add_option("--names", dest="NA",
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


NAME = options.NA.split(",")


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

        if GT == "0":
            A = 0
            B = 1
        else:
            A = 1
            B = 0

        # obtain GTs, note that the REF PL is the second and the ALT Pl is the first in the list, thus we need to switch the order first, how confusing...
        PL = PLi.split(",")[::-1]

        # only consider GT if (1) the read-depth > than MC, (2) the Posterior Likelihood of the GT is > 50 and the PL of the other (non-called) GT is < 30 otherwise mark as ambiguous
        if int(DP) >= int(options.MC) and int(PL[A]) > 50 and int(PL[B]) < 30:
            TestGT.append(GT)

    # print(TestGT)
    # ignore positions with less than MA individuals with an alternative allele
    if TestGT.count("1") < int(options.MA):
        continue

    # loop through all samples
    for i in range(len(header)):
        if header[i] in EX:
            continue
        GT, PLi, DP, AD = pops[i].split(":")

        if GT == "0":
            A = 0
            B = 1
        else:
            A = 1
            B = 0

        # obtain GTs, note that the REF PL is the second and the ALT Pl is the first in the list, thus we need to switch the order first, how confusing...
        PL = PLi.split(",")[:: -1]

        # only consider GT if (1) the read-depth > than MC, (2) the Posterior Likelihood of the GT is > 50 and the PL of the other (non-called) GT is < 30 otherwise mark as ambiguous
        if int(DP) >= int(options.MC) and int(PL[A]) > 50 and int(PL[B]) < 30:
            GTs[header[i]][C] = ALLELE[int(GT)]
        else:
            GTs[header[i]][C] = "N"
            Del[header[i]][C]
    C += 1

# test if length of missing poistions larger than MG. if yes, ignore sample
for k, v in list(GTs.items()):
    if list(v.values()).count("N")/len(v.values()) > float(options.MG):
        del (GTs[k])
        del (Del[k])


# identify positions that contain a missing nucleotide
FinalDel = d(list)
for k, v in Del.items():
    for C in v.keys():
        FinalDel[C].append(k)

# Test if header already printed
TEST = ""
for i in header:

    # ignore samples that were excluded
    if i not in GTs:
        continue

    PL = []
    for j in range(1, C):
        if options.NG:
            if j in FinalDel:
                continue
        PL.append(GTs[i][j])
    if TEST == "":
        print(len(GTs.keys()), len(PL))
        TEST = 1
    print(i+"\t"+"".join(PL))
