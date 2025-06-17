import sys
from collections import defaultdict as d
import re
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun
# edits by  Siyuan Feng

#########################################################   HELP   #########################################################################
usage = """python3 %prog \
      --mpileup data.mpileup \
      --min-cov 10 \
      --max-cov data.cov \
      --min-count 10 \
      --min-freq 0.01 \
      --mis-frac 0.1 \
      --base-quality-threshold 15 \
      --names Kib32,Tam10 \
      --coding 1.8 \
      > output.vcf"""
parser = OptionParser(usage=usage)
helptext = """

H E L P :
_________
"""
group = OptionGroup(parser, helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--mpileup", dest="m", help="A mpileup file")
parser.add_option(
    "--min-cov", dest="minc", help="The minimum coverage threshold: e.g. 10", default=10
)
parser.add_option(
    "--min-count",
    dest="mint",
    help="The minimum number of counts of the alternative allele across all samples pooled",
    default=3,
)
parser.add_option(
    "--min-freq",
    dest="minf",
    help="The minimum Frequency of the alternative allele across all samples pooled",
    default=0.01,
)
parser.add_option(
    "--miss-frac",
    dest="mis",
    help="The minimum Frequency of the alternative allele across all samples pooled",
    default=0.1,
)
parser.add_option(
    "--base-quality-threshold",
    dest="b",
    help="The Base-quality threshold for Qualities encoded in Sanger format (Illumina 1.8 format)",
    default=15,
)
parser.add_option(
    "--names",
    dest="n",
    help="a comma separted list of thenames of all samples in the mpileup file",
)
parser.add_option(
    "--coding", dest="c", help="the Illumina FASTQ quality coding", default=1.8
)
parser.add_option(
    "--allsites",
    dest="AS",
    help="1 if all sites should be report; 0 if only polymorphic sites",
)
parser.add_option_group(group)
(options, args) = parser.parse_args()


################################### functions ######################################


def load_data(x):
    """ import data either from a gzipped or or uncrompessed file or from STDIN"""
    import gzip

    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "r")
    else:
        y = open(x, "r")
    return y


def keywithmaxvalue(d):
    """ This function resturns the key for the maximum value in a dictionary"""
    newhash = d(list)
    for k, v in d.items():
        newhash[v].append(k)
    return newhash[max(newhash.keys())]


def splitter(l, n):
    """ This generator function returns equally sized cunks of an list"""
    # credit: Meric Lieberman, 2012
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i: i + n]


def extract_indel(l, sign):
    """ This function returns an Indel from a sequence string in a pileup"""
    position = l.index(sign)
    numb = ""
    i = 0
    while True:
        if l[position + 1 + i].isdigit():
            numb += l[position + 1 + i]
            i += 1
        else:
            break

    seqlength = int(numb)
    sequence = l[position: position + i + 1 + seqlength]
    indel = sequence.replace(numb, "")

    return sequence, indel


################################## parameters ########################################


data = options.m
minimumcov = int(options.minc)
minimumcount = int(options.mint)
minimumfreq = float(options.minf)
missfrac = float(options.mis)
baseqthreshold = int(options.b)
phred = float(options.c)
names = options.n.split(",") if options.n else None


############################ calculate PHRED cutoff  #############################

# calculate correct PHRED score cutoff: ASCII-pc

if phred >= 1.0 and phred < 1.8:
    pc = 64
else:
    pc = 33


############################ parse MPILEUP ###########################################

# parse mpileup and store alternative alleles:

alleles = d(lambda: d(lambda: d(lambda: d(int))))
for line in load_data(data):
    if len(line.split("\t")) < 2:
        continue

    k = line[:-1].split("\t")
    CHR, POS, REF = k[:3]

    div = list(splitter(k, 3))
    libraries = div[1:]
    # loop through libraries

    for j in range(len(libraries)):
        alleles[j]
        nuc = libraries[j][1]
        qualities = libraries[j][2]

        # test if seq-string is empty
        if nuc == "*":
            continue

        # find and remove read indices and mapping quality string
        nuc = re.sub(r"\^.", r"", nuc)
        nuc = nuc.replace("$", "")

        DelCount = nuc.count("-1")
        allele = d(lambda: 0)
        allele[REF+">Del"] = DelCount

        # find and remove InDels
        while "+" in nuc or "-" in nuc:
            if "+" in nuc:
                insertion, ins = extract_indel(nuc, "+")
                nuc = nuc.replace(insertion, "")
            else:
                deletion, dele = extract_indel(nuc, "-")
                nuc = nuc.replace(deletion, "")

        # test for base quality threshold (if below: ignore nucleotide)
        # print len(nuc),len(qualities)
        nuc = "".join(
            [
                nuc[x]
                for x in range(len(nuc))
                if ord(qualities[x]) - pc >= baseqthreshold
            ]
        )
        nuc = "".join([nuc[x] for x in range(len(nuc)) if nuc[x] != "*"])

        # read all alleles
        for i in range(len(nuc)):

            # ignore single nucleotide deletions
            if nuc[i] == "*" or nuc[i] == "N":
                continue
            # count nucleotides similar to reference base
            if nuc[i] == "," or nuc[i] == ".":
                continue
            # count alternative nucleotides
            allele[REF+">"+nuc[i].upper()] += 1
        for k, v in allele.items():
            FREQ = v / sum(allele.values()) if sum(allele.values()) > 0 else 0
            if FREQ == 0:
                continue
            print(",".join([str(x) for x in [CHR, POS, names[j], k, FREQ]]))
