#!/usr/bin/env python3

from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# from Bio.Alphabet import IUPAC
from collections import Counter
import numpy as np
import gzip
import sys

### GLOBAL VARS
# number of bases at the beginning of the motif before the branch
# for example, the offset for ctaAc is 3
motif_branch_offset = 3
PentamersFasta = "./Misc/BradleyPentamers.fasta"
WindowsFasta = "./Misc/3ssWindows.expanded.bed.fasta"
WindowNamesPassFilter = "BP_calling/MergedGroup/Cola3ssWindowsPassReadFilter.list.txt"
Output = "scratch/out.gz"

if __name__ == "__main__":
    (
        motif_branch_offset,
        PentamersFasta,
        WindowsFasta,
        WindowNamesPassFilter,
        Output,
    ) = sys.argv[1:]
    motif_branch_offset = int(motif_branch_offset)

###

print("Reading branchpoint instances")
instances_u2 = []
instances_u12 = []
for i, record in enumerate(SeqIO.parse(PentamersFasta, "fasta")):
    if i >= 0:
        if record.id.startswith("u2"):
            instances_u2.append(record.seq)
        elif record.id.startswith("u12"):
            instances_u12.append(record.seq)

m = motifs.create(instances_u2)
m_u12 = motifs.create(instances_u12)

print("Reading window sequences to calculate bkgrnd nucleotide freq")
PotentialRegions = ""
for i, record in enumerate(SeqIO.parse(WindowsFasta, "fasta")):
    if i == 1:
        RecordLength = len(record.seq)
    if i < 5000:
        # PotentialRegions.append(record.seq)
        PotentialRegions += record.seq
BaseCounter = Counter(PotentialRegions)
BackgroundNucFreqs = {
    Letter: Count / sum(BaseCounter.values()) for Letter, Count in BaseCounter.items()
}


print("Making position weight matrix")
print(m.counts)
print(m.consensus)
print(m.degenerate_consensus)
# m.weblogo("../output/BradleyWeightedBranchMotif.png")

print("Making position weight matrix for u12")
print(m_u12.counts)
print(m_u12.consensus)
print(m_u12.degenerate_consensus)
# m_u12.weblogo("../output/BradleyWeightedBranchMotif_U12.png")

pwm = m.counts.normalize()
pwm_u12 = m_u12.counts.normalize()
print(pwm.consensus)
print(pwm_u12.consensus)

# pssm = pwm.log_odds(BackgroundNucFreqs)
pssm = pwm.log_odds()
pssm_u12 = pwm_u12.log_odds()

print(pssm)

windowsPassFilter = set()
with open(WindowNamesPassFilter, "rt") as f_windows:
    for line in f_windows:
        windowsPassFilter.add(line.strip())
f_windows.close()

print(pssm.calculate(Seq("CTAACCTAACCTAAC")))

print("Scoring positions and writing output...")
print(
    "u2 motif used for u2 introns, u12 motif used for u12 introns... single output file"
)
BeginningNA = ["NA"] * motif_branch_offset
TailingNA = ["NA"] * (len(m) - motif_branch_offset - 1)
Header = "\t".join(["Pos"] + [str(i) for i in range(1, RecordLength + 1)])

with gzip.open(Output, "wt") as f:
    f.write(Header)
    for i, record in enumerate(SeqIO.parse(WindowsFasta, "fasta")):
        if i >= 0:
            windowName = record.id.split("::")[0]
            if windowName in windowsPassFilter:
                if record.id.split("(")[0].split("_")[-1] == "u12":
                    ScoresOut = (
                        BeginningNA
                        + list(np.around(pssm_u12.calculate(Seq(str(record.seq))), 3))
                        + TailingNA
                    )
                else:
                    ScoresOut = (
                        BeginningNA
                        + list(np.around(pssm.calculate(Seq(str(record.seq))), 3))
                        + TailingNA
                    )
                _ = f.write(
                    "\n"
                    + "\t".join(str(ele) for ele in ([windowName] + ScoresOut[::-1]))
                )
f.close()
