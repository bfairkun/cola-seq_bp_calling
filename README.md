# CoLa-seq_BP_calling

This repository contains code for calling branchpoints (BPs) from cotranscriptional-lariat sequencing (CoLa-seq) data, as described in Yi's unpublished work. This method enriches for lariat derived cDNA molecules that are anchored at a BP at the 5' end (of R1 in paired end sequencing) and terminate (in the case of excised-lariat-derived cDNAs) or cross (in the case of nascent-lariat-intermediate-derived cDNAs) a 3' splice site. However, there is of course some background reads, so to identify true BPs, we implemented a classifier that takes into account the relative peaks of 5' ends from CoLa-seq data, as well as motif and relative distance to 3'ss, in order to classify bases as BPs. The methods will be described in more detail in the publication, and here is a more brief description the approach:

## BP calling approach

1. Identify 3'ss from an annotation source. We define windows of potential BPs as -100 to -5 relative to any annotated 3'ss.
1. Subset reads that terminate at, or cross at 3'ss, suggesting they may be derived from lariats
1. We want quantify peakiness of CoLa-seq 5' ends (relative to background) in the windows upstream of 3'ss. To ensure that our peakiness metrics aren't driven by a single read in a window with no background, we apply a few strategies:
    - Only consider 3'ss windows with >= 10 CoLa-seq 5' ends
    - Add a small pseudocount (0.1) of 5' end coverage along all bases
1. Use kernel density smoothing (`density` function in R, bandwidth=0.5) on the relative coverage in each window. The value of the smoothed kernel density estimate is assumed to be a probability that the base is a BP. Since 5' ends actually mark the base downstream of a branch, this probability function has to be shifted down a base.
1. To integrate positional information, we utilized the empirical distribution of BPs (Pineda et al) relative to 3'ss. Use a smoothed kernel density estimate to minimize noise/peakiness. The value of this smoothed density estimate is assumed to be an independent proababilistic measure that each base is a BP. To integrate that information with the previous step (assuming independence) we can just multiply the two probabilities (or equivalently, sum the log-probabilities)
1. To integrate motif information, we utilized the empirical distribution of BP motifs (Pineda et al) by creating a position weight matrix from Pineda et al. The probability that any particular base is a BP (based on sequence motif information that is) can be calculated by the position weight matrix. As before, we can add the log-odds probability to create a compositie log-probability score.
1. Define a threshold of log-probability scores for which to call bases as BPs. Based on the precision and recall of a high-confidence set of BPs (described in more detail in the manuscriipt), and assuming the same parameters for BP calling as in the manuscript (eg, using 5-base position weight matrix in previous step, using a pseudocount of 0.1 in the step measuring CoLa-seq 5' end coverage, etc), I recommend a threshold of -8.881075, as was used in the manuscript. This threshold yielded a 1% false positive rate on out of the gold standard bpts as described in the manuscript, and yields a similar number and quality of BPs as Pineda et al.
1. For every base that is considered a BP, merge BPs that are within two bases of each other, considering only the highest scoring BP as the true BP.


See the `code/README.md` for more details on how to run the code.
