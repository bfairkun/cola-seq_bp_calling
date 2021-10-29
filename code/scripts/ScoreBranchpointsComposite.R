library(tidyverse)
library(data.table)
library(stringr)
library(knitr)
library(readr)

### This script outputs unscaled posterior probabilites for BP classification.

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 "BP_calling/MergedFromFastq/ColaPutativeBranchFragmentStarts.tab BP_calling/MergedFromFastq/ 0.1 0.5 10 5 100 -8.881075 FALSE", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}


BpToDistanceTrainingData = "Misc/BradleyDistTo3ss.txt.gz"
MotifProbabilityFunction = "BP_calling/MergedGroup/MotifPerBaseScores.txt.gz"

#Params parsed from command line to overwrite the hardcoded params for testing:
args = commandArgs(trailingOnly=TRUE)
ColaPutativeBranchStarts = args[1]
#output dir
DirOut = args[2]
pseudocount = as.numeric(args[3])
cola_data_bandwidth = as.numeric(args[4])
#Minimum cola-seq reads per window filter
MinReadsPerWindow = as.numeric(args[5])
#Distance from 3'ss to downstream edge of the window
Window_downstream_dist_to_3ss = as.numeric(args[6])
Window_upstream_dist_to_3ss = as.numeric(args[7])
# Composite score theshold
Threshold = as.numeric(args[8])
# TRUE or FALSE to output extra stuff to help debug
OutputDebugBedgraphs = as.logical(args[9])

dir.create(FigDir <-DirOut, recursive=T)

## get length pmf
lengths <- read.delim(BpToDistanceTrainingData, col.names = c( "IntronType", "DistToBP" ), sep=" ") %>%
      filter(DistToBP<=Window_upstream_dist_to_3ss)

dx<-density(lengths %>% filter(IntronType=="u2") %>% pull(DistToBP), bw=3)
pmf.length <- as.data.frame(approx(dx$x,dx$y,xout=1:Window_upstream_dist_to_3ss)) %>%
  mutate(PositionScore=case_when(
    is.na(y) ~ log2(min(y, na.rm=T)),
    TRUE~ log2(y)
  )) %>%
  mutate(IntronType="u2") %>%
  dplyr::rename(RelPos=x) %>%
  dplyr::select(IntronType,RelPos,PositionScore)


#same for u12 introns
dx.u12<-density(lengths %>% filter(IntronType=="u12") %>% pull(DistToBP),
            bw=3)
pmf.length.u12 <- as.data.frame(approx(dx.u12$x,dx.u12$y,xout=1:Window_upstream_dist_to_3ss)) %>%
  mutate(PositionScore=case_when(
    is.na(y) ~ log2(min(y, na.rm=T)),
    TRUE~ log2(y)
  )) %>%
  mutate(IntronType="u12") %>%
  dplyr::rename(RelPos=x) %>%
  dplyr::select(IntronType,RelPos,PositionScore)


PositionProbabilities <- bind_rows(pmf.length.u12, pmf.length)

## get motif pmf
MotifScores <- fread(MotifProbabilityFunction, check.names = F, header = T) %>%
  mutate(Intron=str_remove(Pos, "\\([+-]\\)")) %>%
  dplyr::select(-Pos) %>%
  dplyr::select(Intron, everything()) %>%
  gather(-Intron, key="RelPos", value="MotifScore") %>%
  separate(Intron, into=c("gene", "chrom", "windowStart", "windowStop", "strand", "IntronType"), sep="_", remove=F, convert=T) %>%
  mutate(RelPos=case_when(
    strand=="+" ~ as.numeric(RelPos),
    strand=="-" ~ as.numeric(RelPos) -1
  )) %>%
  dplyr::select(region=Intron, RelPos, MotifScore)

## Make a smoothed prbability density function from cola-seq read ends
FragStarts <- fread(ColaPutativeBranchStarts, col.names = c("feature", "Pos", "ReadName", "Extra")) %>%
  separate(feature, into=c("gene", "chrom", "windowStart", "windowStop", "strand", "IntronType"), sep="_", remove=F, convert=T) %>%
  separate(Extra, into=c("MismatchAtEnd", "ThreePrimeFragPos"), sep=";", convert=T) %>%
  mutate(ThreePrimeFragPos=ThreePrimeFragPos-1)

#Filter for 3'ss regions with at least 10 reads. Also convert chromosomal coordinate to position relative to 3'ss
FilteredAndAdjusted <- FragStarts %>%
  add_count(feature) %>%
  filter(n>=MinReadsPerWindow) %>%
  mutate(RelativePos=case_when(
    strand=="-" ~ Pos- windowStart + Window_downstream_dist_to_3ss,
    strand=="+" ~ windowStop - Pos + Window_downstream_dist_to_3ss
  )) %>%
  mutate(ThreeSS_Pos=case_when(
    strand=="-" ~ windowStart - Window_downstream_dist_to_3ss,
    strand=="+" ~ windowStop + Window_downstream_dist_to_3ss
  )) %>%
  mutate(Dist3PrimeEndTo3ss=case_when(
    strand=="-" ~ ThreeSS_Pos - ThreePrimeFragPos,
    strand=="+" ~ ThreePrimeFragPos - ThreeSS_Pos +1
  )) %>%
  mutate(
    LariatType=case_when(
      Dist3PrimeEndTo3ss == 0 ~ "ELI",
      TRUE ~ "Putative LI"
    ))

#Num fragments pass filters that can be used to call branchpoints
print ("Num fragments that pass filters that can be used to call bpts:")
nrow(FilteredAndAdjusted)


#Num LI and ELI in wnidows
table(FilteredAndAdjusted$LariatType)
NumELI.LI.ToSample <- min(table(FilteredAndAdjusted$LariatType))


LI_ELI.dist.hist <- FilteredAndAdjusted %>%
  sample_n(1000000, replace=T) %>%
  ggplot(aes(x=RelativePos, fill=LariatType)) +
  geom_histogram(position="dodge") +
  scale_x_continuous(limits=c(0,100)) +
  scale_x_reverse() +
  xlab("Distance from 5' end to 3'ss") +
  ylab("ReadCount") +
  theme_bw()
# ggsave(filename = paste0(FigDir, "DistanceFrom5PrimeEndToSS.ByReadType.pdf"), LI_ELI.dist.hist,  height=3, width=3)

## define a function to  Get probability function from colaseq data
GetColaseqPrbFunction <- function(FilteredAndAdjusted.df, pseudocount=0, bw=1, shift=1, pseudoprobability=2**-70, from=0, to=100){
  #Add a uniform distirbution of counts, with relative weight determined by pseudocount parameter. pseudocount of 1 is equivalent to a literal pseudocount of 1 added at all positions

FilteredAndAdjusted.pseudocount.added <-
  bind_rows(
    #Uniform distirbution
      (FilteredAndAdjusted.df %>%
      dplyr::select(feature, RelativePos) %>%
      expand(feature, RelativePos) %>%
      mutate(weight=pseudocount)),
    #Actual reads
      (FilteredAndAdjusted.df %>%
      dplyr::select(feature, RelativePos) %>%
      mutate(weight=1))
    ) %>%
  #Normalize weight to sum to 1, as required by density function
  group_by(feature) %>%
  mutate(weightSum = sum(weight)) %>%
  ungroup() %>%
  mutate(RelativeWeights=weight/weightSum) %>%
  dplyr::select(feature, RelativePos, RelativeWeights)

#Reformat data to make it easier to iterate through regions and calculate smoothed densities
group_names <- FilteredAndAdjusted.pseudocount.added %>%
  group_keys(feature) %>% pull(1)
ListOfDf <- FilteredAndAdjusted.pseudocount.added %>%
  group_split(feature) %>%
  set_names(group_names)

#I want the probability density function to be shifted from the cola-seq read ends by 1 to account for that cola-seq reads terminate 1 base downstream of branch

#Calculate points along smoothed kernel density for each region.
densities <- lapply(ListOfDf,function(i) {density(i$RelativePos+shift, weights = i$RelativeWeights, bw=bw)})

#Interpolate to get kernel density function points at discrete bases.
density.interpolations <- lapply(densities, function(i) {approx(i$x, i$y, xout=from:to)})

density.interpolations.df <- bind_rows(density.interpolations, .id = 'region')

#Add a very small (insignificant) amount of probability mass to all points to guarantee there are no NA or 0 values, even if there are no cola-seq reads and the density function sums is rounded to 0.

#Add the small proabability, and pivot the interpolated points to a more readable table format
output.df <- density.interpolations.df %>%
  mutate(Probability=case_when(
    is.na(y) ~ pseudoprobability,
    TRUE ~ y + pseudoprobability
  )) %>%
  dplyr::select(region, RelPos=x, Probability) %>%
  mutate(ColaSeqScore=log2(Probability)) %>%
  dplyr::select(-Probability)

return(output.df)
}

#Use function to get probability function
cola.seq.pmf.PerRegion <- GetColaseqPrbFunction(FilteredAndAdjusted, pseudocount=pseudocount, pseudoprobability = 2**-50, bw=cola_data_bandwidth)

#Only consider introns with motif info and that passed cola-seq filter.
IntronsToAnalyze <- intersect((cola.seq.pmf.PerRegion$region %>% unique()),
                              (MotifScores$region %>% unique()))

#Number of introns
N<-length(IntronsToAnalyze)

#Construct df of positional scores (log2 probabilities)
PositionalScores <- bind_rows(data.frame(region=IntronsToAnalyze),
          data.frame(RelPos=5:100)) %>%
  expand(region, RelPos ) %>%
  drop_na() %>%
  separate(region, into=c("gene", "chrom", "windowStart", "windowStop", "strand", "IntronType"), sep="_", remove=F, convert=T) %>%
  dplyr::select(region, RelPos, IntronType) %>%
  left_join(PositionProbabilities, by=c("IntronType", "RelPos")) %>%
  dplyr::select(-IntronType) %>%
  mutate(RelPos=as.numeric(RelPos))

RelPosRange <- (Window_downstream_dist_to_3ss+1):Window_upstream_dist_to_3ss


# Add log2 probabilities
PerPosScores.df <-
  left_join(
    PositionalScores, cola.seq.pmf.PerRegion, by=c("region", "RelPos")) %>%
  left_join(
    MotifScores,by=c("region", "RelPos")) %>%
  # left_join(
  #   (cola.seq.pmf.PerRegion.ELI %>% dplyr::rename(ColaSeq.ELI.only.score=ColaSeqScore)),
  #    by=c("region", "RelPos")) %>%
  # left_join(
  #   (cola.seq.pmf.PerRegion.LI %>% dplyr::rename(ColaSeq.LI.only.score=ColaSeqScore)),
  #    by=c("region", "RelPos")) %>%
  mutate(CompositeScore=PositionScore+MotifScore+ColaSeqScore,
         NoMotifScore=PositionScore+ColaSeqScore,
         NoPositionScore=MotifScore+ColaSeqScore,
         NoColaScore=PositionScore+MotifScore) %>%
  filter(RelPos %in% RelPosRange) %>%
  separate(region, into=c("gene", "chrom", "windowStart", "windowStop", "strand"), sep = '_', remove=F, convert=T) %>%
  mutate(RelPos = as.numeric(RelPos)) %>%
  mutate(start=case_when(
    strand=="+" ~ windowStop + Window_downstream_dist_to_3ss - RelPos,
    strand=="-" ~ windowStart - Window_downstream_dist_to_3ss + RelPos
  )) %>%
  mutate(stop=start+1) %>%
  mutate(start=as.integer(start), stop=as.integer(stop)) %>%
  left_join((
             FilteredAndAdjusted %>%
                 dplyr::select(region=feature, RelPos=RelativePos) %>%
                 count(region, RelPos, name="ColaReadEndCount")
         ), by=c("region", "RelPos")
  )


## Make some plots of random example introns with a true branchpoint
set.seed(10)
IntronsToPlot <- PerPosScores.df %>% distinct(region) %>%
  dplyr::select(region) %>%
  sample_n(6) %>% pull(region)

ScoresPlot <- PerPosScores.df %>%
  filter(region %in% IntronsToPlot) %>%
  mutate(Label=paste0(gene, ":", strand, ":", chrom, ":", windowStart, "-", windowStop)) %>%
  gather(key="ScoreType", value="log2ProbabilityScore", MotifScore, PositionScore, ColaSeqScore, CompositeScore) %>%
  mutate(ScoreType=recode(ScoreType, MotifScore="Sequence motif", PositionScore="Distance to 3'ss", ColaSeqScore="CoLa-seq 5' ends", CompositeScore="Composite")) %>%
  ggplot(aes(x=as.numeric(RelPos), y=as.numeric(log2ProbabilityScore), color=ScoreType)) +
    geom_line() +
    geom_hline(yintercept=Threshold, linetype="dashed") +
    scale_color_manual(values = c("Composite" = "#F8766D","Sequence motif" = "#A3A500", "CoLa-seq 5' ends" = "#00BF7D", "Distance to 3'ss"="#00B0F6") ) +
    # scale_y_continuous(limits=c(-15,5)) +
    # scale_x_reverse() +
    ylab("log2(Rescaled probability)") +
    facet_wrap(~Label, scales="free_x", ncol=3) +
    theme_bw() +
    theme(legend.position="bottom")
ggsave(filename=paste0(FigDir, "BranchpointScorer6RandomExamples.pdf"), ScoresPlot, width=5, height=3)


#Num regions total
print("Number of regions that were considered for BP calling after passing minimum 5' end count threshold")
PerPosScores.df %>% distinct(region) %>%  nrow()

# Filter for positions  above threshold and write out
PerPosScores.df %>%
  filter(CompositeScore > Threshold) %>%
  dplyr::select(chrom, start, stop, region, CompositeScore, strand) %>%
  arrange(chrom, start) %>%
  write_delim(path=paste0(FigDir, "CalledBranchesUnmerged.bed.gz"), col_names=F, delim = '\t')

#Write out full results for all positions per region
PerPosScores.df %>%
  write_delim(path=paste0(FigDir, "BranchpointCallingAll.tab.gz"), col_names=T, delim = '\t')

# Debugging bedgraph files
# If you change the window region or motif scoring function, you will
# definitely want to carefully check that you did everything correctly - that
# ctaAc for example is a high scoring bpt motif. These bedgraph files can help
# debug
if (OutputDebugBedgraphs == T){
    ScoresToOutput <- c("MotifScore", "ColaSeqScore", "PositionScore", "CompositeScore")
    PerPosScores.df %>%
      filter(chrom=="chr21") %>%
      mutate(ColaSeqRaw=2**ColaSeqScore) %>%
      drop_na() %>%
      dplyr::select(c("region", "chrom", "strand", "start", "stop", ScoresToOutput)) %>%
      arrange(chrom, start) %>%
      distinct(chrom, start, .keep_all=T) %>%
      gather(key = 'ScoreType', value="log2ProbabilityScore", ScoresToOutput) %>%
      dplyr::select(chrom, start, stop, log2ProbabilityScore, ScoreType) %>%
      group_by(ScoreType) %>%
      group_walk(~ write_delim(.x, paste0(FigDir, "chr21.debug.", .y$ScoreType, "_Scores.tab.bedgraph"), delim = '\t', col_names = F))
}
