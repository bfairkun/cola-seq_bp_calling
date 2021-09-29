library(tidyverse)
library(data.table)
library(edgeR)

DataIn_f <- "Misc/ColaPutativeBranchFragmentStarts.demultiplexed.tab.gz"
DataIn <- fread(DataIn_f,  col.names=c("Region", "Count", "NearestBranchpointPos", "DistToNearestBranchpoint"), sep='\t') %>%
  separate(Region, into=c("Region", "sample", "FragStart", "FragEnd"), sep=";", convert=T) %>%
  separate(Region, into=c("gene", "chr", "WindowStart", "WindowEnd", "Strand", "IntronType"), sep="_", convert=T, remove = F)  %>%
  mutate(SpliceAcceptor=case_when(
    Strand=="+" ~ WindowEnd + 5,
    Strand=="-" ~ WindowStart -4
  )) %>%
  mutate(ReadType=case_when(
    FragEnd-SpliceAcceptor==0 ~ "ELI",
    TRUE ~ "putativeNLI"
  )) %>%
  arrange(ReadType) %>% #Order with ELI reads first
  distinct(FragStart, FragEnd, sample, .keep_all = T) #When same fragment appears on more than one line because it is in multiple regions, keep the first to avoid double counting

#Check number of distinct closest bpts
DataIn %>%
  distinct(chr, NearestBranchpointPos) %>% nrow()

#Verify that no reads are double counted because they are in multiple regions.
DataIn %>% filter(sample=="YZ112") %>%
  group_by(chr, FragStart, FragEnd) %>%
  summarize(NumRegions = n_distinct(Region)) %>%
  ungroup() %>%
  count(NumRegions)

#Check that there is pileup of ReadEnds at Splice acceptor. Do for both strands separately in case there are silly bugs
DataIn %>%
  mutate(ReadEndToSpliceAcceptorDist=FragEnd-SpliceAcceptor) %>%
  filter(abs(ReadEndToSpliceAcceptorDist)<10) %>%
ggplot(aes(x=ReadEndToSpliceAcceptorDist, fill=ReadType )) +
    geom_bar() +
    facet_wrap(~Strand) +
    theme_bw()

#Check that there is pileup of ReadStarts at 1bp downstream of bpts
DataIn %>%
  filter(abs(DistToNearestBranchpoint)<10) %>%
  ggplot(aes(x=DistToNearestBranchpoint, fill=ReadType)) +
  geom_bar() +
  facet_wrap(~Strand) +
  theme_bw()

#Ok, those sanity checks are good. Now let's make the count table(s), in the format of leafcutter count.numers files
#The "introns" will be branchpoints, and the "clusters" will be branchpoint regions, merging overlapping branchpoint regions
#First read in file of overlapping branchpoint regions
BptRegionClusters <- read_delim("../output/202012_Stuff/BranchpointRegionsMerged.bed.gz", col_names = c("Chrom", "MergedRegionStart", "MergedRegionStop", "Strand", "OriginalRegions"), delim = '\t')
MaxRegionsInMergedRegion <- BptRegionClusters %>%
  mutate(RegionCount = str_count(OriginalRegions,",") +1 ) %>%
  count(RegionCount) %>% pull(RegionCount) %>% max()
BptRegionKey <- BptRegionClusters %>%
  unite(MergedRegion, Chrom, MergedRegionStart, MergedRegionStop, Strand, sep="_") %>%
  separate(OriginalRegions, into = paste0("OriginalRegion", 1:MaxRegionsInMergedRegion), sep = ",", extra="drop", fill="right") %>%
  mutate(Cluster=paste0("clu_", 1:n())) %>%
  gather(key="RegionNumber", value="OriginalRegion", -MergedRegion, -Cluster) %>%
  filter(!is.na(OriginalRegion)) %>%
  dplyr::select(Region=OriginalRegion, MergedRegion, Cluster)
write_delim(BptRegionKey, "../output/202012_Stuff/Bpt_Cluster_Descriptions.txt.gz")

#Next make count table with just NLI reads, meaning filtering for reads that start at 0 or +1 from branchpoint and summing them.
NLI.Bpt.Table <- DataIn %>%
  filter(ReadType=="putativeNLI") %>%
  # head(1000) %>%
  dplyr::select(Region, NearestBranchpointPos, DistToNearestBranchpoint, Count, sample) %>%
  left_join(BptRegionKey, by="Region") %>%
  filter(DistToNearestBranchpoint %in% c(0,1)) %>%
  group_by(NearestBranchpointPos, sample, MergedRegion, Cluster) %>%
  summarise(MergedCount = sum(Count)) %>%
  ungroup() %>%
  spread(key=sample, value=MergedCount, fill=0) %>%
  separate(MergedRegion, into=c("Chr", "x", "y", "Strand"), sep = "_") %>%
  mutate(RowName=paste(Chr, NearestBranchpointPos, NearestBranchpointPos+1, paste0(Cluster,"_",Strand), sep=":")) %>%
  distinct(Chr, NearestBranchpointPos, .keep_all = T) %>%
  dplyr::select(-c("NearestBranchpointPos", "Cluster", "Chr", "Strand", "x", "y")) %>%
  dplyr::select(RowName, everything()) %>%
  arrange(RowName) %>%
  column_to_rownames("RowName")
nrow(NLI.Bpt.Table)
F_out <- "../output/202012_Stuff/BptCountTables.NLI.txt"
write.table(NLI.Bpt.Table, file = F_out, row.names = T, col.names = T, sep='\t', quote=F)
system(paste("gzip -f", F_out))

#Now Count table from all ELI reads, using all ReadStarts as branchpoints, regardless of if near called bpt
ELI.Bpt.Table <- DataIn %>%
  filter(ReadType=="ELI") %>%
  left_join(BptRegionKey, by="Region") %>%
  dplyr::select(MergedRegion, FragStart, Count, sample, ReadType, Cluster) %>%
  group_by(FragStart, sample, MergedRegion, Cluster) %>%
  summarise(MergedCount = sum(Count)) %>%
  ungroup() %>%
  spread(key=sample, value=MergedCount, fill=0) %>%
  separate(MergedRegion, into=c("Chr", "x", "y", "Strand"), sep = "_") %>%
  mutate(RowName=paste(Chr, FragStart, FragStart+1, paste0(Cluster,"_",Strand), sep=":")) %>%
  dplyr::select(-c("FragStart", "Cluster", "Chr", "Strand", "x", "y")) %>%
  dplyr::select(RowName, everything()) %>%
  arrange(RowName) %>%
  column_to_rownames("RowName")
nrow(ELI.Bpt.Table)
F_out <- "../output/202012_Stuff/BptCountTables.ELI.txt"
write.table(ELI.Bpt.Table, file = F_out, row.names = T, col.names = T, sep='\t', quote=F)
system(paste("gzip -f", F_out))

# Now a count table from NLI+ELI reads, with bpt filtering for reads that start at 0 or +1 from bpt
NLI.ELI.Bpt.Table <- DataIn %>%
  dplyr::select(Region, NearestBranchpointPos, DistToNearestBranchpoint, Count, sample) %>%
  left_join(BptRegionKey, by="Region") %>%
  filter(DistToNearestBranchpoint %in% c(0,1)) %>%
  group_by(NearestBranchpointPos, sample, MergedRegion, Cluster) %>%
  summarise(MergedCount = sum(Count)) %>%
  ungroup() %>%
  spread(key=sample, value=MergedCount, fill=0) %>%
  separate(MergedRegion, into=c("Chr", "x", "y", "Strand"), sep = "_") %>%
  mutate(RowName=paste(Chr, NearestBranchpointPos, NearestBranchpointPos+1, paste0(Cluster,"_",Strand), sep=":")) %>%
  distinct(Chr, NearestBranchpointPos, .keep_all = T) %>%
  dplyr::select(-c("NearestBranchpointPos", "Cluster", "Chr", "Strand", "x", "y")) %>%
  dplyr::select(RowName, everything()) %>%
  arrange(RowName) %>%
  column_to_rownames("RowName")
nrow(NLI.ELI.Bpt.Table)
F_out <- "../output/202012_Stuff/BptCountTables.NLI_ELI.txt"
write.table(NLI.ELI.Bpt.Table, file = F_out, row.names = T, col.names = T, sep='\t', quote=F)
system(paste("gzip -f", F_out))
