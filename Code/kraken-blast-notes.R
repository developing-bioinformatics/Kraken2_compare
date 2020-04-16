kreport = read.table('/projectnb/ct-shbioinf/vrapela/Kraken2_compare/kreport4.tab', sep='\t')
# then dplyr filter for column 4 == 'G'
# then ggplot

#in kout, filter out the U,keep only classified
#third column is taxon ID - can compare against blast
#5643:12  means 12 kmer matched to taxon number 5643, takes into account LCA while blast doesnt
#can compare read numbers from kout and blasthits and see if they match to same LCA
#fourth column, length of sequence

#can run this on blasthits to get one row per read and can then compare lca name with kraken
spcount = prepare %>%
  group_by(QueryID) %>%
  slice(1)

#or run lca 2, filter for things that map and compare to kout and see which is classifying more


#number of unique reads

length(unique(blasthits$QueryID))
sum(kreport$V3)

#can use this to get stuff directly from repository
#this cas core file needed for downloadsra for amplicons
source('raw.githubusercontent.com/developing-bioinformatics/eDNA_BLAST/master/R/core.R')
download_sra(amplicon='ITS2', sample = 'control', dir.out = 'data') #going to data folder
