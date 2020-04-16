library(Biostrings)
library(taxonomizr)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringi)
library(forcats)
library(rBLAST)
library(multidplyr)
library(cowplot)

#devtools::install_github("tidyverse/multidplyr")
#devtools::install_github("mhahsler/rBLAST")

# Download SRA File:
srr=c('SRR11043468')
system(paste('fastq-dump', srr, sep=' '))

#set path to kraken db
kdb='/scratch/vrapela/kdb'
#path to kraken2 executable
kraken2='/projectnb/ct-shbioinf/kraken/kraken2'
report_file='kreport_psba3.tab'
SRR_file=list.files(pattern=srr)
out_file='kout_psba3.txt'
run_kraken=paste(kraken2, '--db', kdb, '--memory-mapping --report', report_file, SRR_file, ">", out_file, sep=' ')
system(run_kraken)


# Read taxonomy database
taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")


# read fastq
dna = readFastq(paste(srr, '.fastq',sep=''))
reads = sread(dna)
qscores = quality(dna) 

# plot readlength
widths = as.data.frame(reads@ranges@width)
(widthplot <- ggplot(widths) +
    geom_histogram(aes(x=reads@ranges@width), binwidth = 10) + 
    theme_linedraw() + 
    xlab('Read Length (bp)') +
    xlim(0,2000) +
    ggtitle('Read length distribution for 550bp amplicon'))
ggsave(widthplot, file='readlengths_psba3.png')

# plot qscores
numqscores = as(qscores, "matrix") # converts to numeric scores automatically
avgscores = apply(numqscores, 1, mean, na.rm=TRUE) #apply the function mean() across 
avgscores = as.data.frame(avgscores)
(qscores = ggplot(avgscores) +
    geom_histogram(aes(x=avgscores), binwidth=0.2) +
    theme_linedraw() +
    xlab('Quality Score') +
    ggtitle('Per Read Average Quality'))
ggsave(qscores, file='quality_psba3.png')


## blast
bl <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa")
cl <- predict(bl, reads, BLAST_args = '-num_threads 8 -evalue 1e-50')
accid = as.character(cl$SubjectID) # accession IDs of BLAST hits

# Plot results

#takes accession number and gets the taxonomic ID
ids<-accessionToTaxa(accid, '/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')
#taxlist displays the taxonomic names from each ID #
taxlist=getTaxonomy(ids, taxaNodes, taxaNames)

cltax=cbind(cl,taxlist)
write.csv(cltax, file='blasthits_psba3.csv')



kreport = read.table('/projectnb/ct-shbioinf/vrapela/Kraken2_compare/kreport_psba3.tab', sep='\t')
blasthits = read.csv('blasthits_psba3.csv')
kout = read.table('/projectnb/ct-shbioinf/vrapela/Kraken2_compare/kout_psba3.txt', sep='\t')

#blast first 

lca2 = function(x) {
  require(dplyr)
  taxnames = c('superkingdom', 'phylum', 'order', 'family', 'genus', 'species')
  x = x %>% filter(!is.na(superkingdom)) %>% filter(superkingdom != 'Bacteria') # need to deal with this more generally
  shortnames = apply(x[,taxnames], 2, unique)
  countshnames = sapply(shortnames, length)
  numcount = countshnames==1
  lastuni = tail(names(shortnames[numcount==T]), n=1)
  nombre = as.data.frame(x[1,which(colnames(x) == lastuni)])
  newtax <- as.list(ifelse(countshnames==1,shortnames,NA))
  
  ret = x %>% 
    mutate(last_common = as.character(nombre[[1]])) %>%
    mutate(level_lca = lastuni) %>%
    mutate(superkingdom = newtax$superkingdom) %>%
    mutate(phylum = newtax$phylum) %>%
    mutate(class = newtax$class) %>%
    mutate(order = newtax$order) %>%
    mutate(family = newtax$family) %>%
    mutate(genus = newtax$genus) %>%
    mutate(species = newtax$species)
  return(ret)
}


cluster <- new_cluster(8)
cluster_library(cluster, 'dplyr')
cluster_copy(cluster, 'lca2')

#split groups across multiple CPU cores
prepare = blasthits %>%
  group_by(QueryID) %>%
  partition(cluster) %>%  #split groups into cluster units
  do({lca2(.)}) %>%
  collect()

#count reads matching species:
spcount = prepare %>%
  group_by(QueryID) %>%
  slice(1) %>%
  group_by(species) %>%
  summarize(count=n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(species))


#count reads matching to genera:
gencount = prepare %>%
  group_by(QueryID) %>%
  slice(1) %>%
  group_by(genus) %>%
  summarize(count=n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(genus))

#count reads matching to families
famcount = prepare %>%
  group_by(QueryID) %>%
  slice(1) %>%
  group_by(family) %>%
  summarize(count=n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(family))


(p1 = ggplot(spcount) +
    geom_col(aes(x=fct_reorder(species, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Species')
)

(p2 = ggplot(gencount) +
    geom_col(aes(x=fct_reorder(genus, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Genus')
)

(p3 = ggplot(famcount) +
    geom_col(aes(x=fct_reorder(family, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Family')
)

gr = plot_grid(p1, p2, p3, ncol=1, nrow=3)
#3 plots in 1 file 

ggsave(gr, file='Blast_grid_taxonomy_psba3.png', height = 9, width = 5, dpi=500)


#now for kraken 


kraken_gencount = kreport %>% filter(V4 == 'G') %>% filter(V2>=5) %>% arrange(desc(V2))
(kpgen1 = ggplot(kraken_gencount) +
    geom_col(aes(x=fct_reorder(V6, V2, .desc=T), y=V2)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Genus') + ylab('Count'))

kraken_famcount = kreport %>% filter(V4 == 'F') %>% filter(V2>=5) %>% arrange(desc(V2))
(kpfam1 = ggplot(kraken_famcount) +
    geom_col(aes(x=fct_reorder(V6, V2, .desc=T), y=V2)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Family') + ylab('Count'))

kraken_spcount = kreport %>% filter(V4 == 'S') %>% filter(V2>=5) %>% arrange(desc(V2))
(kpsp1 = ggplot(kraken_spcount) +
    geom_col(aes(x=fct_reorder(V6, V2, .desc=T), y=V2)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Species') + ylab('Count'))

kgr = plot_grid(kpsp1, kpgen1, kpfam1, ncol=1, nrow=3)

ggsave(kgr, file='Kraken_grid_taxonomy_psba3.png', height = 15, width = 7, dpi=500)

#kraken has alot more than blast
#can search names of results on inaturalist webiste
#is it real, is it in north America


#number of unique reads

length(unique(blasthits$QueryID))
sum(kreport$V3)
kout %>% filter(V1=='C') %>% nrow()

#looks like kraken classifies more than blast can find