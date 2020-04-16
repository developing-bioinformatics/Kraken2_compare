library(ggplot2)
library(dplyr)
library(rBLAST)
library(ShortRead)
library(forcats)
library(taxonomizr)

srr=c('SRR11043497')
system(paste('fastq-dump', srr, sep=' '))
data=readFastq('.', pattern=srr)
reads = sread(data, id=id(dna))
bl <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa")
cl <- predict(bl, reads, BLAST_args = '-num_threads 12 -evalue 1e-100')
taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")
accid = as.character(cl$SubjectID)
ids<-accessionToTaxa(accid,'/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')
taxlist=getTaxonomy(ids, taxaNodes, taxaNames)
cltax=cbind(cl,taxlist)

blasthits=cltax
head(blasthits)


lca = function(x) {
  taxnames = c('superkingdom', 'phylum', 'order', 'family', 'genus', 'species')
  shortnames = apply(x, 2, unique)
  countshnames = sapply(shortnames, length)
  lastuni = tail(names(countshnames[countshnames==1]), n=1)
  nombre = as.data.frame(x[1,which(colnames(x) == lastuni)])
  ret = x %>% 
    mutate(last_common = as.character(nombre[[1]])) %>%
    mutate(level_lca = lastuni)
  return(ret)
}

listhits = blasthits[1:10000,] %>%
  group_by(QueryID) %>%
  group_modify(~ lca(.x))
#splits table into groups and applies lca
unique(listhits$last_common)

allhits = blasthits[1:1000,] %>% #run on first 500K rows
  group_by(QueryID) %>%
  group_modify(~ lca(.x)) %>% #use lca function
  slice(1) %>% #keep only first row in each group (one per read)
  group_by(last_common) %>% # group by LCA taxon
  summarize(lca_count = n()) %>% # count num reads (rows) for each LCA
  arrange(desc(lca_count)) # sort descending
allhits
allhits$last_common

(lca_plot_hw_3_25= ggplot(allhits %>% filter(!is.na(last_common)) ) +#gets rid of NA 
    geom_col(aes(x=last_common, y=lca_count)) + #bar plot but can define height, x is lca nd y is how often appear
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)))

ggsave(lca_plot_hw_3_25, file='lca_plot_hw_3_25.png', height=8, width=8, dpi=700)

#second part of homework
library(multidplyr)

cluster <- new_cluster(8)
cluster_library(cluster, 'dplyr')
cluster_copy(cluster, 'lca')

#run analysis on multiple r sessions to go faster
#split groups across multiple CPU cores
prepare = blasthits %>%
  group_by(QueryID) %>%
  partition(cluster) %>%  #split groups into cluster units
  do({lca(.)}) %>% # Use do to apply function across groups
  collect() %>% #bring it all back together
  slice(1) %>% #keep only first row in each group (one per read)
  group_by(last_common) %>% # group by LCA taxon
  summarize(lca_count = n()) %>% # count num reads (rows) for each LCA
  arrange(desc(lca_count))  %>% # sort descending
  filter(!is.na(last_common)) 
print(prepare)

#get table with lca and number 
#change order of how it shows up in graph so it looks better
#now shows from most to least

(lca_plot = ggplot(prepare)  +
    geom_col(aes(x=fct_reorder(last_common, lca_count, .desc=T), y=lca_count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90))
)
