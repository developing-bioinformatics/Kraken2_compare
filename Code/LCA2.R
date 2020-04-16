library(dplyr)
library(multidplyr)
library(ggplot2)
library(forcats)
install.packages('multidplyr') #not avaliable for R 3.6.2
devtools::install_github("tidyverse/multidplyr")


#load function
lca = function(x) {
  require(dplyr)
  taxnames = c('superkingdom', 'phylum', 'order', 'family', 'genus', 'species')
  shortnames = apply(x[,taxnames], 2, unique)
  countshnames = sapply(shortnames, length)
  lastuni = tail(names(countshnames[countshnames==1]), n=1)
  nombre = as.data.frame(x[1,which(colnames(x) == lastuni)])
  ret = x %>% 
    mutate(last_common = as.character(nombre[[1]])) %>%
    mutate(level_lca = lastuni)
  return(ret)
}

blasthits=read.table('blasthits2.tab')

#parallel processing with multidplyr
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

ggsave(lca_plot, file='lca_plot_class-3-25.png', height=8, width=8, dpi=700)

#results now have mix of taxonomic levels
#now separate by taxonomic level 
#add hits of species to genus hit as well

lca2 = function(x) {
  require(dplyr)
  taxnames = c('superkingdom', 'phylum', 'order', 'family', 'genus', 'species')
  shortnames = apply(x[,taxnames], 2, unique)
  countshnames = sapply(shortnames, length)
  numcount = countshnames==1 #which tx level only have 1 unique name
  lastuni = tail(names(shortnames[numcount==T]), n=1) # pull out what that unique name is 
  nombre = as.data.frame(x[1,which(colnames(x) == lastuni)])
  newtax <- as.list(ifelse(countshnames==1,shortnames,NA)) #checks if only 1 uniqe name, then keep the name that is in the shortnames, if it is anything else, then assign NA
  
  ret = x %>% 
    mutate(last_common = as.character(nombre[[1]])) %>%
    mutate(level_lca = lastuni) %>%
    mutate(superkingdom = newtax$superkingdom) %>% #recode all tax columns using the names we are confident in from new tax
    mutate(phylum = newtax$phylum) %>%
    mutate(class = newtax$class) %>%
    mutate(order = newtax$order) %>%
    mutate(family = newtax$family) %>%
    mutate(genus = newtax$genus) %>%
    mutate(species = newtax$species)
  return(ret)
}

#test function
test = blasthits[500:1000,] %>% 
  group_by(QueryID) %>%
  group_modify(~ lca2(.))

unique(test[,'species'])
#proves function worked
#only species we got back was calitriches stagnalis

#run new function on cluster
cluster <- new_cluster(8)
cluster_library(cluster, 'dplyr')
cluster_copy(cluster, 'lca2')


#split groups across multiple CPU cores
prepare = blasthits %>%
  group_by(QueryID) %>%
  partition(cluster) %>%  #split groups into cluster units
  do({lca2(.)}) %>%
  collect()

#which species are we confident in and how many match to each
#so on and so forth 

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


#plot for genus and species
(p1 = ggplot(gencount) +
    geom_col(aes(x=fct_reorder(genus, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Genus')
)

ggsave(p1, file='gencount_3-25.png', height=8, width=8, dpi=700)


(p2 = ggplot(famcount) +
    geom_col(aes(x=fct_reorder(family, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Family')
)

ggsave(p2, file='famcount_3-25.png', height=8, width=8, dpi=700)
