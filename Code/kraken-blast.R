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
srr=c('SRR11043488')
system(paste('fastq-dump', srr, sep=' '))

#set path to kraken db
kdb='/scratch/vrapela/kdb'
#path to kraken2 executable
kraken2='/projectnb/ct-shbioinf/kraken/kraken2'
report_file='kreport4.tab'
SRR_file=list.files(pattern=srr)
out_file='kout4.txt'
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
ggsave(widthplot, file='readlengths.png')

# plot qscores
numqscores = as(qscores, "matrix") # converts to numeric scores automatically
avgscores = apply(numqscores, 1, mean, na.rm=TRUE) #apply the function mean() across 
avgscores = as.data.frame(avgscores)
(qscores = ggplot(avgscores) +
    geom_histogram(aes(x=avgscores), binwidth=0.2) +
    theme_linedraw() +
    xlab('Quality Score') +
    ggtitle('Per Read Average Quality'))
ggsave(qscores, file='quality.png')


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
write.csv(cltax, file='blasthits.csv')
