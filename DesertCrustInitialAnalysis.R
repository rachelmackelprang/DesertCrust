#Rachel Mackelprang
#Analyze KEGG count data from desert crusts


#load libraries
library(vegan)
library("phyloseq")
library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(phyloseq)

#Import metadata
meta<-read_tsv(file="~/Dropbox/DesertSoils/DesertSoilMetadata.txt")

#Remove variables we don't need. We dont' need latitude because those data are captured by sample month. Sites were sampled at different months
meta_trimmed<-select(des_soil_meta, 'Seq Project Name', 'Collection Month', 'Sample Isolated From', 'Latitude', 'Grouping')

#change anything of type 'chr' to a factor. Make the first column into row names
meta_trimmed<-meta_trimmed %>% as.data.frame() %>% mutate_if(sapply(meta_trimmed, is.character), as.factor) %>% column_to_rownames('Seq Project Name')

#import KEGG data. Counts and proportions
kegg_counts<-read_tsv(file="~/Dropbox/DesertSoils/desert_soil_nomoss_keggcount_all.txt")
#kegg_counts<-read.table(file="~/Dropbox/DesertSoils/desert_soil_nomoss_keggcount_all.txt", header=TRUE, row.names=1, check.names=F)
#des_kegg_prop<-read.table(file="~/Dropbox/DesertSoils/desert_soil_nomoss_keggprop_all.txt", header=TRUE, row.names=1, check.names=F)


#create 'OTU' table using phyloseq
kegg_counts<-kegg_counts %>% as.data.frame() %>% column_to_rownames("KEGG")
kegg_phyloseq<-otu_table(kegg_counts, taxa_are_rows = TRUE)

physeq<-phyloseq(kegg_phyloseq, meta_trimmed)

#prune low abundance genes (observed more than 100 times in at least 10% of samples)
physeq.f<-filter_taxa(physeq, function(x) sum(x<=100) > (0.01 * length(x)), prune=TRUE)


#extract kegg table from phyloseq object and change rows to columns and vice versa
propData<-as.data.frame(t(otu_table(physeq.f.ra)))
