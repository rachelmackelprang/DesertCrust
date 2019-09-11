#Rachel Mackelprang
#ID KEGG genes differing between samples


#load libraries
library(vegan)
library("phyloseq")
library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(phyloseq)
library(edgeR)
library(biomformat)

#load edgeR wraper function from Dan Knights
source("~/Dropbox/GitHub/DesertCrust/wrap.edgeR.r")

#Import metadata
meta<-read_tsv(file="~/Dropbox/DesertSoils/DesertSoilMetadata.txt")
#Remove variables we don't need. We dont' need site because those data are captured by sample month. Sites were sampled at different months
meta_trimmed<-select(meta, 'SeqProjectName', 'CollectionMonth', 'SoilType', 'Grouping', 'CrustCategory', 'MossAssociated', 'CrustLayer')

#change anything of type 'chr' to a factor. Make the first column into row names
meta_trimmed<-meta_trimmed %>% as.data.frame() %>% mutate_if(sapply(meta_trimmed, is.character), as.factor) %>% column_to_rownames('SeqProjectName')


#import data
kegg_counts<-read_tsv(file="~/Dropbox/DesertSoils/desert_soil_nomoss_keggcount_all.txt")
kegg_counts<-kegg_counts %>% as.data.frame() %>% column_to_rownames("KEGG")
OTU<-otu_table(kegg_counts, taxa_are_rows = TRUE)
physeq<-phyloseq(OTU, SDATA)

#prune low abundance genes (keep only those observed more than 100 times in at least 10% of samples)
physeq.f<-filter_taxa(physeq, function(x) sum(x<=100) > (0.01 * length(x)), prune=TRUE)

#extract kegg counts from phyloseq object and change rows to columns and vice versa
kegg_filtered_t<-as.data.frame(t(otu_table(physeq)))


#edgeR
edgeR_result<-glm.edgeR(x=meta_trimmed$SoilType, Y=kegg_filtered_t)
topTags(edgeR_result)

