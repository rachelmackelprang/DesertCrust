#Differential abundance analysis use edgeR
#This identifies genes that are significantly different between sample types

#load libraries
library(tidyverse)
library(phyloseq)
library(edgeR)
source("~/Dropbox/PermafrostProjects/Permafrost_Global/Stats_exploration/edgeR.wrapper.r")


meta<-read_tsv(file="~/Dropbox/DesertSoils/DesertSoilMetadata.txt")
#Remove variables we don't need. We dont' need site because those data are captured by sample month. Sites were sampled at different months
meta_trimmed<-select(meta, 'SeqProjectName', 'CollectionMonth', 'SoilType', 'Grouping', 'CrustCategory', 'MossAssociated', 'CrustLayer')

#change anything of type 'chr' to a factor. Make the first column into row names
meta_trimmed<-meta_trimmed %>% as.data.frame() %>% mutate_if(sapply(meta_trimmed, is.character), as.factor) %>% column_to_rownames('SeqProjectName')
SDATA<-sample_data(meta_trimmed)

#import KEGG count data. 
kegg_counts<-read_tsv(file="~/Dropbox/DesertSoils/desert_soil_nomoss_keggcount_all.txt")


#create 'OTU' table using phyloseq
kegg_counts<-kegg_counts %>% as.data.frame() %>% column_to_rownames("KEGG")
OTU<-otu_table(kegg_counts, taxa_are_rows = TRUE)

physeq<-phyloseq(OTU, SDATA)

#prune low abundance genes (keep only those observed more than 100 times in at least 10% of samples)
physeq.f<-filter_taxa(physeq, function(x) sum(x>100) > (0.01 * length(x)), prune=TRUE)

#Need count matrix for edgeR. Extract counts from phyloseq object
crust_counts<-as.data.frame(t(otu_table(physeq.f)))

#It is important to ensure your kegg table and metadata table samples IDs are lined up correctly. 
#Find the overlapping samples and get just the overlapping samples.
common.ids<-intersect(rownames(meta_trimmed), rownames(crust_counts))
kegg_com<-crust_counts[common.ids,]
map_com<-meta_trimmed[common.ids,]

#This tests the effect of crust category (crust vs hypolithic) while controlling for the effect of moss 
CvsH_mosscontrolled_edgeR<-glm.edgeR(x=map_com$CrustCategory, Y=kegg_com, covariates = map_com$MossAssociated)
CvsH_mosscontrolled_edgeR_table<-topTags(CvsH_mosscontrolled_edgeR, n=Inf)
write.table(CvsH_mosscontrolled_edgeR_table, file="crust_vs_hypo_mosscontrolled.csv", sep=",")

#This tests the effect of moss (moss associated vs non moss associated) while controlling for the effect of crust category 
MvsNM_crustcontrolled_edgeR<-glm.edgeR(x=map_com$MossAssociated, Y=kegg_com, covariates=map_com$CrustCategory)
MvsNM_crustcontrolled_edgeR_table<-topTags(MvsNM_crustcontrolled_edgeR, n=Inf)
write.table(MvsNM_crustcontrolled_edgeR_table, file="moss_vs_nomoss_crustcontrolled.csv", sep=",")

topTags(CvsH_mosscontrolled_edgeR)


