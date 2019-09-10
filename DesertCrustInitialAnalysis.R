#Rachel Mackelprang
#Analyze KEGG count data from desert crusts


#load libraries
library(vegan)
library("phyloseq")
library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(phyloseq)

#some ggplot2 theming
theme_set(theme_bw())


#Import metadata
meta<-read_tsv(file="~/Dropbox/DesertSoils/DesertSoilMetadata.txt")

#Remove variables we don't need. We dont' need latitude because those data are captured by sample month. Sites were sampled at different months
meta_trimmed<-select(des_soil_meta, 'Seq Project Name', 'Collection Month', 'Sample Isolated From', 'Latitude', 'Grouping')

#change anything of type 'chr' to a factor. Make the first column into row names
meta_trimmed<-meta_trimmed %>% as.data.frame() %>% mutate_if(sapply(meta_trimmed, is.character), as.factor) %>% column_to_rownames('Seq Project Name')
SDATA<-sample_data(meta_trimmed)

#import KEGG count data. 
kegg_counts<-read_tsv(file="~/Dropbox/DesertSoils/desert_soil_nomoss_keggcount_all.txt")


#create 'OTU' table using phyloseq
kegg_counts<-kegg_counts %>% as.data.frame() %>% column_to_rownames("KEGG")
OTU<-otu_table(kegg_counts, taxa_are_rows = TRUE)

physeq<-phyloseq(OTU, SDATA)

#prune low abundance genes (observed more than 100 times in at least 10% of samples)
physeq.f<-filter_taxa(physeq, function(x) sum(x<=100) > (0.01 * length(x)), prune=TRUE)

#get relative abundance
physeq.f.ra<-transform_sample_counts(physeq.f, function(x) x*100/sum(x))

#extract kegg relative abundance table from phyloseq object and change rows to columns and vice versa
propData<-as.data.frame(t(otu_table(physeq.f.ra)))

#compute alpha diversity stats (note that alpha diversity in phyloseq uses count data rather than relative abundance)
alpha<-estimate_richness(physeq.f,  measures=c("Shannon", "InvSimpson"))

#plot alpha diversity figure
p = plot_richness(physeq.f, x = "Sample.Isolated.From", color = "Collection.Month", measures = c("Shannon", "InvSimpson"))
p + scale_color_discrete(name="Collection Month") + xlab("")


#export alpha diversity figure to pdf
pdf("~/Dropbox/GitHub/DesertCrust/alphadiversity_figure.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "Sample.Isolated.From", color = "Collection.Month", measures = c("Shannon", "InvSimpson"))
p + scale_color_discrete(name="Collection Month") + xlab("") + scale_color_discrete(name="Collection Month") + xlab("")
dev.off()

#make alpha diveristy boxplot figure
p = plot_richness(physeq.f, x = "Sample.Isolated.From", measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()

#export alpha diversity boxplot figure to pdf
pdf("~/Dropbox/GitHub/DesertCrust/alphadiversity_boxplotfigure.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "Sample.Isolated.From", measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()
dev.off()

#anovas and alpha diversity
alpha_for_anova<-cbind(sample_data(physeq.f), alpha)
shannon.anova<-aov(Shannon ~ Sample.Isolated.From, alpha_for_anova)
invsimp.anova<-aov(InvSimpson ~ Sample.Isolated.From, alpha_for_anova)

#summary(shannon.anova)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Sample.Isolated.From  8 0.2546 0.03183   4.597 0.000468 ***
#Residuals            41 0.2839 0.00692                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#summary(invsimp.anova)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#Sample.Isolated.From  8 501065   62633   4.445 0.000616 ***
#  Residuals            41 577733   14091                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


