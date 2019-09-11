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
physeq.f<-filter_taxa(physeq, function(x) sum(x<=100) > (0.01 * length(x)), prune=TRUE)

#get relative abundance
physeq.f.ra<-transform_sample_counts(physeq.f, function(x) x*100/sum(x))

#extract kegg relative abundance table from phyloseq object and change rows to columns and vice versa
propData<-as.data.frame(t(otu_table(physeq.f.ra)))

#############################
#######Beta diversity########
#############################

propData_bray<-vegdist(propData, "bray")
adonis(propData_bray ~ CollectionMonth, data = meta_trimmed)
#Collection month (ie sample site) not significant


#adonis(propData_bray ~  CrustLayer+CrustCategory + MossAssociated +CollectionMonth + SoilType, data=meta_trimmed)
#Call:
#  adonis(formula = propData_bray ~ CrustLayer + CrustCategory +      MossAssociated + CollectionMonth + SoilType, data = meta_trimmed) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#CrustLayer       2   0.14331 0.071657  2.2497 0.06196  0.047 *  
#  CrustCategory    1   0.04569 0.045691  1.4345 0.01975  0.215    
#MossAssociated   1   0.02447 0.024473  0.7684 0.01058  0.499    
#CollectionMonth  1   0.04688 0.046877  1.4717 0.02027  0.218    
#SoilType         4   0.46002 0.115006  3.6107 0.19889  0.001 ***
#  Residuals       50   1.59256 0.031851         0.68854           
#Total           59   2.31294                  1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#interactions between collection month and other variables? 
adonis(propData_bray ~ CollectionMonth*CrustLayer + CollectionMonth*CrustCategory + CollectionMonth * MossAssociated + CollectionMonth*SoilType,data = meta_trimmed)
#Call:
#  adonis(formula = propData_bray ~ CollectionMonth * CrustLayer +      CollectionMonth * CrustCategory + CollectionMonth * MossAssociated +      CollectionMonth * SoilType, data = meta_trimmed) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#CollectionMonth                 1   0.04697 0.046974  1.7710 0.02031  0.132    
#CrustLayer                      2   0.14476 0.072379  2.7288 0.06259  0.013 *  
#  CrustCategory                   1   0.04500 0.045002  1.6967 0.01946  0.147    
#MossAssociated                  1   0.02362 0.023621  0.8905 0.01021  0.444    
#SoilType                        4   0.46002 0.115006  4.3360 0.19889  0.001 ***
#  CollectionMonth:CrustLayer      2   0.04657 0.023284  0.8778 0.02013  0.482    
#CollectionMonth:CrustCategory   1   0.11987 0.119874  4.5195 0.05183  0.003 ** 
#  CollectionMonth:MossAssociated  1   0.15760 0.157596  5.9417 0.06814  0.001 ***
#  CollectionMonth:SoilType        4   0.15453 0.038632  1.4565 0.06681  0.143    
#Residuals                      42   1.11400 0.026524         0.48164           
#Total                          59   2.31294                  1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(propData_bray ~ CrustCategory*MossAssociated*CrustLayer, data = meta_trimmed)
#Call:
#  adonis(formula = propData_bray ~ CrustCategory * MossAssociated *      CrustLayer, data = meta_trimmed) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#CrustCategory                            2   0.15643 0.078214  2.4336 0.06763  0.026 *  
#  MossAssociated                           1   0.02446 0.024456  0.7610 0.01057  0.535    
#CrustLayer                               1   0.03259 0.032594  1.0141 0.01409  0.366    
#CrustCategory:MossAssociated             1   0.27435 0.274345  8.5361 0.11861  0.001 ***
#  CrustCategory:CrustLayer                 1   0.08805 0.088049  2.7396 0.03807  0.045 *  
#  MossAssociated:CrustLayer                1   0.07998 0.079981  2.4886 0.03458  0.052 .  
#CrustCategory:MossAssociated:CrustLayer  1   0.01799 0.017988  0.5597 0.00778  0.656    
#Residuals                               51   1.63910 0.032139         0.70866           
#Total                                   59   2.31294                  1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#ordinations

ord_bray.pcoa<-ordinate(physeq.f.ra, "PCoA", "bray")
ord_bray.nmds<-ordinate(physeq.f.ra, "NMDS", "bray")
p=plot_ordination(physeq.f.ra, ord_bray.pcoa, color="SoilType", shape="CollectionMonth")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))

p=plot_ordination(physeq.f.ra, ord_bray.pcoa, color="CrustCategory", shape="CollectionMonth")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))



pdf("~/Dropbox/GitHub/DesertCrust/pcoa_month_type.pdf", useDingbats = F)
p=plot_ordination(physeq.f.ra, ord_bray.pcoa, color="SoilType", shape="CollectionMonth")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))
dev.off()

pdf("~/Dropbox/GitHub/DesertCrust/nmds_month_type.pdf", useDingbats = F)
p=plot_ordination(physeq.f.ra, ord_bray.nmds, color="SoilType", shape="CollectionMonth")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))
dev.off()

pdf("~/Dropbox/GitHub/DesertCrust/pcoa_month_crustcategory.pdf", useDingbats = F)
p=plot_ordination(physeq.f.ra, ord_bray.pcoa, color="CrustCategory", shape="CollectionMonth")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))
dev.off()

pdf("~/Dropbox/GitHub/DesertCrust/pcoa_moss_crustcategory.pdf", useDingbats = F)
p=plot_ordination(physeq.f.ra, ord_bray.pcoa, color="CrustCategory", shape="MossAssociated")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))
dev.off()

pdf("~/Dropbox/GitHub/DesertCrust/nmds_moss_crustcategory.pdf", useDingbats = F)
p=plot_ordination(physeq.f.ra, ord_bray.nmds, color="CrustCategory", shape="MossAssociated")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))
dev.off()

#######################
####Alpha Diversity####
#######################

#compute alpha diversity stats (note that alpha diversity in phyloseq uses count data rather than relative abundance)
alpha<-estimate_richness(physeq.f,  measures=c("Shannon", "InvSimpson"))

#plot alpha diversity figure (by soil type)
p = plot_richness(physeq.f, x = "SoilType", color = "CollectionMonth", measures = c("Shannon", "InvSimpson"))
p + scale_color_discrete(name="Collection Month") + xlab("")


#export alpha diversity figure to pdf (by soil type)
pdf("~/Dropbox/GitHub/DesertCrust/alphadiversity_figure.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "SoilType", color = "CollectionMonth", measures = c("Shannon", "InvSimpson"))
p + scale_color_discrete(name="Collection Month") + xlab("") + scale_color_discrete(name="Collection Month") + xlab("")
dev.off()

#make alpha diveristy boxplot figure (by soil type)
p = plot_richness(physeq.f, x = "SoilType", measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()

#export alpha diversity boxplot figure to pdf (by soil type)
pdf("~/Dropbox/GitHub/DesertCrust/alphadiversity_boxplotfigure.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "SoilType", measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()
dev.off()

#anovas and alpha diversity (by soil type)
alpha_for_anova<-cbind(sample_data(physeq.f), alpha)
shannon.anova<-aov(Shannon ~ SoilType, alpha_for_anova)
invsimp.anova<-aov(InvSimpson ~ SoilType, alpha_for_anova)

summary(shannon.anova)
#Df Sum Sq Mean Sq F value Pr(>F)  
#SoilType     8 0.5382 0.06728   2.719 0.0142 *
#  Residuals   51 1.2619 0.02474                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(invsimp.anova)
#Df Sum Sq Mean Sq F value Pr(>F)  
#SoilType     8 0.5382 0.06728   2.719 0.0142 *
#  Residuals   51 1.2619 0.02474                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#post-hoc tests (by soil type). ID pairwise groups that are different. Output gives difference in means, confidence levels, and adjusted p-values for all possible pairs.
tukey.shannon<-TukeyHSD(shannon.anova)
tukey.InvSimpson<-TukeyHSD(invsimp.anova)

write.table(tukey.shannon$SoilType, file="tukey.shannon.csv", sep=",")
write.table(tukey.InvSimpson$SoilType, file="tukey.invsimp.csv", sep=",")


#alpha diversity by crust category, moss associated, and crust layer
shannon_crust.anova<-aov(Shannon ~ CrustCategory, alpha_for_anova) #p = 0.00425
invsimp_crust.anova<-aov(InvSimpson ~ CrustCategory, alpha_for_anova) #p=0.00397


shannon_moss.anova<-aov(Shannon ~ MossAssociated, alpha_for_anova) # not significant
invsimp_moss.anova<-aov(InvSimpson ~ MossAssociated, alpha_for_anova) # not sig

shannon_layer.anova<-aov(Shannon ~ CrustLayer, alpha_for_anova) #not sig
invsimp_layer.anova<-aov(InvSimpson ~ CrustLayer, alpha_for_anova) # not significant

pdf("~/Dropbox/GitHub/DesertCrust/alphadiversity_crustcat.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "CrustCategory",  measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()
dev.off()

pdf("~/Dropbox/GitHub/DesertCrust/alphadiversity_moss.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "MossAssociated",  measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()
dev.off()

pdf("~/Dropbox/GitHub/DesertCrust/alphadiversity_layer.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "CrustLayer",  measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()
dev.off()

#posthoc tests for crust category only b/c it was the only one of the three that was sig
tukey.shannon.crust<-TukeyHSD(shannon_crust.anova) #hypol-biological were sig different. p=0.0044981
tukey.InvSimpson.crust<-TukeyHSD(invsimp_crust.anova)#hypol-biological were sig different. p = 0.0053616



