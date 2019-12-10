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
physeq.f<-filter_taxa(physeq, function(x) sum(x>100) > (0.01 * length(x)), prune=TRUE)

#function(x) sum(x>100) > (0.1 * length(x))

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
  #adonis(formula = propData_bray ~ CrustLayer + CrustCategory +      MossAssociated + CollectionMonth + SoilType, data = meta_trimmed) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)   
#CrustLayer       2  0.013645 0.0068226  2.9456 0.07992  0.015 * 
#  CrustCategory    1  0.004362 0.0043619  1.8832 0.02555  0.126   
#MossAssociated   1  0.001554 0.0015543  0.6711 0.00910  0.577   
#CollectionMonth  1  0.003460 0.0034599  1.4937 0.02026  0.199   
#SoilType         4  0.031902 0.0079754  3.4432 0.18685  0.002 **
#  Residuals       50  0.115813 0.0023163         0.67832          
#Total           59  0.170736                   1.00000          
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#interactions between collection month and other variables? 
#adonis(propData_bray ~ CollectionMonth*CrustLayer + CollectionMonth*CrustCategory + CollectionMonth * MossAssociated + CollectionMonth*SoilType,data = meta_trimmed)
#Call:
#  adonis(formula = propData_bray ~ CollectionMonth * CrustLayer +      CollectionMonth * CrustCategory + CollectionMonth * MossAssociated +      CollectionMonth * SoilType, data = meta_trimmed) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
#CollectionMonth                 1  0.002923 0.0029234  1.3929 0.01712  0.216    
#CrustLayer                      2  0.014253 0.0071263  3.3954 0.08348  0.013 *  
#  CrustCategory                   1  0.004433 0.0044331  2.1122 0.02596  0.102    
#MossAssociated                  1  0.001412 0.0014123  0.6729 0.00827  0.586    
#SoilType                        4  0.031902 0.0079754  3.8000 0.18685  0.001 ***
#  CollectionMonth:CrustLayer      2  0.002342 0.0011712  0.5581 0.01372  0.683    
#CollectionMonth:CrustCategory   1  0.009617 0.0096170  4.5822 0.05633  0.009 ** 
#  CollectionMonth:MossAssociated  1  0.007119 0.0071190  3.3920 0.04170  0.025 *  
#  CollectionMonth:SoilType        4  0.008586 0.0021464  1.0227 0.05029  0.397    
#Residuals                      42  0.088149 0.0020988         0.51629           
#Total                          59  0.170736                   1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(propData_bray ~ CrustCategory*MossAssociated*CrustLayer, data = meta_trimmed)
#Call:
#  adonis(formula = propData_bray ~ CrustCategory * MossAssociated *      CrustLayer, data = meta_trimmed) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
#CrustCategory                            2  0.015036 0.0075179  3.2104 0.08806  0.014 *  
#  MossAssociated                           1  0.001553 0.0015527  0.6630 0.00909  0.546    
#CrustLayer                               1  0.002973 0.0029731  1.2696 0.01741  0.265    
#CrustCategory:MossAssociated             1  0.018856 0.0188563  8.0522 0.11044  0.001 ***
#  CrustCategory:CrustLayer                 1  0.004315 0.0043154  1.8428 0.02528  0.127    
#MossAssociated:CrustLayer                1  0.006930 0.0069305  2.9595 0.04059  0.034 *  
#  CrustCategory:MossAssociated:CrustLayer  1  0.001642 0.0016420  0.7012 0.00962  0.550    
#Residuals                               51  0.119430 0.0023418         0.69950           
#Total                                   59  0.170736                   1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#ordinations

ord_bray.pcoa<-ordinate(physeq.f.ra, "PCoA", "bray")
ord_bray.nmds<-ordinate(physeq.f.ra, "NMDS", "bray") # no convergance for the NMDS. Don't use these data
p=plot_ordination(physeq.f.ra, ord_bray.pcoa, color="CrustCategory", shape="CollectionMonth")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))

p=plot_ordination(physeq.f.ra, ord_bray.pcoa, color="CrustCategory", shape="CollectionMonth")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2")) #don't use these data. no convergance



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
shannon_crust.anova<-aov(Shannon ~ CrustCategory, alpha_for_anova) # not significant
invsimp_crust.anova<-aov(InvSimpson ~ CrustCategory, alpha_for_anova) #p=0.00942


shannon_moss.anova<-aov(Shannon ~ MossAssociated, alpha_for_anova) # not significant
invsimp_moss.anova<-aov(InvSimpson ~ MossAssociated, alpha_for_anova) # not sig

shannon_layer.anova<-aov(Shannon ~ CrustLayer, alpha_for_anova) #p=5.79e-01
invsimp_layer.anova<-aov(InvSimpson ~ CrustLayer, alpha_for_anova) # p=0.00469

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
tukey.shannon.crust<-TukeyHSD(shannon_crust.anova) #not significan
tukey.InvSimpson.crust<-TukeyHSD(invsimp_crust.anova)#hypol-biological were sig different. p = 0.0073955



