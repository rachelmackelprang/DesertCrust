#Rachel Mackelprang
#Analyze KEGG count data from desert crusts


#load libraries
library(vegan)
library("phyloseq")
library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(phyloseq)
library(biomformat)


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

#get relative abundance
physeq.f.ra<-transform_sample_counts(physeq.f, function(x) x*100/sum(x))



#############################
#######Diversity Stats########
#############################

#extract kegg relative abundance table from phyloseq object and change rows to columns and vice versa. Get bray-curtis dissimilarities and run a PERMANOVA. 
propData<-as.data.frame(t(otu_table(physeq.f.ra)))
propData_bray<-vegdist(propData, "bray")

adonis(formula = propData_bray ~ CollectionMonth*CrustCategory*MossAssociated*CrustLayer, data=meta_trimmed)

#Call:
#adonis(formula = propData_bray ~ CollectionMonth * CrustCategory *      MossAssociated * CrustLayer, data = meta_trimmed) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
#CollectionMonth                                          1  0.002923 0.0029234  1.3929 0.01712  0.206    
#CrustCategory                                            2  0.016593 0.0082963  3.9529 0.09718  0.004 ** 
#MossAssociated                                           1  0.001398 0.0013983  0.6662 0.00819  0.565    
#CrustLayer                                               1  0.002107 0.0021070  1.0039 0.01234  0.404    
#CollectionMonth:CrustCategory                            2  0.011309 0.0056543  2.6941 0.06623  0.032 *  
#CollectionMonth:MossAssociated                           1  0.002392 0.0023915  1.1395 0.01401  0.315    
#CrustCategory:MossAssociated                             1  0.019497 0.0194970  9.2897 0.11419  0.001 ***
# CollectionMonth:CrustLayer                               1  0.001743 0.0017430  0.8305 0.01021  0.461    
#CrustCategory:CrustLayer                                 1  0.004917 0.0049166  2.3426 0.02880  0.079 .  
#MossAssociated:CrustLayer                                1  0.009472 0.0094717  4.5130 0.05548  0.015 *  
#CollectionMonth:CrustCategory:MossAssociated             1  0.000325 0.0003251  0.1549 0.00190  0.980    
#CollectionMonth:CrustCategory:CrustLayer                 1  0.002562 0.0025623  1.2209 0.01501  0.277    
#CollectionMonth:MossAssociated:CrustLayer                1  0.002885 0.0028852  1.3747 0.01690  0.210    
#CrustCategory:MossAssociated:CrustLayer                  1  0.002000 0.0019995  0.9527 0.01171  0.372    
#CollectionMonth:CrustCategory:MossAssociated:CrustLayer  1  0.002465 0.0024652  1.1746 0.01444  0.299    
#Residuals                                               42  0.088149 0.0020988         0.51629           
#Total                                                   59  0.170736                   1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




#ordinations

ord_bray.pcoa<-ordinate(physeq.f.ra, "PCoA", "bray")
ord_bray.nmds<-ordinate(physeq.f.ra, "NMDS", "bray", trymax=200)

pdf("~/Dropbox/GitHub/DesertCrust/pcoa_month_crustcategory.pdf", useDingbats = F)
p=plot_ordination(physeq.f.ra, ord_bray.pcoa, color="CrustCategory", shape="CollectionMonth")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))
dev.off()

pdf("~/Dropbox/GitHub/DesertCrust/nmds_month_crustcategory.pdf", useDingbats = F)
p=plot_ordination(physeq.f.ra, ord_bray.nmds, color="CrustCategory", shape="CollectionMonth")
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


#subset hypolithic
hypo_physeq<-subset_samples(physeq.f.ra, CrustCategory == "Hypolithic")
hypo_ord_bray_pcoa<-ordinate(hypo_physeq, "PCoA", "bray")
pdf("~/Dropbox/GitHub/DesertCrust/pcoa_hypoOnly_moss.pdf", useDingbats = F)
p=plot_ordination(hypo_physeq, hypo_ord_bray_pcoa, color="MossAssociated")
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))
dev.off()

#subset biological
bio_physeq<-subset_samples(physeq.f.ra, CrustCategory == "Biological")
bio_ord_bray_pcoa<-ordinate(bio_physeq, "PCoA", "bray")
pdf("~/Dropbox/GitHub/DesertCrust/pcoa_bioOnly_moss.pdf", useDingbats = F)
p=plot_ordination(bio_physeq, bio_ord_bray_pcoa, color="MossAssociated")
dev.off()



#######################
####Alpha Diversity####
#######################

#compute alpha diversity stats (note that alpha diversity in phyloseq uses count data rather than relative abundance)
alpha<-estimate_richness(physeq.f,  measures=c("Shannon", "InvSimpson"))

#plot alpha diversity (by moss associated)
pdf(file="MossAssoc_alpha_dotplot.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "MossAssociated", color = "CrustCategory", measures = c("Shannon", "InvSimpson"))
p + scale_color_discrete(name="Crust Category") + xlab("")
dev.off()

#make alpha diveristy boxplot figure (by moss associated)
pdf(file="MossAssoc_alpha_boxplot.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "MossAssociated", measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()
dev.off()

#plot alpha diversity (by crust category)
pdf(file="CrustCat_alpha_dotplot.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "CrustCategory", color = "MossAssociated", measures = c("Shannon", "InvSimpson"))
p + scale_color_discrete(name="Moss Associated") + xlab("")
dev.off()

#make alpha diveristy boxplot figure (by crust category)
pdf(file="CrustCat_alpha_boxplot.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "CrustCategory", measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()
dev.off()

#plot alpha diversity (by crust layer)
pdf(file="Crustlayer_alpha_dotplot.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "CrustLayer", color = "CrustCategory", measures = c("Shannon", "InvSimpson"))
p + scale_color_discrete(name="Crust Category") + xlab("")
dev.off()

#make alpha diveristy boxplot figure (by crust layer)
pdf(file="CrustLayer_alpha_boxplot.pdf", useDingbats = F)
p = plot_richness(physeq.f, x = "CrustLayer", measures = c("Shannon", "InvSimpson"))
p + xlab("") + geom_boxplot()
dev.off()

#anovas and alpha diversity (by CrustCategory)
alpha_for_anova<-cbind(sample_data(physeq.f), alpha)
shannon.anova<-aov(Shannon ~ CrustCategory, alpha_for_anova)
invsimp.anova<-aov(InvSimpson ~ CrustCategory, alpha_for_anova)

summary(shannon.anova)
#Df  Sum Sq  Mean Sq F value Pr(>F)
#CrustCategory  2 0.00322 0.001612   2.041  0.139
#Residuals     57 0.04503 0.000790        

summary(invsimp.anova)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#CrustCategory  2  17021    8510   5.069 0.00942 **
#Residuals     57  95701    1679                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#post-hoc tests. ID pairwise groups that are different. Output gives difference in means, confidence levels, and adjusted p-values for all possible pairs.
#tukey.shannon<-TukeyHSD(shannon.anova)
tukey.InvSimpson<-TukeyHSD(invsimp.anova)

#write.table(tukey.shannon$CrustCategory, file="tukey.shannon.csv", sep=",")
write.table(tukey.InvSimpson$CrustCategory, file="tukey.invsimp.crustcat.csv", sep=",")

shannon_moss.anova<-aov(Shannon ~ MossAssociated, alpha_for_anova) # not significant
invsimp_moss.anova<-aov(InvSimpson ~ MossAssociated, alpha_for_anova) #not significant

shannon_layer.anova<-aov(Shannon ~ CrustLayer, alpha_for_anova) #p=5.79e-01
invsimp_layer.anova<-aov(InvSimpson ~ CrustLayer, alpha_for_anova) # p=0.00469

shannon_type.anova<-aov(Shannon~SoilType, alpha_for_anova) #p=7.29e-10
invsimp_type.anova<-aov(InvSimpson ~ SoilType, alpha_for_anova) #p=1.29e-09


#post-hoc tests. ID pairwise groups that are different. Output gives difference in means, confidence levels, and adjusted p-values for all possible pairs.
tukey.InvSimpson<-TukeyHSD(invsimp_layer.anova)
tukey.shannon<-TukeyHSD(shannon_layer.anova)
tukey.type.InvSimpson<-TukeyHSD(invsimp_type.anova)
tukey.type.shannon<-TukeyHSD(shannon_type.anova)

write.table(tukey.InvSimpson$CrustLayer, file="tukey.invsimp.crustlayer.csv", sep=",")
write.table(tukey.shannon$CrustLayer, file="tukey.shannon.crustlayer.csv", sep=",")
write.table(tukey.type.InvSimpson$SoilType, file="tukey.type.invsimp.csv", sep=",")
write.table(tukey.type.shannon$SoilType, file="tukey.type.shannon.csv", sep=",")



