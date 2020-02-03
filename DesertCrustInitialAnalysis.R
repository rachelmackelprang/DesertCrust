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
p+geom_point(size=3, alpha=0.75) + theme(plot.title=element_text(size=16, face="bold, vjust =2"))
dev.off()



################################################
#####Differential abundance analysis############
################################################

#Need count matrix for edgeR. Extract from phyloseq object
kegg_counts<-as.data.frame(t(otu_table(physeq.f)))

#It is important to ensure your kegg table and metadata table samples IDs are lined up correctly. 
#Find the overlapping samples and get just the overlapping samples.
kcommon.ids<-intersect(rownames(meta_trimmed), rownames(kegg_counts))
kegg_com<-kegg_counts[kcommon.ids,]
kmap_com<-meta_trimmed[kcommon.ids,]

#use a negative binomial distribution to ID differentially abundant genes. 
#Load edgeR wrapper

library(edgeR)
source("~/Dropbox/PermafrostProjects/Permafrost_Global/Stats_exploration/edgeR.wrapper.r")

crust_edgeR<-glm.edgeR(x=kmap_com$CrustCategory, Y=kegg_com)

#get top results
topTags(crust_edgeR)
#logFC      logCPM       LR       PValue          FDR
#K16023 -3.3120311 -0.91161650 53.16360 3.068932e-13 1.785812e-09
#K11081 -2.0085715 -0.98457446 40.47257 1.993957e-10 5.801419e-07
#K15672 -1.8704153  1.58393743 30.97057 2.619702e-08 5.081349e-05
#K14743 -1.0581863  3.13486237 22.25206 2.391022e-06 3.478339e-03
#K16106 -0.9282737  0.51528524 20.94865 4.717606e-06 4.441608e-03
#K14339 -0.6382925  4.21764275 20.84739 4.973684e-06 4.441608e-03
#K13669 -0.9534547  0.92015053 20.71020 5.343059e-06 4.441608e-03
#K03822 -0.7544529  4.85275008 20.41466 6.235031e-06 4.535206e-03
#K12428 -1.2209972  1.12455337 20.16179 7.116028e-06 4.600908e-03
#K00622 -1.0666923 -0.02473427 19.94313 7.977993e-06 4.642394e-03

moss_edgeR<-glm.edgeR(x=kmap_com$MossAssociated, Y=kegg_com)
topTags(moss_edgeR)

#Coefficient:  xY 
#logFC   logCPM       LR       PValue          FDR
#K08966 -2.3692693 3.337318 65.18864 6.806067e-16 3.960451e-12
#K19694 -1.7839014 5.510622 62.05084 3.347042e-15 7.345485e-12
#K13730 -1.6089053 7.480156 61.80766 3.786983e-15 7.345485e-12
#K13472 -1.6443936 3.811362 60.86518 6.112030e-15 8.891476e-12
#K11751 -0.9221005 5.023282 58.70328 1.833355e-14 2.133658e-11
#K19144 -1.9349088 3.030475 57.08036 4.183629e-14 4.057423e-11
#K01481 -1.4284620 4.496845 54.90192 1.266970e-13 1.053214e-10
#K02723 -2.1232992 1.105740 54.15105 1.856544e-13 1.252573e-10
#K12152 -1.9547038 2.296937 54.06739 1.937301e-13 1.252573e-10
#K05941 -1.0397122 3.071106 52.14732 5.148900e-13 2.996145e-10



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


#post-hoc tests. ID pairwise groups that are different. Output gives difference in means, confidence levels, and adjusted p-values for all possible pairs.
tukey.InvSimpson<-TukeyHSD(invsimp_layer.anova)
tukey.shannon<-TukeyHSD(shannon_layer.anova)

write.table(tukey.InvSimpson$CrustLayer, file="tukey.invsimp.crustlayer.csv", sep=",")
write.table(tukey.shannon$CrustLayer, file="tukey.shannon.crustlayer.csv", sep=",")




