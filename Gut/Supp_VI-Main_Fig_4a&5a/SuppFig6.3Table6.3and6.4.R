library('vegan')
library('phyloseq')
library('stringr')
library('microshades')
library('speedyseq')
library('forcats')
library('cowplot')
library('knitr')
library('ggplot2')
library('PERMANOVA')
library('pairwiseAdonis')
library('ggsignif')
library('tidyverse')
library('ggpubr')
library('betapart')
library('stats')
library('hypervolume')
library('ape')
library("ggplot2")
library("tidyverse")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(1)

#Read in CSV file for OTU table (table.csv for ASVs, tablelevel6.csv for genera)
data <- import_biom('HybridLizardGutsGenus.biom')

#Find OTU table in biom file
OTU_biom<-otu_table(data)

#Find the sample names in the biom file
biom_label<-colnames(OTU_biom)

#Read in the readcounts in each sample (To get this, open the txt file and then save the file as a csv)
read_processing<-read.csv('Read_Processing_Summary.csv', row.names=1,header= TRUE)

#Find the sample names in the readcount file (in this case they were numbers so I converted them to strings to match the biom and csv sample labels)
read_summary_label<-as.character(read_processing$customer_label)

#Find the reads per sample in the readcount file
read_summary_readcount<-read_processing$seqs.after_size_filtration

#For every column in the biom OTU table...
for (k in 1:length(biom_label)){
  #...find the corresponding label in the readcount file...
  find_read_summary_row<-which(read_summary_label==biom_label[k])
  #...find the corresponding number of reads from that sample, multiply the column by the total readcount and round to the nearest integer
  OTU_biom[,k]<-round(OTU_biom[,k]*read_summary_readcount[find_read_summary_row])
}

#Remove the non-bacteria/archaea (make sure the None;Other;etc. row is the FIRST row in the table... it usually is)
badTaxa = c(rownames(OTU_biom)[1])
allTaxa = taxa_names(OTU_biom)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
OTU_biom = prune_taxa(allTaxa, OTU_biom)

#Remove samples with reads below a threshold (here it is 10000) none are lower than 10000 in this dataset
#OTU_biom<-OTU_biom[,-which(colSums(OTU_biom)<10000)]

#Create a taxonomy table (the zymo biom file doesn't seem to include this...)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
#Extract the taxonomic name at each taxonomic level
for (k in 1:length(OTU_biom[,1])){
  temp<-strsplit(rownames(OTU_biom)[k],split=';')
  domain<-c(domain,substr(temp[[1]][1],4,50))
  if (temp[[1]][2]=="p__NA"){temp[[1]][2]<-paste0('unclassified_',substr(temp[[1]][1],4,50))}
  else {temp[[1]][2]<-substr(temp[[1]][2],4,50)}
  phylum<-c(phylum,temp[[1]][2])
  if (temp[[1]][3]=="c__NA"){temp[[1]][3]<-paste0('unclassified_',str_replace(temp[[1]][2],'unclassified_',''))}
  else {temp[[1]][3]<-substr(temp[[1]][3],4,50)}
  class<-c(class,temp[[1]][3])
  if (temp[[1]][4]=="o__NA"){temp[[1]][4]<-paste0('unclassified_',str_replace(temp[[1]][3],'unclassified_',''))}
  else {temp[[1]][4]<-substr(temp[[1]][4],4,50)}
  order<-c(order,temp[[1]][4])
  if (temp[[1]][5]=="f__NA"){temp[[1]][5]<-paste0('unclassified_',str_replace(temp[[1]][4],'unclassified_',''))}
  else {temp[[1]][5]<-substr(temp[[1]][5],4,50)}
  family<-c(family,temp[[1]][5])
  if (temp[[1]][6]=="g__NA"){temp[[1]][6]<-paste0('unclassified_',str_replace(temp[[1]][5],'unclassified_',''))}
  else {temp[[1]][6]<-substr(temp[[1]][6],4,50)}
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}


#bind the taxonomic level names together into a table/matrix
taxmat = cbind(domain,phylum,class,order,family,genus,species)
#The rownames are the rownames from the OTU_biom table
rownames(taxmat) <- rownames(OTU_biom)
#The column names are the taxonomic levels
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
#Convert the taxonomy table to the format required for a phyloseq object
TAX<-tax_table(taxmat)

#Read in metadata
meta_csv<-read.csv('Whiptail_Gut_Metadata.csv', row.names=1,header= TRUE)

#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)


#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq<-phyloseq(OTU_biom,TAX,SAMP)

#Rarefy the OTU table
rZOTU<-rarefy_even_depth(Z_physeq, rngseed=5) #set.seed(5)
#careful! this resets the external seed

type<-sample_data(rZOTU)$speciesxsite

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate a variety of different BETA diversity metrics using the VEGAN package

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
euclidean<-vegdist(t(otu_table(rZOTU)),method='euclidean',upper=TRUE,diag=TRUE,binary=FALSE)
jaccard<-vegdist(t(otu_table(rZOTU)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(rZOTU)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform 2D PCA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Run PCA analysis
pca <- prcomp(t(otu_table(rZOTU)))

#Find variation explained by axes
variance <- (pca$sdev)^2
varPercent <- variance/sum(variance) * 100

#Find loadings on axes
loadings <- pca$rotation

#Find coordinates for plotting (identical to PCoA on Euclidean distances)
scores <- pca$x
scores_2<-scores[,1:2]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform 2D PCoA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Perform PcoA
#pcoa_euclidean<-pcoa(euclidean)$vectors[,1:2]
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]

#Perform nmds
nmds_euclidean<-metaMDS(euclidean,dimensions=2)
nmds_jaccard<-metaMDS(jaccard, dimensions=2)
nmds_bray<-metaMDS(bray,dimensions=2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scores_2 <- data.frame(scores_2) #convert to dataframe for ggplot2
scores_2$spec_site <- as.factor(type) #pull grouping variable from existing vector
pca_hull <- scores_2 %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(PC1, PC2))
pca_centroid <- scores_2 %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_PC1 = mean(PC1),    
  mean_PC2 = mean(PC2),)

#plot PCA
PCA <- ggplot() + geom_point(data=scores_2, aes(x=PC1, y=PC2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=pca_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PC1") + ylab("PC2") + geom_polygon(data = pca_hull, alpha = 0.2, aes(x=PC1, y=PC2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('A.') + annotate("text", x=8000, y=-5000, label= "Euclidean", size=13)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the Jaccard PCoA (equivalent to PCA)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcoa_jaccard <- data.frame(pcoa_jaccard) #convert to dataframe for ggplot2
pcoa_jaccard$spec_site <- as.factor(type) #pull grouping variable from existing vector
pcoa_jacc_hull <- pcoa_jaccard %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(Axis.1, Axis.2))
pcoa_jacc_centroid <- pcoa_jaccard %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_PC1 = mean(Axis.1),    
  mean_PC2 = mean(Axis.2),)

#plot Jaccard PCoA
PCoA_jacc <- ggplot() + geom_point(data=pcoa_jaccard, aes(x=Axis.1, y=Axis.2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=pcoa_jacc_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PCoA1") + ylab("PCoA2") + geom_polygon(data = pca_hull, alpha = 0.2, aes(Axis.1, y=Axis.2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) + 
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('B.') + annotate("text", x=0.15, y=-0.31, label= "Jaccard", size=13)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the Bray PCoA (equivalent to PCA)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcoa_bray <- data.frame(pcoa_bray) #convert to dataframe for ggplot2
pcoa_bray$spec_site <- as.factor(type) #pull grouping variable from existing vector
pcoa_bray_hull <- pcoa_bray %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(Axis.1, Axis.2))
pcoa_bray_centroid <- pcoa_bray %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_PC1 = mean(Axis.1),    
  mean_PC2 = mean(Axis.2),)

#plot Bray-Curtis PCoA
PCoA_bray <- ggplot() + geom_point(data=pcoa_bray, aes(x=Axis.1, y=Axis.2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=pcoa_bray_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PCoA1") + ylab("PCoA2") + geom_polygon(data = pca_hull, alpha = 0.2, aes(x=Axis.1, y=Axis.2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) + 
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('C.') + annotate("text", x=0.40, y=-0.5, label= "Bray-Curtis", size=13)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Plot Euclidean NMDS

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmds_euclidean_scores <- data.frame(scores(nmds_euclidean, display="sites")) #convert to NMDS scores to dataframe for ggplot2
nmds_euclidean_scores$spec_site <- as.factor(type) #pull grouping variable from existing vector
euclidean_hull <- nmds_euclidean_scores %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(NMDS1, NMDS2))
euclidean_centroid <- nmds_euclidean_scores %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_NMDS1 = mean(NMDS1),    
  mean_NMDS2 = mean(NMDS2),)

#plot Euclidean NMDS
NMDS_euclidean <- ggplot() + geom_point(data=nmds_euclidean_scores, aes(x=NMDS1, y=NMDS2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=euclidean_centroid, aes(x=mean_NMDS1,y=mean_NMDS2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("NMDS1") + ylab("NMDS2") + geom_polygon(data = euclidean_hull, alpha = 0.2, aes(x=NMDS1, y=NMDS2,  fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('D.') + annotate("text", x=5000, y=-8000, label= "Euclidean", size=13)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Plot Jaccard NMDS

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmds_jaccard_scores <- data.frame(scores(nmds_jaccard, display="sites")) #convert to NMDS scores to dataframe for ggplot2
nmds_jaccard_scores$spec_site <- as.factor(type) #pull grouping variable from existing vector
jaccard_hull <- nmds_jaccard_scores %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(NMDS1, NMDS2))
jaccard_centroid <- nmds_jaccard_scores %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_NMDS1 = mean(NMDS1),    
  mean_NMDS2 = mean(NMDS2),)

#plot Jaccard NMDS
NMDS_jaccard <- ggplot() + geom_point(data=nmds_jaccard_scores, aes(x=NMDS1, y=NMDS2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=jaccard_centroid, aes(x=mean_NMDS1,y=mean_NMDS2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("NMDS1") + ylab("NMDS2") + geom_polygon(data = jaccard_hull, alpha = 0.2, aes(x=NMDS1, y=NMDS2,  fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('E.') + annotate("text", x=0.25, y=-.4, label= "Jaccard", size=13)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Plot Bray-Curtis NMDS

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmds_bray_scores <- data.frame(scores(nmds_bray, display="sites")) #convert to NMDS scores to dataframe for ggplot2
nmds_bray_scores$spec_site <- as.factor(type) #pull grouping variable from existing vector
bray_hull <- nmds_bray_scores %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(NMDS1, NMDS2))
bray_centroid <- nmds_bray_scores %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_NMDS1 = mean(NMDS1),    
  mean_NMDS2 = mean(NMDS2),)

#plot Bray-Curtis NMDS
NMDS_bray <- ggplot() + geom_point(data=nmds_bray_scores, aes(x=NMDS1, y=NMDS2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=bray_centroid, aes(x=mean_NMDS1,y=mean_NMDS2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("NMDS1") + ylab("NMDS2") + geom_polygon(data = bray_hull, alpha = 0.2, aes(x=NMDS1, y=NMDS2,  fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('F.') + annotate("text", x=-.4, y=.35, label= "Bray-Curtis", size=13)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#make panel plot (Fig. 6.3)

library("gridExtra")
library("grid")

panel_list <- list(PCA, PCoA_jacc, PCoA_bray, NMDS_euclidean, NMDS_jaccard, NMDS_bray)

panel_layout <- rbind(c(1,2,3),
                      c(4,5,6))

grid_panel <- grid.arrange(grobs = panel_list, layout_matrix = panel_layout)

grid.draw(grid_panel) # interactive device

ggsave("Figure_SI6.3.png", grid_panel, height = 20, width = 30)


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Table 6.3 information

#Percent Variation/Cumulative Percent Variation for each PCoA
pv_euclidean<-pcoa(euclidean)$values[,2]*100
cs_euclidean<-cumsum(pv_euclidean)
pv_jaccard<-pcoa(jaccard)$values[,2]*100
cs_jaccard<-cumsum(pv_jaccard)
pv_bray<-pcoa(bray)$values[,2]*100
cs_bray<-cumsum(pv_bray)

round2 <- function(x) { #rounding function to two decimal places
  format(round((cumsum(x)),digits=2),nsmall=2)
}

#create table 6.3
tbl_6.3 <- data.frame(
  "metric" = c("Euclidean", "Jaccard", "Bray-Curtis"),
  "percent1" = c(round2(pv_euclidean[1]), round2(pv_jaccard[1]), round2(pv_bray[1])),
  "cumulative1" = c(round2(cs_euclidean[1]), round2(cs_jaccard[1]), round2(cs_bray[1])),
  "percent2" = c(round2(pv_euclidean[2]), round2(pv_jaccard[2]), round2(pv_bray[2])),
  "cumulative2" = c(round2(cs_euclidean[2]), round2(cs_jaccard[2]), round2(cs_bray[2])),
  "percent3" = c(round2(pv_euclidean[3]), round2(pv_jaccard[3]), round2(pv_bray[3])),
  "cumulative3" = c(round2(cs_euclidean[3]), round2(cs_jaccard[3]), round2(cs_bray[3])),
  "percent4" = c(round2(pv_euclidean[4]), round2(pv_jaccard[4]), round2(pv_bray[4])),
  "cumulative4" = c(round2(cs_euclidean[4]), round2(cs_jaccard[4]), round2(cs_bray[4])),
  "percent5" = c(round2(pv_euclidean[5]), round2(pv_jaccard[5]), round2(pv_bray[5])),
  "cumulative5" = c(round2(cs_euclidean[5]), round2(cs_jaccard[5]), round2(cs_bray[5])))

write.csv(tbl_6.3,'Table_6.3.csv')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Table 6.4 information

#Find five largest loadings for PCA1
loads1<-loadings[,1]*loadings[,1]*100
toploadings1<-sort(loadings[,1]*loadings[,1],decreasing=TRUE)[1:5]*100

toploads1<-c()
#Find the largest loads in the list and print them with their signs
toploads1<-c(toploads1,loadings[which(loads1==toploadings1[1]),1])
toploads1<-c(toploads1,loadings[which(loads1==toploadings1[2]),1])
toploads1<-c(toploads1,loadings[which(loads1==toploadings1[3]),1])
toploads1<-c(toploads1,loadings[which(loads1==toploadings1[4]),1])
toploads1<-c(toploads1,loadings[which(loads1==toploadings1[5]),1])

#Make a list of seq names for top five largest loadings
topseqs1<-names(toploadings1)

#Find five largest loadings for PCA2
loads2<-loadings[,2]*loadings[,2]*100
toploadings2<-sort(loadings[,2]*loadings[,2],decreasing=TRUE)[1:5]*100

toploads2<-c()
#Find the largest loads in the list and print them with their signs
toploads2<-c(toploads2,loadings[which(loads2==toploadings2[1]),2])
toploads2<-c(toploads2,loadings[which(loads2==toploadings2[2]),2])
toploads2<-c(toploads2,loadings[which(loads2==toploadings2[3]),2])
toploads2<-c(toploads2,loadings[which(loads2==toploadings2[4]),2])
toploads2<-c(toploads2,loadings[which(loads2==toploadings2[5]),2])

#Make a list of seq names for top five largest loadings
topseqs2<-names(toploadings2)

toploadings<-c(toploadings1,toploadings2)
toploads<-c(toploads1,toploads2)
topnames<-c(topseqs1,topseqs2)

tops<-data.frame(topnames,toploads,toploadings)

write.csv(tops,'Table_6.4.csv')

