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
data <- import_biom('HybridLizardSkins.biom')

#Find ASV table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

#Read in metadata
meta_csv<-read.csv('Whiptail_Skin_Metadata.csv', row.names=1,header= TRUE)

#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)

#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-multi2di(read.tree('rooted_tree_out/tree.nwk'))

#Rename the tree tips to match the ASV table names
cut_names<-c()
for (k in 1:length(taxa_names(tree_file))){
  #Depending on the Qiime2 output, you may need to split the names with a space or with an underscore
  #cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],' ')[[1]][1])
  cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],'_')[[1]][1])
}
taxa_names(tree_file)<-cut_names

#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq<-phyloseq(OTU_biom,TAX,SAMP,tree_file)

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
ZNC_physeq = prune_taxa(assigned_taxa_seqs, Z_physeq)

#Define interesting treatment groups
type_full<-sample_data(ZNC_physeq)$speciesxsite

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm')


#Rarefy the OTU table specifically for this bootstrap sample
rtempZOTU<-rarefy_even_depth(ZNC_physeq,verbose=FALSE,rngseed = sample(1000,1))

#Remove taxa that drop out of the rarefied OTU table (this makes UniFrac run more smoothly... and you don't need them anyhow)
rZOTU<-prune_taxa(row.names(tax_table(rtempZOTU)),rtempZOTU)

#Find the treatment group for each of your samples
type<-sample_data(rZOTU)$speciesxsite

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm')

#Convert names of groups to numbers
no_type<-rep(1,length(type))
for (k in 1:length(types)){
  no_type[which(type==types[k])]<-k
}

#Define the color scheme for plotting
colvec<-c("red", "magenta","lightslateblue", "blue")



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

#                 Calculate a variety of different BETA diversity metrics using the VEGAN package

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
euclidean<-vegdist(t(otu_table(rZOTU)),method='euclidean',upper=TRUE,diag=TRUE,binary=FALSE)
jaccard<-vegdist(t(otu_table(rZOTU)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(rZOTU)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(rZOTU,weighted=FALSE)
wunifrac<-UniFrac(rZOTU,weighted=TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform 2D PCoA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Perform PcoA
#pcoa_euclidean<-pcoa(euclidean)$vectors[,1:2]
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]
pcoa_unifrac<-pcoa(unifrac)$vectors[,1:2]
pcoa_wunifrac<-pcoa(wunifrac)$vectors[,1:2]


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
  ggtitle('A.') + annotate("text", x=30000, y=20000, label= "Euclidean", size=13)

#Main Fig. 5C
ggplot() + geom_point(data=scores_2, aes(x=PC1, y=PC2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=pca_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PC1") + ylab("PC2") + geom_polygon(data = pca_hull, alpha = 0.2, aes(x=PC1, y=PC2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank())
  

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
  geom_point(data=pcoa_jacc_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site), colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PCoA1") + ylab("PCoA2") + geom_polygon(data = pcoa_jacc_hull, alpha = 0.2, aes(x=Axis.1, y=Axis.2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) + 
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('B.') + annotate("text", x=0.05, y=.2, label= "Jaccard", size=13)

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
  geom_point(data=pcoa_bray_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site), colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PCoA1") + ylab("PCoA2") + geom_polygon(data = pcoa_bray_hull, alpha = 0.2, aes(x=Axis.1, y=Axis.2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) + 
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('C.') + annotate("text", x=-0.40, y=-0.2, label= "Bray-Curtis", size=13)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the UniFrac PCoA (equivalent to PCA)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcoa_unifrac <- data.frame(pcoa_unifrac) #convert to dataframe for ggplot2
pcoa_unifrac$spec_site <- as.factor(type) #pull grouping variable from existing vector
pcoa_unifrac_hull <- pcoa_unifrac %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(Axis.1, Axis.2))
pcoa_unifrac_centroid <- pcoa_unifrac %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_PC1 = mean(Axis.1),    
  mean_PC2 = mean(Axis.2),)

#plot Unifrac PCoA
PCoA_unifrac <- ggplot() + geom_point(data=pcoa_unifrac, aes(x=Axis.1, y=Axis.2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=pcoa_unifrac_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site), colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PCoA1") + ylab("PCoA2") + geom_polygon(data = pcoa_unifrac_hull, alpha = 0.2, aes(x=Axis.1, y=Axis.2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) + 
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('D.') + annotate("text", x=-0.35, y=0.1, label= "Unifrac", size=13)

#Main Fig. 4C
ggplot() + geom_point(data=pcoa_unifrac, aes(x=Axis.1, y=Axis.2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=pcoa_unifrac_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site), colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PCoA1") + ylab("PCoA2") + geom_polygon(data = pcoa_unifrac_hull, alpha = 0.2, aes(x=Axis.1, y=Axis.2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) + 
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) 
 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the weighted Unifrac PCoA (equivalent to PCA)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pcoa_wunifrac <- data.frame(pcoa_wunifrac) #convert to dataframe for ggplot2
pcoa_wunifrac$spec_site <- as.factor(type) #pull grouping variable from existing vector
pcoa_wunifrac_hull <- pcoa_wunifrac %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(Axis.1, Axis.2))
pcoa_wunifrac_centroid <- pcoa_wunifrac %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_PC1 = mean(Axis.1),    
  mean_PC2 = mean(Axis.2),)

#plot Weighted Unifrac PCoA
PCoA_wunifrac <- ggplot() + geom_point(data=pcoa_wunifrac, aes(x=Axis.1, y=Axis.2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=pcoa_wunifrac_centroid, aes(x=mean_PC1,y=mean_PC2, fill=spec_site), colour="black", pch=21, stroke=1.5, size=7) +
  xlab("PCoA1") + ylab("PCoA2") + geom_polygon(data = pcoa_wunifrac_hull, alpha = 0.2, aes(x=Axis.1, y=Axis.2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) + 
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('E.') + annotate("text", x=-0.065, y=0.045, label= "Weighted Unifrac", size=13)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#make panel plot (Fig. 6.4)

library("gridExtra")
library("grid")

panel_list <- list(PCA, PCoA_jacc, PCoA_bray, PCoA_unifrac, PCoA_wunifrac)
  
panel_layout <- rbind(c(1,2,3),
                   c(4,5,NA))
                   
grid_panel <- grid.arrange(grobs = panel_list, layout_matrix = panel_layout)

grid.draw(grid_panel) # interactive device

ggsave("Figure_SI6.4.png", grid_panel, height = 20, width = 30)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Table 6.5 information

#Percent Variation/Cumulative Percent Variation for each PCoA
pv_euclidean<-pcoa(euclidean)$values[,2]*100
cs_euclidean<-cumsum(pv_euclidean)
pv_jaccard<-pcoa(jaccard)$values[,2]*100
cs_jaccard<-cumsum(pv_jaccard)
pv_bray<-pcoa(bray)$values[,2]*100
cs_bray<-cumsum(pv_bray)
pv_unifrac<-pcoa(unifrac)$values[,2]*100
cs_unifrac<-cumsum(pv_unifrac)
pv_wunifrac<-pcoa(wunifrac)$values[,2]*100
cs_wunifrac<-cumsum(pv_wunifrac)

round2 <- function(x) { #rounding function to two decimal places
  format(round((cumsum(x)),digits=2),nsmall=2)
}

#create table 6.5
tbl_6.5 <- data.frame(
  "metric" = c("Euclidean", "Jaccard", "Bray-Curtis", "Unifrac", "W-Unifrac"),
  "percent1" = c(round2(pv_euclidean[1]), round2(pv_jaccard[1]), round2(pv_bray[1]), round2(pv_unifrac[1]), round2(pv_wunifrac[1])),
  "cumulative1" = c(round2(cs_euclidean[1]), round2(cs_jaccard[1]), round2(cs_bray[1]), round2(cs_unifrac[1]), round2(cs_wunifrac[1])),
  "percent2" = c(round2(pv_euclidean[2]), round2(pv_jaccard[2]), round2(pv_bray[2]), round2(pv_unifrac[2]), round2(pv_wunifrac[2])),
  "cumulative2" = c(round2(cs_euclidean[2]), round2(cs_jaccard[2]), round2(cs_bray[2]), round2(cs_unifrac[2]), round2(cs_wunifrac[2])),
  "percent3" = c(round2(pv_euclidean[3]), round2(pv_jaccard[3]), round2(pv_bray[3]), round2(pv_unifrac[3]), round2(pv_wunifrac[3])),
  "cumulative3" = c(round2(cs_euclidean[3]), round2(cs_jaccard[3]), round2(cs_bray[3]), round2(cs_unifrac[3]), round2(cs_wunifrac[3])),
  "percent4" = c(round2(pv_euclidean[4]), round2(pv_jaccard[4]), round2(pv_bray[4]), round2(pv_unifrac[4]), round2(pv_wunifrac[4])),
  "cumulative4" = c(round2(cs_euclidean[4]), round2(cs_jaccard[4]), round2(cs_bray[4]), round2(cs_unifrac[4]), round2(cs_wunifrac[4])),
  "percent5" = c(round2(pv_euclidean[5]), round2(pv_jaccard[5]), round2(pv_bray[5]), round2(pv_unifrac[5]), round2(pv_wunifrac[5])),
  "cumulative5" = c(round2(cs_euclidean[5]), round2(cs_jaccard[5]), round2(cs_bray[5]), round2(cs_unifrac[5]), round2(cs_wunifrac[5])))

write.csv(tbl_6.5,'Table_6.5.csv')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Table 6.6 information

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ASV_tax_assignments.csv',header = TRUE)
assigned_taxa_seqs<-assigned_taxa[,1]


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

nametop1<-c()
#Find taxon names for top five largest loadings
nametop1<-c(nametop1,assigned_taxa[which(assigned_taxa[,1]==topseqs1[1]),2])
nametop1<-c(nametop1,assigned_taxa[which(assigned_taxa[,1]==topseqs1[2]),2])
nametop1<-c(nametop1,assigned_taxa[which(assigned_taxa[,1]==topseqs1[3]),2])
nametop1<-c(nametop1,assigned_taxa[which(assigned_taxa[,1]==topseqs1[4]),2])
nametop1<-c(nametop1,assigned_taxa[which(assigned_taxa[,1]==topseqs1[5]),2])

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

nametop2<-c()
#Find taxon names for top five largest loadings
nametop2<-c(nametop2,assigned_taxa[which(assigned_taxa[,1]==topseqs2[1]),2])
nametop2<-c(nametop2,assigned_taxa[which(assigned_taxa[,1]==topseqs2[2]),2])
nametop2<-c(nametop2,assigned_taxa[which(assigned_taxa[,1]==topseqs2[3]),2])
nametop2<-c(nametop2,assigned_taxa[which(assigned_taxa[,1]==topseqs2[4]),2])
nametop2<-c(nametop2,assigned_taxa[which(assigned_taxa[,1]==topseqs2[5]),2])

toploadings<-c(toploadings1,toploadings2)
topnames<-c(nametop1,nametop2)
toploads<-c(toploads1,toploads2)

tops<-data.frame(topnames,toploads,toploadings)

write.csv(tops,'Table_6.6.csv')

