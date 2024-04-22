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
data <- import_biom('HybridLizardGuts.biom')

#Find ASV table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

#Read in metadata
meta_csv<-read.csv('Whiptail_Gut_Metadata.csv', row.names=1,header= TRUE)

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

#                 Perform 2D nmds

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Perform nmds
nmds_euclidean<-metaMDS(euclidean,dimensions=2)
nmds_jaccard<-metaMDS(jaccard, dimensions=2)
nmds_bray<-metaMDS(bray,dimensions=2)
nmds_unifrac<-metaMDS(unifrac,dimensions=2)
nmds_wunifrac<-metaMDS(wunifrac,dimensions=2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the Euclidean NMDS

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
  xlab("NMDS1") + ylab("NMDS2") + geom_polygon(data = euclidean_hull, alpha = 0.2, aes(x=NMDS1, y=NMDS2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('A.') + annotate("text", x=-7000, y=4000, label= "Euclidean", size=13)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the Jaccard NMDS

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
  xlab("NMDS1") + ylab("NMDS2") + geom_polygon(data = jaccard_hull, alpha = 0.2, aes(x=NMDS1, y=NMDS2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('B.') + annotate("text", x=-0.17, y=0.35, label= "Jaccard", size=13)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the Bray-Curtis NMDS

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
  xlab("NMDS1") + ylab("NMDS2") + geom_polygon(data = bray_hull, alpha = 0.2, aes(x=NMDS1, y=NMDS2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('C.') + annotate("text", x=0.4, y=0.35, label= "Bray-Curtis", size=13)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the Unifrac NMDS

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmds_unifrac_scores <- data.frame(scores(nmds_unifrac, display="sites")) #convert to NMDS scores to dataframe for ggplot2
nmds_unifrac_scores$spec_site <- as.factor(type) #pull grouping variable from existing vector
unifrac_hull <- nmds_unifrac_scores %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(NMDS1, NMDS2))
unifrac_centroid <- nmds_unifrac_scores %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_NMDS1 = mean(NMDS1),    
  mean_NMDS2 = mean(NMDS2),)

#plot Unifrac NMDS
NMDS_unifrac <- ggplot() + geom_point(data=nmds_unifrac_scores, aes(x=NMDS1, y=NMDS2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=unifrac_centroid, aes(x=mean_NMDS1,y=mean_NMDS2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("NMDS1") + ylab("NMDS2") + geom_polygon(data = unifrac_hull, alpha = 0.2, aes(x=NMDS1, y=NMDS2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('D.') + annotate("text", x=0.25, y=0.25, label= "Unifrac", size=13)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the weighted Unifrac NMDS

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmds_wunifrac_scores <- data.frame(scores(nmds_wunifrac, display="sites")) #convert to NMDS scores to dataframe for ggplot2
nmds_wunifrac_scores$spec_site <- as.factor(type) #pull grouping variable from existing vector
wunifrac_hull <- nmds_wunifrac_scores %>% group_by(spec_site) %>% #calculate convex hulls
  slice(chull(NMDS1, NMDS2))
wunifrac_centroid <- nmds_wunifrac_scores %>% group_by(spec_site) %>% summarise( #calculate group centroid
  mean_NMDS1 = mean(NMDS1),    
  mean_NMDS2 = mean(NMDS2),)

#plot Weighted Unifrac NMDS
NMDS_wunifrac <- ggplot() + geom_point(data=nmds_wunifrac_scores, aes(x=NMDS1, y=NMDS2, fill=spec_site), colour="black", pch=21, size=3) + theme_bw(40) +
  geom_point(data=wunifrac_centroid, aes(x=mean_NMDS1,y=mean_NMDS2, fill=spec_site),  colour="black", pch=21, stroke=1.5, size=7) +
  xlab("NMDS1") + ylab("NMDS2") + geom_polygon(data = wunifrac_hull, alpha = 0.2, aes(x=NMDS1, y=NMDS2, fill = spec_site), colour = "black") + 
  scale_colour_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                      values= c("red", "magenta", "lightslateblue", "blue")) +  
  scale_fill_manual(name = "Species-Site", breaks = c("SBluG_ino", "SBluG_neo", "SNW_neo", "SNW_marm"),
                    values= c("red", "magenta", "lightslateblue", "blue")) + theme(legend.position="none", panel.grid=element_blank()) +
  ggtitle('E.') + annotate("text", x=0.3, y=0.4, label= "Weighted Unifrac", size=13)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#make panel plot (Fig. 6.2)

library("gridExtra")
library("grid")

panel_list <- list(NMDS_euclidean, NMDS_jaccard, NMDS_bray, NMDS_unifrac, NMDS_wunifrac)

panel_layout <- rbind(c(1,2,3),
                      c(4,5,NA))

grid_panel <- grid.arrange(grobs = panel_list, layout_matrix = panel_layout)

grid.draw(grid_panel) # interactive device

ggsave("Figure_SI6.2.png", grid_panel, height = 20, width = 30)

