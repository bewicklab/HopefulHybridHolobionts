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
library('HybridMicrobiomes')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(10)


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


#Find the treatment group for each of your samples
type<-sample_data(ZNC_physeq)$speciesxsite

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Perform Bootstrap on pooled hybrid population from both locations

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Convert names of groups to numbers with 1= parent 1, 2 = hybrid and 3 = parent 2
no_type<-rep(1,length(type))
no_type[which(type==types[1])]<-1
no_type[which(type==types[2])]<-2
no_type[which(type==types[3])]<-2
no_type[which(type==types[4])]<-3
full_type<-no_type

#Perform Bootstrap
bootstraps<-FourHbootstrap(ZNC_physeq,full_type,0.5, 1000,12,rarefy_each_step = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Perform Bootstrap on hybrid population from SNW (parents from SNW and SBluG)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Consider only hybrids at SNW
ZNC_physeq_SNW = prune_samples(sample_data(ZNC_physeq)$speciesxsite %in% c('SBluG_ino','SNW_neo','SNW_marm'), ZNC_physeq)

#Find the treatment group for each of your samples
type_SNW<-sample_data(ZNC_physeq_SNW)$speciesxsite

#Convert names of groups to numbers, 1= parent 1, 2= hybrid, 3= parent 2
no_type_SNW<-rep(1,length(type_SNW))
no_type_SNW[which(type_SNW==types[1])]<-1
no_type_SNW[which(type_SNW==types[2])]<-2
no_type_SNW[which(type_SNW==types[3])]<-2
no_type_SNW[which(type_SNW==types[4])]<-3
full_type_SNW<-no_type_SNW

#Perform bootstrap
bootstraps_SNW<-FourHbootstrap(ZNC_physeq_SNW,full_type_SNW,0.5, 1000,12, rarefy_each_step = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Perform Bootstrap on hybrid population from SBluG (parents from SNW and SBluG)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Consider only hybrids at SBluG
ZNC_physeq_SBluG = prune_samples(sample_data(ZNC_physeq)$speciesxsite %in% c('SBluG_ino','SBluG_neo','SNW_marm'), ZNC_physeq)

#Find the treatment group for each of your samples
type_SBluG<-sample_data(ZNC_physeq_SBluG)$speciesxsite

#Convert names of groups to numbers, 1 = parent 1, 2 = hybrid, 3= parent 2
no_type_SBluG<-rep(1,length(type_SBluG))
no_type_SBluG[which(type_SBluG==types[1])]<-1
no_type_SBluG[which(type_SBluG==types[2])]<-2
no_type_SBluG[which(type_SBluG==types[3])]<-2
no_type_SBluG[which(type_SBluG==types[4])]<-3
full_type_SBluG<-no_type_SBluG

#Perform bootstrap
bootstraps_SBluG<-FourHbootstrap(ZNC_physeq_SBluG,full_type_SBluG,0.5, 1000,12, rarefy_each_step = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot results, including bootstraps, bootstrap centroid and null plane
    #for Supp. Video 3.3

    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#pooled hybrid populations
FourHquaternary(bootstraps,col='purple')
FourHnullplane(bootstraps,col='purple')

#hybrids from SNW
FourHquaternary(bootstraps_SNW,col='lightslateblue',addplot=TRUE)
FourHnullplane(bootstraps_SNW,col='lightslateblue')

#hybrids from SBluG
FourHquaternary(bootstraps_SBluG,col='magenta',addplot=TRUE)
FourHnullplane(bootstraps_SBluG,col='magenta')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Find Centroids and Information on Location relative to Null Plane

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printout<-c()

#Find the centroids/FourH metric
printout<-rbind(printout,FourHcentroid(bootstraps))
printout<-rbind(printout,FourHcentroid(bootstraps_SNW))
printout<-rbind(printout,FourHcentroid(bootstraps_SBluG))

printout2<-c()

Fn<-FourHnullplaneD(bootstraps)
printout2<-rbind(printout2,c(Fn$mean_intersection_preference,Fn$p_value_intersection_preference))
FnSNW<-FourHnullplaneD(bootstraps_SNW)
printout2<-rbind(printout2,c(FnSNW$mean_intersection_preference,FnSNW$p_value_intersection_preference))
FnSBluG<-FourHnullplaneD(bootstraps_SBluG)
printout2<-rbind(printout2,c(FnSBluG$mean_intersection_preference,FnSBluG$p_value_intersection_preference))

write.csv(printout,'SkinASV4H.csv') #Table S3.1 rows 4-6 & Table S3.2 sigma rows 4-6
write.csv(printout2,'SkinASVnull.csv') #Table S3.2 rows 4-6 I-Inull & p

