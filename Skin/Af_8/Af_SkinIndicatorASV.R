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
library('indicspecies')
library('data.table')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
assigned_taxa<-read.csv('ASV_tax_assignments.csv',header = FALSE)
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
ZNC_physeq = prune_taxa(assigned_taxa_seqs, Z_physeq)

#Remove samples below a certain number of reads
highreads<-which(colSums(otu_table(ZNC_physeq))>10000)
ZNC_physeq<-prune_samples(sample_names(ZNC_physeq) %in% c(sample_names(ZNC_physeq)[highreads]),ZNC_physeq)

#Rarefy the OTU table
rZOTU<-rarefy_even_depth(ZNC_physeq,rngseed = 1)

type_species<-sample_data(rZOTU)$species

type_hybrid<-sample_data(rZOTU)$hybrid


abundant_taxa<-rownames(otu_table(rZOTU))[which(apply(otu_table(rZOTU),MARGIN=1,FUN=max)/colSums(otu_table(rZOTU))[1]>0.01)]


set.seed(1)

indval_hybrid = multipatt(t(otu_table(rZOTU)), type_hybrid, control = how(nperm=50000), duleg=TRUE)
summary(indval_hybrid,indvalcomp=TRUE)
#extract table of stats
indisp.sign<-as.data.table(indval_hybrid$sign, keep.rownames=TRUE)
indisp.sign<-indisp.sign[which(indisp.sign$rn %in% abundant_taxa),]
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
df_hybrid<-data.frame(indisp.sign[p.value<=0.05, ])
genus_names_hybrid<-c()
for (j in 1:length(df_hybrid$rn)){
  genus_names_hybrid<-c(genus_names_hybrid,assigned_taxa[which(assigned_taxa[,1]==df_hybrid[j,1]),2])
}
df_hybrid<-data.frame(genus_names_hybrid,df_hybrid)
write.csv(df_hybrid,'hybrid_ASV_indicators.csv')

set.seed(1)
indval_species = multipatt(t(otu_table(rZOTU)), type_species, control = how(nperm=50000), duleg=TRUE)
summary(indval_species,indvalcomp=TRUE)
#extract table of stats
indisp_species.sign<-as.data.table(indval_species$sign, keep.rownames=TRUE)
indisp_species.sign<-indisp_species.sign[which(indisp_species.sign$rn %in% abundant_taxa),]
#add adjusted p-value
indisp_species.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
df_species<-data.frame(indisp_species.sign[p.value<=0.05, ])
genus_names_species<-c()
for (j in 1:length(df_species$rn)){
  genus_names_species<-c(genus_names_species,assigned_taxa[which(assigned_taxa[,1]==df_species[j,1]),2])
}
df_species<-data.frame(genus_names_species,df_species)
write.csv(df_species,'species_ASV_indicators.csv')






