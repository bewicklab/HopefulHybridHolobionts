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
library('ape')
library('picante')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Read in CSV file for OTU table (table.csv for ASVs, tablelevel6.csv for genera)
data <- import_biom('HybridLizardGuts.biom')

#Find OTU table in biom file
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

#Remove samples below a certain number of reads (I found that the prune_taxa command was somehow adding back in the low read samples... so I'm removing them here...)
highreads<-which(colSums(otu_table(ZNC_physeq))>10000)
ZNC_physeq<-prune_samples(sample_names(ZNC_physeq) %in% c(sample_names(ZNC_physeq)[highreads]),ZNC_physeq)

#Rarefy the dataset
rZOTU<-rarefy_even_depth(ZNC_physeq, rngseed = 1)
#Define interesting treatment groups
type_full<-sample_data(rZOTU)$speciesxsite

core<-0.5

#Sharon's code...
#Are microbes 100% unique or just core unique on each population???

#Make phyloseq objects for the various groupings
marmZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite=='SNW_marm',rZOTU)
inoZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite=='SBluG_ino',rZOTU)
parentZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite %in% c('SBluG_ino','SNW_marm'),rZOTU)
neoSNWZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite=='SNW_neo',rZOTU)
neoSBluGZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite=='SBluG_neo',rZOTU)
neoZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite %in% c('SBluG_neo','SNW_neo'),rZOTU)

#Find core taxa from each phyloseq object
core_marm<-which(rowSums(sign(otu_table(marmZOTU)))>=core*length(colnames(otu_table(marmZOTU))))
core_marm_names<-rownames(otu_table(marmZOTU))[core_marm]
core_ino<-which(rowSums(sign(otu_table(inoZOTU)))>=core*length(colnames(otu_table(inoZOTU))))
core_ino_names<-rownames(otu_table(inoZOTU))[core_ino]
core_neoSNW<-which(rowSums(sign(otu_table(neoSNWZOTU)))>=core*length(colnames(otu_table(neoSNWZOTU))))
core_neoSNW_names<-rownames(otu_table(neoSNWZOTU))[core_neoSNW]
core_neoSBluG<-which(rowSums(sign(otu_table(neoSBluGZOTU)))>=core*length(colnames(otu_table(neoSBluGZOTU))))
core_neoSBluG_names<-rownames(otu_table(neoSBluGZOTU))[core_neoSBluG]
core_neo<-which(rowSums(sign(otu_table(neoZOTU)))>=core*length(colnames(otu_table(neoZOTU))))
core_neo_names<-rownames(otu_table(neoZOTU))[core_neo]

in_marm<-which(rowSums(sign(otu_table(marmZOTU)))>=1)
in_marm_names<-rownames(otu_table(marmZOTU))[in_marm]
in_ino<-which(rowSums(sign(otu_table(inoZOTU)))>=1)
in_ino_names<-rownames(otu_table(inoZOTU))[in_ino]
in_neoSNW<-which(rowSums(sign(otu_table(neoSNWZOTU)))>=1)
in_neoSNW_names<-rownames(otu_table(neoSNWZOTU))[in_neoSNW]
in_neoSBluG<-which(rowSums(sign(otu_table(neoSBluGZOTU)))>=1)
in_neoSBluG_names<-rownames(otu_table(neoSBluGZOTU))[in_neoSBluG]
in_neo<-which(rowSums(sign(otu_table(neoZOTU)))>=1)
in_neo_names<-rownames(otu_table(neoZOTU))[in_neo]


#Taxa that are part of both progenitor cores (progenitor cores calculated separately), but not on neos
core_progenitor_not_on_neo<-setdiff(intersect(core_marm_names,core_ino_names),in_neo_names)
print('Core progenitor not on hybrid:')
for (k in 1:length(core_progenitor_not_on_neo)){
  print(tax_table(marmZOTU)[which(rownames(tax_table(marmZOTU))==core_progenitor_not_on_neo[k])])
}
print('')

#Taxa that are part of the marm core, but not on neos
core_marm_not_on_neo<-setdiff(core_marm_names,in_neo_names)
print('Core marm not on hybrid:')
for (k in 1:length(core_marm_not_on_neo)){
  print(tax_table(marmZOTU)[which(rownames(tax_table(marmZOTU))==core_marm_not_on_neo[k])])
}
print('')

#Check to see if taxa that are part of the marm core, are on inos (though not part of core), but not on neos
core_marm_not_on_neo_or_ino<-setdiff(intersect(core_marm_names,in_ino_names),in_neo_names)
print('Core marm not on hybrid or ino:')
for (k in 1:length(core_marm_not_on_neo_or_ino)){
  print(tax_table(marmZOTU)[which(rownames(tax_table(marmZOTU))==core_marm_not_on_neo_or_ino[k])])
}
print('')

#Taxa that are part of the ino core, but not on neos
print('Core ino not on hybrid:')
core_ino_not_on_neo<-setdiff(core_ino_names,in_neo_names)
for (k in 1:length(core_ino_not_on_neo)){
  print(tax_table(inoZOTU)[which(rownames(tax_table(inoZOTU))==core_ino_not_on_neo[k])])
}

print('Core hybrid not on marms:')
core_neo_not_on_marm<-setdiff(core_neo_names,in_marm_names)
for (k in 1:length(core_neo_not_on_marm)){
  print(tax_table(neoZOTU)[which(rownames(tax_table(neoZOTU))==core_neo_not_on_marm[k])])
}#Nothing

print('Core hybrid not on inos:')
core_neo_not_on_ino<-setdiff(core_neo_names,in_ino_names)
for (k in 1:length(core_neo_not_on_ino)){
  print(tax_table(neoZOTU)[which(rownames(tax_table(neoZOTU))==core_neo_not_on_ino[k])])
}#Nothing



