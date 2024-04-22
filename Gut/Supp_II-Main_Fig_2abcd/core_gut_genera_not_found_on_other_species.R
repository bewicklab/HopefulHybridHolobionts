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
ZNC_physeq<-phyloseq(OTU_biom,TAX,SAMP)

#Rarefy the dataset
rZOTU<-rarefy_even_depth(ZNC_physeq, rngseed = 10)

#Define interesting treatment groups
type_full<-sample_data(ZNC_physeq)$speciesxsite


core<-0.5

#Sharon's code...
#Are microbes 100% unique or just core unique on each population???

#Make phyloseq objects for the various groupings
marmZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite=='SNW_marm',rZOTU)
inoZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite=='SBluG_ino',rZOTU)
neoSNWZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite=='SNW_neo',rZOTU)
neoSBluGZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite=='SBluG_neo',rZOTU)
neoZOTU<-prune_samples(sample_data(rZOTU)$speciesxsite %in% c('SBluG_neo','SNW_neo'),rZOTU)

#Find core taxa from each phyloseq object
core_marm<-which(rowSums(sign(otu_table(marmZOTU)))>=core*length(colnames(otu_table(marmZOTU))))
core_marm_names<-tax_table(marmZOTU)[core_marm,6]
core_ino<-which(rowSums(sign(otu_table(inoZOTU)))>=core*length(colnames(otu_table(inoZOTU))))
core_ino_names<-tax_table(inoZOTU)[core_ino,6]
core_neoSNW<-which(rowSums(sign(otu_table(neoSNWZOTU)))>=core*length(colnames(otu_table(neoSNWZOTU))))
core_neoSNW_names<-tax_table(neoSNWZOTU)[core_neoSNW,6]
core_neoSBluG<-which(rowSums(sign(otu_table(neoSBluGZOTU)))>=core*length(colnames(otu_table(neoSBluGZOTU))))
core_neoSBluG_names<-tax_table(neoSBluGZOTU)[core_neoSBluG,6]
core_neo<-which(rowSums(sign(otu_table(neoZOTU)))>=core*length(colnames(otu_table(neoZOTU))))
core_neo_names<-tax_table(neoZOTU)[core_neo,6]

in_marm<-which(rowSums(sign(otu_table(marmZOTU)))>=1)
in_marm_names<-tax_table(marmZOTU)[in_marm,6]
in_ino<-which(rowSums(sign(otu_table(inoZOTU)))>=1)
in_ino_names<-tax_table(inoZOTU)[in_ino,6]
in_neoSNW<-which(rowSums(sign(otu_table(neoSNWZOTU)))>=1)
in_neoSNW_names<-tax_table(neoSNWZOTU)[in_neoSNW,6]
in_neoSBluG<-which(rowSums(sign(otu_table(neoSBluGZOTU)))>=1)
in_neoSBluG_names<-otu_table(neoSBluGZOTU)[in_neoSBluG,6]
in_neo<-which(rowSums(sign(otu_table(neoZOTU)))>=1)
in_neo_names<-tax_table(neoZOTU)[in_neo,6]

#Taxa that are part of both progenitor cores (progenitor cores calculated separately), but not on neos
core_progenitor_not_on_neo<-setdiff(intersect(core_marm_names,core_ino_names),in_neo_names)
print('Core progenitor not on hybrid:')
print(core_progenitor_not_on_neo)
print('')

#Taxa that are part of the marm core, but not on neos
core_marm_not_on_neo<-setdiff(core_marm_names,in_neo_names)
print('Core marm not on hybrid:')
print(core_marm_not_on_neo)
print('')

#Check to see if taxa that are part of the marm core, are on inos (though not part of core), but not on neos
core_marm_not_on_neo<-setdiff(intersect(core_marm_names,in_ino_names),in_neo_names)
print('Core marm not on hybrid:')
print(core_marm_not_on_neo)
print('')

#Taxa that are part of the ino core, but not on neos
print('Core ino not on hybrid:')
core_ino_not_on_neo<-setdiff(core_ino_names,in_neo_names)
print(core_ino_not_on_neo)
print('')

#Taxa that are part of the neo core, but not on marms
print('Core hybrid not on marms:')
core_neo_not_on_marm<-setdiff(core_neo_names,in_marm_names)
print(core_neo_not_on_marm)
print('')

#Taxa that are part of the neo core, but not on inos
print('Core hybrid not on inos:')
core_neo_not_on_ino<-setdiff(core_neo_names,in_ino_names)
print(core_neo_not_on_ino)
print('')

#Taxa that are part of the marm core, but not on inos
print('Core marms not on inos:')
core_marms_not_on_ino<-setdiff(core_marm_names,in_ino_names)
print(core_marms_not_on_ino)
print('')

#Taxa that are part of the ino core, but not on marms
print('Core inos not on marms:')
core_inos_not_on_marm<-setdiff(core_ino_names,in_marm_names)
print(core_inos_not_on_marm)
print('')






