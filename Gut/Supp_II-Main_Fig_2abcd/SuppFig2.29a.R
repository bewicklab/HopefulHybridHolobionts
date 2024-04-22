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
library("VennDiagram")
library('scales')
library('ape')
library('picante')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#Remove samples below a certain number of reads
highreads<-which(colSums(otu_table(ZNC_physeq))>10000)
ZNC_physeq<-prune_samples(sample_names(ZNC_physeq) %in% c(sample_names(ZNC_physeq)[highreads]),ZNC_physeq)

#Rarefy the OTU table
rZOTU<-rarefy_even_depth(ZNC_physeq,rngseed = 1)

type_full<-sample_data(rZOTU)$speciesxsite

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Find the core microbiome for each group

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Define the core fraction
core<-0.5

#The treatment group lables
types<-unique(type_full)

#Find the samples corresponding to each treatment group
pickgroups<-c()
for (k in 1:length(types)){
  #Find the animals from each group
  pick<-which(sample_data(rZOTU)$speciesxsite==types[k])
  pickgroups[[types[k]]]<-pick
}


#Initialize a vector to store the names of the core microbial taxa from each treatment group
corenames<-c()
allnames<-c()
#For each treatment group...
for (k in 1:length(types)){
  #Make a new phyloseq object with just the animals from that treatment group
  rZOTUtemp<-prune_samples(sample_names(rZOTU)[pickgroups[[types[k]]]],rZOTU)
  #Sum up the presences of each microbial taxon in the subsetted phyloseq object... count how many microbial taxa are present on as many or more individuals than the defined core fraction
  #Make a list of the rownumbers associated with microbial taxa that satisfy the core requirement
  find_core<-which(rowSums(sign(otu_table(rZOTUtemp)))>=core*length(pickgroups[[types[k]]]))
  #Find the names of the microbial taxa that satisfy the core requirement
  find_names<-rownames(otu_table(rZOTUtemp))[find_core]
  #Make a list of the rownumbers associated with microbial taxa that are present on at least one individual
  find_all<-which(rowSums(sign(otu_table(rZOTUtemp)))>=1)
  #Find the names of the microbial taxa that satisfy the core requirement
  find_names_all<-rownames(otu_table(rZOTUtemp))[find_all]

  #Initialize a vector of the microbial taxon names at the defined level of interest
  list_names<-c()
  list_names_all<-c()

  #For each taxon in the list...
  for (j in 1:length(find_names)){
      list_names<-c(list_names,find_names[j])
  }
  #For each taxon in the list...
  for (j in 1:length(find_names_all)){
    list_names_all<-c(list_names_all,find_names_all[j])
  }


  #Add the list of core microbial names to a 'list', where each entry is a list of microbial core taxa names for a particular treatment group
  corenames[[types[k]]]<-list_names
  allnames[[types[k]]]<-list_names_all

}

#Sort your core microbiome names in the order you want them to appear in the Venn Diagram
corenames_sort<-corenames[c("SBluG_ino" , "SNW_marm" , "SBluG_neo" , "SNW_neo")]
allnames_sort<-allnames[c("SBluG_ino" , "SNW_marm" , "SBluG_neo" , "SNW_neo")]

#Draw the venn diagram
vd <- venn.diagram(
  x = corenames_sort,
  #The names you want to put by each 'circle' in the Venn diagram - should be in same order as your lists from line 188
  category.names = c("SBluG_ino" , "SNW_marm" , "SBluG_neo" , "SNW_neo"),
  #The colors you want to use for your treatment groups - again, should be in same order as your lists from line 188
  fill=c(alpha("firebrick1", .3), alpha("blue", .3), alpha("magenta", .3), alpha("lightslateblue", .3)),
  #The colors you want to use for the lines around your treatment groups
  col=c("firebrick1", "blue", "magenta", "lightslateblue"),
  #Size of the numbers in the Venn diagram
  cex = 2.5, fontface = "bold", fontfamily = "sans", #numbers
  #Size of the names you use for each circle in your Venn diagram
  cat.cex = 2.5, cat.fontface = "bold", cat.default.pos = "outer", #set names
  #cat.pos shows the position of each circle in the Venn diagram... in degrees... with 0 at 12 o'clock - again this should be in the same order as your lists from line 188
  #cat.dist shows the distance of each circle in the Venn diagram from the 'imaginary' circle of their centroids - again this should be in the same order as your lists from line 188
  cat.pos = c(-32, 32, -30, 30), cat.dist = c(0.28, 0.28, 0.18, 0.18), cat.fontfamily = "sans",
  #increase margin to fit the whole figure and circle labels in
  filename = NULL,margin=0.2
)
grid.draw(vd)

in1234<-intersect(intersect(corenames_sort[[2]],corenames_sort[[1]]),intersect(corenames_sort[[3]],corenames_sort[[4]]))
in123<-setdiff(intersect(intersect(corenames_sort[[2]],corenames_sort[[1]]),corenames_sort[[3]]),corenames_sort[[4]])
in124<-setdiff(intersect(intersect(corenames_sort[[2]],corenames_sort[[1]]),corenames_sort[[4]]),corenames_sort[[3]])
in134<-setdiff(intersect(intersect(corenames_sort[[3]],corenames_sort[[1]]),corenames_sort[[4]]),corenames_sort[[2]])
in234<-setdiff(intersect(intersect(corenames_sort[[3]],corenames_sort[[2]]),corenames_sort[[4]]),corenames_sort[[1]])
in12<-setdiff(setdiff(intersect(corenames_sort[[1]],corenames_sort[[2]]),corenames_sort[[3]]),corenames_sort[[4]])
in13<-setdiff(setdiff(intersect(corenames_sort[[1]],corenames_sort[[3]]),corenames_sort[[2]]),corenames_sort[[4]])
in14<-setdiff(setdiff(intersect(corenames_sort[[1]],corenames_sort[[4]]),corenames_sort[[3]]),corenames_sort[[2]])
in23<-setdiff(setdiff(intersect(corenames_sort[[3]],corenames_sort[[2]]),corenames_sort[[1]]),corenames_sort[[4]])
in24<-setdiff(setdiff(intersect(corenames_sort[[4]],corenames_sort[[2]]),corenames_sort[[3]]),corenames_sort[[1]])
in34<-setdiff(setdiff(intersect(corenames_sort[[3]],corenames_sort[[4]]),corenames_sort[[1]]),corenames_sort[[2]])
in1<-setdiff(setdiff(setdiff(corenames_sort[[1]],corenames_sort[[2]]),corenames_sort[[3]]),corenames_sort[[4]])
in2<-setdiff(setdiff(setdiff(corenames_sort[[2]],corenames_sort[[1]]),corenames_sort[[3]]),corenames_sort[[4]])
in3<-setdiff(setdiff(setdiff(corenames_sort[[3]],corenames_sort[[2]]),corenames_sort[[1]]),corenames_sort[[4]])
in4<-setdiff(setdiff(setdiff(corenames_sort[[4]],corenames_sort[[2]]),corenames_sort[[3]]),corenames_sort[[1]])

#There are no taxa unique to SBluG_ino
#checker1<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in1),which(sample_data(rZOTU)$speciesxsite %in% c(types[3],types[4]))])
#Find the taxa unique to SNW_marm and pull out abundances on either hybrid array
checker2<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in2),which(sample_data(rZOTU)$speciesxsite %in% c(types[3],types[4]))])
#See if it's detected on SBluG_ino
checker22<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in2),which(sample_data(rZOTU)$speciesxsite %in% c(types[3],types[4], types[1]))])
#Find the taxa unique to parents and pull out abundances on either hybrid array
checker12<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in12),which(sample_data(rZOTU)$speciesxsite %in% c(types[3],types[4]))])

#Find the taxa with zero abundance on either hybrid population
#seqs1<-which(rowSums(checker1)==0)
seqs2<-which(rowSums(checker2)==0)
seqs22<-which(rowSums(checker22)==0)
seqs12<-which(rowSums(checker12)==0)

#Find the names of those taxa
assigned_taxa[which(assigned_taxa[,1] %in% rownames(checker2[seqs2])),2]
assigned_taxa[which(assigned_taxa[,1] %in% rownames(checker2[seqs22])),2]
assigned_taxa[which(assigned_taxa[,1] %in% rownames(checker12[seqs12])),2]

#Find the taxa unique to the hybrid populations and pull out abundances on either parent poulation
checker3<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in3),which(sample_data(rZOTU)$speciesxsite %in% c(types[1],types[2]))])
checker4<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in4),which(sample_data(rZOTU)$speciesxsite %in% c(types[1],types[2]))])
checker34<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in34),which(sample_data(rZOTU)$speciesxsite %in% c(types[1],types[2]))])

#Find the taxa with zero abundance on either parent population
seqs3<-which(rowSums(checker3)==0)
seqs4<-which(rowSums(checker4)==0)
seqs34<-which(rowSums(checker34)==0)

#None are not found on the parents

for (k in 1:length(in1234)){
  findme<-which(assigned_taxa[,1]==in1234[k])
  in1234[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in123)){
  findme<-which(assigned_taxa[,1]==in123[k])
  in123[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in124)){
  findme<-which(assigned_taxa[,1]==in124[k])
  in124[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in134)){
  findme<-which(assigned_taxa[,1]==in134[k])
  in134[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in234)){
  findme<-which(assigned_taxa[,1]==in234[k])
  in234[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in12)){
  findme<-which(assigned_taxa[,1]==in12[k])
  in12[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in13)){
  findme<-which(assigned_taxa[,1]==in13[k])
  in13[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in14)){
  findme<-which(assigned_taxa[,1]==in14[k])
  in14[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in23)){
  findme<-which(assigned_taxa[,1]==in23[k])
  in23[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in24)){
  findme<-which(assigned_taxa[,1]==in24[k])
  in24[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in34)){
  findme<-which(assigned_taxa[,1]==in34[k])
  in34[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in1)){
  findme<-which(assigned_taxa[,1]==in1[k])
  in1[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in2)){
  findme<-which(assigned_taxa[,1]==in2[k])
  in2[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in3)){
  findme<-which(assigned_taxa[,1]==in3[k])
  in3[k]<-assigned_taxa[findme,2]
}

for (k in 1:length(in4)){
  findme<-which(assigned_taxa[,1]==in4[k])
  in4[k]<-assigned_taxa[findme,2]
}










maxer<-max(max(max(max(max(max(max(max(max(max(max(max(max(max(length(in1234),length(in123)),length(in124)),length(in134)),length(in234)),length(in12)),length(in13)),length(in14)),length(in23)),length(in24)),length(in34)),length(in1)),length(in2)),length(in3)),length(in4))

while (length(in1234)<maxer){
  in1234<-c(in1234,'')
}

while (length(in123)<maxer){
  in123<-c(in123,'')
}

while (length(in124)<maxer){
  in124<-c(in124,'')
}

while (length(in134)<maxer){
  in134<-c(in134,'')
}

while (length(in234)<maxer){
  in234<-c(in234,'')
}

while (length(in12)<maxer){
  in12<-c(in12,'')
}

while (length(in13)<maxer){
  in13<-c(in13,'')
}

while (length(in14)<maxer){
  in14<-c(in14,'')
}

while (length(in23)<maxer){
  in23<-c(in23,'')
}

while (length(in24)<maxer){
  in24<-c(in24,'')
}

while (length(in34)<maxer){
  in34<-c(in34,'')
}

while (length(in1)<maxer){
  in1<-c(in1,'')
}

while (length(in2)<maxer){
  in2<-c(in2,'')
}

while (length(in3)<maxer){
  in3<-c(in3,'')
}

while (length(in4)<maxer){
  in4<-c(in4,'')
}


in1234<-c('SBluG_ino+SNW_marm+SBluG_neo+SNW_neo',in1234)
in123<-c('SBluG_ino+SNW_marm+SBluG_neo',in123)
in124<-c('SBluG_ino+SNW_marm+SNW_neo',in124)
in134<-c('SBluG_ino+SBluG_neo+SNW_neo',in134)
in234<-c('SNW_marm+SBluG_neo+SNW_neo',in234)
in12<-c('SBluG_ino+SNW_marm',in12)
in13<-c('SBluG_ino+SBluG_neo',in13)
in14<-c('SBluG_ino+SNW_neo',in14)
in23<-c('SNW_marm+SBluG_neo',in23)
in24<-c('SNW_marm++SNW_neo',in24)
in34<-c('SBluG_neo+SNW_neo',in34)
in1<-c('SBluG_ino',in1)
in2<-c('SNW_marm',in2)
in3<-c('SBluG_neo',in3)
in4<-c('SNW_neo',in4)

sharing<-data.frame(in1234,in123,in124,in134,in234,in12,in13,in14,in23,in24,in34,in1,in2,in3,in4)

write.csv(sharing,'ASV_gut_Venn.csv') #Fig. 2.29a



