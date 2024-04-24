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
data <- import_biom('HybridLizardSkinsGenus.biom')

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
meta_csv<-read.csv('Whiptail_Skin_Metadata.csv', row.names=1,header= TRUE)

#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)


#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq<-phyloseq(OTU_biom,TAX,SAMP)

#Rarefy the OTU table
rZOTU<-rarefy_even_depth(Z_physeq, rngseed=5)

type<-sample_data(rZOTU)$speciesxsite


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Find the core microbiome for each group

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Define the core fraction
core<-0.5

#The treatment group lables
types<-unique(type)

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

checker1<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in1),which(sample_data(rZOTU)$speciesxsite %in% c(types[3],types[4]))])
checker2<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in2),which(sample_data(rZOTU)$speciesxsite %in% c(types[3],types[4]))])
checker12<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in12),which(sample_data(rZOTU)$speciesxsite %in% c(types[3],types[4]))])

seqs1<-which(rowSums(checker1)==0)
seqs2<-which(rowSums(checker2)==0)
seqs12<-which(rowSums(checker12)==0)

checker3<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in3),which(sample_data(rZOTU)$speciesxsite %in% c(types[1],types[2]))])
checker4<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in4),which(sample_data(rZOTU)$speciesxsite %in% c(types[1],types[2]))])
checker34<-(otu_table(rZOTU)[which(rownames(otu_table(rZOTU)) %in% in34),which(sample_data(rZOTU)$speciesxsite %in% c(types[1],types[2]))])

seqs3<-which(rowSums(checker3)==0)
seqs4<-which(rowSums(checker4)==0)
seqs34<-which(rowSums(checker34)==0)


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

write.csv(sharing,'Genus_skin_Venn.csv')



