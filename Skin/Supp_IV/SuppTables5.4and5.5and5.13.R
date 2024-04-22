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
ZNC_physeq<-phyloseq(OTU_biom,TAX,SAMP)

set.seed(3)

#The number from each treatment group you want to select to create a 'balanced' design
sample_size<-15

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm')

#Make a list of the individual animals you're going to randomly choose
picklist<-c()
pickgroups<-matrix(0,ncol=sample_size,nrow=length(types))
for (k in 1:length(types)){
  #Randomly select the appropriate number of animals from each group
  pick<-sample(which(sample_data(ZNC_physeq)$speciesxsite==types[k]),sample_size,replace=FALSE)
  #Add those animals to your list
  picklist<-c(picklist,pick)
  pickgroups[k,]<-pick

}

#Make a phyloseq object of only the individual animals you selected
ZB_physeq<-prune_samples(sample_names(ZNC_physeq)[picklist],ZNC_physeq)

#Rarefy the OTU table specifically for this bootstrap sample
rtempZOTU<-rarefy_even_depth(ZB_physeq,verbose=FALSE,rngseed = sample(1000,1))

#Remove taxa that drop out of the rarefied OTU table (this makes UniFrac run more smoothly... and you don't need them anyhow)
rZOTU<-prune_taxa(row.names(tax_table(rtempZOTU)),rtempZOTU)


#Define interesting treatment groups
type<-sample_data(rZOTU)$speciesxsite

#Rarefy the OTU table specifically for this bootstrap sample
#rtempZOTU<-rarefy_even_depth(ZNC_physeq,verbose=FALSE,rngseed = sample(1000,1))

#Remove taxa that drop out of the rarefied OTU table (this makes UniFrac run more smoothly... and you don't need them anyhow)
#rZOTU<-prune_taxa(row.names(tax_table(rtempZOTU)),rtempZOTU)

#Find the treatment group for each of your samples
#type<-sample_data(rZOTU)$speciesxsite

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
#types<-c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate a variety of different BETA diversity metrics using the VEGAN package

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
euclidean<-vegdist(t(otu_table(rZOTU)),method='euclidean',upper=TRUE,diag=TRUE,binary=FALSE)
jaccard<-vegdist(t(otu_table(rZOTU)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(rZOTU)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)

#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
euclidean_matrix<-as.matrix(euclidean,labels=TRUE)
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Euclidean Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_euclidean<-list(Data=t(otu_table(rZOTU)),D=euclidean_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_euclidean=PERMANOVA(inputpermanova_euclidean, factor(type),nperm=10000)

#If PERMANOVA by treatment group is significant, figure out which groups differ from which other groups
if (permanova_by_group_euclidean$pvalue<0.05 && length(unique(type))>2){

  #Perform a post-hoc pairwise test between each pair of treatment groups
  posthoc_permanova_euclidean<-pairwise.adonis(euclidean_matrix,as.factor(type),p.adjust.m='BH')

  #Write the p-values (corrected for multiple comparisons) from the pairwise.adonis function as a concise data.frame, rather than a list of comparisons
  #(This can be useful for quickly summarizing/interpreting results)
  a1<-c()     #This will be a list of the first treatment group for each pairwise comparison
  a2<-c()     #This will be a list of the second treatment group for each pairwise comparison
  #For each pairwise comparison outputted by the pairwise.adonis function...
  for (k in 1:length(posthoc_permanova_euclidean[,1])){
    #Find the first and second treatment groups being compared
    ff<-str_split(posthoc_permanova_euclidean[,1][k],' vs ')
    a1<-c(a1,ff[[1]][1])
    a2<-c(a2,ff[[1]][2])
  }

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  permanova_matrix_euclidean<-matrix(data='',nrow=length(types),ncol=length(types))
  for (k in 1:length(posthoc_permanova_euclidean[,1])){
    #For each comparison in the list outputted from the pairwise.adonis matrix...
    hit1<-which(types==a1[k])
    hit2<-which(types==a2[k])
    #Record the adjusted p-value in the matrix
    permanova_matrix_euclidean[max(hit1,hit2),min(hit1,hit2)]<-round(posthoc_permanova_euclidean$p.adjusted[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(permanova_matrix_euclidean)<-types
  colnames(permanova_matrix_euclidean)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  permanova_pvalues_euclidean<-data.frame(permanova_matrix_euclidean[2:4,1:3])
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Jaccard Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_jaccard<-list(Data=t(otu_table(rZOTU)),D=jaccard_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_jaccard=PERMANOVA(inputpermanova_jaccard, factor(type),nperm=10000)

#If PERMANOVA by treatment group is significant, figure out which groups differ from which other groups
if (permanova_by_group_jaccard$pvalue<0.05 && length(unique(type))>2){

  #Perform a post-hoc pairwise test between each pair of treatment groups
  posthoc_permanova_jaccard<-pairwise.adonis(jaccard_matrix,as.factor(type),p.adjust.m='BH')

  #Write the p-values (corrected for multiple comparisons) from the pairwise.adonis function as a concise data.frame, rather than a list of comparisons
  #(This can be useful for quickly summarizing/interpreting results)
  a1<-c()     #This will be a list of the first treatment group for each pairwise comparison
  a2<-c()     #This will be a list of the second treatment group for each pairwise comparison
  #For each pairwise comparison outputted by the pairwise.adonis function...
  for (k in 1:length(posthoc_permanova_jaccard[,1])){
    #Find the first and second treatment groups being compared
    ff<-str_split(posthoc_permanova_jaccard[,1][k],' vs ')
    a1<-c(a1,ff[[1]][1])
    a2<-c(a2,ff[[1]][2])
  }

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  permanova_matrix_jaccard<-matrix(data='',nrow=length(types),ncol=length(types))
  for (k in 1:length(posthoc_permanova_jaccard[,1])){
    #For each comparison in the list outputted from the pairwise.adonis matrix...
    hit1<-which(types==a1[k])
    hit2<-which(types==a2[k])
    #Record the adjusted p-value in the matrix
    permanova_matrix_jaccard[max(hit1,hit2),min(hit1,hit2)]<-round(posthoc_permanova_jaccard$p.adjusted[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(permanova_matrix_jaccard)<-types
  colnames(permanova_matrix_jaccard)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  permanova_pvalues_jaccard<-data.frame(permanova_matrix_jaccard[2:4,1:3])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Bray Curtis Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_bray<-list(Data=t(otu_table(rZOTU)),D=bray_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_bray=PERMANOVA(inputpermanova_bray, factor(type),nperm=10000)

#If PERMANOVA by treatment group is significant, figure out which groups differ from which other groups
if (permanova_by_group_bray$pvalue<0.05 && length(unique(type))>2){

  #Perform a post-hoc pairwise test between each pair of treatment groups
  posthoc_permanova_bray<-pairwise.adonis(bray_matrix,as.factor(type),p.adjust.m='BH')

  #Write the p-values (corrected for multiple comparisons) from the pairwise.adonis function as a concise data.frame, rather than a list of comparisons
  #(This can be useful for quickly summarizing/interpreting results)
  a1<-c()     #This will be a list of the first treatment group for each pairwise comparison
  a2<-c()     #This will be a list of the second treatment group for each pairwise comparison
  #For each pairwise comparison outputted by the pairwise.adonis function...
  for (k in 1:length(posthoc_permanova_bray[,1])){
    #Find the first and second treatment groups being compared
    ff<-str_split(posthoc_permanova_bray[,1][k],' vs ')
    a1<-c(a1,ff[[1]][1])
    a2<-c(a2,ff[[1]][2])
  }

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  permanova_matrix_bray<-matrix(data='',nrow=length(types),ncol=length(types))
  for (k in 1:length(posthoc_permanova_bray[,1])){
    #For each comparison in the list outputted from the pairwise.adonis matrix...
    hit1<-which(types==a1[k])
    hit2<-which(types==a2[k])
    #Record the adjusted p-value in the matrix
    permanova_matrix_bray[max(hit1,hit2),min(hit1,hit2)]<-round(posthoc_permanova_bray$p.adjusted[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(permanova_matrix_bray)<-types
  colnames(permanova_matrix_bray)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  permanova_pvalues_bray<-data.frame(permanova_matrix_bray[2:4,1:3])
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test for Homogeneity of Variances (Table 5.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Test for homogeneity of variances for euclidean
bd<-betadisper(euclidean,as.factor(type))
dispersion_anova_euclidean<-anova(bd)

#Test for homogeneity of variances for jaccard
bd<-betadisper(jaccard,as.factor(type))
dispersion_anova_jaccard<-anova(bd)

#Test for homogeneity of variances for bray-curtis
bd<-betadisper(bray,as.factor(type))
dispersion_anova_bray<-anova(bd)


#Table 5.13 (use adjusted p-values to generate Table 5.4 manually)
posthocs<-rbind(posthoc_permanova_euclidean,posthoc_permanova_jaccard,posthoc_permanova_bray)
write.csv(posthocs,'posthocs_genera.csv')



