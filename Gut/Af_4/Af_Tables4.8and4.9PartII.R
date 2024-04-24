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
type<-sample_data(rZOTU)$species

types<-unique(type)

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


}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform ANOSIM Tests on Euclidean Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Perform anosim comparing treatment groups
anosim_by_group_euclidean<-anosim(euclidean_matrix, type, permutations = 10000)

#If anosim by treatment group is significant, figure out which groups differ from which other groups
#This is done by performing anosim individually on each set of pairs
if (anosim_by_group_euclidean$signif<0.05){

  anosim_pairwise_ps<-c()   #This will be a list of pairwise p-values
  anosim_pairwise_ts<-c()   #This will be a list of pairwise test statistics
  first_type_ps<-c()        #This is a list of the first treatment group for each comparison
  second_type_ps<-c()       #This is a list of the second treatment group for each comparison

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  anosim_matrix_euclidean<-matrix(data='',nrow=length(types),ncol=length(types))

  kv<-c()   #This will be a list of the row number for each p-value in the list of comparisons
  jv<-c()   #This will be a list of the column number for each p-value in the list of comparisons

  #For each treatment group
  for (k in 1:length(types)){
    #Compared to every other UNIQUE treatment group
    for (j in 1:k){
      if (k!=j){
        #Find the distances associated with the two focal treatment groups being compared
        temp_matrix<-euclidean_matrix[c(which(type==types[k]),which(type==types[j])),c(which(type==types[k]),which(type==types[j]))]
        #Define the treatment group for each microbial community in the paired down distance matrix from the line above
        temp_type<-c(rep(types[k],length(which(type==types[k]))),rep(types[j],length(which(type==types[j]))))
        #Perform anosim on the paired treatment groups
        anosim_temp<-anosim(temp_matrix,temp_type,permutation = 10000)
        #Add the p-value for the anosim pairwise comparison to the list of anosim pairwise p-values
        anosim_pairwise_ps<-c(anosim_pairwise_ps,anosim_temp$signif)
        #Add the test statistic for the anosim pairwise comparison to the list of anosim pairwise test statistics
        anosim_pairwise_ts<-c(anosim_pairwise_ts,anosim_temp$statistic)
        #Record the type of the first treatment group for this comparison
        first_type_ps<-c(first_type_ps,types[k])
        #Record the type of the second treatment group for this comparison
        second_type_ps<-c(second_type_ps,types[j])
        #Record the row number for this comparison (facilitates entering the corrected p-values into a matrix in the line below)
        kv<-c(kv,k)
        #Record the column number for this comparison
        jv<-c(jv,j)
      }
    }
  }

  #Correct the anosim p-values for multiple comaprisons
  corrected_pairwise_ps<-p.adjust(anosim_pairwise_ps,method='BH',n=length(anosim_pairwise_ps))

  #Enter the corrected p-values into a matrix using the lists of the row and column numbers for each p-value in the corrected p-value list
  for (k in 1:length(corrected_pairwise_ps)){
    anosim_matrix_euclidean[kv[k],jv[k]]<-round(corrected_pairwise_ps[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(anosim_matrix_euclidean)<-types
  colnames(anosim_matrix_euclidean)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  anosim_pvalues_euclidean<-data.frame(anosim_matrix_euclidean)

  #Mark the various significance levels of the pairwise anosim comparisons using the 'star' system
  #<0.001: '***', <0.01:'**', <0.05:'*',<0.1:'.'
  anosim_temp<-c()
  compare_temp<-c()
  for (k in 1:length(corrected_pairwise_ps)){
    if (corrected_pairwise_ps[k]>0.1){
      significance<-''
    }
    else if (corrected_pairwise_ps[k]>0.05){
      significance<-'.'
    }
    else if (corrected_pairwise_ps[k]>0.01){
      significance<-'*'
    }
    else if (corrected_pairwise_ps[k]>0.001){
      significance<-'**'
    }
    else {
      significance<-'***'
    }

    #Make a list of the pairwise comparisons for each p-value and each test statistic
    compare_temp<-c(compare_temp,paste0(first_type_ps[k],'+',second_type_ps[k]))
    #add the p-values, test statistics and significance for each of the pairwise comparisons into a master list
    anosim_temp<-rbind(anosim_temp,c(anosim_pairwise_ts[k],anosim_pairwise_ps[k],corrected_pairwise_ps[k],significance))
  }

  #Turn the master list of information into a dataframe
  posthoc_anosim_euclidean<-data.frame(anosim_temp)

  #Assign the column names of the above dataframe based on the information they contain
  colnames(posthoc_anosim_euclidean)<-c('ANOSIM statistic','p-value','corrected p-value','significance')
  #Assign the row names of the above dataframe based on the pairwise comparison being made in each row
  rownames(posthoc_anosim_euclidean)<-compare_temp
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform ANOSIM Tests on Jaccard Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Perform anosim comparing treatment groups
anosim_by_group_jaccard<-anosim(jaccard_matrix, type, permutations = 10000)

#If anosim by treatment group is significant, figure out which groups differ from which other groups
#This is done by performing anosim individually on each set of pairs
if (anosim_by_group_jaccard$signif<0.05){

  anosim_pairwise_ps<-c()   #This will be a list of pairwise p-values
  anosim_pairwise_ts<-c()   #This will be a list of pairwise test statistics
  first_type_ps<-c()        #This is a list of the first treatment group for each comparison
  second_type_ps<-c()       #This is a list of the second treatment group for each comparison

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  anosim_matrix_jaccard<-matrix(data='',nrow=length(types),ncol=length(types))

  kv<-c()   #This will be a list of the row number for each p-value in the list of comparisons
  jv<-c()   #This will be a list of the column number for each p-value in the list of comparisons

  #For each treatment group
  for (k in 1:length(types)){
    #Compared to every other UNIQUE treatment group
    for (j in 1:k){
      if (k!=j){
        #Find the distances associated with the two focal treatment groups being compared
        temp_matrix<-jaccard_matrix[c(which(type==types[k]),which(type==types[j])),c(which(type==types[k]),which(type==types[j]))]
        #Define the treatment group for each microbial community in the paired down distance matrix from the line above
        temp_type<-c(rep(types[k],length(which(type==types[k]))),rep(types[j],length(which(type==types[j]))))
        #Perform anosim on the paired treatment groups
        anosim_temp<-anosim(temp_matrix,temp_type,permutation = 10000)
        #Add the p-value for the anosim pairwise comparison to the list of anosim pairwise p-values
        anosim_pairwise_ps<-c(anosim_pairwise_ps,anosim_temp$signif)
        #Add the test statistic for the anosim pairwise comparison to the list of anosim pairwise test statistics
        anosim_pairwise_ts<-c(anosim_pairwise_ts,anosim_temp$statistic)
        #Record the type of the first treatment group for this comparison
        first_type_ps<-c(first_type_ps,types[k])
        #Record the type of the second treatment group for this comparison
        second_type_ps<-c(second_type_ps,types[j])
        #Record the row number for this comparison (facilitates entering the corrected p-values into a matrix in the line below)
        kv<-c(kv,k)
        #Record the column number for this comparison
        jv<-c(jv,j)
      }
    }
  }

  #Correct the anosim p-values for multiple comaprisons
  corrected_pairwise_ps<-p.adjust(anosim_pairwise_ps,method='BH',n=length(anosim_pairwise_ps))

  #Enter the corrected p-values into a matrix using the lists of the row and column numbers for each p-value in the corrected p-value list
  for (k in 1:length(corrected_pairwise_ps)){
    anosim_matrix_jaccard[kv[k],jv[k]]<-round(corrected_pairwise_ps[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(anosim_matrix_jaccard)<-types
  colnames(anosim_matrix_jaccard)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  anosim_pvalues_jaccard<-data.frame(anosim_matrix_jaccard)

  #Mark the various significance levels of the pairwise anosim comparisons using the 'star' system
  #<0.001: '***', <0.01:'**', <0.05:'*',<0.1:'.'
  anosim_temp<-c()
  compare_temp<-c()
  for (k in 1:length(corrected_pairwise_ps)){
    if (corrected_pairwise_ps[k]>0.1){
      significance<-''
    }
    else if (corrected_pairwise_ps[k]>0.05){
      significance<-'.'
    }
    else if (corrected_pairwise_ps[k]>0.01){
      significance<-'*'
    }
    else if (corrected_pairwise_ps[k]>0.001){
      significance<-'**'
    }
    else {
      significance<-'***'
    }

    #Make a list of the pairwise comparisons for each p-value and each test statistic
    compare_temp<-c(compare_temp,paste0(first_type_ps[k],'+',second_type_ps[k]))
    #add the p-values, test statistics and significance for each of the pairwise comparisons into a master list
    anosim_temp<-rbind(anosim_temp,c(anosim_pairwise_ts[k],anosim_pairwise_ps[k],corrected_pairwise_ps[k],significance))
  }

  #Turn the master list of information into a dataframe
  posthoc_anosim_jaccard<-data.frame(anosim_temp)

  #Assign the column names of the above dataframe based on the information they contain
  colnames(posthoc_anosim_jaccard)<-c('ANOSIM statistic','p-value','corrected p-value','significance')
  #Assign the row names of the above dataframe based on the pairwise comparison being made in each row
  rownames(posthoc_anosim_jaccard)<-compare_temp
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform ANOSIM Tests on Bray Curtis Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Perform anosim comparing treatment groups
anosim_by_group_bray<-anosim(bray_matrix, type, permutations = 10000)

#If anosim by treatment group is significant, figure out which groups differ from which other groups
#This is done by performing anosim individually on each set of pairs
if (anosim_by_group_bray$signif<0.05){

  anosim_pairwise_ps<-c()   #This will be a list of pairwise p-values
  anosim_pairwise_ts<-c()   #This will be a list of pairwise test statistics
  first_type_ps<-c()        #This is a list of the first treatment group for each comparison
  second_type_ps<-c()       #This is a list of the second treatment group for each comparison

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  anosim_matrix_bray<-matrix(data='',nrow=length(types),ncol=length(types))

  kv<-c()   #This will be a list of the row number for each p-value in the list of comparisons
  jv<-c()   #This will be a list of the column number for each p-value in the list of comparisons

  #For each treatment group
  for (k in 1:length(types)){
    #Compared to every other UNIQUE treatment group
    for (j in 1:k){
      if (k!=j){
        #Find the distances associated with the two focal treatment groups being compared
        temp_matrix<-bray_matrix[c(which(type==types[k]),which(type==types[j])),c(which(type==types[k]),which(type==types[j]))]
        #Define the treatment group for each microbial community in the paired down distance matrix from the line above
        temp_type<-c(rep(types[k],length(which(type==types[k]))),rep(types[j],length(which(type==types[j]))))
        #Perform anosim on the paired treatment groups
        anosim_temp<-anosim(temp_matrix,temp_type,permutation = 10000)
        #Add the p-value for the anosim pairwise comparison to the list of anosim pairwise p-values
        anosim_pairwise_ps<-c(anosim_pairwise_ps,anosim_temp$signif)
        #Add the test statistic for the anosim pairwise comparison to the list of anosim pairwise test statistics
        anosim_pairwise_ts<-c(anosim_pairwise_ts,anosim_temp$statistic)
        #Record the type of the first treatment group for this comparison
        first_type_ps<-c(first_type_ps,types[k])
        #Record the type of the second treatment group for this comparison
        second_type_ps<-c(second_type_ps,types[j])
        #Record the row number for this comparison (facilitates entering the corrected p-values into a matrix in the line below)
        kv<-c(kv,k)
        #Record the column number for this comparison
        jv<-c(jv,j)
      }
    }
  }

  #Correct the anosim p-values for multiple comaprisons
  corrected_pairwise_ps<-p.adjust(anosim_pairwise_ps,method='BH',n=length(anosim_pairwise_ps))

  #Enter the corrected p-values into a matrix using the lists of the row and column numbers for each p-value in the corrected p-value list
  for (k in 1:length(corrected_pairwise_ps)){
    anosim_matrix_bray[kv[k],jv[k]]<-round(corrected_pairwise_ps[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(anosim_matrix_bray)<-types
  colnames(anosim_matrix_bray)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  anosim_pvalues_bray<-data.frame(anosim_matrix_bray)

  #Mark the various significance levels of the pairwise anosim comparisons using the 'star' system
  #<0.001: '***', <0.01:'**', <0.05:'*',<0.1:'.'
  anosim_temp<-c()
  compare_temp<-c()
  for (k in 1:length(corrected_pairwise_ps)){
    if (corrected_pairwise_ps[k]>0.1){
      significance<-''
    }
    else if (corrected_pairwise_ps[k]>0.05){
      significance<-'.'
    }
    else if (corrected_pairwise_ps[k]>0.01){
      significance<-'*'
    }
    else if (corrected_pairwise_ps[k]>0.001){
      significance<-'**'
    }
    else {
      significance<-'***'
    }

    #Make a list of the pairwise comparisons for each p-value and each test statistic
    compare_temp<-c(compare_temp,paste0(first_type_ps[k],'+',second_type_ps[k]))
    #add the p-values, test statistics and significance for each of the pairwise comparisons into a master list
    anosim_temp<-rbind(anosim_temp,c(anosim_pairwise_ts[k],anosim_pairwise_ps[k],corrected_pairwise_ps[k],significance))
  }

  #Turn the master list of information into a dataframe
  posthoc_anosim_bray<-data.frame(anosim_temp)

  #Assign the column names of the above dataframe based on the information they contain
  colnames(posthoc_anosim_bray)<-c('ANOSIM statistic','p-value','corrected p-value','significance')
  #Assign the row names of the above dataframe based on the pairwise comparison being made in each row
  rownames(posthoc_anosim_bray)<-compare_temp
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test for Homogeneity of Variances (another measure of beta diversity)

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




