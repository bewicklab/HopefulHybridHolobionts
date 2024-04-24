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
unifrac<-UniFrac(rZOTU,weighted=FALSE)
wunifrac<-UniFrac(rZOTU,weighted=TRUE)

#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
euclidean_matrix<-as.matrix(euclidean,labels=TRUE)
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)

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

#                 Perform PERMANOVA Tests on Unifrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_unifrac<-list(Data=t(otu_table(rZOTU)),D=unifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_unifrac=PERMANOVA(inputpermanova_unifrac, factor(type),nperm=10000)

#If PERMANOVA by treatment group is significant, figure out which groups differ from which other groups
if (permanova_by_group_unifrac$pvalue<0.05 && length(unique(type))>2){

  #Perform a post-hoc pairwise test between each pair of treatment groups
  posthoc_permanova_unifrac<-pairwise.adonis(unifrac_matrix,as.factor(type),p.adjust.m='BH')


}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on weighted Unifrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_wunifrac<-list(Data=t(otu_table(rZOTU)),D=wunifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_wunifrac=PERMANOVA(inputpermanova_wunifrac, factor(type),nperm=10000)

#If PERMANOVA by treatment group is significant, figure out which groups differ from which other groups
if (permanova_by_group_wunifrac$pvalue<0.05 && length(unique(type))>2){

  #Perform a post-hoc pairwise test between each pair of treatment groups
  posthoc_permanova_wunifrac<-pairwise.adonis(wunifrac_matrix,as.factor(type),p.adjust.m='BH')


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

#                 Perform ANOSIM Tests on unifrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Perform anosim comparing treatment groups
anosim_by_group_unifrac<-anosim(unifrac_matrix, type, permutations = 10000)

#If anosim by treatment group is significant, figure out which groups differ from which other groups
#This is done by performing anosim individually on each set of pairs
if (anosim_by_group_unifrac$signif<0.05){

  anosim_pairwise_ps<-c()   #This will be a list of pairwise p-values
  anosim_pairwise_ts<-c()   #This will be a list of pairwise test statistics
  first_type_ps<-c()        #This is a list of the first treatment group for each comparison
  second_type_ps<-c()       #This is a list of the second treatment group for each comparison

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  anosim_matrix_unifrac<-matrix(data='',nrow=length(types),ncol=length(types))

  kv<-c()   #This will be a list of the row number for each p-value in the list of comparisons
  jv<-c()   #This will be a list of the column number for each p-value in the list of comparisons

  #For each treatment group
  for (k in 1:length(types)){
    #Compared to every other UNIQUE treatment group
    for (j in 1:k){
      if (k!=j){
        #Find the distances associated with the two focal treatment groups being compared
        temp_matrix<-unifrac_matrix[c(which(type==types[k]),which(type==types[j])),c(which(type==types[k]),which(type==types[j]))]
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
    anosim_matrix_unifrac[kv[k],jv[k]]<-round(corrected_pairwise_ps[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(anosim_matrix_unifrac)<-types
  colnames(anosim_matrix_unifrac)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  anosim_pvalues_unifrac<-data.frame(anosim_matrix_unifrac)

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
  posthoc_anosim_unifrac<-data.frame(anosim_temp)

  #Assign the column names of the above dataframe based on the information they contain
  colnames(posthoc_anosim_unifrac)<-c('ANOSIM statistic','p-value','corrected p-value','significance')
  #Assign the row names of the above dataframe based on the pairwise comparison being made in each row
  rownames(posthoc_anosim_unifrac)<-compare_temp
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform ANOSIM Tests on weighted unifrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Perform anosim comparing treatment groups
anosim_by_group_wunifrac<-anosim(wunifrac_matrix, type, permutations = 10000)

#If anosim by treatment group is significant, figure out which groups differ from which other groups
#This is done by performing anosim individually on each set of pairs
if (anosim_by_group_wunifrac$signif<0.05){

  anosim_pairwise_ps<-c()   #This will be a list of pairwise p-values
  anosim_pairwise_ts<-c()   #This will be a list of pairwise test statistics
  first_type_ps<-c()        #This is a list of the first treatment group for each comparison
  second_type_ps<-c()       #This is a list of the second treatment group for each comparison

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  anosim_matrix_wunifrac<-matrix(data='',nrow=length(types),ncol=length(types))

  kv<-c()   #This will be a list of the row number for each p-value in the list of comparisons
  jv<-c()   #This will be a list of the column number for each p-value in the list of comparisons

  #For each treatment group
  for (k in 1:length(types)){
    #Compared to every other UNIQUE treatment group
    for (j in 1:k){
      if (k!=j){
        #Find the distances associated with the two focal treatment groups being compared
        temp_matrix<-wunifrac_matrix[c(which(type==types[k]),which(type==types[j])),c(which(type==types[k]),which(type==types[j]))]
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
    anosim_matrix_wunifrac[kv[k],jv[k]]<-round(corrected_pairwise_ps[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(anosim_matrix_wunifrac)<-types
  colnames(anosim_matrix_wunifrac)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  anosim_pvalues_wunifrac<-data.frame(anosim_matrix_wunifrac)

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
  posthoc_anosim_wunifrac<-data.frame(anosim_temp)

  #Assign the column names of the above dataframe based on the information they contain
  colnames(posthoc_anosim_wunifrac)<-c('ANOSIM statistic','p-value','corrected p-value','significance')
  #Assign the row names of the above dataframe based on the pairwise comparison being made in each row
  rownames(posthoc_anosim_wunifrac)<-compare_temp
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

#Test for homogeneity of variances for unifrac
bd<-betadisper(unifrac,as.factor(type))
dispersion_anova_unifrac<-anova(bd)

#Test for homogeneity of variances for weighted unifrac
bd<-betadisper(wunifrac,as.factor(type))
dispersion_anova_wunifrac<-anova(bd)



