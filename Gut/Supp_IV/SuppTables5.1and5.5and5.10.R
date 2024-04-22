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

  #Write the p-values (corrected for multiple comparisons) from the pairwise.adonis function as a concise data.frame, rather than a list of comparisons
  #(This can be useful for quickly summarizing/interpreting results)
  a1<-c()     #This will be a list of the first treatment group for each pairwise comparison
  a2<-c()     #This will be a list of the second treatment group for each pairwise comparison
  #For each pairwise comparison outputted by the pairwise.adonis function...
  for (k in 1:length(posthoc_permanova_unifrac[,1])){
    #Find the first and second treatment groups being compared
    ff<-str_split(posthoc_permanova_unifrac[,1][k],' vs ')
    a1<-c(a1,ff[[1]][1])
    a2<-c(a2,ff[[1]][2])
  }

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  permanova_matrix_unifrac<-matrix(data='',nrow=length(types),ncol=length(types))
  for (k in 1:length(posthoc_permanova_unifrac[,1])){
    #For each comparison in the list outputted from the pairwise.adonis matrix...
    hit1<-which(types==a1[k])
    hit2<-which(types==a2[k])
    #Record the adjusted p-value in the matrix
    permanova_matrix_unifrac[max(hit1,hit2),min(hit1,hit2)]<-round(posthoc_permanova_unifrac$p.adjusted[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(permanova_matrix_unifrac)<-types
  colnames(permanova_matrix_unifrac)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  permanova_pvalues_unifrac<-data.frame(permanova_matrix_unifrac[2:4,1:3])
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

  #Write the p-values (corrected for multiple comparisons) from the pairwise.adonis function as a concise data.frame, rather than a list of comparisons
  #(This can be useful for quickly summarizing/interpreting results)
  a1<-c()     #This will be a list of the first treatment group for each pairwise comparison
  a2<-c()     #This will be a list of the second treatment group for each pairwise comparison
  #For each pairwise comparison outputted by the pairwise.adonis function...
  for (k in 1:length(posthoc_permanova_wunifrac[,1])){
    #Find the first and second treatment groups being compared
    ff<-str_split(posthoc_permanova_wunifrac[,1][k],' vs ')
    a1<-c(a1,ff[[1]][1])
    a2<-c(a2,ff[[1]][2])
  }

  #Define a matrix where the rows and columns are the treatment groups and the matrix elements are the p-values for the rowxcolumn comparison
  permanova_matrix_wunifrac<-matrix(data='',nrow=length(types),ncol=length(types))
  for (k in 1:length(posthoc_permanova_wunifrac[,1])){
    #For each comparison in the list outputted from the pairwise.adonis matrix...
    hit1<-which(types==a1[k])
    hit2<-which(types==a2[k])
    #Record the adjusted p-value in the matrix
    permanova_matrix_wunifrac[max(hit1,hit2),min(hit1,hit2)]<-round(posthoc_permanova_wunifrac$p.adjusted[k],4)
  }
  #The rows and columns of this matrix are the treatment groups
  rownames(permanova_matrix_wunifrac)<-types
  colnames(permanova_matrix_wunifrac)<-types
  #Convert this into a dataframe, removing the first row and last column (this is just a lower matrix, so we don't look at k=j or k>j since these are identical/redundant respectively)
  permanova_pvalues_wunifrac<-data.frame(permanova_matrix_wunifrac[2:4,1:3])
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

#Test for homogeneity of variances for unifrac
bd<-betadisper(unifrac,as.factor(type))
dispersion_anova_unifrac<-anova(bd)

#Test for homogeneity of variances for weighted unifrac
bd<-betadisper(wunifrac,as.factor(type))
dispersion_anova_wunifrac<-anova(bd)

#Table 5.10 (use adjusted p-values to generate Table 5.1 manually)
posthocs<-rbind(posthoc_permanova_euclidean,posthoc_permanova_jaccard,posthoc_permanova_bray,posthoc_permanova_unifrac,posthoc_permanova_wunifrac)
write.csv(posthocs,'posthocs_ASV.csv')


