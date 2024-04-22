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
type<-sample_data(rZOTU)$hybrid

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

#Use adjusted p-values from tests below to generate Table 5.6

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Euclidean Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_euclidean<-list(Data=t(otu_table(rZOTU)),D=euclidean_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_euclidean=PERMANOVA(inputpermanova_euclidean, factor(type),nperm=10000)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Jaccard Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_jaccard<-list(Data=t(otu_table(rZOTU)),D=jaccard_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_jaccard=PERMANOVA(inputpermanova_jaccard, factor(type),nperm=10000)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Bray Curtis Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_bray<-list(Data=t(otu_table(rZOTU)),D=bray_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_bray=PERMANOVA(inputpermanova_bray, factor(type),nperm=10000)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Unifrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_unifrac<-list(Data=t(otu_table(rZOTU)),D=unifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_unifrac=PERMANOVA(inputpermanova_unifrac, factor(type),nperm=10000)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on weighted Unifrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_wunifrac<-list(Data=t(otu_table(rZOTU)),D=wunifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_wunifrac=PERMANOVA(inputpermanova_wunifrac, factor(type),nperm=10000)


permanova_pvalues<-c(permanova_by_group_euclidean$pvalue,permanova_by_group_jaccard$pvalue,permanova_by_group_bray$pvalue,permanova_by_group_unifrac$pvalue,permanova_by_group_wunifrac$pvalue)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test for Homogeneity of Variances (Table 5.7)

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



