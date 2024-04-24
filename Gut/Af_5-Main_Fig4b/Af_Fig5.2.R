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
library('HybridMicrobiomes')

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
Z_physeq<-phyloseq(OTU_biom,TAX,SAMP)

#Determine which sample is which lizard species
species_list<-sample_data(Z_physeq)$species
species_no<-rep(1,length(species_list))
species_no[which(species_list=='neomexicanus')]<-2
species_no[which(species_list=='marmoratus')]<-3

#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq_SNW<-prune_samples(sample_data(Z_physeq)$speciesxsite %in% c('SNW_neo','SNW_marm','SBluG_ino'),Z_physeq)

#Determine which sample is which lizard species
species_list_SNW<-sample_data(Z_physeq_SNW)$species
species_no_SNW<-rep(1,length(species_list_SNW))
species_no_SNW[which(species_list_SNW=='neomexicanus')]<-2
species_no_SNW[which(species_list_SNW=='marmoratus')]<-3


#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq_SBluG<-prune_samples(sample_data(Z_physeq)$speciesxsite %in% c('SBluG_neo','SNW_marm','SBluG_ino'),Z_physeq)

#Determine which sample is which lizard species
species_list_SBluG<-sample_data(Z_physeq_SBluG)$species
species_no_SBluG<-rep(1,length(species_list_SBluG))
species_no_SBluG[which(species_list_SBluG=='neomexicanus')]<-2
species_no_SBluG[which(species_list_SBluG=='marmoratus')]<-3


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 PCoA Euclidean

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


triangle_euclidean_all_2<-TriangleHbootstrap(Z_physeq,species_no,50,12,dimensions=2,seed=1,ordination='PCoA')
triangle_euclidean_all_2_SNW<-TriangleHbootstrap(Z_physeq_SNW,species_no_SNW,50,12,dimensions=2,seed=1,ordination='PCoA')
triangle_euclidean_all_2_SBluG<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,50,12,dimensions=2,seed=1,ordination='PCoA')

triangle_euclidean_all_2_null1<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1,ordination='PCoA')
triangle_euclidean_all_2_null6<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1,null_model=6,ordination='PCoA')
triangle_euclidean_all_2_null7<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1,null_model=7,ordination='PCoA')

# Fig. 5.2a
TriangleHplot(triangle_euclidean_all_2,col='purple')
TriangleHplot(triangle_euclidean_all_2_SNW,col='lightslateblue',addplot=TRUE)
TriangleHplot(triangle_euclidean_all_2_SBluG,col='magenta',addplot=TRUE)
TriangleHplot(triangle_euclidean_all_2_null1,col='mediumpurple1',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_euclidean_all_2_null6,col='red',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_euclidean_all_2_null7,col='blue',addplot=TRUE,centroid_shape='triangle')

triangle_euclidean_all_20<-TriangleHbootstrap(Z_physeq,species_no,50,12,dimensions=20,seed=1,ordination='PCoA')
triangle_euclidean_all_20_SNW<-TriangleHbootstrap(Z_physeq_SNW,species_no_SNW,50,12,dimensions=20,seed=1,ordination='PCoA')
triangle_euclidean_all_20_SBluG<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,50,12,dimensions=20,seed=1,ordination='PCoA')

triangle_euclidean_all_20_null1<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=20,seed=1,ordination='PCoA')
triangle_euclidean_all_20_null6<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=20,seed=1,null_model=6,ordination='PCoA')
triangle_euclidean_all_20_null7<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=20,seed=1,null_model=7,ordination='PCoA')

# Fig. 5.2d
TriangleHplot(triangle_euclidean_all_20,col='purple')
TriangleHplot(triangle_euclidean_all_20_SNW,col='lightslateblue',addplot=TRUE)
TriangleHplot(triangle_euclidean_all_20_SBluG,col='magenta',addplot=TRUE)
TriangleHplot(triangle_euclidean_all_20_null1,col='mediumpurple1',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_euclidean_all_20_null6,col='red',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_euclidean_all_20_null7,col='blue',addplot=TRUE,centroid_shape='triangle')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 PCoA Jaccard

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


triangle_jaccard_all_2<-TriangleHbootstrap(Z_physeq,species_no,50,12,dimensions=2,seed=1,metric='jaccard',ordination='PCoA')
triangle_jaccard_all_2_SNW<-TriangleHbootstrap(Z_physeq_SNW,species_no_SNW,50,12,dimensions=2,seed=1,metric='jaccard',ordination='PCoA')
triangle_jaccard_all_2_SBluG<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,50,12,dimensions=2,seed=1,metric='jaccard',ordination='PCoA')

triangle_jaccard_all_2_null10<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1, null_model = 10,metric='jaccard',ordination='PCoA')
triangle_jaccard_all_2_null6<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1,null_model=6,metric='jaccard',ordination='PCoA')
triangle_jaccard_all_2_null7<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1,null_model=7,metric='jaccard',ordination='PCoA')

# Fig. 5.2b
TriangleHplot(triangle_jaccard_all_2,col='purple')
TriangleHplot(triangle_jaccard_all_2_SNW,col='lightslateblue',addplot=TRUE)
TriangleHplot(triangle_jaccard_all_2_SBluG,col='magenta',addplot=TRUE)
TriangleHplot(triangle_jaccard_all_2_null10,col='mediumpurple1',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_jaccard_all_2_null6,col='red',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_jaccard_all_2_null7,col='blue',addplot=TRUE,centroid_shape='triangle')

triangle_jaccard_all_20<-TriangleHbootstrap(Z_physeq,species_no,50,12,dimensions=20,seed=1,metric='jaccard',ordination='PCoA')
triangle_jaccard_all_20_SNW<-TriangleHbootstrap(Z_physeq_SNW,species_no_SNW,50,12,dimensions=20,seed=1,metric='jaccard',ordination='PCoA')
triangle_jaccard_all_20_SBluG<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,50,12,dimensions=20,seed=1,metric='jaccard',ordination='PCoA')

triangle_jaccard_all_20_null10<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=20,seed=1,null_model=10,metric='jaccard',ordination='PCoA')
triangle_jaccard_all_20_null6<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=20,seed=1,null_model=6,metric='jaccard',ordination='PCoA')
triangle_jaccard_all_20_null7<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=20,seed=1,null_model=7,metric='jaccard',ordination='PCoA')

# Fig. 5.2e
TriangleHplot(triangle_jaccard_all_20,col='purple')
TriangleHplot(triangle_jaccard_all_20_SNW,col='lightslateblue',addplot=TRUE)
TriangleHplot(triangle_jaccard_all_20_SBluG,col='magenta',addplot=TRUE)
TriangleHplot(triangle_jaccard_all_20_null10,col='mediumpurple1',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_jaccard_all_20_null6,col='red',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_jaccard_all_20_null7,col='blue',addplot=TRUE,centroid_shape='triangle')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 PCoA Bray

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


triangle_bray_all_2<-TriangleHbootstrap(Z_physeq,species_no,50,12,dimensions=2,seed=1,metric='bray',ordination='PCoA')
triangle_bray_all_2_SNW<-TriangleHbootstrap(Z_physeq_SNW,species_no_SNW,50,12,dimensions=2,seed=1,metric='bray',ordination='PCoA')
triangle_bray_all_2_SBluG<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,50,12,dimensions=2,seed=1,metric='bray',ordination='PCoA')

triangle_bray_all_2_null1<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1,metric='bray',ordination='PCoA')
triangle_bray_all_2_null6<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1,null_model=6,metric='bray',ordination='PCoA')
triangle_bray_all_2_null7<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=2,seed=1,null_model=7,metric='bray',ordination='PCoA')

# Fig. 5.2c
TriangleHplot(triangle_bray_all_2,col='purple',xlow = -2,xhigh=2,yhigh=3)
TriangleHplot(triangle_bray_all_2_SNW,col='lightslateblue',addplot=TRUE)
TriangleHplot(triangle_bray_all_2_SBluG,col='magenta',addplot=TRUE)
TriangleHplot(triangle_bray_all_2_null1,col='mediumpurple1',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_bray_all_2_null6,col='red',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_bray_all_2_null7,col='blue',addplot=TRUE,centroid_shape='triangle')

triangle_bray_all_20<-TriangleHbootstrap(Z_physeq,species_no,50,12,dimensions=15,seed=2,metric='bray',ordination='PCoA')
triangle_bray_all_20_SNW<-TriangleHbootstrap(Z_physeq_SNW,species_no_SNW,50,12,dimensions=15,seed=2,metric='bray',ordination='PCoA')
triangle_bray_all_20_SBluG<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,50,12,dimensions=15,seed=2,metric='bray',ordination='PCoA')

triangle_bray_all_20_null1<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=15,seed=2,metric='bray',ordination='PCoA')
triangle_bray_all_20_null6<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=15,seed=2,null_model=6,metric='bray',ordination='PCoA')
triangle_bray_all_20_null7<-TriangleHnull(Z_physeq,species_no,50,12,dimensions=15,seed=2,null_model=7,metric='bray',ordination='PCoA')

# Fig. 5.2f
TriangleHplot(triangle_bray_all_20,col='purple')
TriangleHplot(triangle_bray_all_20_SNW,col='lightslateblue',addplot=TRUE)
TriangleHplot(triangle_bray_all_20_SBluG,col='magenta',addplot=TRUE)
TriangleHplot(triangle_bray_all_20_null1,col='mediumpurple1',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_bray_all_20_null6,col='red',addplot=TRUE,centroid_shape='triangle')
TriangleHplot(triangle_bray_all_20_null7,col='blue',addplot=TRUE,centroid_shape='triangle')







