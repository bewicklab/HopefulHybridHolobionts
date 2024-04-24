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
library('ape')
library('HybridMicrobiomes')

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


triangle_euclidean_all_2<-TriangleHbootstrap(Z_physeq,species_no,15,12,dimensions=2,seed=1,ordination='PCoA')

triangle_euclidean_all_2_null1<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,ordination='PCoA')
triangle_euclidean_all_2_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=6, ordination='PCoA')
triangle_euclidean_all_2_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=7, ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_euclidean_all_2[,1:2],triangle_euclidean_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_euclidean_all_2[,1:2],triangle_euclidean_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_euclidean_all_2[,1:2],triangle_euclidean_all_2_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_euclidean_all_2[,1:2],triangle_euclidean_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_euclidean_all_2[,1:2],triangle_euclidean_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_euclidean_all_2[,1:2],triangle_euclidean_all_2_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-c(median(triangle_euclidean_all_2[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_euclidean_all_2[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value)

triangle_euclidean_all_20<-TriangleHbootstrap(Z_physeq,species_no,15,12,dimensions=20,seed=1,ordination='PCoA')

triangle_euclidean_all_20_null1<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=1,ordination='PCoA')
triangle_euclidean_all_20_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=1,null_model=6,ordination='PCoA')
triangle_euclidean_all_20_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=1,null_model=7,ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_euclidean_all_20[,1:2],triangle_euclidean_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_euclidean_all_20[,1:2],triangle_euclidean_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_euclidean_all_20[,1:2],triangle_euclidean_all_20_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_euclidean_all_20[,1:2],triangle_euclidean_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_euclidean_all_20[,1:2],triangle_euclidean_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_euclidean_all_20[,1:2],triangle_euclidean_all_20_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_euclidean_all_20[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_euclidean_all_20[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 PCoA Jaccard

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


triangle_jaccard_all_2<-TriangleHbootstrap(Z_physeq,species_no,15,12,dimensions=2,seed=1,dist='jaccard',ordination='PCoA')

triangle_jaccard_all_2_null10<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=10,dist='jaccard',ordination='PCoA')
triangle_jaccard_all_2_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=6,dist='jaccard',ordination='PCoA')
triangle_jaccard_all_2_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=7,dist='jaccard',ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_jaccard_all_2[,1:2],triangle_jaccard_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_jaccard_all_2[,1:2],triangle_jaccard_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_jaccard_all_2[,1:2],triangle_jaccard_all_2_null10[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_jaccard_all_2[,1:2],triangle_jaccard_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_jaccard_all_2[,1:2],triangle_jaccard_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_jaccard_all_2[,1:2],triangle_jaccard_all_2_null10[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_jaccard_all_2[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_jaccard_all_2[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))

triangle_jaccard_all_20<-TriangleHbootstrap(Z_physeq,species_no,15,12,dimensions=20,seed=1,dist='jaccard',ordination='PCoA')

triangle_jaccard_all_20_null10<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=1,null_model=10,dist='jaccard',ordination='PCoA')
triangle_jaccard_all_20_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=1,null_model=6,dist='jaccard',ordination='PCoA')
triangle_jaccard_all_20_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=1,null_model=7,dist='jaccard',ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_jaccard_all_20[,1:2],triangle_jaccard_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_jaccard_all_20[,1:2],triangle_jaccard_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_jaccard_all_20[,1:2],triangle_jaccard_all_20_null10[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_jaccard_all_20[,1:2],triangle_jaccard_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_jaccard_all_20[,1:2],triangle_jaccard_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_jaccard_all_20[,1:2],triangle_jaccard_all_20_null10[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_jaccard_all_20[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_jaccard_all_20[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 PCoA Bray

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


triangle_bray_all_2<-TriangleHbootstrap(Z_physeq,species_no,15,12,dimensions=2,seed=1,dist='bray',ordination='PCoA')

triangle_bray_all_2_null1<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,dist='bray',ordination='PCoA')
triangle_bray_all_2_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=6,dist='bray',ordination='PCoA')
triangle_bray_all_2_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=7,dist='bray',ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_bray_all_2[,1:2],triangle_bray_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_bray_all_2[,1:2],triangle_bray_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_bray_all_2[,1:2],triangle_bray_all_2_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_bray_all_2[,1:2],triangle_bray_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_bray_all_2[,1:2],triangle_bray_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_bray_all_2[,1:2],triangle_bray_all_2_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_bray_all_2[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_bray_all_2[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))

triangle_bray_all_20<-TriangleHbootstrap(Z_physeq,species_no,15,12,dimensions=15,seed=2,dist='bray',ordination='PCoA')

triangle_bray_all_20_null1<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=15,seed=2,dist='bray',ordination='PCoA')
triangle_bray_all_20_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=15,seed=2,null_model=6,dist='bray',ordination='PCoA')
triangle_bray_all_20_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=15,seed=2,null_model=7,dist='bray',ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_bray_all_20[,1:2],triangle_bray_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_bray_all_20[,1:2],triangle_bray_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_bray_all_20[,1:2],triangle_bray_all_20_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_bray_all_20[,1:2],triangle_bray_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_bray_all_20[,1:2],triangle_bray_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_bray_all_20[,1:2],triangle_bray_all_20_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_bray_all_20[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_bray_all_20[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))




tests_df<-data.frame(tests)
rownames(tests_df)<-c('euclidean_2','euclidean_20','jaccard_2','jaccard_20','bray_2','bray_20')
colnames(tests_df)<-c('medianx','>momx','<dadx','><hybridx','mediany','>momy','>dady','>hybridy')

write.csv(tests_df,'Table5.10.csv')
