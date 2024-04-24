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


triangle_euclidean_all_2<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=2,seed=1,ordination='PCoA')

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

triangle_euclidean_all_20<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=20,seed=1,ordination='PCoA')

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


triangle_jaccard_all_2<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=2,seed=1,dist='jaccard',ordination='PCoA')

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

triangle_jaccard_all_20<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=20,seed=1,dist='jaccard',ordination='PCoA')

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


triangle_bray_all_2<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=2,seed=1,dist='bray',ordination='PCoA')

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

triangle_bray_all_20<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=20,seed=2,dist='bray',ordination='PCoA')

triangle_bray_all_20_null1<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=2,dist='bray',ordination='PCoA')
triangle_bray_all_20_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=2,null_model=6,dist='bray',ordination='PCoA')
triangle_bray_all_20_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=2,null_model=7,dist='bray',ordination='PCoA')

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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 PCoA Unifrac

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


triangle_unifrac_all_2<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=2,seed=1,dist='unifrac',ordination='PCoA')

triangle_unifrac_all_2_null10<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model = 10,dist='unifrac',ordination='PCoA')
triangle_unifrac_all_2_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=6,dist='unifrac',ordination='PCoA')
triangle_unifrac_all_2_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=7,dist='unifrac',ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_unifrac_all_2[,1:2],triangle_unifrac_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_unifrac_all_2[,1:2],triangle_unifrac_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_unifrac_all_2[,1:2],triangle_unifrac_all_2_null10[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_unifrac_all_2[,1:2],triangle_unifrac_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_unifrac_all_2[,1:2],triangle_unifrac_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_unifrac_all_2[,1:2],triangle_unifrac_all_2_null10[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_unifrac_all_2[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_unifrac_all_2[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))

triangle_unifrac_all_20<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=20,seed=2,dist='unifrac',ordination='PCoA')

triangle_unifrac_all_20_null10<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=2,null_model=10,dist='unifrac',ordination='PCoA')
triangle_unifrac_all_20_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=2,null_model=6,dist='unifrac',ordination='PCoA')
triangle_unifrac_all_20_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=20,seed=2,null_model=7,dist='unifrac',ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_unifrac_all_20[,1:2],triangle_unifrac_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_unifrac_all_20[,1:2],triangle_unifrac_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_unifrac_all_20[,1:2],triangle_unifrac_all_20_null10[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_unifrac_all_20[,1:2],triangle_unifrac_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_unifrac_all_20[,1:2],triangle_unifrac_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_unifrac_all_20[,1:2],triangle_unifrac_all_20_null10[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_unifrac_all_20[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_unifrac_all_20[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 PCoA weighted UniFrac

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


triangle_wunifrac_all_2<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=2,seed=1,dist='wunifrac',ordination='PCoA')

triangle_wunifrac_all_2_null1<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,dist='wunifrac',ordination='PCoA')
triangle_wunifrac_all_2_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=6,dist='wunifrac',ordination='PCoA')
triangle_wunifrac_all_2_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=2,seed=1,null_model=7,dist='wunifrac',ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_wunifrac_all_2[,1:2],triangle_wunifrac_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_wunifrac_all_2[,1:2],triangle_wunifrac_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_wunifrac_all_2[,1:2],triangle_wunifrac_all_2_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_wunifrac_all_2[,1:2],triangle_wunifrac_all_2_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_wunifrac_all_2[,1:2],triangle_wunifrac_all_2_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_wunifrac_all_2[,1:2],triangle_wunifrac_all_2_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_wunifrac_all_2[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_wunifrac_all_2[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))

triangle_wunifrac_all_20<-TriangleHbootstrap(Z_physeq_SBluG,species_no_SBluG,15,12,dimensions=15,seed=2,dist='wunifrac',ordination='PCoA')

triangle_wunifrac_all_20_null1<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=15,seed=2,dist='wunifrac',ordination='PCoA')
triangle_wunifrac_all_20_null6<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=15,seed=2,null_model=6,dist='wunifrac',ordination='PCoA')
triangle_wunifrac_all_20_null7<-TriangleHnull(Z_physeq,species_no,15,12,dimensions=15,seed=2,null_model=7,dist='wunifrac',ordination='PCoA')

#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
momxtest<-TriangleHcompare(rbind(triangle_wunifrac_all_20[,1:2],triangle_wunifrac_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly less than progenitor 2 along parental axis of variation
dadxtest<-TriangleHcompare(rbind(triangle_wunifrac_all_20[,1:2],triangle_wunifrac_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)
#Test if hybrid is significantly greater than progenitor 1 along parental axis of variation
hybridxtest<-TriangleHcompare(rbind(triangle_wunifrac_all_20[,1:2],triangle_wunifrac_all_20_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='two.sided',dimension=1)

#Test if hybrid is significantly greater than progenitor 1 along perpendicular axis
momytest<-TriangleHcompare(rbind(triangle_wunifrac_all_20[,1:2],triangle_wunifrac_all_20_null6[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than progenitor 2 along perpendicular axis
dadytest<-TriangleHcompare(rbind(triangle_wunifrac_all_20[,1:2],triangle_wunifrac_all_20_null7[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)
#Test if hybrid is significantly less than mixed progenitor 1 and 2 along perpendicular axis
hybridytest<-TriangleHcompare(rbind(triangle_wunifrac_all_20[,1:2],triangle_wunifrac_all_20_null1[,1:2]),c(rep(1,15),rep(2,15)),comparison = 'MW',hypothesis='greater',dimension=2)

tests<-rbind(tests,c(median(triangle_wunifrac_all_20[,1]),momxtest$p.value,dadxtest$p.value,hybridxtest$p.value,median(triangle_wunifrac_all_20[,2]),momytest$p.value,dadytest$p.value,hybridytest$p.value))

tests_df<-data.frame(tests)
rownames(tests_df)<-c('euclidean_2','euclidean_20','jaccard_2','jaccard_20','bray_2','bray_20','unifrac_2','unifrac_20','wunifrac_2','wunifrac_20')
colnames(tests_df)<-c('medianx','>momx','<dadx','><hybridx','mediany','>momy','>dady','>hybridy')

write.csv(tests_df,'Table5.3.csv')



