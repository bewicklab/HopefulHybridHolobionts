library(vegan)
library(phyloseq)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(microshades)
library(speedyseq)
library(cowplot)
library(stringr)
library(knitr)

#Installation of microshades and speedyseq is a bit different; all the other packages can be installed as usual
#remotes::install_github("KarstensLab/microshades")
#remotes::install_github("mikemc/speedyseq")


#%%%%%%%%%%%%%%%%%%%%%%   GET THE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#Rarefy the OTU table
rZOTU<-rarefy_even_depth(ZNC_physeq, rngseed=5)

#Define interesting treatment groups
type<-sample_data(ZNC_physeq)$speciesxsite


#%%%%%%%%%%%%%%%%%%%%%%   MAKE FANCY BAR GRAPHS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum
summed_phyla<-rowSums(data.frame(otu_table(tax_glom(rZOTU,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list<-data.frame(tax_table(tax_glom(rZOTU,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_phyla_list<-phyla_list[order(-summed_phyla)]

#Prepare the OTU table in the correct format for the plotting program (this must be done with level 6 data including genera!)
mdf_prep <- prep_mdf(rZOTU)
#Tell the function to plot the five most abundant phyla (you could pick something different... )
color_objs_GP <- create_color_dfs(mdf_prep,selected_groups = c('Proteobacteria','Actinobacteria','Bacteroidetes','Chloroflexi','Firmicutes'),cvd=TRUE)
#Extract the OTU table and color choices (again, putting things in the right format for the function)
mdf_GP <- color_objs_GP$mdf
cdf_GP <- color_objs_GP$cdf
#Define the legend
GP_legend <-custom_legend(mdf_GP, cdf_GP)

#Expand the number of Actinobacteria shown
new_groups <- extend_group(mdf_GP, cdf_GP, "Phylum", "Genus", "Actinobacteria", existing_palette = "micro_cvd_orange", new_palette = "micro_orange", n_add = 3)
GP_legend_new <-custom_legend(new_groups$mdf, new_groups$cdf)
#Expand the number of Proteobacteria shown
new_groups2 <- extend_group(new_groups$mdf, new_groups$cdf, "Phylum", "Genus", "Proteobacteria", existing_palette = "micro_cvd_green", new_palette = "micro_green", n_add = 3)
#Expand the number of Bacteroidetes shown
new_groups3 <- extend_group(new_groups2$mdf, new_groups2$cdf, "Phylum", "Genus", "Bacteroidetes", existing_palette = "micro_cvd_blue", new_palette = "micro_blue", n_add = 3)
#Expand the number of Firmicutes shown
new_groups4 <- extend_group(new_groups3$mdf, new_groups3$cdf, "Phylum", "Genus", "Firmicutes", existing_palette = "micro_cvd_purple", new_palette = "micro_purple", n_add = 3)
#Expand the number of Firmicutes shown
new_groups5 <- extend_group(new_groups4$mdf, new_groups4$cdf, "Phylum", "Genus", "Chloroflexi", existing_palette = "micro_cvd_turquoise", new_palette = "micro_brown", n_add = 3)

#Define your new legend with the expanded groups
GP_legend_new <-custom_legend(new_groups5$mdf, new_groups5$cdf,legend_key_size=0.6,legend_text_size = 12)

#Define the plot that you will be making
plot2 <- plot_microshades(new_groups4$mdf, new_groups4$cdf)
#Define all the formatting aspects of the plot
plot_diff2 <- plot2 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 6)) +
  facet_grid(~fct_relevel(speciesxsite,'SBluG_ino','SBluG_neo','SNW_neo','SNW_marm'), scale="free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = 30))+
  theme(axis.text.x = element_text(size= 6)) +
  theme(plot.margin = margin(6,20,6,6))+
  theme(axis.title.x = element_text(size= 20))+
  theme(axis.title.y = element_text(size= 20))

#Make the plot
plot_grid(plot_diff2, GP_legend_new,  rel_widths = c(1, .25))


#%%%%%%%%%%%%%%%%%%%%%%   FIND PHYLUM PERCENTAGES ACROSS ALL SPECIES/POPULATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_phyla_order<-summed_phyla[order(-summed_phyla)]
#Turn the vector of sums into a dataframe
abundant_phyla<-data.frame(summed_phyla_order[1:5]*100/sum(summed_phyla))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_phyla)<-sorted_phyla_list[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND PHYLUM PERCENTAGES ACROSS BOTH A. NEOMEXICANUS POPULATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Prune phyloseq object o just neomexicanus
rZOTU_neo<-prune_samples(sample_data(rZOTU)$species=='neomexicanus',rZOTU)
#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_phyla_neo<-rowSums(data.frame(otu_table(tax_glom(rZOTU_neo,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list_neo<-data.frame(tax_table(tax_glom(rZOTU_neo,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_phyla_list_neo<-phyla_list_neo[order(-summed_phyla_neo)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_phyla_order_neo<-summed_phyla_neo[order(-summed_phyla_neo)]
#Turn the vector of sums into a dataframe
abundant_phyla_neo<-data.frame(summed_phyla_order_neo[1:5]*100/sum(summed_phyla_neo))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_phyla_neo)<-sorted_phyla_list_neo[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND PHYLUM PERCENTAGES ACROSS THE A. NEOMEXICANUS POPULATION FROM SNW   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Prune phyloseq object o just neomexicanus
rZOTU_SNW_neo<-prune_samples(sample_data(rZOTU)$speciesxsite=='SNW_neo',rZOTU)
#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_phyla_SNW_neo<-rowSums(data.frame(otu_table(tax_glom(rZOTU_SNW_neo,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list_SNW_neo<-data.frame(tax_table(tax_glom(rZOTU_SNW_neo,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_phyla_list_SNW_neo<-phyla_list_SNW_neo[order(-summed_phyla_SNW_neo)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_phyla_order_SNW_neo<-summed_phyla_SNW_neo[order(-summed_phyla_SNW_neo)]
#Turn the vector of sums into a dataframe
abundant_phyla_SNW_neo<-data.frame(summed_phyla_order_SNW_neo[1:5]*100/sum(summed_phyla_SNW_neo))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_phyla_SNW_neo)<-sorted_phyla_list_SNW_neo[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND PHYLUM PERCENTAGES ACROSS THE A. NEOMEXICANUS POPULATION FROM SBluG   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Prune phyloseq object o just neomexicanus
rZOTU_SBluG_neo<-prune_samples(sample_data(rZOTU)$speciesxsite=='SBluG_neo',rZOTU)
#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_phyla_SBluG_neo<-rowSums(data.frame(otu_table(tax_glom(rZOTU_SBluG_neo,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list_SBluG_neo<-data.frame(tax_table(tax_glom(rZOTU_SBluG_neo,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_phyla_list_SBluG_neo<-phyla_list_SBluG_neo[order(-summed_phyla_SBluG_neo)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_phyla_order_SBluG_neo<-summed_phyla_SBluG_neo[order(-summed_phyla_SBluG_neo)]
#Turn the vector of sums into a dataframe
abundant_phyla_SBluG_neo<-data.frame(summed_phyla_order_SBluG_neo[1:5]*100/sum(summed_phyla_SBluG_neo))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_phyla_SBluG_neo)<-sorted_phyla_list_SBluG_neo[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND PHYLUM PERCENTAGES ACROSS THE A. MARMORATUS POPULATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Prune phyloseq object o just neomexicanus
rZOTU_marm<-prune_samples(sample_data(rZOTU)$species=='marmoratus',rZOTU)
#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_phyla_marm<-rowSums(data.frame(otu_table(tax_glom(rZOTU_marm,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list_marm<-data.frame(tax_table(tax_glom(rZOTU_marm,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_phyla_list_marm<-phyla_list_marm[order(-summed_phyla_marm)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_phyla_order_marm<-summed_phyla_marm[order(-summed_phyla_marm)]
#Turn the vector of sums into a dataframe
abundant_phyla_marm<-data.frame(summed_phyla_order_marm[1:5]*100/sum(summed_phyla_marm))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_phyla_marm)<-sorted_phyla_list_marm[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND PHYLUM PERCENTAGES ACROSS THE A. INORNATUS POPULATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Prune phyloseq object o just neomexicanus
rZOTU_ino<-prune_samples(sample_data(rZOTU)$species=='inornatus',rZOTU)
#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_phyla_ino<-rowSums(data.frame(otu_table(tax_glom(rZOTU_ino,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list_ino<-data.frame(tax_table(tax_glom(rZOTU_ino,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_phyla_list_ino<-phyla_list_ino[order(-summed_phyla_ino)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_phyla_order_ino<-summed_phyla_ino[order(-summed_phyla_ino)]
#Turn the vector of sums into a dataframe
abundant_phyla_ino<-data.frame(summed_phyla_order_ino[1:5]*100/sum(summed_phyla_ino))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_phyla_ino)<-sorted_phyla_list_ino[1:5]


#%%%%%%%%%%%%%%%%%%%%%%   FIND GENUS PERCENTAGES ACROSS ALL POPULATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_genus<-rowSums(data.frame(otu_table(tax_glom(rZOTU,taxrank="Genus"))))
#Find the names of the different phylum
genus_list<-data.frame(tax_table(tax_glom(rZOTU,taxrank="Genus")))[,6]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_genus_list<-genus_list[order(-summed_genus)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_genus_order<-summed_genus[order(-summed_genus)]
#Turn the vector of sums into a dataframe
abundant_genus<-data.frame(summed_genus_order[1:5]*100/sum(summed_genus))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_genus)<-sorted_genus_list[1:5]


#%%%%%%%%%%%%%%%%%%%%%%   FIND GENUS PERCENTAGES ACROSS BOTH A. NEOMEXICANUS POPULATIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_genus_neo<-rowSums(data.frame(otu_table(tax_glom(rZOTU_neo,taxrank="Genus"))))
#Find the names of the different phylum
genus_list_neo<-data.frame(tax_table(tax_glom(rZOTU_neo,taxrank="Genus")))[,6]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_genus_list_neo<-genus_list_neo[order(-summed_genus_neo)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_genus_order_neo<-summed_genus_neo[order(-summed_genus_neo)]
#Turn the vector of sums into a dataframe
abundant_genus_neo<-data.frame(summed_genus_order_neo[1:5]*100/sum(summed_genus_neo))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_genus_neo)<-sorted_genus_list_neo[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND GENUS PERCENTAGES ACROSS THE A. NEOMEXICANUS POPULATION FROM SNW   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_genus_SNW_neo<-rowSums(data.frame(otu_table(tax_glom(rZOTU_SNW_neo,taxrank="Genus"))))
#Find the names of the different phylum
genus_list_SNW_neo<-data.frame(tax_table(tax_glom(rZOTU_SNW_neo,taxrank="Genus")))[,6]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_genus_list_SNW_neo<-genus_list_SNW_neo[order(-summed_genus_SNW_neo)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_genus_order_SNW_neo<-summed_genus_SNW_neo[order(-summed_genus_SNW_neo)]
#Turn the vector of sums into a dataframe
abundant_genus_SNW_neo<-data.frame(summed_genus_order_SNW_neo[1:5]*100/sum(summed_genus_SNW_neo))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_genus_SNW_neo)<-sorted_genus_list_SNW_neo[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND GENUS PERCENTAGES ACROSS THE A. NEOMEXICANUS POPULATION FROM SBluG   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_genus_SBluG_neo<-rowSums(data.frame(otu_table(tax_glom(rZOTU_SBluG_neo,taxrank="Genus"))))
#Find the names of the different phylum
genus_list_SBluG_neo<-data.frame(tax_table(tax_glom(rZOTU_SBluG_neo,taxrank="Genus")))[,6]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_genus_list_SBluG_neo<-genus_list_SBluG_neo[order(-summed_genus_SBluG_neo)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_genus_order_SBluG_neo<-summed_genus_SBluG_neo[order(-summed_genus_SBluG_neo)]
#Turn the vector of sums into a dataframe
abundant_genus_SBluG_neo<-data.frame(summed_genus_order_SBluG_neo[1:5]*100/sum(summed_genus_SBluG_neo))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_genus_SBluG_neo)<-sorted_genus_list_SBluG_neo[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND GENUS PERCENTAGES ACROSS THE A. MARMORATUS POPULATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_genus_marm<-rowSums(data.frame(otu_table(tax_glom(rZOTU_marm,taxrank="Genus"))))
#Find the names of the different phylum
genus_list_marm<-data.frame(tax_table(tax_glom(rZOTU_marm,taxrank="Genus")))[,6]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_genus_list_marm<-genus_list_marm[order(-summed_genus_marm)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_genus_order_marm<-summed_genus_marm[order(-summed_genus_marm)]
#Turn the vector of sums into a dataframe
abundant_genus_marm<-data.frame(summed_genus_order_marm[1:5]*100/sum(summed_genus_marm))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_genus_marm)<-sorted_genus_list_marm[1:5]

#%%%%%%%%%%%%%%%%%%%%%%   FIND GENUS PERCENTAGES ACROSS THE A. INORNATUS POPULATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum in the neomexicanus phyloseq table
summed_genus_ino<-rowSums(data.frame(otu_table(tax_glom(rZOTU_ino,taxrank="Genus"))))
#Find the names of the different phylum
genus_list_ino<-data.frame(tax_table(tax_glom(rZOTU_ino,taxrank="Genus")))[,6]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant...
sorted_genus_list_ino<-genus_list_ino[order(-summed_genus_ino)]
#Sort the read tallies for each phylum names so that the first is the most abundant, the second the second most abundant...
summed_genus_order_ino<-summed_genus_ino[order(-summed_genus_ino)]
#Turn the vector of sums into a dataframe
abundant_genus_ino<-data.frame(summed_genus_order_ino[1:5]*100/sum(summed_genus_ino))
#Name the rows of the dataframe according to their phylum names
rownames(abundant_genus_ino)<-sorted_genus_list_ino[1:5]

