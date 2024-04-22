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
ZNC_physeq<-phyloseq(OTU_biom,TAX,SAMP)

#Define interesting treatment groups
type_full<-sample_data(ZNC_physeq)$speciesxsite


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Bootstrap core richness

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(4)

#The number from each treatment group you want to pool to assess core richness for each bootstrap
sample_size<-12

#The number of bootstrap samples you want to do
bootstrap_no<-15

#Define the core fraction
core<-0.5

#The treatment group lables
types<-c("SBluG_ino" , "SNW_marm" , "SBluG_neo" , "SNW_neo")

#Define a matrix that will keep track of each treatment group's core richness and core faith's PD for each bootstrap sample
boot_tracker<-data.frame()

#For each bootstrap sample...
for (j in 1:bootstrap_no){
  print(paste('bootstrap',j))   #It was going slow... this way you can see progress

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

  #Make a phyloseq object of only the individual animals you selected for this bootstrap sample
  ZB_physeq<-prune_samples(sample_names(ZNC_physeq)[picklist],ZNC_physeq)

  #Rarefy the OTU table specifically for this bootstrap sample
  rZOTU<-rarefy_even_depth(ZB_physeq,verbose=FALSE,rngseed=sample(1000,1))

  #Initialize a vector of counts of core microbes for each treatment group
  corecount<-c()

  corenames<-c()
  allnames<-c()
  
  #For each treatment group...
  for (k in 1:length(types)){

    #Pull out the microbiomes from that treatment group and put them in their own phyloseq object
    rZOTUtemp<-prune_samples(sample_names(ZNC_physeq)[pickgroups[k,]],rZOTU)
    
    #Make a list of the rownumbers associated with microbial taxa that satisfy the core requirement
    find_core<-which(rowSums(sign(otu_table(rZOTUtemp)))>=core*length(pickgroups[k,]))
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
    for (jj in 1:length(find_names)){
      list_names<-c(list_names,find_names[jj])
    }
    #For each taxon in the list...
    for (jj in 1:length(find_names_all)){
      list_names_all<-c(list_names_all,find_names_all[jj])
    }
    
    
    #Add the list of core microbial names to a 'list', where each entry is a list of microbial core taxa names for a particular treatment group
    corenames[[types[k]]]<-list_names
    allnames[[types[k]]]<-list_names_all
    
    ##########################ASV Richness#########################################################################################################################################
    #Sum up the presences of each microbial taxon in the subsetted phyloseq object... count how many microbial taxa are present on as many or more individuals than the defined core fraction
    corecount<-c(corecount,length(which(rowSums(sign(otu_table(rZOTUtemp)))>=core*sample_size)))
    


  }
  
  #Sort your core microbiome names in the order you want them to appear in the Venn Diagram
  corenames_sort<-corenames[c("SBluG_ino" , "SNW_marm" , "SBluG_neo" , "SNW_neo")]
  allnames_sort<-allnames[c("SBluG_ino" , "SNW_marm" , "SBluG_neo" , "SNW_neo")]
  
  in1<-na.omit(setdiff(setdiff(setdiff(corenames_sort[[1]],corenames_sort[[2]]),corenames_sort[[3]]),corenames_sort[[4]]))
  in2<-na.omit(setdiff(setdiff(setdiff(corenames_sort[[2]],corenames_sort[[1]]),corenames_sort[[3]]),corenames_sort[[4]]))
  in3<-na.omit(setdiff(setdiff(corenames_sort[[3]],corenames_sort[[2]]),corenames_sort[[1]]))
  in4<-na.omit(setdiff(setdiff(corenames_sort[[4]],corenames_sort[[2]]),corenames_sort[[1]]))
  
  #Percentage of taxa unique to SBluG_neo
  uniques<-c(length(in1),length(in2),length(in3),length(in4))
  


  #Record the core counts and Faith's pd measures for each treatment group from this particular bootstrap
  boot_tracker<-rbind(boot_tracker,uniques)

}

#Name the bootstrap columns based on the treatment groups
colnames(boot_tracker)<-c("SBluG_ino" , "SNW_marm" , "SBluG_neo" , "SNW_neo")

typer<-c("SBluG_ino" , "SNW_marm" , "SBluG_neo" , "SNW_neo")

#Turn your bootstraps into four lists...
#One with the treatment group for each measurement
type<-c()
#One with the core richness for each measurement
richness<-c()


for (k in 1:length(typer)){
  #A list of the treatment group for each pooled bootstrap
  type<-c(type,rep(typer[k],bootstrap_no))
  #A list of the core richness for each pooled bootstrap
  richness<-c(richness,boot_tracker[,k])
 
}

#Put your richness, faith's pd and treatment group data into a dataframe
diversity_df<-data.frame(type,richness)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm'))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test whether there are differences in CORE RICHNESS between groups

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness
kw_richness<-kruskal.test(richness ~ type, data = diversity_df)
if (kw_richness$p.value<0.05){
  print('Core Richness')
  pw_richness<-pairwise.wilcox.test(diversity_df$richness, diversity_df$type,p.adjust.method = "BH")
  print(pw_richness)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Make Violin plots of the various diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
p_value_ns<-'no'
pshow<-'star'

##########################Richness Violin Plot#########################################################################################################################################

#Define violin plot
p_richness <- ggplot(diversity_df, aes(x=type, y=richness)) + geom_violin(aes(fill=type)) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_richness<-p_richness+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_richness<-p_richness+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_richness<-p_richness+scale_fill_manual(values=c("red", "magenta", "lightslateblue","blue"))
#Add boxplots inside the violins
p_richness<-p_richness+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_richness<-p_richness+annotate("text", size=6,x=2.7, y=5, label= "Kruskal-Wallis")
p_richness<-p_richness+annotate("text", size=6,x=2.7, y=0, label= paste("p = ",round(kw_richness$p.value,4)))

#If there are significant differences in richness between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_richness$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-65
  y_step<-3
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_richness$p.value))){
    for (j in 1:k){
      if (rownames(pw_richness$p.value)[k]!=colnames(pw_richness$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_richness$p.value[k,j]<0.1 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_richness$p.value)[k])
          group2<-c(group2,colnames(pw_richness$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_richness$p.value[k,j])),6)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_richness<-p_richness+stat_pvalue_manual(stat.test_richness,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_richness #Fig. 2.15f






