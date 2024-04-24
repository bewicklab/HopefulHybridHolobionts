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
data <- import_biom('HybridLizardSkins.biom')

#Find OTU table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

#Read in metadata
meta_csv<-read.csv('Whiptail_Skin_Metadata.csv', row.names=1,header= TRUE)

#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)

#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-multi2di(read.tree('rooted_tree_out/tree.nwk'))

#Rename the tree tips to match the ASV table names
cut_names<-c()
for (k in 1:length(taxa_names(tree_file))){
  #Depending on the Qiime2 output, you may need to split the names with a space or with an underscore
  cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],' ')[[1]][1])
  #cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],'_')[[1]][1])
}
taxa_names(tree_file)<-cut_names


#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq<-phyloseq(OTU_biom,TAX,SAMP,tree_file)


#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
ZNC_physeq = prune_taxa(assigned_taxa_seqs, Z_physeq)

#Make ultrametric tree
ultrametric_tree<-chronoMPL(multi2di(phy_tree(Z_physeq)))
#plot.phylo(ultrametric_tree)

#Remove samples below a certain number of reads (I found that the prune_taxa command was somehow adding back in the low read samples... so I'm removing them here...)
highreads<-which(colSums(otu_table(ZNC_physeq))>10000)
ZNC_physeq<-prune_samples(sample_names(ZNC_physeq) %in% c(sample_names(ZNC_physeq)[highreads]),ZNC_physeq)

#Define interesting treatment groups
type_full<-sample_data(ZNC_physeq)$speciesxsite

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Bootstrap gamma 'diversity'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(2)

#The number from each treatment group you want to pool to assess gamma diversity for each bootstrap
sample_size<-12

#The number of bootstrap samples you want to do
bootstrap_no<-15

#The treatment group lables
types<-unique(type_full)

#Define a matrix that will keep track of each treatment group gamma diversity for each bootstrap sample
boot_tracker<-data.frame()
faith_tracker<-data.frame()


#For each bootstrap sample...
for (j in 1:bootstrap_no){

  #Make a list of the individual animals you're going to randomly choose
  picklist<-c()
  for (k in 1:length(types)){
    #Randomly select the appropriate number of animals from each group
    pick<-sample(which(sample_data(ZNC_physeq)$speciesxsite==types[k]),sample_size,replace=FALSE)
    #Add those animals to your list
    picklist<-c(picklist,pick)
  }

  #Make a phyloseq object of only the individual animals you selected for this bootstrap sample
  ZB_physeq<-prune_samples(sample_names(ZNC_physeq)[picklist],ZNC_physeq)

  #Rarefy the OTU table specifically for this bootstrap sample
  rZOTU<-rarefy_even_depth(ZB_physeq,verbose=FALSE,rngseed=sample(1000,1))

  #Merge the samples by treatment group
  mZOTU<-merge_samples(rZOTU,sample_data(rZOTU)$speciesxsite,fun=sum)

  #Find Faith's pd for each sample
  faithpd<-pd((otu_table(mZOTU)), tree_file, include.root=FALSE)
  faithsv<-faithpd[,1]


  #Find the total number of microbial taxa present across all individual animals within a treatment group (gamma diversity)
  boot_tracker<-rbind(boot_tracker,rowSums(sign(otu_table(mZOTU))))
  faith_tracker<-rbind(faith_tracker,faithsv)


  }

colnames(boot_tracker)<-rownames(otu_table(mZOTU))
dataframe_order_types<-colnames(boot_tracker)


#Now turn your bootstraps into two lists
type<-c()
richness<-c()
faiths<-c()

for (k in 1:length(dataframe_order_types)){
  #A list of the treatment group for each pooled bootstrap
  type<-c(type,rep(dataframe_order_types[k],bootstrap_no))
  #A list of the gamma richness for each pooled bootstrap
  richness<-c(richness,boot_tracker[,k])
  #A list of the gamma faiths for each pooled bootstrap
  faiths<-c(faiths,faith_tracker[,k])
}

#Put your richness and treatment group data into a dataframe
diversity_df<-data.frame(type,richness,faiths)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm'))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test whether there are differences in ALPHA diversity between groups

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness
kw_richness<-kruskal.test(richness ~ type, data = diversity_df)
if (kw_richness$p.value<0.05){
  print('Gamma Richness')
  pw_richness<-pairwise.wilcox.test(diversity_df$richness, diversity_df$type,p.adjust.method = "BH")
  print(pw_richness)
}

#faiths
kw_faiths<-kruskal.test(faiths ~ type, data = diversity_df)
if (kw_faiths$p.value<0.05){
  print('Gamma faiths')
  pw_faiths<-pairwise.wilcox.test(diversity_df$faiths, diversity_df$type,p.adjust.method = "BH")
  print(pw_faiths)
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
p_richness<-p_richness+annotate("text", size=6,x=1.5, y=7000, label= "Kruskal-Wallis")
p_richness<-p_richness+annotate("text", size=6,x=1.5, y=6750, label= paste("p = ",round(kw_richness$p.value,4)))

#If there are significant differences in richness between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_richness$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-8800
  y_step<-250
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
p_richness #SI Fig. 1.9d


##########################Faith's Violin Plot#########################################################################################################################################

#Define violin plot
p_faiths <- ggplot(diversity_df, aes(x=type, y=faiths)) + geom_violin(aes(fill=type)) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_faiths<-p_faiths+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_faiths<-p_faiths+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_faiths<-p_faiths+scale_fill_manual(values=c("red", "magenta", "lightslateblue","blue"))
#Add boxplots inside the violins
p_faiths<-p_faiths+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_faiths<-p_faiths+annotate("text", size=6,x=1.5, y=310, label= "Kruskal-Wallis")
p_faiths<-p_faiths+annotate("text", size=6,x=1.5, y=300, label= paste("p = ",round(kw_faiths$p.value,4)))

#If there are significant differences in faiths between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_faiths$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-360
  y_step<-7
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_faiths$p.value))){
    for (j in 1:k){
      if (rownames(pw_faiths$p.value)[k]!=colnames(pw_faiths$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_faiths$p.value[k,j]<0.1 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_faiths$p.value)[k])
          group2<-c(group2,colnames(pw_faiths$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_faiths$p.value[k,j])),6)
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
  stat.test_faiths<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_faiths<-p_faiths+stat_pvalue_manual(stat.test_faiths,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_faiths #SI Fig. 1.9e & Main Fig. 1f

