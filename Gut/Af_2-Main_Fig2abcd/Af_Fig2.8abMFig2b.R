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
data <- import_biom('HybridLizardGuts.biom')

#Find OTU table in biom file
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

#Remove samples below a certain number of reads (I found that the prune_taxa command was somehow adding back in the low read samples... so I'm removing them here...)
highreads<-which(colSums(otu_table(ZNC_physeq))>10000)
ZNC_physeq<-prune_samples(sample_names(ZNC_physeq) %in% c(sample_names(ZNC_physeq)[highreads]),ZNC_physeq)

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
faith_boot_tracker<-data.frame()    #for Faith's PD

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
  #The tree pared down to the microbial taxa in the rarefied OTU table
  straight_tree<-phy_tree(rZOTU)

  #Initialize a vector of counts of core microbes for each treatment group
  corecount<-c()
  allcount<-c()
  percentcore<-c()

  #Initialize vectors of Faith's PD measures of the core microbes for each treatment group
  Faithcount<-c()
  Faithcountall<-c()
  Faithpercent<-c()

  #For each treatment group...
  for (k in 1:length(types)){

    #Pull out the microbiomes from that treatment group and put them in their own phyloseq object
    rZOTUtemp<-prune_samples(sample_names(ZNC_physeq)[pickgroups[k,]],rZOTU)

    ##########################ASV Richness#########################################################################################################################################
    #Sum up the presences of each microbial taxon in the subsetted phyloseq object... count how many microbial taxa are present on as many or more individuals than the defined core fraction
    corecount<-c(corecount,length(which(rowSums(sign(otu_table(rZOTUtemp)))>=core*sample_size)))
    allcount<-c(corecount,length(which(rowSums(sign(otu_table(rZOTUtemp)))>=1)))
    percentcore<-c(percentcore,100*length(which(rowSums(sign(otu_table(rZOTUtemp)))>=core*sample_size))/length(which(rowSums(sign(otu_table(rZOTUtemp)))>=1)))

    ##########################ASV Faith's PD#########################################################################################################################################

    #Pull out the microbial taxa that are part of the host's 'core'
    coretaxa<-rownames(otu_table(rZOTUtemp))[which(rowSums(sign(otu_table(rZOTUtemp)))>=core*sample_size)]
    alltaxa<-rownames(otu_table(rZOTUtemp))[which(rowSums(sign(otu_table(rZOTUtemp)))>=1)]

    #Make a new phyloseq object pruned to only include what is considered a 'core' microbe
    rZOTUtempcore = prune_taxa(coretaxa, rZOTUtemp)
    rZOTUtempall = prune_taxa(alltaxa, rZOTUtemp)

    #Set the counts of all core microbes to 1 across all samples
    #This seems strange, but it will allow us to consider all 'core' microbes present on any animal as if they were present on all animals
    #We can then calculate a Faith's pd for a single animal that contains ALL the 'core' microbes, even though not all the animals contain all the core microbes unless core = 1
    #We can now use any of the animals with ALL the 'core' microbes (which is now all the animals) to get a Faith's pd measure of the 'core' microbiome itself, which may not be present on any single animal
    onematrix<-sign(t(otu_table(rZOTUtempcore))+1)
    onematrixall<-sign(t(otu_table(rZOTUtempall))+1)

    #Find Faith's pd for the core taxa from this particular bootstrap of this species (non-ultrametric tree)
    faithpd<-pd(onematrix, straight_tree, include.root=FALSE)
    faiths<-faithpd[,1]
    faithpdall<-pd(onematrixall, straight_tree, include.root=FALSE)
    faithsall<-faithpdall[,1]
    #use the Faith's pd from the first animal (remember they are now all the same)
    Faithcount<-c(Faithcount,faiths[1])
    Faithcountall<-c(Faithcountall,faithsall[1])
    Faithpercent<-c(Faithpercent, 100*faiths[1]/faithsall[1])

  }

  #Record the core counts and Faith's pd measures for each treatment group from this particular bootstrap
  boot_tracker<-rbind(boot_tracker,percentcore)
  faith_boot_tracker<-rbind(faith_boot_tracker,Faithpercent)

}

#Name the bootstrap columns based on the treatment groups
colnames(boot_tracker)<-types
colnames(faith_boot_tracker)<-types

#Turn your bootstraps into four lists...
#One with the treatment group for each measurement
type<-c()
#One with the core richness for each measurement
richness<-c()
#One with the faith's pd from the non-ultrametric tree
faith<-c()


for (k in 1:length(types)){
  #A list of the treatment group for each pooled bootstrap
  type<-c(type,rep(types[k],bootstrap_no))
  #A list of the core richness for each pooled bootstrap
  richness<-c(richness,boot_tracker[,k])
  #A list of the faith's pd (non-ultrametric) for each pooled bootstrap
  faith<-c(faith,faith_boot_tracker[,k])

}

#Put your richness, faith's pd and treatment group data into a dataframe
diversity_df<-data.frame(type,richness,faith)
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

#Faith's PD (non-ultrametric tree)
kw_faith<-kruskal.test(faith ~ type, data = diversity_df)
if (kw_faith$p.value<0.05){
  print('Core Faiths PD')
  pw_faith<-pairwise.wilcox.test(diversity_df$faith, diversity_df$type,p.adjust.method = "BH")
  print(pw_faith)
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
p_richness<-p_richness+annotate("text", size=6,x=3, y=0.3, label= "Kruskal-Wallis")
p_richness<-p_richness+annotate("text", size=6,x=3, y=0.2, label= paste("p = ",round(kw_richness$p.value,4)))

#If there are significant differences in richness between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_richness$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-1.1
  y_step<-0.07
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
p_richness #SI Fig. 2.8a


##########################Faith's PD Violin Plot#########################################################################################################################################

#Define violin plot
p_faith <- ggplot(diversity_df, aes(x=type, y=faith)) + geom_violin(aes(fill=type)) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_faith<-p_faith+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_faith<-p_faith+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_faith<-p_faith+scale_fill_manual(values=c("red", "magenta", "lightslateblue","blue"))
#Add boxplots inside the violins
p_faith<-p_faith+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_faith<-p_faith+annotate("text", size=6,x=3.2, y=1.35, label= "Kruskal-Wallis")
p_faith<-p_faith+annotate("text", size=6,x=3.2, y=1, label= paste("p = ",round(kw_faith$p.value,4)))

#If there are significant differences in faith between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_faith$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-5.0
  y_step<-0.3
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_faith$p.value))){
    for (j in 1:k){
      if (rownames(pw_faith$p.value)[k]!=colnames(pw_faith$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_faith$p.value[k,j]<0.1 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_faith$p.value)[k])
          group2<-c(group2,colnames(pw_faith$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_faith$p.value[k,j])),6)
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
  stat.test_faith<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_faith<-p_faith+stat_pvalue_manual(stat.test_faith,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_faith #SI Fig. 2.8b & Main Fig. 2b



