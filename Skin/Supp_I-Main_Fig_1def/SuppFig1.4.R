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
library('picante')
library('SYNCSA')
library('pracma')
library('ggpubr')
library('tidyverse')
library('ggpattern')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make ASV Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Read in CSV file for OTU table (table.csv for ASVs, tablelevel6.csv for genera)
data <- import_biom('HybridLizardSkins.biom')

#Find ASV table in biom file
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

#Remove samples below a certain number of reads (no samples have fewer than 10000 reads)
#highreads<-which(colSums(otu_table(ZNC_physeq))>10000)
#ZNC_physeq<-prune_samples(sample_names(ZNC_physeq) %in% c(sample_names(ZNC_physeq)[highreads]),ZNC_physeq)

#Rarefy the OTU table
rZOTU<-rarefy_even_depth(ZNC_physeq,rngseed = 1)

#Define the catogory as hybrid or non-hybrid
type<-sample_data(rZOTU)$hybrid


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate a variety of different ALPHA diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Find the richness for each sample
richness<-colSums(sign(otu_table(rZOTU)))

#Find the shannon diversity for each sample
shannon<-diversity(otu_table(t(rZOTU)),index='shannon')

#Find the simpson's index for each sample
simpson<-diversity(otu_table(t(rZOTU)),index='simpson')

#Find Faith's pd for each sample
faithpd<-pd(t(otu_table(rZOTU)), tree_file, include.root=FALSE)
faiths<-faithpd[,1]

#The Rao's entropy package requires an ultrametric tree.... make the ultrametric tree and plot it to see how it's different
ultrametric_tree<-chronoMPL(multi2di(phy_tree(Z_physeq)))
plot.phylo(ultrametric_tree)

#Calculate Rao's entropy on the ultrametric tree
raopd<-raoD(comm=t(otu_table(rZOTU)),phy=ultrametric_tree)
raos<-raopd$Dkk

#Put all of your different diversity metrics into a dataframe
diversity_df<-data.frame(type,richness,shannon,simpson,faiths)#, raos)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('asexual','bisexual'))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test whether there are differences in ALPHA diversity between groups

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness
mw_richness<-wilcox.test(diversity_df$richness, sign(diversity_df$type=='asexual'))

#Shannon
mw_shannon<-wilcox.test(diversity_df$shannon, sign(diversity_df$type=='asexual'))

#Simpson
mw_simpson<-wilcox.test(diversity_df$simpson, sign(diversity_df$type=='asexual'))

#Faith's PD
mw_faiths<-wilcox.test(diversity_df$faiths, sign(diversity_df$type=='asexual'))

#Rao's
mw_raos<-wilcox.test(diversity_df$raos, sign(diversity_df$type=='asexual'))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Make Violin plots of the various diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
p_value_ns<-'no'
pshow<-'star'

##########################Richness Violin Plot#########################################################################################################################################

#Define violin plot
p_richness <- ggplot(diversity_df, aes(x=type, y=richness)) + geom_violin_pattern(aes(fill=type,pattern_fill=type),pattern='stripe',pattern_spacing=0.05,pattern_density=0.5) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_richness<-p_richness+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_richness<-p_richness+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))+theme(legend.key.size = unit(1.5, 'cm'))
#Choose the violin colors for each group
p_richness<-p_richness+scale_fill_manual(values=c("magenta", "red"))+ scale_pattern_fill_manual(values = c(asexual='lightslateblue', bisexual='blue'))
#Add boxplots inside the violins
p_richness<-p_richness+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_richness<-p_richness+annotate("text", size=6,x=1, y=740, label= "Mann-Whitney U")
p_richness<-p_richness+annotate("text", size=6,x=1, y=620, label= paste("p = ",round(mw_richness$p.value,4)))

#If there are significant differences in richness between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (mw_richness$p.value<0.05 || p_value_ns=='yes'){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-2500
  #For each pairwise comparison...
  #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
  group1<-'asexual'
  group2<-'bisexual'
  p.adj<-mw_richness$p.value
  #Add the y-position of the bracket and bump the y-location of the bracket up one for next time

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
  }else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_richness<-p_richness+stat_pvalue_manual(stat.test_richness,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_richness #SI Fig. 1.4a

##########################Shannon Diversity Violin Plot#########################################################################################################################################

#Define violin plot
p_shannon <- ggplot(diversity_df, aes(x=type, y=shannon)) + geom_violin_pattern(aes(fill=type,pattern_fill=type),pattern='stripe',pattern_spacing=0.05,pattern_density=0.5) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_shannon<-p_shannon+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_shannon<-p_shannon+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))+theme(legend.key.size = unit(1.5, 'cm'))
#Choose the violin colors for each group
p_shannon<-p_shannon+scale_fill_manual(values=c("magenta", "red"))+scale_pattern_fill_manual(values = c(asexual='lightslateblue', bisexual='blue'))
#Add boxplots inside the violins
p_shannon<-p_shannon+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_shannon<-p_shannon+annotate("text", size=6,x=1.5, y=3.5, label= "Mann-Whitney U")
p_shannon<-p_shannon+annotate("text", size=6,x=1.5, y=3.2, label= paste("p = ",round(mw_shannon$p.value,4)))

#If there are significant differences in shannon between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (mw_shannon$p.value<0.05 || p_value_ns=='yes'){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-7.5
  #For each pairwise comparison...
  #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
  group1<-'asexual'
  group2<-'bisexual'
  p.adj<-mw_shannon$p.value
  #Add the y-position of the bracket and bump the y-location of the bracket up one for next time

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
  }else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_shannon<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_shannon<-p_shannon+stat_pvalue_manual(stat.test_shannon,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_shannon #SI Fig. 1.4b

##########################Simpson's Diversity Violin Plot#########################################################################################################################################

#Define violin plot
p_simpson <- ggplot(diversity_df, aes(x=type, y=simpson)) + geom_violin_pattern(aes(fill=type,pattern_fill=type),pattern='stripe',pattern_spacing=0.05,pattern_density=0.5) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_simpson<-p_simpson+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_simpson<-p_simpson+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))+theme(legend.key.size = unit(1.5, 'cm'))
#Choose the violin colors for each group
p_simpson<-p_simpson+scale_fill_manual(values=c("magenta", "red"))+ scale_pattern_fill_manual(values = c(asexual='lightslateblue', bisexual='blue'))
#Add boxplots inside the violins
p_simpson<-p_simpson+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_simpson<-p_simpson+annotate("text", size=6,x=1.5, y=0.75, label= "Mann-Whitney U")
p_simpson<-p_simpson+annotate("text", size=6,x=1.5, y=0.72, label= paste("p = ",round(mw_simpson$p.value,4)))

#If there are significant differences in simpson between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (mw_simpson$p.value<0.05 || p_value_ns=='yes'){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-7.5
  #For each pairwise comparison...
  #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
  group1<-'asexual'
  group2<-'bisexual'
  p.adj<-mw_simpson$p.value
  #Add the y-position of the bracket and bump the y-location of the bracket up one for next time

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
  }else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_simpson<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_simpson<-p_simpson+stat_pvalue_manual(stat.test_simpson,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_simpson #SI Fig. 1.4c

##########################Faith's PD Violin Plot#########################################################################################################################################

#Define violin plot
p_faiths <- ggplot(diversity_df, aes(x=type, y=faiths)) + geom_violin_pattern(aes(fill=type,pattern_fill=type),pattern='stripe',pattern_spacing=0.05,pattern_density=0.5) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_faiths<-p_faiths+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_faiths<-p_faiths+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))+theme(legend.key.size = unit(1.5, 'cm'))
#Choose the violin colors for each group
p_faiths<-p_faiths+scale_fill_manual(values=c("magenta", "red"))+ scale_pattern_fill_manual(values = c(asexual='lightslateblue', bisexual='blue'))
#Add boxplots inside the violins
p_faiths<-p_faiths+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_faiths<-p_faiths+annotate("text", size=6,x=1, y=80, label= "Mann-Whitney U")
p_faiths<-p_faiths+annotate("text", size=6,x=1, y=73, label= paste("p = ",round(mw_faiths$p.value,4)))

#If there are significant differences in faiths between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (mw_faiths$p.value<0.05 || p_value_ns=='yes'){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-170
  #For each pairwise comparison...
  #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
  group1<-'asexual'
  group2<-'bisexual'
  p.adj<-mw_faiths$p.value
  #Add the y-position of the bracket and bump the y-location of the bracket up one for next time

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
  }else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_faiths<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_faiths<-p_faiths+stat_pvalue_manual(stat.test_faiths,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_faiths #SI Fig. 1.4d


##########################Rao's Entropy Violin Plot#########################################################################################################################################

#Define violin plot
p_raos <- ggplot(diversity_df, aes(x=type, y=raos)) + geom_violin_pattern(aes(fill=type,pattern_fill=type),pattern='stripe',pattern_spacing=0.05,pattern_density=0.5) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_raos<-p_raos+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_raos<-p_raos+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))+theme(legend.key.size = unit(1.5, 'cm'))
#Choose the violin colors for each group
p_raos<-p_raos+scale_fill_manual(values=c("magenta", "red"))+ scale_pattern_fill_manual(values = c(asexual='lightslateblue', bisexual='blue'))
#Add boxplots inside the violins
p_raos<-p_raos+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_raos<-p_raos+annotate("text", size=6,x=1.5, y=0.35, label= "Mann-Whitney U")
p_raos<-p_raos+annotate("text", size=6,x=1.5, y=0.33, label= paste("p = ",round(mw_raos$p.value,4)))

#If there are significant differences in raos between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (mw_raos$p.value<0.05 || p_value_ns=='yes'){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-0.45
  #For each pairwise comparison...
  #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
  group1<-'asexual'
  group2<-'bisexual'
  p.adj<-mw_raos$p.value
  #Add the y-position of the bracket and bump the y-location of the bracket up one for next time

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
  }else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_raos<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_raos<-p_raos+stat_pvalue_manual(stat.test_raos,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_raos #SI Fig. 1.4e



