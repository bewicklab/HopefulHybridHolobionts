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
library('ape')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#Plot the phylogenetic tree - does it look crazy? Some 'bad stuff' may be in there
#plot.phylo(phy_tree(Z_physeq))

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
Z_physeq = prune_taxa(assigned_taxa_seqs, Z_physeq)


#Remove samples below a certain number of reads (I found that the prune_taxa command was somehow adding back in the low read samples... so I'm removing them here...)
highreads<-which(colSums(otu_table(Z_physeq))>10000)
ZNC_physeq<-prune_samples(sample_names(Z_physeq) %in% c(sample_names(Z_physeq)[highreads]),Z_physeq)


#Define interesting treatment groups
type_full<-sample_data(ZNC_physeq)$speciesxsite


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate a variety of different BETA diversity metrics using the BETAPART package

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

community_matrix<-otu_table(ZNC_physeq)
class(community_matrix) <- "matrix"

betapart_jaccard<-beta.multi(sign(t(community_matrix)),index.family='jaccard')
betapart_bray<-beta.multi.abund(t(community_matrix),index.family='bray')
betapart_unifrac<-phylo.beta.multi(sign(t(community_matrix)),phy_tree(ZNC_physeq),index.family='jaccard')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Bootstrap beta diversity

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(3)

#The number from each treatment group you want to pool to assess core richness for each bootstrap
sample_size<-12

#The number of bootstrap samples you want to do
bootstrap_no<-15

#The treatment group lables
types<-unique(type_full)

#Define a matrix that will keep track of each treatment group core richness for each bootstrap sample
boot_tracker_turnover_jaccard<-data.frame()
boot_tracker_nestedness_jaccard<-data.frame()
boot_tracker_overall_jaccard<-data.frame()

#Define a matrix that will keep track of each treatment group core richness for each bootstrap sample
boot_tracker_turnover_jaccard_phy<-data.frame()
boot_tracker_nestedness_jaccard_phy<-data.frame()
boot_tracker_overall_jaccard_phy<-data.frame()


#boot_tracker_turnover_bray<-data.frame()
#boot_tracker_nestedness_bray<-data.frame()
boot_tracker_overall_bray<-data.frame()

#For each bootstrap sample...
for (j in 1:bootstrap_no){
  print(paste('bootstrap',j))

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
  rZOTU<-rarefy_even_depth(ZB_physeq,verbose=FALSE,rngseed = sample(1000,1))




  #Initialize a vector of counts of core microbes for each treatment group
  turnover_jaccard<-c()
  nestedness_jaccard<-c()
  overall_jaccard<-c()
  turnover_jaccard_phy<-c()
  nestedness_jaccard_phy<-c()
  overall_jaccard_phy<-c()
  # turnover_bray<-c()
  # nestedness_bray<-c()
  overall_bray<-c()

  #For each treatment group...
  for (k in 1:length(types)){
    #Pull out the microbiomes from that treatment group and put them in their own phyloseq object
    rZOTUtemp<-prune_samples(sample_names(ZNC_physeq)[pickgroups[k,]],rZOTU)
    #Sum up the presences of each microbial taxon in the subsetted phyloseq object... count how many microbial taxa are present on as many or more individuals than the defined core fraction
    temp_community_matrix<-otu_table(rZOTUtemp)
    temp_tree<-phy_tree(rZOTUtemp)
    class(temp_community_matrix) <- "matrix"
    temp_betapart_jaccard<-beta.multi(sign(t(temp_community_matrix)),index.family='jaccard')
    temp_betapart_jaccard_phy<-phylo.beta.multi(sign(t(temp_community_matrix)),temp_tree,index.family='jaccard')
    temp_betapart_bray<-beta.multi.abund(t(temp_community_matrix),index.family='bray')

    turnover_jaccard<-c(turnover_jaccard,temp_betapart_jaccard$beta.JTU)
    nestedness_jaccard<-c(nestedness_jaccard,temp_betapart_jaccard$beta.JNE)
    overall_jaccard<-c(overall_jaccard,temp_betapart_jaccard$beta.JAC)
    turnover_jaccard_phy<-c(turnover_jaccard_phy,temp_betapart_jaccard_phy$phylo.beta.JTU)
    nestedness_jaccard_phy<-c(nestedness_jaccard_phy,temp_betapart_jaccard_phy$phylo.beta.JNE)
    overall_jaccard_phy<-c(overall_jaccard_phy,temp_betapart_jaccard_phy$phylo.beta.JAC)

    # turnover_bray<-c(turnover_bray,temp_betapart_bray$beta.BRAY.BAL)
    # nestedness_bray<-c(nestedness_bray,temp_betapart_bray$beta.BRAY.GRA)
    overall_bray<-c(overall_bray,temp_betapart_bray$beta.BRAY)
  }

  #Record the core counts for each treatment group from this particular bootstrap
  boot_tracker_turnover_jaccard<-rbind(boot_tracker_turnover_jaccard,turnover_jaccard)
  boot_tracker_nestedness_jaccard<-rbind(boot_tracker_nestedness_jaccard,nestedness_jaccard)
  boot_tracker_overall_jaccard<-rbind(boot_tracker_overall_jaccard,overall_jaccard)
  boot_tracker_turnover_jaccard_phy<-rbind(boot_tracker_turnover_jaccard_phy,turnover_jaccard_phy)
  boot_tracker_nestedness_jaccard_phy<-rbind(boot_tracker_nestedness_jaccard_phy,nestedness_jaccard_phy)
  boot_tracker_overall_jaccard_phy<-rbind(boot_tracker_overall_jaccard_phy,overall_jaccard_phy)

  # boot_tracker_turnover_bray<-rbind(boot_tracker_turnover_bray,turnover_bray)
  # boot_tracker_nestedness_bray<-rbind(boot_tracker_nestedness_bray,nestedness_bray)
  boot_tracker_overall_bray<-rbind(boot_tracker_overall_bray,overall_bray)

  }

#Name the bootstrap columns based on the treatment groups
colnames(boot_tracker_turnover_jaccard)<-types
#colnames(boot_tracker_turnover_bray)<-types
colnames(boot_tracker_nestedness_jaccard)<-types
# colnames(boot_tracker_nestedness_bray)<-types
 colnames(boot_tracker_overall_jaccard)<-types
 colnames(boot_tracker_turnover_jaccard_phy)<-types
 #colnames(boot_tracker_turnover_bray)<-types
 colnames(boot_tracker_nestedness_jaccard_phy)<-types
 # colnames(boot_tracker_nestedness_bray)<-types
 colnames(boot_tracker_overall_jaccard_phy)<-types
 colnames(boot_tracker_overall_bray)<-types

#Turn your bootstraps into two lists...
#One with the treatment group for each measurement
type<-c()
#And one with the core richness for each measurement
jaccard_turnover<-c()
jaccard_turnover_phy<-c()
#bray_turnover<-c()
jaccard_nestedness<-c()
jaccard_nestedness_phy<-c()
#bray_nestedness<-c()
jaccard_overall<-c()
jaccard_overall_phy<-c()
bray_overall<-c()

for (k in 1:length(types)){
  #A list of the treatment group for each pooled bootstrap
  type<-c(type,rep(types[k],bootstrap_no))
  #A list of the core richness for each pooled bootstrap
  jaccard_turnover<-c(jaccard_turnover,boot_tracker_turnover_jaccard[,k])
  jaccard_nestedness<-c(jaccard_nestedness,boot_tracker_nestedness_jaccard[,k])
  jaccard_overall<-c(jaccard_overall,boot_tracker_overall_jaccard[,k])
  jaccard_turnover_phy<-c(jaccard_turnover_phy,boot_tracker_turnover_jaccard_phy[,k])
  jaccard_nestedness_phy<-c(jaccard_nestedness_phy,boot_tracker_nestedness_jaccard_phy[,k])
  jaccard_overall_phy<-c(jaccard_overall_phy,boot_tracker_overall_jaccard_phy[,k])

#  bray_turnover<-c(bray_turnover,boot_tracker_turnover_bray[,k])
#  bray_nestedness<-c(bray_nestedness,boot_tracker_nestedness_bray[,k])
  bray_overall<-c(bray_overall,boot_tracker_overall_bray[,k])
}

#Put your richness and treatment group data into a dataframe
diversity_df<-data.frame(type,jaccard_turnover,jaccard_nestedness,jaccard_overall,bray_overall,jaccard_turnover_phy,jaccard_nestedness_phy,jaccard_overall_phy)#,bray_turnover,bray_nestedness)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm'))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test whether there are differences in BETA DIVERSITY between groups

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Jaccard Turnover
kw_jaccard_turnover<-kruskal.test(jaccard_turnover ~ type, data = diversity_df)
if (kw_jaccard_turnover$p.value<0.05){
  print('Jaccard Turnover')
  pw_jaccard_turnover<-pairwise.wilcox.test(diversity_df$jaccard_turnover, diversity_df$type,p.adjust.method = "BH")
  print(pw_jaccard_turnover)
}

#Jaccard Nestedness
kw_jaccard_nestedness<-kruskal.test(jaccard_nestedness ~ type, data = diversity_df)
if (kw_jaccard_nestedness$p.value<0.05){
  print('Jaccard Nestedness')
  pw_jaccard_nestedness<-pairwise.wilcox.test(diversity_df$jaccard_nestedness, diversity_df$type,p.adjust.method = "BH")
  print(pw_jaccard_nestedness)
}

#Jaccard Overall
kw_jaccard_overall<-kruskal.test(jaccard_overall ~ type, data = diversity_df)
if (kw_jaccard_overall$p.value<0.05){
  print('Jaccard Overall')
  pw_jaccard_overall<-pairwise.wilcox.test(diversity_df$jaccard_overall, diversity_df$type,p.adjust.method = "BH")
  print(pw_jaccard_overall)
}

#Unifrac Turnover
kw_jaccard_turnover_phy<-kruskal.test(jaccard_turnover_phy ~ type, data = diversity_df)
if (kw_jaccard_turnover_phy$p.value<0.05){
  print('Unifrac Turnover')
  pw_jaccard_turnover_phy<-pairwise.wilcox.test(diversity_df$jaccard_turnover_phy, diversity_df$type,p.adjust.method = "BH")
  print(pw_jaccard_turnover_phy)
}

#Unifrac Nestedness
kw_jaccard_nestedness_phy<-kruskal.test(jaccard_nestedness_phy ~ type, data = diversity_df)
if (kw_jaccard_nestedness_phy$p.value<0.05){
  print('Unifrac Nestedness')
  pw_jaccard_nestedness_phy<-pairwise.wilcox.test(diversity_df$jaccard_nestedness_phy, diversity_df$type,p.adjust.method = "BH")
  print(pw_jaccard_nestedness_phy)
}

#Unifrac Overall
kw_jaccard_overall_phy<-kruskal.test(jaccard_overall_phy ~ type, data = diversity_df)
if (kw_jaccard_overall_phy$p.value<0.05){
  print('Unifrac Overall')
  pw_jaccard_overall_phy<-pairwise.wilcox.test(diversity_df$jaccard_overall_phy, diversity_df$type,p.adjust.method = "BH")
  print(pw_jaccard_overall_phy)
}



#Bray Overall
kw_bray_overall<-kruskal.test(bray_overall ~ type, data = diversity_df)
if (kw_bray_overall$p.value<0.05){
  print('Core Richness')
  pw_bray_overall<-pairwise.wilcox.test(diversity_df$bray_overall, diversity_df$type,p.adjust.method = "BH")
  print(pw_bray_overall)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Make Violin plots of the various diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
p_value_ns<-'no'
pshow<-'star'


##########################Jaccard Overall Violin Plot#########################################################################################################################################

#Define violin plot
p_jaccard_overall <- ggplot(diversity_df, aes(x=type, y=jaccard_overall)) + geom_violin(aes(fill=type)) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_jaccard_overall<-p_jaccard_overall+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_jaccard_overall<-p_jaccard_overall+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_jaccard_overall<-p_jaccard_overall+scale_fill_manual(values=c("red", "magenta", "lightslateblue","blue"))
#Add boxplots inside the violins
p_jaccard_overall<-p_jaccard_overall+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_jaccard_overall<-p_jaccard_overall+annotate("text", size=6,x=3, y=0.93, label= "Kruskal-Wallis")
p_jaccard_overall<-p_jaccard_overall+annotate("text", size=6,x=3, y=0.927, label= paste("p = ",round(kw_jaccard_overall$p.value,4)))

#If there are significant differences in jaccard_overall between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_jaccard_overall$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-0.94
  y_step<-0.0025
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_jaccard_overall$p.value))){
    for (j in 1:k){
      if (rownames(pw_jaccard_overall$p.value)[k]!=colnames(pw_jaccard_overall$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_jaccard_overall$p.value[k,j]<0.05 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_jaccard_overall$p.value)[k])
          group2<-c(group2,colnames(pw_jaccard_overall$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_jaccard_overall$p.value[k,j])),6)
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
  stat.test_jaccard_overall<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_jaccard_overall<-p_jaccard_overall+stat_pvalue_manual(stat.test_jaccard_overall,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_jaccard_overall #SI Fig. 1.8a



##########################Unifrac Overall Violin Plot#########################################################################################################################################

#Define violin plot
p_jaccard_overall_phy <- ggplot(diversity_df, aes(x=type, y=jaccard_overall_phy)) + geom_violin(aes(fill=type)) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_jaccard_overall_phy<-p_jaccard_overall_phy+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_jaccard_overall_phy<-p_jaccard_overall_phy+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_jaccard_overall_phy<-p_jaccard_overall_phy+scale_fill_manual(values=c("red", "magenta", "lightslateblue","blue"))
#Add boxplots inside the violins
p_jaccard_overall_phy<-p_jaccard_overall_phy+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
#p_jaccard_overall_phy<-p_jaccard_overall_phy+annotate("text", size=6,x=2.5, y=0.835, label= "Kruskal-Wallis")
#p_jaccard_overall_phy<-p_jaccard_overall_phy+annotate("text", size=6,x=2.5, y=0.831, label= paste("p = ",round(kw_jaccard_overall_phy$p.value,4)))

#If there are significant differences in jaccard_overall_phy between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_jaccard_overall_phy$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-0.845
  y_step<-0.003
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_jaccard_overall_phy$p.value))){
    for (j in 1:k){
      if (rownames(pw_jaccard_overall_phy$p.value)[k]!=colnames(pw_jaccard_overall_phy$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_jaccard_overall_phy$p.value[k,j]<0.05 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_jaccard_overall_phy$p.value)[k])
          group2<-c(group2,colnames(pw_jaccard_overall_phy$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_jaccard_overall_phy$p.value[k,j])),6)
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
  stat.test_jaccard_overall_phy<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_jaccard_overall_phy<-p_jaccard_overall_phy+stat_pvalue_manual(stat.test_jaccard_overall_phy,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_jaccard_overall_phy #SI Fig. 1.8c & Main Fig. 1e


##########################Bray-Curtis Overall Violin Plot#########################################################################################################################################

#Define violin plot
p_bray_overall <- ggplot(diversity_df, aes(x=type, y=bray_overall)) + geom_violin(aes(fill=type)) + theme(axis.title.x = element_blank())
#Choose the size of font for the axes titles and labels
p_bray_overall<-p_bray_overall+theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_bray_overall<-p_bray_overall+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_bray_overall<-p_bray_overall+scale_fill_manual(values=c("red", "magenta", "lightslateblue","blue"))
#Add boxplots inside the violins
p_bray_overall<-p_bray_overall+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_bray_overall<-p_bray_overall+annotate("text", size=6,x=2.5, y=0.9, label= "Kruskal-Wallis")
p_bray_overall<-p_bray_overall+annotate("text", size=6,x=2.5, y=0.89, label= paste("p = ",round(kw_bray_overall$p.value,4)))

#If there are significant differences in bray_overall between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_bray_overall$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-0.925
  y_step<-0.008
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_bray_overall$p.value))){
    for (j in 1:k){
      if (rownames(pw_bray_overall$p.value)[k]!=colnames(pw_bray_overall$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_bray_overall$p.value[k,j]<0.05 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_bray_overall$p.value)[k])
          group2<-c(group2,colnames(pw_bray_overall$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_bray_overall$p.value[k,j])),6)
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
  stat.test_bray_overall<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_bray_overall<-p_bray_overall+stat_pvalue_manual(stat.test_bray_overall,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_bray_overall #SI Fig. 1.8b
