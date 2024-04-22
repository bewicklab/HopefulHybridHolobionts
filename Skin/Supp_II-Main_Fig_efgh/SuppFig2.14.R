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


#The number from each treatment group you want to pool to assess core richness for each bootstrap
sample_size<-12

#The number of bootstrap samples you want to do
bootstrap_no<-15

richness_SBluG_ino<-c()
richness_SBluG_neo<-c()
richness_SNW_neo<-c()
richness_SNW_marm<-c()

mylist_rich=list()

for (kk in 1:sample_size){
  
  set.seed(4)
  
  #Define the core fraction
  core<-kk/sample_size
  
  #Define the core fraction
  #core<-0.5
  
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
    allcount<-c()
    percentcore<-c()
    
    #For each treatment group...
    for (k in 1:length(types)){
      
      #Pull out the microbiomes from that treatment group and put them in their own phyloseq object
      rZOTUtemp<-prune_samples(sample_names(ZNC_physeq)[pickgroups[k,]],rZOTU)
      
      ##########################ASV Richness#########################################################################################################################################
      #Sum up the presences of each microbial taxon in the subsetted phyloseq object... count how many microbial taxa are present on as many or more individuals than the defined core fraction
      corecount<-c(corecount,length(which(rowSums(sign(otu_table(rZOTUtemp)))>=core*sample_size)))
      allcount<-c(corecount,length(which(rowSums(sign(otu_table(rZOTUtemp)))>=1)))
      percentcore<-c(percentcore,100*length(which(rowSums(sign(otu_table(rZOTUtemp)))>=core*sample_size))/length(which(rowSums(sign(otu_table(rZOTUtemp)))>=1)))
      
      
    }
    
    #Record the core counts and Faith's pd measures for each treatment group from this particular bootstrap
    boot_tracker<-rbind(boot_tracker,percentcore)
    
  }
  
  #Name the bootstrap columns based on the treatment groups
  colnames(boot_tracker)<-types
  
  richness_SBluG_ino<-c(richness_SBluG_ino,median(boot_tracker[,1]))
  richness_SBluG_neo<-c(richness_SBluG_neo,median(boot_tracker[,3]))
  richness_SNW_neo<-c(richness_SNW_neo,median(boot_tracker[,4]))
  richness_SNW_marm<-c(richness_SNW_marm,median(boot_tracker[,2]))
  
  
  
  
  #Turn your bootstraps into four lists...
  #One with the treatment group for each measurement
  type<-c()
  #One with the core richness for each measurement
  richness<-c()
  
  
  for (k in 1:length(types)){
    #A list of the treatment group for each pooled bootstrap
    type<-c(type,rep(types[k],bootstrap_no))
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
  if (is.na(kw_richness$p.value)){kw_richness$p.value<-1}
  if (kw_richness$p.value<0.05 && abs(max(diversity_df$richness)-min(diversity_df$richness))>1E-10){
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
  
  #If there are significant differences in richness between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
  if (kw_richness$p.value<0.05 && abs(max(diversity_df$richness)-min(diversity_df$richness))>1E-10){
    #Names of the groups you're comparing
    group1<-c()
    group2<-c()
    #p-value for the pairwise comparisons
    p.adj<-c()
    #locations of the p-value brackets
    ypos<-c()
    #  new_y<-22
    #  y_step<-2
    new_y<-1.1*max(diversity_df[,2])
    y_step<-0.1*(max(diversity_df[,2])-min(diversity_df[,2]))
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
    
  }
  
  p_richness <- NULL
  mylist_rich[[kk]] <- local({
    kk <- kk  
    
    #Define violin plot
    p_richness <- ggplot(diversity_df, aes(x=type, y=richness)) + geom_violin(aes(fill=type)) + theme(axis.title.x = element_blank())
    #Choose the size of font for the axes titles and labels
    p_richness<-p_richness+theme(axis.title = element_text(size = 20))+theme(axis.text.y = element_text(size = 15))+theme(axis.text.x = element_text(size = 12))
    #Choose the size of font for the legend title and lables
    p_richness<-p_richness+theme(legend.position="none")
    #p_richness<-p_richness+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
    #Choose the violin colors for each group
    p_richness<-p_richness+scale_fill_manual(values=c("red", "magenta", "lightslateblue","blue"))
    #Add boxplots inside the violins
    p_richness<-p_richness+geom_boxplot(aes(fill=type),width=0.1)
    #Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
    p_richness<-p_richness+annotate("text", -Inf, Inf, hjust = 0, vjust = 1, size=6, label= paste('richness:',paste0(round(core*100,2),'%'))) #x=3, y=0.9*min(diversity_df[,2]),
    #p_richness<-p_richness+annotate("text", size=6,x=3, y=5, label= paste("p = ",round(kw_richness$p.value,4)))
    if (kw_richness$p.value<0.05 && abs(max(diversity_df$richness)-min(diversity_df$richness))>1E-10){
      #Add the pairwise comparisons to your plot
      p_richness<-p_richness+stat_pvalue_manual(stat.test_richness,label=pdisplay,y.position=ypos,size=5)
    } 
    #Make your plot
    plot(p_richness)
    
  })
  
  #save richness plots to list
  myplots_richness <- append(mylist_rich, p_richness)
  
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Plot Results

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


############### Richness of the core microbiome ####################

library("gridExtra")
library("grid")

kaxis<-(1:1:sample_size)
rich_line_core <- ggplot() + theme_classic(16) + ylab("Core Richness % (Skin Genera)") + xlab("Minimum Hosts") +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) + expand_limits(y=0) + #setting x-axis interval and y-axis minimum
  geom_line(aes(x=kaxis, y=richness_SBluG_ino), colour='red') + geom_point(aes(x=kaxis, y=richness_SBluG_ino), size=3.5, colour='red') +
  geom_line(aes(x=kaxis, y=richness_SBluG_neo), colour='magenta') + geom_point(aes(x=kaxis, y=richness_SBluG_neo), size=3.5, colour='magenta') +
  geom_line(aes(x=kaxis, y=richness_SNW_neo), colour='lightslateblue') + geom_point(aes(x=kaxis, y=richness_SNW_neo), size=3.5, colour='lightslateblue') +
  geom_line(aes(x=kaxis, y=richness_SNW_marm), colour='blue') + geom_point(aes(x=kaxis, y=richness_SNW_marm), size=3.5, colour='blue') +
  theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))

rich_log_line_core<-ggplot() + theme_classic(16) + ylab("Log Core Richness % (Skin Genera)") + xlab("Minimum Hosts") +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) + expand_limits(y=0) + #setting x-axis interval and y-axis minimum
  geom_line(aes(x=if (length(richness_SBluG_ino[log(richness_SBluG_ino)>=0])==12) (kaxis) else head(kaxis,-(12-length(richness_SBluG_ino[log(richness_SBluG_ino)>=0]))), y=log(richness_SBluG_ino[log(richness_SBluG_ino)>=0])), colour='red') + 
  geom_point(aes(x=if (length(richness_SBluG_ino[log(richness_SBluG_ino)>=0])==12) (kaxis) else head(kaxis,-(12-length(richness_SBluG_ino[log(richness_SBluG_ino)>=0]))), y=log(richness_SBluG_ino[log(richness_SBluG_ino)>=0])), size=3.5, colour='red') +
  geom_line(aes(x=if (length(richness_SBluG_neo[log(richness_SBluG_neo)>=0])==12) (kaxis) else head(kaxis,-(12-length(richness_SBluG_neo[log(richness_SBluG_neo)>=0]))), y=log(richness_SBluG_neo[log(richness_SBluG_neo)>=0])), colour='magenta') + 
  geom_point(aes(x=if (length(richness_SBluG_neo[log(richness_SBluG_neo)>=0])==12) (kaxis) else head(kaxis,-(12-length(richness_SBluG_neo[log(richness_SBluG_neo)>=0]))), y=log(richness_SBluG_neo[log(richness_SBluG_neo)>=0])), size=3.5, colour='magenta') +
  geom_line(aes(x=if (length(richness_SNW_neo[log(richness_SNW_neo)>=0])==12) (kaxis) else head(kaxis,-(12-length(richness_SNW_neo[log(richness_SNW_neo)>=0]))), y=log(richness_SNW_neo[log(richness_SNW_neo)>=0])), colour='lightslateblue') + 
  geom_point(aes(x=if (length(richness_SNW_neo[log(richness_SNW_neo)>=0])==12) (kaxis) else head(kaxis,-(12-length(richness_SNW_neo[log(richness_SNW_neo)>=0]))), y=log(richness_SNW_neo[log(richness_SNW_neo)>=0])), size=3.5, colour='lightslateblue') +
  geom_line(aes(x=if (length(richness_SNW_marm[log(richness_SNW_marm)>=0])==12) (kaxis) else head(kaxis,-(12-length(richness_SNW_marm[log(richness_SNW_marm)>=0]))), y=log(richness_SNW_marm[log(richness_SNW_marm)>=0])), colour='blue') + 
  geom_point(aes(x=if (length(richness_SNW_marm[log(richness_SNW_marm)>=0])==12) (kaxis) else head(kaxis,-(12-length(richness_SNW_marm[log(richness_SNW_marm)>=0]))), y=log(richness_SNW_marm[log(richness_SNW_marm)>=0])), size=3.5, colour='blue') +
  theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))


gs_rich <- list(rich_line_core, rich_log_line_core, myplots_richness[[1]], myplots_richness[[2]], myplots_richness[[3]], myplots_richness[[4]], 
                myplots_richness[[5]], myplots_richness[[6]], myplots_richness[[7]], myplots_richness[[8]], 
                myplots_richness[[9]], myplots_richness[[10]], myplots_richness[[11]], myplots_richness[[12]])

lay_rich <- rbind(c(1,1,2,2),
                  c(1,1,2,2),
                  c(1,1,2,2),
                  c(3,4,5,6),
                  c(7,8,9,10),
                  c(11,12,13,14))

grid_rich <- grid.arrange(grobs = gs_rich, layout_matrix = lay_rich)

grid.draw(grid_rich) # interactive device

ggsave("Figure_SI2.14.png", grid_rich, height = 20, width = 20)


