require(tidyverse)
require(ggtree)
require(randomcoloR)
require(rlist)
require(tidytree)
require(ape)
require(phytools)
require(ggnewscale)
require(ggpubr)
require(cowplot)
require(ggtext)
require(circlize)
require(flipPlots)


##############################################################################################
#
# Phylogenetic Analysis per Host
#
##############################################################################################

# first draw the basic tree
core = read.tree(file = "Fastree_sept_rooted.tree")

# Now get host and localisation data
meta = read.csv("dataset/dataset_sample_2.csv",header =T)
rownames(meta)= as.character(meta$sample)
meta %>% group_by(phylogroup) %>% count(sort =T) %>% ungroup %>%  mutate(percentage = round(n/sum(n)*100,1)) %>% ggtexttable(rows = NULL)

metaHost = meta %>% mutate(Host = as.character(Host)) %>% select(Host)

colorHost = c("Human"="#000000",
              "Cattle"="#8DD3C7",
              "Chicken"="#e8eb34",
              "Pig" = "#BEBADA",
              "Dog" = "#FB8072",
              "Cat" ="#BC80BD",
              "Duck" = "#FDB462",
              "Goose" = "#B3DE69",
              "Horse" = "#FCCDE5",
              "Turkey" = "#D9D9D9",
              "Water" = "#80B1D3",
              "Animal" = "#b14624")
names(colorHost)

p = ggtree(core,layout="circular",size = 0.2,branch.length = 0.01)
#p = ggtree(core,size = 0.2,branch.length = 0.01)
pg = as.data.frame(meta$phylogroup)
colnames(pg)="phylogroup"
rownames(pg) = rownames(meta)
pg$phylogroup= as.character(pg$phylogroup)
pg$node = rownames(pg)
p$data = p$data %>% left_join(pg,by=c("label"="node"))

phyancestors = p$data$phylogroup %>% unique
names(phyancestors)=phyancestors
phyancestors = phyancestors[1:7]

for(i in names(phyancestors) %>% unique)
{
  p$data %>% filter(phylogroup==i) %>% select(label) %>% as.vector -> F_labels
  phyancestors[i] = findMRCA(core,F_labels$label)
}
dataphy = phyancestors %>% as.data.frame %>% mutate(phylogroup = names(phyancestors))
names(dataphy) = c("id","phylogroup")
dataphy$id =as.integer(dataphy$id)
dataphy

colphy = RColorBrewer::brewer.pal(n=7,name="Set2")
dataphy$color = colphy

dataphy %>% glimpse

test = randomColor(count=7)
test

p = p + geom_highlight(data = dataphy,mapping=aes(node = id,fill = color))  +
  scale_fill_brewer(palette = "Set2") + guides(fill=F)
p = p + new_scale_fill()

metaHost = metaHost %>% mutate(Host2= ifelse(Host=="Human","Human",ifelse(Host=="Water","Water","Animal")))

p = gheatmap(p,metaHost[,c("Host","Host2"), drop=F], width=0.2,font.size=8, colnames = F,
             hjust = 0, offset = 0.0) +  scale_fill_manual(name="Host",breaks = names(colorHost),values=colorHost) +
  scale_color_manual(values=test) + guides(color = guide_legend(override.aes = list(size=3))) +
  theme(legend.box = "vertical",legend.position = "top")



plot_phylogeny = p + geom_cladelab(data = dataphy,mapping=aes(node =id,label =phylogroup,color=color),
                  offset=0.13,offset.text=0.03,show.legend=F,fontsize=3) + scale_color_brewer(palette="Set2")  +
                  guides(color=F) + theme(legend.key.size = unit(0.4, 'cm')) +
  theme(
    plot.margin=unit(c(0, 0, 0, 0), units="line"),
    legend.position="top",
    legend.margin=unit(0, "lines"))
plot_phylogeny

ggsave("phylogeny_Host.pdf")

#####################################################################################
# Add genetic distance between Hosts
#####################################################################################

colorHost = c("Human"="#000000",
              "Cattle"="#8DD3C7",
              "Chicken"="#e8eb34",
              "Pig" = "#BEBADA",
              "Dog" = "#FB8072",
              "Cat" ="#BC80BD",
              "Duck" = "#FDB462",
              "Goose" = "#B3DE69",
              "Horse" = "#FCCDE5",
              "Turkey" = "#D9D9D9",
              "Water" = "#80B1D3")

## ADD ONE FOR INTERHOSTS DATA
colorHost_inter = append(colorHost,c("Inter"="white"))

# first draw the basic tree
core = read.tree(file = "dataset/Fastree_sept_rooted.tree")
## leave some space for inlay

## Metadata parsing
metadata = read.csv("dataset/Metadata.csv")
metadata = metadata %>%
  mutate(Host = ifelse(Host =="Poultry","Chicken",as.character(Host)))
meta = metadata %>%  select(sample = True_indiv_name,Host,Numero_foyer,Fokontany)
host = meta %>% select(sample,Host)
rownames(host) = meta$sample

geo = meta %>% select(sample,Numero_foyer,Fokontany)
geo
rownames(geo) = meta$sample
host %>% View
geo %>%  View
## Get pairwise distances from tree
distTips = cophenetic.phylo(core)
distTips %>% dim
# Deal with matrix to keep only upper diag values (avoid duplicated values)
# diag = T to also remove values on the diagonal (distance between same samples)
lowdiag = lower.tri(distTips,diag=T)
distTips[lowdiag]=NA
distTips = distTips %>% as.table %>% as.data.frame %>% filter(!is.na(Freq))

distTips


# it checks out, we have the correct number of unique pairs I think)
distTips %>% dim
((510*510-510)/2)
# Now make a graph with color according to type (within / between)

distTips %>%  glimpse

######### Dist_tips according to host

distTipsHost = distTips %>% left_join(host,by=c("Var1"="sample")) %>%
  mutate(host1 = Host) %>% select(-Host) %>%
  left_join(host,by=c("Var2"="sample")) %>%
  mutate(host2 = Host) %>% select(-Host) %>%
  mutate(type=ifelse(host1==host2,host1,"Between Host")) %>%
  select(-host1,-host2) %>% mutate(InterIntra = ifelse(type=="Between Host","Between Host","Within Host"))

distancePlot = distTipsHost %>%  mutate(type = fct_infreq(type)) %>%
  ggplot() +
  geom_violin(aes(x=type,y=Freq,fill=type,),
              color="black",show.legend= F) +
  geom_boxplot(aes(x=type,y=Freq),width=.15,outlier.shape=NA) +
  geom_point(data = distTipsHost %>% filter(type == "Horse"),
             aes(x=type,y=Freq,fill=type),pch=21,color="black",show.legend=F,size=3) +
  scale_fill_manual(values = colorHost_inter) +
  geom_text(data = distTipsHost %>% group_by(type) %>% summarize(n=n()),aes(x=type,y=-0.02,label=n)) +
  theme_pubr() + xlab("Host pair") + ylab("Genetic distances") + theme(axis.text.x = element_text(angle = 45,hjust = 1))

distancePlot


globalDistancePlotHost  = distTipsHost %>% ggplot() + geom_density(aes(x=Freq,color=InterIntra)) +
  theme_pubr()  + xlab("Genetic distances")  + scale_color_discrete(name="")
globalDistancePlotHost




# ggdraw() + draw_plot(plot_phylogeny,x = 0,y = 0,width = 0.7,height = 1) +
#   draw_plot(distancePlot,x = 0.7,y = 0,width = 0.3,height = 0.5)  +
#   draw_plot(globalDistancePlotHost ,x = 0.7,y = 0.5,width = 0.3,height = 0.5)
#

plot_phylogeny
ggsave("Phylogenie_host.svg")
distancePlot
ggsave("DistanceHost.svg",width = 7)
globalDistancePlotHost
ggsave("DistanceHostglob.svg")


###### Figure 3 (Phylogenie + host distance)  is then assembled based on this three files



###### Now for the supplementary figures and analysis related to household:



### First get the distances

distTips_foyer = distTips %>% left_join(geo,by=c("Var1"="sample")) %>%
  mutate(Foyer1 = Numero_foyer) %>% select(-Numero_foyer) %>%
  left_join(geo,by=c("Var2"="sample")) %>%
  mutate(Foyer2 = Numero_foyer) %>% select(-Numero_foyer) %>%
  mutate(type=ifelse(Foyer1==Foyer2,Foyer1,"Between Foyer")) %>%
  mutate(InterIntra = ifelse(type=="Between Foyer","Between Foyer","Within Foyer")) %>%
  mutate(IntraFoyer = ifelse(type=="Between Foyer","Between Household",Foyer1))



p = ggtree(core,layout="circular",size = 0.2)

meta = read.csv("dataset/Metadata.csv",header =T)
meta %>% glimpse
metad = meta %>% select(sample= True_indiv_name,Numero_foyer)
rowsNames = metad$sample

metad = metad %>% select(-sample)
rownames(metad)= rowsNames

metad = metad %>%
  mutate(Household = as.factor(Numero_foyer))
metad %>% head
metad = metad %>% select(-Numero_foyer)
metad %>% glimpse
metad$Household %>% unique
pal = randomcoloR::distinctColorPalette(k=70)
pal
metad %>% row.names

p1 = gheatmap(p,metad, width=0.2, font.size=8, colnames = F,
              hjust = 0) + scale_fill_manual(name="Household",values = pal) +
  guides(fill="none")
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm")) +
  theme_tree() + theme(legend.position = NA)
p1
ggsave("Phylogenie_household.svg")

p2 = distTips_foyer %>%
  ggplot() +
  geom_boxplot(aes(x=IntraFoyer,y=Freq,fill=IntraFoyer),show.legend =F,outlier.shape=NA) +
  scale_fill_manual(values = pal) +
  theme_pubr() + theme(axis.text.x = element_text(angle = 90)) + xlab("") + ylab("Paiwise genetic distances")
p2
ggsave("Distance_genetique_foyers.svg",width = 12)
#### Will have to be edited in inkscape to keep only "Between Household"
p3  = distTips_foyer %>%
  mutate(InterIntra=fct_recode(InterIntra,"Between households"="Between Foyer",
                               "Within households" = "Within Foyer")) %>%
  ggplot() +
  geom_density(aes(x=Freq,color=InterIntra),alpha = 0.8,size = 0.8) +
  theme_pubr()  + xlab("Pairwise genetic distances")  +
  scale_color_brewer(name="",type="qual", palette=1) + theme_pubr()
p3
ggsave("GloblaDistance_household.svg")
p4  = distTips_foyer %>%
  mutate(InterIntra=fct_recode(InterIntra,"Between households"="Between Foyer",
                               "Within households" = "Within Foyer")) %>%
  ggplot() +
  geom_violin(aes(y=Freq,x=InterIntra,fill=InterIntra)) +
  geom_boxplot(aes(y=Freq,x=InterIntra),outlier.shape = NA,width =0.1) +
  theme_pubr()  + ylab("Pairwise genetic distances")  +
  scale_fill_brewer(name="",type="qual", palette=1) + xlab("") + theme(axis.text.x=element_blank())
p4
ggsave("ViolinPlot_househould.svg")

# to assemble in order to produce figure

############################################################################################

###### Now for the supplementary figures and analysis related to fokontany:

distTipsGeoFokontany = distTips %>% left_join(geo,by=c("Var1"="sample")) %>%
  mutate(Fokontany1 = Fokontany) %>% select(-Fokontany) %>%
  left_join(geo,by=c("Var2"="sample")) %>%
  mutate(Fokontany2 = Fokontany) %>% select(-Fokontany) %>%
  mutate(type=ifelse(Fokontany1==Fokontany2,Fokontany1,"Between Fokontany")) %>%
  select(-Fokontany1,-Fokontany2) %>% mutate(InterIntra = ifelse(type=="Between Fokontany","Between Fokontany","Within Fokontany"))

pal=randomcoloR::distinctColorPalette(k=2)

distTipsGeoFokontany %>%
  ggplot() +
  geom_violin(aes(x=InterIntra,y=Freq,fill=InterIntra)) +
  geom_boxplot(aes(x=InterIntra,y=Freq),show.legend =F,outlier.shape=NA,width=0.05) +
  scale_fill_manual(name=NULL,values = pal) +
  theme_pubr() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + xlab("") + ylab("Paiwise genetic distances")
ggsave("figures/Violin_DistanceFokontanyGlobal.svg")


# distTipsGeoFokontany %>% glimpse
# distTipsGeoFokontany %>%
#   mutate(InterIntra=fct_recode(InterIntra,"Between households"="Between Foyer",
#                                "Within households" = "Within Foyer")) %>%
#   ggplot() +
#   geom_violin(aes(y=Freq,x=InterIntra,fill=InterIntra)) +
#   geom_boxplot(aes(y=Freq,x=InterIntra),outlier.shape = NA,width =0.1) +
#   theme_pubr()  + ylab("Pairwise genetic distances")  +
#   scale_fill_brewer(name="",type="qual", palette=1) + xlab("") + theme(axis.text.x=element_blank())



globalDistancePlotGeoFokontany  = distTipsGeoFokontany %>% ggplot() + geom_density(aes(x=Freq,color=InterIntra),alpha = 0.8,size = 0.8) +
  theme_pubr()  + xlab("Genetic distances")  + scale_color_brewer(name="",type="qual",
                                                                     palette=1)
globalDistancePlotGeoFokontany
ggsave("GlobaldistanceFokontany.svg")

p = ggtree(core,layout="circular",size = 0.2)
meta = read.csv("dataset/Metadata.csv",header =T)
meta %>% glimpse
metad = meta %>% select(sample= True_indiv_name,Fokontany)
rowsNames = metad$sample
metad = metad %>% select(-sample)
rownames(metad)= rowsNames
metad = metad %>%
  mutate(Fokontany = as.factor(Fokontany))
metad %>% head
metad$Fokontany %>% unique
pal = randomcoloR::distinctColorPalette(k=10)
gheatmap(p,metad, width=0.2, font.size=8, colnames = F,
         hjust = 0) + scale_fill_manual(name="Fokontany",values = pal) +
  theme(legend.direction = "vertical",legend.position ="left") +
  theme_tree()
ggsave("Phylogeny_fokontany.svg",width=6)

distTipsGeoFokontany %>%
  ggplot() +
  geom_violin(aes(x=type,y=Freq,fill=type),show.legend =F) +
  geom_boxplot(aes(x=type,y=Freq),show.legend =F,outlier.shape=NA,width = 0.1) +
  scale_fill_manual(values = pal) +
  theme_pubr() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + xlab("") + ylab("Paiwise genetic distances")
ggsave("Violin_foko_distance.svg")

##### to assemble in order to produce figure sup



#################################################################################
################ Now let's do the heatmaps (genes and ST by Host)
#################################################################################


###### Fist resistance genes
data  = readxl::read_xlsx("dataset/ResistanceGenesWithMetadataWideFormat.xlsx")
data = data %>%  mutate(Sample = True_indiv_name)
data %>% select(Sample) %>% unique %>% count
######### get list of ecoli samples #######
GoodSamples = read.table("dataset/Sample_510.txt")
samples = GoodSamples %>% select(V2)
colnames(samples)="id"
######## save the good ones
data = data %>% filter(True_indiv_name %in% samples$id)
df = data %>% select(Sample,Host,Numero_foyer,Fokontany,mdf.A.:blaTEM.32) %>%
  gather("Resistance.gene","presence",mdf.A.:blaTEM.32)
df$Abtype = ifelse(substring(df$Resistance.gene,1,6)=="blaCTX" , "CTX" ,
                   ifelse(substring(df$Resistance.gene,1,6)=="blaSHV","SHV",
                          ifelse(substring(df$Resistance.gene,1,6)=="blaTEM","TEM","Others")))
# Prepare dataset for matrix, change sample and gene size here
gene_counts = df %>%  filter(presence ==1 ) %>%
  group_by(Host,Resistance.gene) %>% count
gene_counts = gene_counts %>% pivot_wider(names_from=Host,values_from=n)
gene_counts = gene_counts %>%  pivot_longer(cols=2:12,names_to="Host",values_to="n")
hg = df %>% filter(presence ==1 ) %>% group_by(Resistance.gene) %>% count %>% select(Resistance.gene)
counts = df[row.names(df %>% select(Sample) %>% unique),"Host"] %>% table %>% as.data.frame
counts %>% glimpse
counts = rbind(counts,data.frame("Host"="All",Freq = sum(counts$Freq)))
colnames(counts) = c("Host","counts")
hs = counts

formatgene = function(genename){
  tmp  = genename %>% gsub("\\.","-", . ) %>%  gsub("bla","",.)
  res = paste("<i>bla</i><sub>",tmp,"</sub>",sep="")
  return(res)
}
######## gyrb & parC

readxl::read_xlsx("dataset/parC_QRDR_AA.xlsx") -> parC
splitsample = function(test){
  return(str_split(test,".scfd")[[1]][1])
}
parC %>% rowwise %>% mutate(sample = splitsample(isolate)) %>% select(sample,phenotype) -> parC
readxl::read_xlsx("dataset/Metadata.xlsx") -> meta
meta %>%
  left_join(parC,by = c("Reads_names"="sample")) %>%
  select(sample = Reads_names,Host,phenotype) -> parC

parC = parC %>% mutate(mut = ifelse(phenotype=="WT" | is.na(phenotype),0,1)) %>%
  group_by(Host) %>% summarize(Resistance.gene = "parC",n = sum(mut))
parC = parC[c(2,1,3)]
parC
parC = parC %>% mutate(Host = ifelse(Host=="Poultry","Chicken",Host)) %>% droplevels()
readxl::read_xlsx("dataset/gyrA_QRDR_AA.xlsx") -> gyrB
splitsample = function(test){
  return(str_split(test,".scfd")[[1]][1])
}
gyrB %>% rowwise %>% mutate(sample = splitsample(isolate)) %>% select(sample,phenotype) -> gyrB
readxl::read_xlsx("dataset/Metadata.xlsx") -> meta
meta %>%
  left_join(gyrB,by = c("Reads_names"="sample")) %>%
  select(sample = Reads_names,Host,phenotype) -> gyrB
gyrB= gyrB %>% mutate(mut = ifelse(phenotype=="WT" | is.na(phenotype),0,1)) %>%
  group_by(Host) %>% summarize(Resistance.gene = "gyrB",n = sum(mut))
gyrB = gyrB[c(2,1,3)]
gyrB = gyrB %>% mutate(Host = ifelse(Host=="Poultry","Chicken",Host)) %>% droplevels()
gene_counts = rbind(as.data.frame(gene_counts),as.data.frame(gyrB),as.data.frame(parC))
gene_counts = gene_counts %>% filter(Host!="Poultry")
gene_counts = gene_counts %>% mutate(n=ifelse(n==0,NA,n)) %>% droplevels()

all = gene_counts %>% group_by(Resistance.gene) %>% summarize(Host="All",n=sum(n,na.rm=T))

rbind(as.data.frame(all),as.data.frame(gene_counts)) %>%
  mutate(type = ifelse(Host == "All","All","Host")) %>%
  filter(grepl("bla",Resistance.gene) | Resistance.gene=="parC" | Resistance.gene=="gyrB") %>% left_join(counts,by="Host") %>%
  mutate(mutype = ifelse(grepl("bla",Resistance.gene),"betalactamase","other")) %>%
  mutate(freq = n/counts*100) %>% ungroup %>% filter(!is.na(n)) %>%
  mutate(Host = paste(Host,"  <br><span style='font-size:9pt'>(N=",counts,")</span>",sep="") ) %>%
  mutate(Host = fct_reorder(Host,counts,unique,.desc = T)) %>%
  mutate(Resistance.gene = formatgene(Resistance.gene)) %>%
  mutate(Resistance.gene = fct_reorder(Resistance.gene,n,sum,.desc=F)) %>%
  ggplot() + geom_tile(aes(x=Host,y=Resistance.gene,fill=freq)) +
  geom_text(aes(x=Host, y=Resistance.gene, label=round(freq,1))) +
  theme_bw() + scale_fill_viridis_c(name="Percentage",direction=-1) +
  ylab("Resistance genes") + xlab("Hosts") +
  theme_minimal() + theme(panel.grid=element_blank()) +
  theme(panel.grid = element_blank(),text= element_text(family = "Helvetica",size=14)) +
  theme(axis.text.y = element_markdown(hjust = 0),axis.text.x = element_markdown()) +
  facet_grid(rows = vars(mutype),cols = vars(type),scales = "free",space="free") + xlab(NULL)
ggsave("Figure_heatmap_prevalence_gene_hote.svg",height = 7,width=10)
ggsave("Figure_heatmap_prevalence_gene_hote.pdf",height = 7,width=10)
ggsave("Figure_heatmap_prevalence_gene_hote.png",height = 7,width=10)

######### No association genes / host
gene_counts %>%
  filter(grepl("bla",Resistance.gene) |
           Resistance.gene=="parC" |
           Resistance.gene=="gyrB") %>%
  pivot_wider(values_fill = 0,values_from = n,names_from=Resistance.gene) -> genedf
genedf %>% select(-Host) %>% as.matrix -> genematrix
genematrix[is.na(genematrix)] = 0
rownames(genematrix) = genedf$Host
chisq.test(genematrix,simulate.p.value = T)


################### Now for the ST Heatmap #######################

ST = read.csv("dataset/ST_data_final.txt",sep="\t",h=F,stringsAsFactors = F) %>% select(V1,V2)
colnames(ST) = c("sample","ST")
samples$id2 = gsub("_.+","",samples$id)
ST %>% filter(sample %in% samples$id2) -> ST


data %>% select(True_indiv_name,Host) %>% mutate(True_indiv_name =gsub("_.+","", True_indiv_name)) %>%
  left_join(ST,by=c("True_indiv_name"="sample")) %>% group_by(ST,Host) %>% count -> ST_counts
ST_counts$ST %>% unique

ST_counts %>% group_by(ST) %>% summarize(Host = "All",n=sum(n)) -> All_count
All_count
ST_counts = rbind(as.data.frame(ST_counts),as.data.frame(All_count))
ST_counts %>% group_by(Host) %>% summarize(total = sum(n)) %>% ungroup -> ST_host_counts

dev.off()
ST_counts %>%
  mutate(type = ifelse(Host =="All","All","Host")) %>%
  #filter(Host!="Horse") %>%
  left_join(ST_host_counts,by="Host") %>% ungroup %>%
  mutate(Host = paste(Host,"  <br><span style='font-size:9pt'>(N=",total,")</span>",sep="") ) %>%
  mutate(Host = fct_reorder(factor(Host),n,sum,.desc=T)) %>%
  mutate(ST = fct_reorder(factor(ST),n,sum,.desc = F)) %>%
  mutate(ST = fct_recode(ST,Unknown="0")) %>%
  mutate(ST = fct_lump_n(ST,25)) %>%
  mutate(ST = fct_relevel(ST, "Other", after = 0)) %>%
  mutate(freq = n/total*100) %>%
  ggplot() + geom_tile(aes(x=Host,y=ST,fill=freq)) +
  geom_text(aes(x=Host, y=ST, label=round(freq,1)),check_overlap = TRUE) +
  theme_bw() + scale_fill_viridis_c(option = "plasma",name="Percentage",direction=-1) +
  ylab("Sequence Type") + xlab("Hosts") +
  theme_minimal() + theme(panel.grid=element_blank()) +
  theme(panel.grid = element_blank(),text= element_text(family = "Helvetica",size=14)) +
  theme(axis.text.y = element_markdown(hjust = 0),axis.text.x = element_markdown()) +
  facet_grid(cols = vars(type),scales = "free",space="free") + xlab(NULL)
ggsave("Figure_ST_frequency_host.pdf",height = 7,width=10)
ggsave("Figure_ST_frequency_host.svg",height = 7,width=10)
ggsave("Figure_ST_frequency_host.png",height = 7,width=10)

#################################################################################

# NO association between ST and Host

ST_counts %>% filter(Host!="All") %>%
  pivot_wider(values_fill=0 ,values_from = n,names_from = Host) %>%
  as.data.frame -> ST_dataframe
ST_matrix = as.matrix(ST_dataframe)
ST_matrix[,-1] -> ST_matrix
rownames(ST_matrix) = ST_dataframe$ST
chisq.test(ST_matrix,simulate.p.value = T)



##########################################################################
#       Now the pieChart of clusters !
##########################################################################

clust = read.table("dataset/cluster_phydelity_k2_sol0_Fastree_sept_rooted.txt",header=T)
clust %>%  glimpse
clust %>% select(CLUSTER) %>% unique %>% dim

# Add infos in SNP number
distanceSNP = read.csv2(file="dataset/Distances_in_cluster_phydelity_norecomb.txt",sep=" ",h=T)
colnames(distanceSNP) = c("Cluster","TAXA1","TAXA2","SNP")
distanceSNP %>%  glimpse

distanceSNP %>%
  group_by(Cluster) %>%
  summarise(meanSNP = mean(SNP),maxSNP=max(SNP))

clust = clust %>% left_join(distanceSNP %>%
                              group_by(Cluster) %>%
                              summarise(meanSNP = mean(SNP),maxSNP=max(SNP)) ,by=c("CLUSTER"="Cluster"))
# filter cluster with a #SNP > 20 in any transmission
clust_toremove = distanceSNP %>% filter(SNP>20) %>% select(Cluster) %>% unique
clust = clust %>% filter(!CLUSTER %in% clust_toremove$Cluster )

meta = readxl::read_xlsx("dataset/Metadata.xlsx")
df = meta %>%  select(True_indiv_name,Host,Numero_foyer) %>%
  left_join(clust,by = c("True_indiv_name" = "TAXA"))

df %>% select(CLUSTER) %>% unique %>% dim

colorHost = c("Human"="#000000",
              "Cattle"="#8DD3C7",
              "Poultry"="#e8eb34",
              "Pig" = "#BEBADA",
              "Dog" = "#FB8072",
              "Cat" ="#BC80BD",
              "Duck" = "#FDB462",
              "Goose" = "#B3DE69",
              "Horse" = "#FCCDE5",
              "Chicken" = "#D9D9D9",
              "Water" = "#80B1D3")

######### get list of ecoli samples #######
GoodSamples = read.table("dataset/Sample_510.txt")
samples = GoodSamples %>% select(V2)
colnames(samples)="id"
samples$id2 = gsub("_.+","",samples$id)

######## get list of bad samples
df %>% filter(!True_indiv_name %in% samples$id2)
######## save the good ones
df = df %>% filter(True_indiv_name %in% samples$id)
df %>% glimpse
df %>% filter(CLUSTER!="NA") %>% select(Numero_foyer) %>% unique %>% dim
df %>% write.csv("dataset/CompositionClusterSup20.csv")

df[which(df$Host=="Poultry"),]$Host="Chicken"


left_join(df %>% group_by(Host) %>% count() %>% mutate(all = n)  %>% select(-n),
          df %>% filter(CLUSTER!="NA") %>% group_by(Host) %>% count() %>% mutate(cluster = n)  %>% select(-n)) %>%
  gather(key="type",value = "count",-"Host" ) %>%
  ggplot() + geom_bar(aes(x=Host,y=count,fill = Host),stat ="identity") + facet_grid(type ~ .) +
  scale_fill_manual(values = colorHost) + theme_pubr()
ggsave("GlobalCompositionCluster.svg",width =6)

###### shamelessly stolen from some internet page to make nice looking pie chart (yuck) in R
cp <- coord_polar("y")
cp$is_free <- function() TRUE


meanSNPs = df %>% filter(CLUSTER!="NA") %>% select(CLUSTER,meanSNP) %>% unique
labmeanSNPS = as.character(signif(meanSNPs$meanSNP,digits = 3))
names(labmeanSNPS) = meanSNPs$CLUSTER
library(extrafont)
font_import()
df %>% filter(CLUSTER !="NA") %>%
  mutate(CLUSTER = as.factor(CLUSTER)) %>% select(CLUSTER,Numero_foyer) %>% unique %>%
  group_by(CLUSTER) %>%  count() %>% transmute(nfoyer = n) -> clusterfoyer
pie  = function(nf)
{
df %>%
  filter(CLUSTER !="NA") %>%
  mutate(CLUSTER = as.factor(CLUSTER)) %>%
  group_by(CLUSTER,Host,meanSNP) %>% count %>%
  ungroup %>% group_by(CLUSTER) %>% mutate(size = sqrt(sum(n))) %>% ungroup %>%
  mutate(CLUSTER = fct_reorder(CLUSTER,meanSNP,mean)) %>% left_join(clusterfoyer,by="CLUSTER") %>%
  filter(nfoyer ==nf) %>%
  ggplot(aes(group=CLUSTER,x=size/2, y=n,width=size, fill=Host))+
  geom_bar(width = 1, stat = "identity")+
  scale_fill_manual(values = colorHost) + guides(fill=F) +
  facet_wrap(~ CLUSTER,ncol = 13,drop=T,scales="free_y",
             labeller = labeller(CLUSTER = labmeanSNPS)) +
  cp +
  theme_void() +
  theme(strip.background = element_blank(),
        #strip.text.x = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.ticks=element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1,legend.position = "top",
        text = element_text(family = "Times"))
}

library(egg)
library(grid)
p1_fixed <- set_panel_size(pie(1),
                          width  = unit(0.7, "cm"),
                          height = unit(0.7, "cm"))
p2_fixed <- set_panel_size(pie(2),
                           width  = unit(0.7, "cm"),
                           height = unit(0.7, "cm"))
p3_fixed <- set_panel_size(pie(3),
                           width  = unit(0.7, "cm"),
                           height = unit(0.7, "cm"))
p4_fixed <- set_panel_size(pie(4),
                           width  = unit(0.7, "cm"),
                           height = unit(0.7, "cm"))
p5_fixed <- set_panel_size(pie(5),
                           width  = unit(0.7, "cm"),
                           height = unit(0.7, "cm"))

grid.newpage()
g = grid.arrange(grobs = list(p1_fixed,p2_fixed,p3_fixed,p4_fixed,p5_fixed),ncol=1,nrow=5)
ggsave("PieCharts.svg",g)
#### OK this one is tricky, best to edit it in inkscape to produce the figures but this a good basis


require(EMT)

allpop = df %>% group_by(Host) %>%
  count %>%
  mutate(freq = n/sum(.$n)) %>% select(-n)


multitest = function(clustid){
  testpop = allpop %>% left_join(df %>%
                                   filter(CLUSTER==clustid) %>%
                                   select(Host) %>%
                                   group_by(Host) %>%
                                   count,by="Host",) %>%
    mutate(n=ifelse(is.na(n),0,n))
  testpopstat  = multinomial.test(testpop$n,testpop$freq)
  foyerclust = df %>%  filter(CLUSTER==clustid) %>% select(Numero_foyer)
  foyerpop = df %>% filter(Numero_foyer %in% foyerclust$Numero_foyer) %>% group_by(Host) %>%
    count %>%
    mutate(freq = n/sum(.$n)) %>% select(-n)
  testfoyer = foyerpop %>% left_join(df %>%
                                       filter(CLUSTER==clustid) %>%
                                       select(Host) %>%
                                       group_by(Host) %>%
                                       count,by="Host",) %>%
    mutate(n=ifelse(is.na(n),0,n))
  testfoyerstat =  multinomial.test(testfoyer$n,testfoyer$freq)
  return(c(pop = testpopstat$p.value,foyer = testfoyerstat$p.value))
}

pvalcluster = df %>% filter(CLUSTER!="NA") %>%  group_by(CLUSTER) %>% summarize(pvalpop = multitest(CLUSTER)["pop"],
                                                                                pvalfoy = multitest(CLUSTER)["foyer"])

pvalcluster %>% filter(pvalpop < 0.05)

pvalcluster %>% dim
annotation=df %>% filter(CLUSTER!="NA") %>%
  group_by(CLUSTER,meanSNP,maxSNP) %>%
  summarize(nHousehold = length(unique(Numero_foyer)),
            nHost = length(Host),
            Households = paste0(Numero_foyer,collapse="",sep=","),
            Samples = paste0(True_indiv_name,collapse="",sep=","),
            HostType = paste0(Host,collapse="",sep=",") ) %>%
  left_join(pvalcluster)

write.csv2(annotation,"tables/table_sup_Clusters.csv")



##########################################################################
#       Now just the finisher figure : The Sankey Diagram !
##########################################################################


####### Here we will combine plasmids data with Host data


ST_data = read.csv2("dataset/ST_data_final.txt",sep="\t",header=F)
colnames(ST_data) = c("Sample","ST")

ST_data
stpla = read.csv2("dataset/ST_plasm.csv")
stpla$Sample = gsub("_.+","",stpla$Sample)

stpla = stpla %>% select(Sample,ST2,Host,Resgenes,groupeinc) %>%
  na.omit %>% mutate(ST2= as.factor(ST2)) %>% left_join(ST_data,by="Sample")

tt = stpla %>% select(Sample,ST,Host,Resgenes,groupeinc) %>%
  na.omit %>% mutate(ST= as.factor(ST)) %>%  filter(ST!= "Unknown")

df = stpla %>% select(ST,Host,Resgenes,groupeinc) %>%
  na.omit %>% mutate(ST= as.factor(ST)) %>%
  mutate(ST = fct_lump(ST,n=30)) %>%
  group_by_all() %>%  count
dim(df)

stpla %>% select(ST,Host,Resgenes,groupeinc) %>%
  na.omit %>% mutate(ST= as.factor(ST)) %>%
  mutate(ST = fct_lump(ST,n=30)) %>% group_by(Resgenes,groupeinc) %>% count(sort=T) %>%
  write.table("Association_Resgenes_incgroup.csv",sep=",",row.names = F)

df[,c(2,1,4,3)]
df$Host = ifelse(df$Host=="Poultry","Chicken",df$Host)
df$ST = as.factor(ifelse(df$ST==0,"Unknown",df$ST))
sankey = SankeyDiagram(max.categories = 31,data = df[,c(2,1,4,3)],weights = df$n,link.color = "Source",label.show.percentages = T,label.show.varname=F)
sankey

# install.packages("webshot")
# webshot::install_phantomjs()

htmlwidgets::saveWidget(sankey, file ="sankey1.html") # save html widget
webshot::webshot("sankey1.html",file="sankey1.pdf")


############# Chi.square test for composition ###########
TODO

dataset = read_csv2("dataset/Table_info_samples.csv")
dataset %>% select(Host,ST) %>% table %>% chisq.test(simulate.p.value = TRUE)
dataset %>% select(Host,phylogroup) %>% table %>% chisq.test(simulate.p.value = TRUE)
dataset %>% select(Fokontany,phylogroup) %>% table %>% chisq.test(simulate.p.value = TRUE)
dataset %>% select(Fokontany,ST) %>% table %>% chisq.test(simulate.p.value = TRUE)
dataset %>% select(Household,phylogroup) %>% table %>% chisq.test(simulate.p.value = TRUE)
dataset %>% select(Household,ST) %>% table %>% chisq.test(simulate.p.value = TRUE)






