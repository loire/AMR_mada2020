require(tidyverse)
require(ggtree)
require(randomcoloR)
require(rlist)
require(tidytree)
require(ape)
require(geosphere)

####### DEFINE COLOR PALETTE FOR ALL GRAPHS WITH HOSTS

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
core = read.tree(file = "Fastree_sept_rooted.tree")
## leave some space for inlay
ggtree(core,layout="fan",open.angle = 90)


## Metadata parsing
metadata = read.csv("Metadata.csv")
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
  geom_point(data = distTipsHost %>% filter(type == "Horse"),
             aes(x=type,y=Freq,fill=type),pch=21,color="black",show.legend=F,size=3) +
  scale_fill_manual(values = colorHost_inter) +
  geom_text(data = distTipsHost %>% group_by(type) %>% summarize(n=n()),aes(x=type,y=-0.02,label=n)) +
  theme_minimal() + xlab("Host pair") + ylab("Genetic distances")

distancePlot
ggsave("distance_host.pdf",height = 3,width = 7)

globalDistancePlotHost  = distTipsHost %>% ggplot() + geom_density(aes(x=Freq,color=InterIntra)) +
  theme_minimal()  + xlab("Genetic distances")  + scale_color_discrete(name="")
globalDistancePlotHost

############# Now with Geography
host %>% glimpse
geo %>% glimpse
distTips_foyer = distTips %>% left_join(geo,by=c("Var1"="sample")) %>%
  mutate(Foyer1 = Numero_foyer) %>% select(-Numero_foyer) %>%
  left_join(geo,by=c("Var2"="sample")) %>%
  mutate(Foyer2 = Numero_foyer) %>% select(-Numero_foyer) %>%
  mutate(type=ifelse(Foyer1==Foyer2,Foyer1,"Between Foyer")) %>%
  mutate(InterIntra = ifelse(type=="Between Foyer","Between Foyer","Within Foyer")) %>%
  mutate(IntraFoyer = ifelse(type=="Between Foyer","Between Household",Foyer1))
pal
distTips_foyer %>%
  ggplot() +
  geom_boxplot(aes(x=IntraFoyer,y=Freq,fill=IntraFoyer),show.legend =F) +
  scale_fill_manual(values = pal) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + xlab("") + ylab("Genetic distance")


metadata %>% filter(Numero_foyer==31)

distTips_foyer %>% ggplot() +
  geom_violin(aes(x=InterIntra,y=Freq,fill=InterIntra),
              draw_quantiles = c(0.5),trim=T,color="black",scale="area") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90))



globalDistancePlotFoyer  = distTips_foyer %>% ggplot() +
  geom_density(aes(x=Freq,color=InterIntra)) +
  theme_minimal()  + xlab("Genetic distances")  +
  scale_color_brewer(name="",type="qual", palette=1)

globalDistancePlotFoyer
###### Now with Fokontany

distTipsGeoFokontany = distTips %>% left_join(geo,by=c("Var1"="sample")) %>%
  mutate(Fokontany1 = Fokontany) %>% select(-Fokontany) %>%
  left_join(geo,by=c("Var2"="sample")) %>%
  mutate(Fokontany2 = Fokontany) %>% select(-Fokontany) %>%
  mutate(type=ifelse(Fokontany1==Fokontany2,Fokontany1,"Between Fokontany")) %>%
  select(-Fokontany1,-Fokontany2) %>% mutate(InterIntra = ifelse(type=="Between Fokontany","Between Fokontany","Within Fokontany"))

globalDistancePlotGeoFokontany  = distTipsGeoFokontany %>% ggplot() + geom_density(aes(x=Freq,color=InterIntra)) +
  theme_minimal()  + xlab("Genetic distances")  + scale_color_brewer(name="",type="qual",
                                                                     palette=1)
globalDistancePlotGeoFokontany


#### Now with geographic distances
metadata = read.csv2("Metadata.csv",sep=",")
meta = metadata %>%
  select(sample = True_indiv_name,Host,Numero_foyer,Fokontany,Point_GPS_S,Point_GPS_E)
meta %>% glimpse
meta$Point_GPS_E = as.numeric(as.character(meta$Point_GPS_E))
meta$Point_GPS_S = as.numeric(as.character(meta$Point_GPS_S))

data1 = left_join(distTips,meta,by=c("Var1"="sample"))
data2 = left_join(data1,meta,by=c("Var2"="sample"))
data2 %>% glimpse

getDistMeter = function(a,b,c,d){
  return(as.numeric(distm(x=c(a,b),y=c(c,d))))
}
getDistByhand = function(a,b,c,d){
  return(sqrt((b-a)^2+(d-c)^2))
}

getDistMeter(10,23,14,14)
getDistByhand(10,23,14,14)

data2 = data2 %>% rowwise %>%
  mutate(meterDist = getDistMeter(-Point_GPS_S.x,Point_GPS_E.x,-Point_GPS_S.y,Point_GPS_E.y))
data2 %>% glimpse

ggplot(data2) + geom_point(aes(x=meterDist,y=Freq),size=1) + xlab("Distance (meter)") + ylab("Distance (genetic)")
ggsave("Distance_Geo_vs_genetic.pdf")
test = lm(data2$Freq ~ data2$meterDist)
summary(test)

require(cowplot)
plot_grid(globalDistancePlotFoyer,globalDistancePlotGeoFokontany,nrow=2)
ggsave("distance_geo.pdf",height = 6,width = 10)




distancePlot = distTipsHost %>%  mutate(type = fct_infreq(type)) %>%
  ggplot() +
  geom_violin(aes(x=type,y=Freq,fill=type,),
              ,color="black",show.legend= F) +
  geom_point(data = distTipsHost %>% filter(type == "Horse"),
             aes(x=type,y=Freq,fill=type),pch=21,color="black",show.legend=F,size=3) +
  scale_fill_manual(values = colorHost_inter) +
  geom_text(data = distTipsHost %>% group_by(type) %>% summarize(n=n()),aes(x=type,y=-0.02,label=n)) +
  theme_minimal() + xlab("Host pair") + ylab("Genetic distances")

distancePlot
ggsave("distance_host.pdf",height = 3,width = 7)

globalDistancePlotHost  = distTipsHost %>% ggplot() + geom_density(aes(x=Freq,color=InterIntra)) +
  theme_minimal()  + xlab("Genetic distances")  + scale_color_discrete(name="")
globalDistancePlotHost


require(cowplot)
plot_grid(p1,p2,nrow=2)
ggsave("distance_host.pdf",height = 6,width = 10)
##########################################################################
#
# distances = read.csv("paired_distance.txt",sep="")
# distances %>% glimpse

ggsave("tree.pdf",dpi = 300,height=15,width=20 )


### Reminder of how  to get the order of leafs
tmp = fortify(core)
d = subset(tmp,isTip)
order = with(d, label[order(y, decreasing=T)])


# ####### Adding old cluster to the tree
#
# d2 =  read.csv2("cluster0.01",sep=",")
# tmp = d2 %>%  filter(sequencesperCluster>=2) %>%  filter(clustername!=0)
# tmp$clustername %>% as.factor %>% levels
# tmp %>% glimpse
#
# ###### Now new version
# curated_clust = read.csv2("clusters_curated.txt",sep = "\t")
# curated_clust %>%  glimpse
# seqpercur = curated_clust %>% group_by(clustername) %>% count
# tmp2  = left_join(curated_clust,seqpercur,by="clustername")
# tmp2 %>%  glimpse

##### And phydelity data:
phyclust = read.csv2("../ClusterDescription/cluster_phydelity_k2_sol0_Fastree_sept_rooted.txt",sep="\t")
phyclust %>%  glimpse
colnames(phyclust) = c("clustername","leafname")
seqperphyclust = phyclust %>%   group_by(clustername) %>% count
phyclust = left_join(phyclust,seqperphyclust,by="clustername")

phyclust






require(cowplot)


p3 = phyclust %>% group_by(clustername) %>%  count %>%  ggplot() + geom_histogram(aes(n),binwidth = 1) + xlim(c(0,12))
p1 = ggplot(data = seqperphyclust[-1,]) + geom_histogram(aes(n),binwidth = 1) + xlim(c(0,12))
p2 =  ggplot(data = seqperphyclust) + geom_histogram(aes(n),binwidth = 1) + xlim(c(0,12))

plot_grid(p3,p1,p2,ncol=1)

nodesclu = c()
for (i in phyclust$clustername %>% as.factor %>% levels){
  tt = phyclust %>% filter(clustername == i) %>% select(leafname) %>% droplevels
  if(length(tt$leafname)>0){
    nodesclu=c(nodesclu,do.call(MRCA,list.prepend(as.list(as.character(tt$leafname)),core)))
  }
}
nodesclu


#
# p + geom_cladelabel(node=541, label=NA,extend = 1,
#                 color="red2", align=TRUE) +
#   geom_cladelabel(node=812, label=NA,extend = 1,
#                   color="green", align=TRUE)



p = ggtree(core,layout="circular",branch.length = 0.01)
tmp_tree = p + geom_cladelabel(node=541, label=NA,offset = 0.1,
                    color=col, align=TRUE)
p
for(i in nodesclu){
  col = randomColor(1)
  p = p + geom_cladelabel(node=i, label=NA,offset = 1,
                          color=col, align=TRUE)
}
p

# COOoooooool

# Now get host and localisation data
meta = read.csv("Metadata.csv",header =T)
meta %>% glimpse
metad = meta %>% select(sample= True_indiv_name,Host)
rowsNames = metad$sample

metad = metad %>% select(-sample)
rownames(metad)= rowsNames

metad = metad %>%
  mutate(Host = as.factor(Host)) %>%
  mutate(Host = fct_recode(Host,Chicken = "Poultry"))


test_metad = metad
test_metad$sample = rownames(metad)
test_metad = left_join(test_metad,phyclust,by=c("sample" = "leafname")) %>%
  select(-sample,-n)
rownames(test_metad) = rowsNames

clusters_name = test_metad %>% filter(clustername!="NA") %>%
  select(clustername) %>% unique
length(clusters_name$clustername)

clusters_name

colsclust = randomColor(119,hue="red")
names(colsclust) = clusters_name$clustername
colsclust["NA"]="white"
# save file
#write.csv(metad,file = "metasimple.csv",row.names = F)

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

#### Plot with fokontany
p = ggtree(core,layout="circular",size = 0.2)

meta = read.csv("Metadata.csv",header =T)
meta %>% glimpse
metad = meta %>% select(sample= True_indiv_name,Fokontany)
rowsNames = metad$sample

metad = metad %>% select(-sample)
rownames(metad)= rowsNames

metad = metad %>%
  mutate(Fokontany = as.factor(Fokontany))
metad %>% head

metad$Fokontany %>% unique
pal = randomcoloR::distinctColorPalette(k=9)
gheatmap(p,metad, width=0.2, font.size=8, colnames = F,
         hjust = 0) + scale_fill_manual(name="Fokontany",values = pal) +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm")) +
  theme_tree() + theme(legend.position = "top")



#### Plot with household
p = ggtree(core,layout="circular",size = 0.2)

meta = read.csv("Metadata.csv",header =T)
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
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm")) +
  theme_tree() + theme(legend.position = "right")
p1

p2 = distTips_foyer %>%
  ggplot() +
  geom_boxplot(aes(x=IntraFoyer,y=Freq,fill=IntraFoyer),show.legend =F) +
  scale_fill_manual(values = pal) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + xlab("") + ylab("Paiwise genetic distances")
p1
ggsave("Household_tree.pdf",scale = 4)
p2

p3  = distTips_foyer %>%
  mutate(InterIntra=fct_recode(InterIntra,"Between households"="Between Foyer",
                               "Within households" = "Within Foyer")) %>%
  ggplot() +
  geom_density(aes(x=Freq,fill=InterIntra),color="black",alpha = 0.5,size = 0.5) +
  theme_minimal()  + xlab("Pairwise genetic distances")  +
  scale_fill_brewer(name="",type="qual", palette=1,)

p3
p4  = distTips_foyer %>%
  mutate(InterIntra=fct_recode(InterIntra,"Between households"="Between Foyer",
                               "Within households" = "Within Foyer")) %>%
  ggplot() +
  geom_violin(aes(y=Freq,x=InterIntra,fill=InterIntra)) +
  theme_minimal()  + ylab("Pairwise genetic distances")  +
  scale_fill_brewer(name="",type="qual", palette=1) + xlab("") + theme(axis.text.x=element_blank())
p4

test = cowplot::plot_grid(p2,NULL,p3,p4,nrow = 1,rel_widths = c(1.7,0.3,1.2,1.2),
                          labels=c("B","","C","D"))
cowplot::plot_grid(p1,test,nrow=2,rel_heights = c(1.5,1),labels=c("A",""))
ggsave("Household_tree.pdf",scale = 3)












######## TEST
p = ggtree(core,layout="rectangular",size = 0.2) +
  geom_tiplab(align=T,linetype="dotted",linesize=0.1) + geom_treescale()
p
# for(i in nodesclu){
#   col = randomColor(1,hue="red")
#   tmp_tree = tmp_tree + geom_hilight(node=i, fill=col, extend = 0.08)
# }
#tmp_tree = p + geom_cladelabel(node=541, label=NA,offset = 0.1,barsize = 4,
#                              color=col, align=TRUE)




gheatmap(p,test_metad, width=0.2, font.size=8, colnames = F,
      hjust = 0) +
  scale_fill_manual(values=colorHost) +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm")) +
  theme_tree() + theme(legend.position = "top")


ggsave("test.pdf",width=20,height = 24,units="cm")


p = ggtree(core,layout="circular",size = 0.3) + geom_tiplab()
p
ggsave("test.pdf",height=49,width=49)
BV85_S113
H39_S46
ancestor = do.call(MRCA,list.prepend(as.list(as.character(c("CN14_S198","H283_S147"))),core))
require(treeio)
ggtree(sub,layout="rectangular") + geom_treescale() + geom_tiplab(align =T )

sub = tree_subset(core,ancestor,levels_back = 0)

subf= ggtree(sub,layout ="rectangular") + geom_tiplab(align = T)
gheatmap(subf,test_metad, width=0.2, font.size=8, colnames = F,
         hjust = 0) +
  scale_fill_manual(values=colorHost) +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm")) +
  theme_tree() + theme(legend.position = "top")


subtip = subf$data %>% filter(isTip==T) %>% select(label) %>% as.vector
phyclust = phyclust %>% filter(as.character(leafname) %in% subtip$label)


nodesclu = c()
for (i in phyclust$clustername %>% as.factor %>% levels){
  tt = phyclust %>% filter(clustername == i) %>% select(leafname) %>% droplevels
  if(length(tt$leafname)>0){
    nodesclu=c(nodesclu,do.call(MRCA,list.prepend(as.list(as.character(tt$leafname)),sub)))
  }
}
nodesclu

for(i in nodesclu){
   col = randomColor(1,hue="red")
   subf = subf + geom_hilight(node=i, fill=col, extend = 0.08)
}
gheatmap(subf,test_metad, width=0.2, font.size=8, colnames = F,
         hjust = 0) +
  scale_fill_manual(values=colorHost) +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm")) +
  theme_tree() + theme(legend.position = NULL)
ggsave("treezoom.pdf",height=10,width=13,unit="cm")

############

gheatmap(subf,metad, width=1, font.size=8,
         colnames_position= "top" , colnames_angle = 90,
         colnames_offset_y = 0, hjust = 0) +
  scale_fill_manual(values=colorHost) +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm"))


# Get resistance genes data
Res = read.csv("RES_finder_matrix.tsv",sep=" ",header=T,row.names = 1,stringsAsFactors = F)
blaRes = Res %>% dplyr::select(matches("bla"))

blaRes %>% rowSums %>% max


# save file
write.csv(blaRes,file = "blaRES.csv")

# Keep genes present at least 20 times in the dataset
rowSums(blaRes[,colSums(blaRes)<20])
HblaRes  = blaRes[,colSums(blaRes) >= 20]
HblaRes$other = ifelse(rowSums(blaRes[,colSums(blaRes)<20])==0,0,1)
# Change factor levels
HblaRes[HblaRes==0] = "absent"
HblaRes[HblaRes==1] = "present"
rownames(HblaRes)
HblaRes %>% glimpse


# Now join res and meta dataset
HblaRes$sample = rownames(HblaRes)
HblaRes %>%  glimpse
metad %>%  group_by(Host) %>%  count

# refactor Host as human or Other
metad = metad %>% mutate(Host = ifelse(Host=="Human","Human","Other"))

# actual join
tmp2 = left_join(metad,HblaRes,by = "sample")
tmp2 %>%  glimpse
rownames(tmp2) = tmp2$sample
rownames(tmp2)
alldat = tmp2
alldat = tmp2 %>%  select(-sample) %>%  transmute_all(.funs = as.character)
alldat %>%  glimpse

item = c("Human","Other","present","absent")
item

rownames(alldat) = rownames(tmp2)
rownames(alldat)


require(RColorBrewer)



# dEal with colors
colHost = RColorBrewer::brewer.pal(11,"Set3")
colFok = RColorBrewer::brewer.pal(9,"Pastel1")
item
itemcol = c("red" ,"#B3DE69","black","white")
names(itemcol) = item


gheatmap(p,alldat, width=1, font.size=8,
         colnames_position= "top" , colnames_angle = 90,
         colnames_offset_y = 0, hjust = 0) +
  scale_fill_manual(values=itemcol,breaks =item) +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm"))

ggsave("tree_curated_clust.pdf",width=20,height = 45)



gheatmap(p2,alldat, width=0.7, font.size=6,
         colnames_position= "top" ,offset = 0, colnames_angle = 70,
         colnames_offset_y = 0, hjust = 0) +
  scale_fill_manual(values=itemcol,breaks =item) +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm")) +scale_y_continuous(expand = c(0.05,0))

ggsave("tree_circular_curated_clust.pdf",width=20,height = 45)

#### Circular map with only Hosts
hostd = meta %>% select(sample = True_indiv_name,Host)
rownames(hostd) = hostd$sample
hostd = hostd %>%  mutate(Host = as.character(Host)) %>%  select(-sample)
rownames(hostd) = meta$True_indiv_name
item = hostd$Host %>%  as.factor %>%  levels
itemcol = randomColor(length(item))
names(itemcol) = item

gheatmap(p2,hostd, width=0.1, font.size=6,
         colnames_position= "top" ,offset = 0, colnames = F,colnames_angle = 70,
         colnames_offset_y = 0, hjust = 0) +
  scale_fill_manual(values=itemcol,breaks =item) +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm"))
ggsave("tree_host_circular_curated.pdf",width=10,height = 10)



### Circular map with only Fokontany

foyd = meta %>% select(sample = True_indiv_name,localisation = Fokontany)
rownames(foyd) = hostd$sample
foyd = foyd %>%  mutate(localisation = as.character(localisation)) %>%  select(-sample)
rownames(foyd) = meta$True_indiv_name
item = foyd$localisation %>%  as.factor %>%  levels
itemcol = randomColor(length(item))
names(itemcol) = item

gheatmap(p2,foyd, width=0.1, font.size=6,
         colnames_position= "top" ,offset = 0, colnames = F,colnames_angle = 70,
         colnames_offset_y = 0, hjust = 0) +
  scale_fill_brewer(type="qual",palette="Set3") +
  theme(plot.margin=unit(c(3,1,1.5,1.2),"cm"))
ggsave("tree_loc_circular_curated.pdf",width=10,height = 10)


## Cicular map with host and Fokontany
itemcol


