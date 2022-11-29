library(tidyverse)
require(circlize)
require(flipPlots)

stpla = read.csv2("dataset/ST_plasm.csv")

tt = stpla %>% select(Sample,ST2,Host,Resgenes,groupeinc) %>%
  na.omit %>% mutate(ST2= as.factor(ST2)) %>%  filter(ST2!= "Unknown")
tt$Sample %>% unique %>% length


df = stpla %>% select(ST2,Host,Resgenes,groupeinc) %>%
  na.omit %>% mutate(ST2= as.factor(ST2)) %>%
  mutate(ST2 = fct_lump(ST2,n=30)) %>%
  group_by_all() %>%  count
dim(df)

stpla %>% select(ST2,Host,Resgenes,groupeinc) %>%
  na.omit %>% mutate(ST2= as.factor(ST2)) %>%
  mutate(ST2 = fct_lump(ST2,n=30)) %>% group_by(Resgenes,groupeinc) %>% count(sort=T) %>%
  write.table("Association_Resgenes_incgroup.csv",sep=",",row.names = F)


df[,c(2,1,4,3)]

sankey = SankeyDiagram(max.categories = 31,data = df[,c(2,1,4,3)],weights = df$n,link.color = "Source",label.show.percentages = T,label.show.varname=F)
sankey

df2 = df %>% filter(ST2 != "Other") %>% filter(ST2!= "Unknown") %>% filter(groupeinc!="other-rep") %>% droplevels
sankey2 = SankeyDiagram(max.categories = 31,data = df2[,c(2,4,3)],weights = df2$n,link.color = "Source",label.show.counts = F,label.show.varname=F)
sankey2

df  %>% filter(ST2!= "Unknown")



install.packages("webshot")
webshot::install_phantomjs()





plasmid_data = read.csv2("plasmides/infos_plasmides.tsv",sep="\t",header=T)
plasmid_data %>% glimpse

ST_data = read.csv2("ST_data.txt",sep="\t",header=F)
ST_data = ST_data %>% select(V1,V3)
ST_data %>%  glimpse
colnames(ST_data) = c("Sample","ST")

plasmid_ST = left_join(plasmid_data,ST_data,by="Sample")

meta = readxl::read_xlsx("Metadata.xlsx")
meta %>%  glimpse

plasmid_ST_meta = left_join(meta,plasmid_ST,by =c("True_indiv_name"= "Sample"))
plasmid_ST_meta %>%  glimpse

plasmid_ST_meta  %>% mutate(ST = ifelse(ST=="-","novel",ST)) %>%
  group_by(Host,ST) %>%
  count %>% na.omit() %>% ungroup %>% filter(ST=="131")


plasmid_ST_meta  %>% mutate(ST = ifelse(ST=="-","novel",ST)) %>%
  group_by(Host,ST) %>%  select(Host,ST) %>%
  count %>% na.omit() %>% ungroup  %>% filter(n>1) %>%
  chordDiagram()

resistance = read.csv2("Resistances/RES_finder_matrix.tsv",sep=" ")
resistance %>%  glimpse
bla = resistance %>%  select(matches("bla"))
bla$Sample = resistance$sample
bla %>%  glimpse


pdf("blaGenes_repartation.pdf")
bla %>% gather(blagene, presence,1:23) %>% filter(presence==1) %>%
  left_join(meta,by= c("Sample" = "True_indiv_name")) %>%
  select(Host,blagene) %>% group_by(Host,blagene) %>%  count %>%  ungroup %>%
  na.omit()  %>%
  chordDiagram()
dev.off()


require(flipPlots)

test  =plasmid_ST_meta %>% separate_rows(Resgenes,sep=";") %>%
  select(Host,ST,file_id,Resgenes) %>%
  mutate(ST = ifelse(ST=="-","novel",ST)) %>%  group_by_all() %>%
  filter(grepl("bla",Resgenes)) %>% count
test$Host = as.factor(test$Host)
test$ST = as.factor(test$ST)
test$Resgenes = as.factor(test$Resgenes)
pdf("Sankey_alldata.pdf")
SankeyDiagram(data = test[,-5],weights=test$n,link.color="Target",label.show.counts = T)
dev.off()


plasmid_ST_meta =  plasmid_ST_meta %>% mutate(sample = True_indiv_name)

cluster = read.csv2("Phylogenie/cluster0.01",sep=",")
clusterd = cluster %>%  filter(sequencesperCluster>6) %>%  select(sample = leafname,clustername)

test2 = left_join(clusterd,plasmid_ST_meta,by="sample") %>% separate_rows(Resgenes,sep=";") %>%
  select(Host,clustername,file_id,Resgenes) %>% group_by_all() %>%
  filter(grepl("bla",Resgenes)) %>% count %>% ungroup %>% transmute_all(.funs = as.factor) %>% mutate(n=as.numeric(n))
test2 %>%  glimpse
SankeyDiagram(max.categories = 30,data = test2[,-5],weights = test2$n,link.color = "Source",label.show.counts = F,label.show.varname=F)



htmlwidgets::saveWidget(sankey, file ="sankey1.html") # save html widget
webshot::webshot("sankey1.html",file="sankey1.pdf")
htmlwidgets::saveWidget(sankey2, file ="sankey2.html") # save html widget
webshot::webshot("sankey2.html",file="sankey2.pdf")






