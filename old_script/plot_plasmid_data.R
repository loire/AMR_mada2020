require(tidyverse)






datamob = read_delim("plasmid_mob.tsv",delim="\t",col_names=F)
datamob %>% glimpse
colnames(datamob)[1]="plasmid"
colnames(datamob)[2]="sample"
colnames(datamob)[7]="replicase"
colnames(datamob)[9]="relaxase"
datamob = datamob %>% select(plasmid,sample,replicase,relaxase) %>% 
  mutate(replicase=as.factor(replicase)) %>%  
  mutate(relaxase=as.factor(relaxase)) 

dataRes = read_delim("Resistance_in_plasmid.tsv",delim="\t",col_names=F)
colnames(dataRes)[1]="plasmid"
colnames(dataRes)[2]="resgenes"

jdf = left_join(datamob,dataRes) %>%  
  mutate(SampleType = ifelse(substring(sample,1,1) %in% c("E","H"),
                                        substring(sample,1,1),
                                        substring(sample,1,2))) %>% 
  mutate(SampleType = as.factor(SampleType)) 
jdf$SampleType = fct_collapse(jdf$SampleType, PO=c("PG","PH","PP"))


  

bladf = jdf %>% separate_rows(replicase,sep=",") %>%
  separate_rows(relaxase,sep=",") %>%  
  separate_rows(resgenes,sep=",") %>% 
  mutate(resgenes=as.factor(resgenes)) %>% 
  mutate(relaxase = as.factor(relaxase)) %>% 
  mutate(replicase = as.factor(replicase)) %>% 
  mutate(Abtype = ifelse(substring(resgenes,1,6)=="blaCTX" , 
                       "CTX" ,  ifelse(substring(resgenes,1,6)=="blaSHV",
                                       "SHV", ifelse(substring(resgenes,1,6)=="blaTEM",
                                                     "TEM",
                                                     "Others")))) %>%
  filter(Abtype %in% c("CTX","TEM","SHV")) %>%  select(-Abtype)

  


bladf2 = bladf %>%  gather(key="feature",value="gene",3:5) %>%  mutate(feature="feature") %>% 
  select(-feature) %>% distinct
bladf2 %>%  dim
  
blamat = table(bladf2$plasmid,bladf2$gene,useNA="no") %>%  as.matrix
blamat %>%  dim
mjdf = jdf %>% filter(plasmid %in% rownames(blamat)) 

require(ggdendro)
require(factoextra)
require(dendextend)


dd = blamat %>% 
  dist(method="binary") %>%  
  hclust %>%  
  as.dendrogram %>% dendro_data

hc = blamat %>% 
  dist(method="binary") %>%  
  hclust 
library(gplots)  
cutree(hc,h=0) %>% as.data.frame %>% group_by(as.factor(.)) %>% count

ct = cutree(hc,h=0)
ct[ct==1]
ct[ct==2]


dd$labels$reservoir= ifelse(substring(dd$labels$label,1,1) %in% c("E","H"),
                            substring(dd$labels$label,1,1),
                            substring(dd$labels$label,1,2))
dd$labels$reservoir = ifelse(substring(dd$labels$label,1,2) %in% c("PG","PH","PP"),"PO",dd$labels$reservoir)

ggplot(segment(dd)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  theme_bw() + 
  geom_rect(data=label(dd),
            aes(xmin=x-0.5,xmax=x+0.5, ymin=-0.01,ymax=y-0.05, fill=reservoir)) +
  theme(axis.ticks.x =element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),) + ylab("distance") + 
  ggtitle("Clustering hierarchique des plasmides (n=313) circulant au sein des reservoirs") + 
  scale_fill_brewer(palette='Set3')
ggsave("clustering_plasmides.png",width=10)  




subblamat = blamat
rownames(subblamat) = ifelse(substring(rownames(subblamat),1,1) %in% c("E","H"),
                          substring(rownames(subblamat),1,1),
                          substring(rownames(subblamat),1,2))
rownames(subblamat) = ifelse(substring(rownames(subblamat),1,2) %in% c("PG","PH","PP"),"PO",rownames(subblamat))

subblamat = subblamat[which(rownames(subblamat) %in% c("H","PO","PC","CH","CN")),]

sdd = subblamat %>% 
  dist(method="binary") %>%  
  hclust %>%  
  as.dendrogram %>% dendro_data

dim(subblamat)
ggplot(segment(sdd)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  theme_bw() + 
  geom_rect(data=label(sdd),
            aes(xmin=x-0.5,xmax=x+0.5, ymin=-0.01,ymax=y-0.05, fill=label)) +
  theme(axis.ticks.x =element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),) + ylab("distance") + ggtitle("Clustering hierarchique des plasmides (n=255) circulant au sein des reservoirs principaux") + 
  scale_fill_brewer(palette='Set1')

ggsave("clustering_plasmides_abundant.png",width=10)





bladf %>% group_by(SampleType) %>% summary
 
blamat %>% as.data.frame %>% mutate(reservoir=ifelse(substring(Var1,1,1) %in% c("E","H"),
                                                     substring(Var1,1,1),
                                                     substring(Var1,1,2))) %>% 
  mutate(reservoir = ifelse(substring(Var1,1,2) %in% c("PG","PH","PP"),"PO",reservoir)) %>% 
  group_by(reservoir,Var2) %>% summarize(n=sum(Freq)/n()) %>% spread(Var2,value=n) %>% t


labels(dd) = mjdf$SampleType
dd %>% labels

require(circlize)











bladf %>% filter(resgenes %in% c("blaTEM-1B","blaCTX-M-15")) %>%
  na.omit %>% select(replicase) %>% as.vector %>%  levels


bladf %>% filter(substring(resgenes,1,6)=="blaCTX") %>% na.omit %>%
  group_by(SampleType,resgenes,replicase) %>% count %>%  ungroup %>% 
  mutate(replicase = fct_collapse(replicase,
                                  IncX=c("IncX1"),
                                  colRNAI=c("ColRNAI_rep_cluster_1987","ColRNAI_rep_cluster_1857","ColRNAI_rep_cluster_1291"),
                                  IncI=c("IncI1","IncI2"),
                                  IncQ = c("IncQ1"),
                                  IncF = c("IncFIA","IncFIIA","IncFII","IncFIB"))) %>% 
  mutate(replicase = fct_reorder(replicase,n,.fun=mean)) %>% 
  ggplot() + geom_bar(aes(x=SampleType,y=n,fill=replicase),stat="identity",position="fill") +
  facet_wrap(~ resgenes) +  scale_fill_brewer(name="Groupe d'incompatibilité",type="qual",palette ="Set3",direction =-1)  + theme_bw() +xlab("reservoir") + ylab("proportion")
ggsave("Association_replicase_CTX.png")
bladf %>% filter(substring(resgenes,1,6)=="blaTEM") %>% na.omit %>%
  group_by(SampleType,resgenes,replicase) %>% count %>%  ungroup %>% 
  mutate(replicase = fct_collapse(replicase,
                                  IncX=c("IncX1"),
                                  colRNAI=c("ColRNAI_rep_cluster_1987","ColRNAI_rep_cluster_1857","ColRNAI_rep_cluster_1291"),
                                  IncI=c("IncI1","IncI2"),
                                  IncQ = c("IncQ1"),
                                  IncF = c("IncFIA","IncFIIA","IncFII","IncFIB"))) %>% 
  mutate(replicase = fct_reorder(replicase,n,.fun=mean)) %>% 
  ggplot() + geom_bar(aes(x=SampleType,y=n,fill=replicase),stat="identity",position="fill") +
  facet_wrap(~ resgenes) + theme_bw()+xlab("reservoir") + ylab("proportion") +
  scale_fill_brewer(name="Groupe d'incompatibilité",type="qual",palette ="Set3",direction =-1) 
ggsave("Association_replicase_TEM.png")


bladf %>% filter(substring(resgenes,1,6)=="blaSHV") %>% na.omit %>%
  group_by(SampleType,resgenes,replicase) %>% count %>%  ungroup %>% 
  mutate(replicase = fct_collapse(replicase,
                                  IncX=c("IncX1"),
                                  colRNAI=c("ColRNAI_rep_cluster_1987","ColRNAI_rep_cluster_1857","ColRNAI_rep_cluster_1291"),
                                  IncI=c("IncI1","IncI2"),
                                  IncQ = c("IncQ1"),
                                  IncF = c("IncFIA","IncFIIA","IncFII","IncFIB"))) %>% 
  mutate(replicase = fct_reorder(replicase,n,.fun=mean)) %>% 
  ggplot() + geom_bar(aes(x=SampleType,y=n,fill=replicase),stat="identity",position="fill") +
  facet_wrap(~ resgenes) + theme_bw()+xlab("reservoir") + ylab("proportion") +
  scale_fill_brewer(name="Groupe d'incompatibilité",type="qual",palette ="Set3",direction =-1) 

ggsave("Association_replicase_TEM.png")




jdf %>% filter(grepl("IncF",replicase) & grepl("IncY",replicase) & grepl("blaTEM-1B",resgenes)& grepl("blaCTX-M-15",resgenes)) 
jdf %>% filter(grep("blaTEM-1B",resgenes))

grepl("blaTEM-1B",jdf$resgenes)












bladf %>% filter(resgenes %in% c("blaTEM-1B","blaCTX-M-15")) %>% na.omit %>%
  group_by(SampleType,resgenes,replicase) %>% count %>%  ungroup %>% 
  mutate(replicase = fct_collapse(replicase,
                                  IncX=c("IncX1"),
                                  colRNAI=c("ColRNAI_rep_cluster_1987","ColRNAI_rep_cluster_1857","ColRNAI_rep_cluster_1291"),
                                  IncI=c("IncI1","IncI2"),
                                  IncQ = c("IncQ1"),
                                  IncF = c("IncFIA","IncFIIA","IncFII","IncFIB"))) %>% 
  mutate(replicase = fct_reorder(replicase,n,.fun=mean)) %>% 
  ggplot() + geom_bar(aes(x=SampleType,y=n,fill=replicase),stat="identity",position="fill") +
  facet_wrap(~ resgenes) + scale_fill_brewer(name="Groupe d'incompatibilité",type="qual",palette ="Set3",direction =-1) + theme_bw() +xlab("reservoir") + ylab("proportion")
ggsave("Association_replicase_CTXM15-TEM1.png")


bladf %>% filter(resgenes %in% c("blaTEM-1B","blaCTX-M-15")) %>% na.omit %>%
  group_by(SampleType,resgenes,relaxase) %>% count %>%  ungroup %>% 
  mutate(relaxase = fct_reorder(relaxase,n,.fun=mean)) %>% 
  ggplot() + geom_bar(aes(x=SampleType,y=n,fill=relaxase),stat="identity",position="fill") +
  facet_wrap(~ resgenes) + scale_fill_brewer(name="MOB type",type="qual",palette ="Set1",direction =-1) + theme_bw() +xlab("reservoir") + ylab("proportion")
ggsave("Association_MOB_CTXM15-TEM1.png")





bladf %>% filter(resgenes %in% c("blaTEM-1B","blaCTX-M-15")) %>% na.omit %>%
    group_by(SampleType,resgenes,relaxase) %>% count %>%  ungroup %>% 
  ggplot() + geom_bar(aes(x=SampleType,y=n,fill=relaxase),stat="identity",position="fill") +
  facet_wrap(~ resgenes)

png()
for (i in levels(bladf$SampleType)) {
bladf %>%  filter(SampleType==i) %>% group_by(replicase,relaxase,resgenes) %>% 
  count %>% na.omit() %>% ungroup %>%  select(replicase,resgenes,n) %>% filter(n>1) %>%   
  chordDiagram() 
}
bladf %>%  filter(SampleType=="PO") %>% group_by(replicase,relaxase,resgenes) %>% 
  count %>% na.omit() %>% ungroup %>%  select(replicase,resgenes,n) %>% filter(n>0) %>%   
  head


bladf %>%  filter(SampleType=="BV") %>% group_by(replicase,relaxase,resgenes) %>% 
  count %>% na.omit() %>% ungroup %>%  select(replicase,resgenes,n) %>% filter(n>0) %>%   
  chordDiagram() 

levels(bladf$SampleType)


bladf 
require(devtools)
install_github("Displayr/flipPlots")
require(flipPlots)

test = bladf %>%  filter(SampleType=="PO") %>% group_by(replicase,relaxase,resgenes) %>% 
  count %>% na.omit() %>% ungroup %>%  select(resgenes,replicase,n) %>% filter(n>0)
SankeyDiagram(data = test[,-3],weights=test$n,link.color="Source")

test = bladf %>% group_by(replicase,relaxase,resgenes) %>% 
  count %>% na.omit() %>% ungroup %>%  select(resgenes,replicase,n) %>% filter(n>0)
SankeyDiagram(data = test[,-3],max.categories = 25,weights=test$n,link.color="Source")











