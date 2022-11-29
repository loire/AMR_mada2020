require(tidyverse)

clust = read.table("dataset/cluster_phydelity_k2_sol0_Fastree_sept_rooted.txt",header=T)
clust %>%  glimpse
clust %>% select(CLUSTER) %>% unique %>% dim

# Add infos in SNP number
distanceSNP = read.csv2(file="dataset/Distances_in_cluster_phydelity_norecomb.txt",sep=" ",h=T)
colnames(distanceSNP) = c("Cluster","TAXA1","TAXA2","SNP")
distanceSNP %>%  glimpse

require(tidyverse)

distanceSNP %>%
  group_by(Cluster) %>%
  summarise(meanSNP = mean(SNP),maxSNP=max(SNP))

clust = clust %>% left_join(distanceSNP %>%
                      group_by(Cluster) %>%
                      summarise(meanSNP = mean(SNP),maxSNP=max(SNP)) ,by=c("CLUSTER"="Cluster"))



clust
distanceSNP$SNP > 20
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
              "Turkey" = "#D9D9D9",
              "Water" = "#80B1D3")




######### get list of ecoli samples #######
GoodSamples = read.table("ataset/Sample_510.txt")
samples = GoodSamples %>% select(V2)
colnames(samples)="id"
######## get list of bad samples
df %>% filter(!True_indiv_name %in% samples$id)
######## save the good ones
df = df %>% filter(True_indiv_name %in% samples$id)

df %>% glimpse

df %>% filter(CLUSTER!="NA") %>% select(Numero_foyer) %>% unique %>% dim


df %>% write.csv("dataset/CompositionClusterSup20.csv")


left_join(df %>% group_by(Host) %>% count() %>% mutate(all = n)  %>% select(-n),
  df %>% filter(CLUSTER!="NA") %>% group_by(Host) %>% count() %>% mutate(cluster = n)  %>% select(-n)) %>%
  gather(key="type",value = "count",-"Host" ) %>%
  ggplot() + geom_bar(aes(x=Host,y=count,fill = Host),stat ="identity") + facet_grid(type ~ .) +
  scale_fill_manual(values = colorHost)




# df %>%
#   filter(CLUSTER !="NA") %>%
#   group_by(CLUSTER) %>% summarise(diversity = length(unique(Host)),count = length(Host)) %>%
#   mutate(stat = diversity/count)
#
# df %>% group_by(CLUSTER,Host) %>%  count %>% filter(CLUSTER!="NA") %>%
#   ggplot() + geom_bar(aes(x=as.factor(CLUSTER),y=n,fill=Host),position="stack", stat="identity") +
#   scale_fill_manual(values = colorHost)

# df %>%
#   filter(CLUSTER !="NA") %>%
#   mutate(CLUSTER = as.factor(CLUSTER)) %>%
#   group_by(CLUSTER,Host) %>% count %>%
#   ggplot(aes(x=CLUSTER, y=n, fill=Host))+
#   geom_bar(width = 1, stat = "identity")+
#   scale_fill_manual(values = colorHost) +
#   coord_polar("y", start=0) +
#   facet_wrap(~CLUSTER,ncol = 10) +
#   theme(strip.background = element_blank(),
#         strip.text.x = element_blank())

cp <- coord_polar("y")
cp$is_free <- function() TRUE

df %>% glimpse

meanSNPs = df %>% filter(CLUSTER!="NA") %>% select(CLUSTER,meanSNP) %>% unique
meanSNPs

labmeanSNPS = as.character(signif(meanSNPs$meanSNP,digits = 3))
names(labmeanSNPS) = meanSNPs$CLUSTER
labmeanSNPS

library(extrafont)
font_import()


df %>%
  filter(CLUSTER !="NA") %>%
  mutate(CLUSTER = as.factor(CLUSTER)) %>%
  group_by(CLUSTER,Host,meanSNP) %>% count %>%
  ungroup %>% group_by(CLUSTER) %>% mutate(size = sqrt(sum(n))) %>% ungroup %>%
  mutate(CLUSTER = fct_reorder(CLUSTER,meanSNP,mean)) %>%
  ggplot(aes(group=CLUSTER,x=size/2, y=n,width=size, fill=Host))+
  geom_bar(width = 1, stat = "identity")+
  scale_fill_manual(values = colorHost) +
  facet_wrap(~CLUSTER,ncol = 13,drop=T,scales="free_y",
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
        text = element_text(family = "Times")
        )

ggsave("clusterplot.pdf",dpi=600)
ggsave("Figure4_B.svg",width=7,height = 8)

df

####### Multinomial test for cluster proba ######
##### pop distribution: All samples  ############
##### foyer distribution: Samples count in foyers of individuals in clusters ####

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
  group_by(CLUSTER) %>%
  summarize(nfoyer = length(unique(Numero_foyer)),size = length(Host),indivs = paste0(Host,collapse="",sep=",") ) %>%
  left_join(pvalcluster)

annotation %>% arrange(arrange(pvalfoy)) %>%  write_csv2("Clusters_composition_and_pvals.csv")

annotationtmp = annotation %>% filter(CLUSTER%in%c(916,798))

df %>% filter(CLUSTER != "NA")

df

contact = read.csv2("../../OddRatio/metadonnees_typologiesoins.csv",sep=",")
contact %>% select(human = Numéro.ID.individu,animal = Numéro.ID.animaux.espèce.1) -> contact_simple
contact_simple %>% separate_rows(animal,sep=";") %>%
  group_by(human) %>% summarize(animals_contact = paste(animal,collapse=";")) %>% filter(animals_contact!=";;") -> contact_all


df  %>%  mutate(sample = str_split_fixed(True_indiv_name,pattern ="_",n=2)[,1]) -> df
df %>% mutate(sample = gsub("-redeposit","",sample)) -> df


testanimal = "PG165;PG166;PG167;PG168;PG169"
testcluster = "CH89;E25;PG166;PG167"

contact_in_cluster = function(cluster_str,animals_str){
  cluster_sample = str_split(cluster_str,";",simplify = F)[[1]]
  animal_contact = str_split(animals_str,";",simplify = F)[[1]]
  return(length(intersect(cluster_sample, animal_contact)))
}

df

left_join(df,contact_all,by=c("sample"="human")) %>%
  filter(CLUSTER!="NA") %>%
  group_by(CLUSTER) %>%
  summarize(maxSNP = max(maxSNP),foyer_number=length(unique(Numero_foyer)),
            samples = paste(sample,collapse=";"),
            sample,animals_contact) %>% filter(grepl("^H",sample)) %>%
  mutate(human_in_cluster=sample) %>% select(-sample) %>%
  mutate(contact_number_in_cluster= contact_in_cluster(samples,animals_contact)) -> dfout

write_excel_csv(dfout,file = "Cluster_with_contact_info.csv")

dfout %>% filter(contact_number_in_cluster>0)





cluster %>% filter(TAXA == "H117")





p = df %>%
  left_join(annotation) %>%
  filter(CLUSTER !="NA")  %>%
  mutate(CLUSTER = as.factor(CLUSTER)) %>%
  group_by(CLUSTER,Host,pvalfoy) %>% count %>%
  ungroup %>% group_by(CLUSTER) %>% mutate(size = sqrt(sum(n))) %>% ungroup %>%
  mutate(CLUSTER = fct_reorder(CLUSTER,pvalfoy,mean)) %>%
  ggplot() +
  geom_bar(aes(group=CLUSTER,x=size/2, y=n, fill=Host,width=size),width=1, stat = "identity")+
  scale_fill_manual(values = colorHost) +
  facet_wrap(~CLUSTER,ncol = 13,drop=T,scales="free_y") +
  cp +
  theme_void() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1,legend.position = "top")
p

annotation %>% arrange(pvalfoy)
annotation
p  +
  geom_text(annotation,mapping =aes(label=paste("#Foyers = ",nfoyer),
                                       x=-Inf,y=-Inf,hjust = -1,vjust = -1)) +
  geom_text(annotation,mapping =aes(label=paste("p = ",pvalfoy),
                                       x=-Inf,y=-Inf,hjust = -1.5,vjust = -4))




distanceSNP = read.csv2(file="Distances_in_cluster_phydelity.txt",sep=" ")



#### OLD STUFF

# library(devtools)
# devtools::install_github("briatte/ggnet")
# require(network)
# require(sna)
# require(ggnet)
# require(scales)
# edges = distance %>%   %>% select("from" = sample1,"to"= sample2,"weigth" = SNP)
# nodes = df  %>% select("id"= True_indiv_name,Host,Numero_foyer)
# net = network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE,names.eval = "weights")
# read.csv2("Res",sep=" ")




