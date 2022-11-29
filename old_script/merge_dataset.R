require(tidyverse)
samples_list = read.csv(file = "dataset/Sample_510.txt",sep=" ",h=F)
samples_list = samples_list$V2
samples_list %>% length
phylogroups = read.csv(file="dataset/Phylogroups.txt",sep="\t",h=F)
ST = read.csv(file = "dataset/ST_data_updated_10_2022.txt",sep="\t",h=F)

colnames(ST) = c("sample","ST")

phylogroups %>% select(V1,V5) %>% mutate(sample = sub(".assembly.fasta","",V1)) %>%
  select(sample,phylogroup = V5) %>%  filter(sample %in% samples_list) -> dataset
dataset = unique(dataset)

dataset = left_join(dataset,ST,by="sample")
dataset = dataset %>% mutate(ST=ifelse(ST==0,"unknown",ST))

res = read.csv("dataset/RES_finder_matrix.tsv",sep=" ")

dataset = left_join(dataset,res,by="sample")

host = read.csv("dataset/host.csv")
host %>% mutate(Host=ifelse(Host=="Poultry","Chicken",as.character(Host))) -> host

dataset = left_join(host,dataset,by="sample")
dataset = dataset %>% filter(sample %in% samples_list)

dataset$Host %>% unique

meta = readxl::read_excel("dataset/Metadata.xlsx")
meta %>% select(sample = True_indiv_name,Household = Numero_foyer,Fokontany) -> othermeta

dataset = left_join(othermeta,dataset,by="sample") %>% filter(sample %in% samples_list)

dataset$Host %>% unique


trad = meta %>% select(Reads_names,True_indiv_name)
gyrb = read.csv("dataset/gyrA_QRDR_AA.csv",sep="\t")
parc = read.csv("dataset/parC_QRDR_AA.csv",sep="\t")

gyrb = gyrb %>% rowwise %>% mutate(Reads_names = str_split_fixed(isolate,"\\.",2)[1])
gyrb = gyrb %>% mutate(gyrb = ifelse(phenotype=="WT",0,1))
gyrb = gyrb %>% left_join(trad,by="Reads_names")

parc = parc %>% rowwise %>% mutate(Reads_names = str_split_fixed(isolate,"\\.",2)[1])
parc = parc %>% mutate(parc = ifelse(phenotype=="WT",0,1))
parc = parc %>% left_join(trad,by="Reads_names")

pointmutdata = left_join(parc,gyrb,by="Reads_names") %>% select(sample = True_indiv_name.x,parc,gyrb)

left_join(dataset,pointmutdata,by="sample") -> dataset

write.csv(dataset,file = "dataset_sample.csv")

