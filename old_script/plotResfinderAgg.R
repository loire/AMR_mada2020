
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
hg
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
ggsave("Figure5.svg",height = 7,width=10)


gene_counts %>%
  filter(grepl("bla",Resistance.gene) |
           Resistance.gene=="parC" |
           Resistance.gene=="gyrB") %>%
  pivot_wider(values_fill = 0,values_from = n,names_from=Resistance.gene) -> genedf
genedf %>% select(-Host) %>% as.matrix -> genematrix
genematrix[is.na(genematrix)] = 0
rownames(genematrix) = genedf$Host
chisq.test(genematrix,simulate.p.value = T)

genematrix

####################### ST ############################

ST = read.csv("dataset/ST_data_final.txt",sep="\t",h=F,stringsAsFactors = F) %>% select(V1,V2)
colnames(ST) = c("sample","ST")
ST
ST %>% filter(sample %in% samples$id) -> ST

data %>% select(True_indiv_name,Host) %>%
  left_join(ST,by=c("True_indiv_name"="sample")) %>% group_by(ST,Host) %>% count -> ST_counts

ST_counts %>% group_by(ST) %>% summarize(Host = "All",n=sum(n)) -> All_count
All_count
ST_counts = rbind(as.data.frame(ST_counts),as.data.frame(All_count))
ST_counts %>% group_by(Host) %>% summarize(total = sum(n)) %>% ungroup -> ST_host_counts
ST_host_counts

ST_counts
require(ggpubr)
require(ggtext)
ST_counts %>%
  mutate(type = ifelse(Host =="All","All","Host")) %>%
  #filter(Host!="Horse") %>%
  left_join(ST_host_counts,by="Host") %>% ungroup %>%
  mutate(Host = paste(Host,"  <br><span style='font-size:9pt'>(N=",total,")</span>",sep="") ) %>%
  mutate(Host = fct_reorder(factor(Host),n,sum,.desc=T)) %>%
  mutate(ST = fct_reorder(factor(ST),n,sum,.desc = F)) %>%
  mutate(ST = fct_recode(ST,Unkown="0")) %>%
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
ggsave("ST_frequency_host.pdf",height = 7,width=10)
ggsave("Figure2_A.svg",height = 7,width=10)


ST_counts %>% filter(Host!="All") %>%
  pivot_wider(values_fill=0 ,values_from = n,names_from = Host) %>%
  as.data.frame -> ST_dataframe
ST_matrix = as.matrix(ST_dataframe)
ST_matrix[,-1] -> ST_matrix
rownames(ST_matrix) = ST_dataframe$ST
chisq.test(ST_matrix,simulate.p.value = T)

ST_matrix

require(tidyverse)
###################################
Phylogroup = read.csv("Phylogroup.txt",sep=" ",h=F,stringsAsFactors = F) %>% select(V1,V2)
colnames(Phylogroup) = c("sample","Phylogroup")
Phylogroup %>% filter(sample %in% GoodSamples$V2) -> Phylogroup

data %>% select(True_indiv_name,Host) %>%
  left_join(Phylogroup,by=c("True_indiv_name"="sample")) %>%
  group_by(Phylogroup,Host) %>% count -> Phylogroup_counts
Phylogroup_counts %>% group_by(Phylogroup) %>% summarize(Host = "All",n=sum(n)) -> All_count
All_count
Phylogroup_counts = rbind(as.data.frame(Phylogroup_counts),as.data.frame(All_count))
Phylogroup_counts %>% group_by(Host) %>% summarize(total = sum(n)) %>% ungroup -> Phylogroup_host_counts
Phylogroup_host_counts

Phylogroup_host_counts

Phylogroup_counts %>% filter(Host!="All") %>%
  pivot_wider(names_from = "Host",values_from="n",values_fill= 0) -> PhylogroupMat
PhylogroupMat %>% select(-Phylogroup) %>% chisq.test(simulate.p.value = T)

Phylogroup_counts %>%
  mutate(type = ifelse(Host =="All","All","Host")) %>%
  #filter(Host!="Horse") %>%
  left_join(Phylogroup_host_counts,by="Host") %>% ungroup %>%
  mutate(Host = paste(Host,"  <br><span style='font-size:9pt'>(N=",total,")</span>",sep="") ) %>%
  mutate(Host = fct_reorder(factor(Host),n,sum,.desc=T)) %>%
  mutate(Phylogroup = fct_reorder(factor(Phylogroup),n,sum,.desc = F)) %>%
  mutate(Phylogroup = fct_lump_n(Phylogroup,25)) %>%
  mutate(freq = n/total*100) %>%
  ggplot() + geom_tile(aes(x=Host,y=Phylogroup,fill=freq)) +
  geom_text(aes(x=Host, y=Phylogroup, label=round(freq,1)),check_overlap = TRUE) +
  theme_bw() + scale_fill_viridis_c(option = "plasma",name="Percentage",direction=-1) +
  ylab("Sequence Type") + xlab("Hosts") +
  theme_minimal() + theme(panel.grid=element_blank()) +
  theme(panel.grid = element_blank(),text= element_text(family = "Helvetica",size=14)) +
  theme(axis.text.y = element_markdown(hjust = 0),axis.text.x = element_markdown()) +
  facet_grid(cols = vars(type),scales = "free",space="free") + xlab(NULL)
ggsave("Phylogroup_frequency_host.pdf",height = 7,width=10)
ggsave("Figure2_B.svg",height = 7,width=10)


####### Plasmids and genes
plas = read.csv("infos_plasmides.tsv",sep="\t")
plas %>%  glimpse

plas = data %>% select(Sample,Host) %>% left_join(plas,by = "Sample")

plas %>% separate_rows(Resgenes,sep=";") %>%
  group_by(rep_type.s.,Resgenes) %>%
  summarize(n = n()) %>% arrange(desc(n))
plas %>% xlsx::write.xlsx(file = "Plasmides_and_resistance.xlsx")

ST %>%  left_join(plas, by = c("sample"="Sample")) %>% xlsx::write.xlsx(file = "ST_Plasmides_and_resistance.xlsx")
