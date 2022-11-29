require(tidyverse)
ST =read.csv2("ST_data_updated.txt",sep=" ",h=F)
met = readxl::read_xlsx("Metadata.xlsx")
met %>%  glimpse
colnames(ST) = c("Sample","ST")
data = met %>% select(True_indiv_name,Host) %>% left_join(ST,by = c("True_indiv_name"="Sample"))

hostlevel=data %>% group_by(Host) %>% summarize(count = n()) %>% arrange(desc(count)) %>% select(Host) %>% as.vector()
hostlevel = c("All",hostlevel$Host)

heat1 = data %>%
  group_by(Host,ST) %>%
  summarize(Host,STcount = n()) %>% unique
heat1 = heat1[,c(2,1,3)]

heat2 = data %>%
  group_by(ST) %>%
  summarize(Host="All",STcount = n())
above5 = heat2 %>% filter(STcount > 5) %>% select(ST)

STlevel = heat2 %>% filter(STcount > 5) %>% filter(ST!="0") %>% arrange(STcount) %>% select(ST)
STlevel

dataf = rbind(as.data.frame(heat2),as.data.frame(heat1)) %>%
  filter(ST %in% above5$ST) %>%
  filter(ST !="0")
#dataf = dataf %>% complete(expand(dataf,ST,Host),fill=list(STcount=0))

dataf %>%
  mutate(type = ifelse(Host=="All","All","Host")) %>%
  mutate(Host = factor(Host,levels = hostlevel) ) %>%
  mutate(ST = factor(ST,levels = STlevel$ST)) %>%
  ggplot() + geom_tile(aes(x=Host,y=ST,fill=STcount)) +
  geom_text(aes(x=Host,y=ST,label=STcount)) +
  scale_fill_viridis_c(name = "count",direction = -1,option="C") +
  theme_minimal() + theme(panel.grid=element_blank()) +
  facet_grid(cols = vars(type),scales = "free",space="free") + xlab(NULL) + ylab("SequenceType")

+
  facet_grid(rows=vars(species),cols = vars(type),scales = "free",space="free") +
  xlab(NULL)
ggsave("Figure2_A.svg",width=8,height=6)



