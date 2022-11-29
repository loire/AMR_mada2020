require(readxl)
Phylo = read.csv("dataset/dataset_sample_2.csv",header =T)
Phylo = Phylo %>% select(-(1:5))
Phylo %>% glimpse

mapping = readxl::read_xlsx("dataset/mapping_stats.xlsx")

mapping %>% select(-1) -> Mapping
colnames(Mapping)[1]="Sample"
Mapping %>% left_join(Phylo,by = c("Sample"="sample")) -> Info_sample
write_excel_csv2(Info_sample,file = "Table_info_samples.csv")


