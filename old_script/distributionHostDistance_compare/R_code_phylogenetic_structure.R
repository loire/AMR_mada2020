#############################################################################################
## dirty piece of R code to investigate the amount of structure within a phylogenetic tree ##
#############################################################################################

## required packages

library(ape)
library(adephylo)
library(phangorn)
library(lattice)
library(network)
library(vioplot)
library(dplyr)
# install.packages("remotes")
library(remotes)
# remotes::install_github("prmac/slatkin.maddison")
library(slatkin.maddison)
library(ggplot2)
## some helpful functions

# trim a string and extract the xth element after a sep

trim_names <- function(X,sep,pos,fix){
  strsplit(X,split=sep,fixed=fix)[[1]][pos]
}

# list all descendant tips of each internal node in a tree

get_descent <- function(arbre,nb=FALSE){
  out <- list()
  all <- mrca(arbre)
  N<- arbre$Nnode
  tip <- length(arbre$tip.label)
  for (i in 1:N)
  {
    if(!nb) out[[i]] <- rownames(all)[apply(all==(tip+i),1,sum)>0]
    if(nb) out[[i]] <- match(rownames(all)[apply(all==(tip+i),1,sum)>0],arbre$tip.label)
  }
  out
}

##########################################################################################

## load trees

# toy tree with strong genetic structure

tree_toy_s <- read.tree("Xf_new_names")
factor_s=as.character(sapply(tree_toy_s$tip.label,trim_names,sep="_",pos=3,fix=TRUE))
plot.phylo(tree_toy_s,show.tip.label = T,tip.color = as.color(factor_s),cex=0.7)

# toy tree with random genetic structure

tree_toy_ns <- tree_toy_s
new_names=sample(tree_toy_ns$tip.label, length(tree_toy_ns$tip.label))
tree_toy_ns$tip.label=new_names
factor_ns=as.character(sapply(tree_toy_ns$tip.label,trim_names,sep="_",pos=3,fix=TRUE))
plot.phylo(tree_toy_ns,show.tip.label = T,tip.color = as.color(factor_ns),cex=0.7)

# amr tree

tree_amr <- read.tree("Fastree_sept_rooted_new_names.tree")
factor_amr=as.character(sapply(tree_amr$tip.label,trim_names,sep="_",pos=3,fix=TRUE))
plot.phylo(tree_amr, type="fan", show.tip.label = T,tip.color = as.color(factor_amr),cex=0.7)



##########################################################################################

# Test structure using the Association index statistic, as described in oi:10.1016/j.meegid.2007.08.001 

#tree=tree_toy_s 
#tree=tree_toy_ns 
tree=tree_amr

pos_factor=3
# get for each internal node, the list of descending tips
data=get_descent(tree)
# compute AI on real tree
calculation_AI=NULL
for (i in 1:length(data)){
summary=as.character(sapply(data[[i]],trim_names,sep="_",pos=pos_factor,fix=TRUE))
freq=max(table(summary))/sum(table(summary))
nbre_tip_minus_1=length(data[[i]])-1
AI_node=(1-freq)/(2^nbre_tip_minus_1)
calculation_AI=append(calculation_AI,AI_node)
}
AI_data=sum(calculation_AI)

# make trees with randomized tips (random structure)

reps=100
rtree_list <- list()
for(i in 1:reps) {
	temp <- tree
	randtips <- sample(temp$tip.label, length(temp$tip.label))
	temp$tip.label <- randtips
	rtree_list[[i]] <- temp
}
# calculate AI on each randomized trees
list_AI_randomized=NULL
for(i in 1:length(rtree_list)) {
print(i)
data=get_descent(rtree_list[[i]])
calculation_AI_randomized=NULL
for (j in 1:length(data)){
	summary=as.character(sapply(data[[j]],trim_names,sep="_",pos=pos_factor,fix=TRUE))
	freq=max(table(summary))/sum(table(summary))
	nbre_tip_minus_1=length(data[[j]])-1
	AI_node_randomized=(1-freq)/(2^nbre_tip_minus_1)
	calculation_AI_randomized=append(calculation_AI_randomized,AI_node_randomized)
}
list_AI_randomized=append(list_AI_randomized,sum(calculation_AI_randomized))
}
# compare AI randomized and real values
hist(list_AI_randomized,xlim=c(min(list_AI_randomized)-10,max(list_AI_randomized)+1))
abline(v = AI_data, col="red", lwd=3, lty=2)
# compute pvalue
num_smaller_than_obs <- length(which(list_AI_randomized <= AI_data))
pval <- num_smaller_than_obs / reps
pval

# if pval<0.05 then closely related strains are more likely 
# to share the same trait values (here geographical origin) 
# than one would expect by chance alone.

##########################################################################################
# Test structure using the Slatkin_madison index 

# tree=tree_toy_s 
# tree=tree_toy_ns 
tree=tree_amr

plot(tree,cex=0.5)
trait <- as.factor(as.character(sapply(tree$tip.label,trim_names,sep="_",pos=3,fix=TRUE)))
sm_test(trait,tree,rep=999)

# if pval<0.05 then closely related strains are more likely 
# to share the same trait values (here geographical origin) 
# than one would expect by chance alone.

##########################################################################################

# compute and analyze pairwise distance between tips

# tree=tree_toy_s 
# tree=tree_toy_ns 
tree=tree_amr

# compute all pairwise distances
dist_tips=distTips(tree, tips = "all", method = "patristic", useC = TRUE)
dist1=as.matrix(dist_tips)
dist1[upper.tri(dist1,diag = TRUE)] <- NA
# convert the matrix into a dataframe
lignes=rep(row.names(dist1),dim(dist1)[1])
values=c(dist1)
taga=NULL
colonnes=NULL
for (i in 1:dim(dist1)[2]){
	taga=rep(colnames(dist1)[i],dim(dist1)[1])
	colonnes=append(colonnes,taga)
}
matrice_linearise=cbind(lignes,colonnes,values)
matrice_lineariseb=as.data.frame(matrice_linearise,stringsAsFactors=F)
matrice_lineariseb$values=as.numeric(matrice_lineariseb$values)
all_values=na.omit(matrice_lineariseb)


all_values$factor_lignes=as.character(sapply(all_values$lignes,trim_names,sep="_",pos=3,fix=TRUE))
all_values$factor_colonnes=as.character(sapply(all_values$colonnes,trim_names,sep="_",pos=3,fix=TRUE))
all_values$type <- ifelse(all_values$factor_lignes == all_values$factor_colonnes, "within", "between")
all_values$type_bis <- ifelse(all_values$type=="within", all_values$factor_colonnes, paste(all_values$factor_colonnes,all_values$factor_lignes,sep="-"))


all_values %>% filter(grepl("PG360",lignes))  %>% filter(grepl("CT13",colonnes))
distTips %>% filter(grepl("CT13",Var1))

distTips$Var1 %>% unique %>% length
distTips$Var2 %>% unique %>% as.character %>% sort == distTips$Var1 %>% unique %>% as.character %>% sort

all_values %>%  dim
distTips %>%  glimpse
PG360_F67_Poultry CT13_F12_Cat

# plot within vs between hosts

# densityplot(~ values, group = type, data = all_values, auto.key = T, plot.points = FALSE)
# boxplot(values ~ type, data=all_values)
# vioplot(values ~ type, data=all_values)

plot_density <- ggplot(data = all_values, aes(x = values, color=type)) +	geom_density()
plot_density
plot_boxplot <- ggplot(all_values, aes(x = type, y = values, fill = type)) +
								geom_boxplot() +
								scale_x_discrete(name = "Comparisons") +
								scale_y_continuous(name = "Pairwise phylogenetic distance")
plot_violon <-  ggplot(all_values, aes(x = type, y = values, fill = type)) +
								geom_violin() +
								scale_x_discrete(name = "Comparisons") +
								scale_y_continuous(name = "Pairwise phylogenetic distance")

wilcox.test(all_values[all_values$type == "within", "values"],all_values[all_values$type == "between", "values"])

# if pval significant, then difference between distributions

# plot all distances vs between only distances vs within each factor

all=select(all_values, values)
all$type="all"
within=select(all_values[all_values$type == "within",],values,type_bis)
names(within)[2]="type"
between=select(all_values[all_values$type == "between",],values,type)

plot_values=rbind(all,within,between)

plot_density <- ggplot(data = plot_values, aes(x = values, color=type)) +	geom_density()

plot_boxplot <- ggplot(plot_values, aes(x = type, y = values, fill = type)) +
	geom_boxplot() +
	scale_x_discrete(name = "Comparisons") +
	scale_y_continuous(name = "Pairwise phylogenetic distance")

plot_boxplot <- ggplot(within, aes(x = type, y = values, fill = type)) +
	geom_boxplot() +
	scale_x_discrete(name = "Comparisons") +
	scale_y_continuous(name = "Pairwise phylogenetic distance")

plot_violon <-  ggplot(plot_values, aes(x = type, y = values, fill = type)) +
	geom_violin() +
	scale_x_discrete(name = "Comparisons") +
	scale_y_continuous(name = "Pairwise phylogenetic distance")

plot_violon


tree_amr$tip.label %>% View

