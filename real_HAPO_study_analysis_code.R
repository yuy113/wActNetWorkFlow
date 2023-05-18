

#full names of metabolites in 

dat.full.names.subnetwork.e<-read.csv(file = 'name_map_subnetwork_e.csv')

head(dat.full.names.subnetwork.e)


dat.full.names.subnetwork.e[,9]

full.names.subnetwork.e<-dat.full.names.subnetwork.e[,1]

met.fasting.fullname.final

met.fasting.shortname.final<-names(met.fasting.fullname.final)

names(met.fasting.shortname.final)<-met.fasting.fullname.final

short.names.subnetwork.e<-met.fasting.shortname.final[full.names.subnetwork.e]


pval.subnetwork.e.fasting<-pval.met.fasting32[short.names.subnetwork.e]

coeff.subnetwork.e.fasting<-coeff.met.fasting32[short.names.subnetwork.e]



dat.fasting.subnetwork.e<-data.frame(p_value=pval.subnetwork.e.fasting,estimate_coeff=coeff.subnetwork.e.fasting,metabolites_name=full.names.subnetwork.e,
                                     short_name=short.names.subnetwork.e,est_upper_ci=)
dat.fasting.subnetwork.e$class_label<-dat.full.names.subnetwork.e[,9]

save(dat.fasting.subnetwork.e, file = "data_fasting_sub_e.RData")
#corresponding regression coefficients and P-values, node scores
load("data_fasting_sub_e.RData")
dim(dat.fasting.subnetwork.e)

dat.fasting.subnetwork.e$est_upper_ci<-coeff.met.fasting32.subnetwork.e[3,]

dat.fasting.subnetwork.e$est_lower_ci<-coeff.met.fasting32.subnetwork.e[2,]


head(dat.fasting.subnetwork.e)

save(dat.fasting.subnetwork.e, file = "data_fasting_sub_e.RData")

dat.fasting.subnetwork.e$est_ci<-paste("(",as.character(format(dat.fasting.subnetwork.e$est_lower_ci,digits=1)),",",as.character(format(dat.fasting.subnetwork.e$est_upper_ci,digits=1)),")",sep="")

? as.character(dat.fasting.subnetwork.e$est_lower_ci[1])


rownames(dat.fasting.subnetwork.e)
library(xtable)
rownames(dat.fasting.subnetwork.e)<-NULL
xtable(dat.fasting.subnetwork.e[,c(3,2,1,5,8)], digits=8,auto=T)

? xtable
##############################################################
###without age adjusted
pvalue.met32<-function(i){
  
  
  
  dat.analy.fasting_temp2<-data.frame(met=metobolites.final.fasting[,i],cov.other)
  lm1<-lm(bmi_ogtt~met+age_ogtt,data=dat.analy.fasting_temp2)
  summary(lm1)$coefficients[2,4]
  
  
  
}

coeff.met32<-function(i){
  

  
  dat.analy.fasting_temp2<-data.frame(met=metobolites.final.fasting[,i],cov.other)
  lm1<-lm(bmi_ogtt~met+age_ogtt,data=dat.analy.fasting_temp2)
  est<-summary(lm1)$coefficients[2,1]
  
  est_lowci<-est+qnorm(0.025)*summary(lm1)$coefficients[2,2]
  est_upperci<-est+qnorm(0.975)*summary(lm1)$coefficients[2,2]
  
  c(est,est_lowci,est_upperci)
}


pval.met.fasting32<-sapply(1:ncol(metobolites.final.fasting),pvalue.met32)

coeff.met.fasting32<-sapply(1:ncol(metobolites.final.fasting),coeff.met32)

colnames(coeff.met.fasting32)<-colnames(metobolites.final.fasting1)


coeff.met.fasting32.subnetwork.e<-coeff.met.fasting32[,rownames(dat.fasting.subnetwork.e)
]


dim(coeff.met.fasting32.subnetwork.e)



names(pval.met.fasting32)<-colnames(metobolites.final.fasting1)

coeff.met.fasting32.subnetwork.e[3,]




pathway_met_fasting <- read.csv(file = 'name_map_network_linking_pathway.csv')
head(pathway_met_fasting)
dim(pathway_met_fasting)


dat.met.names<-read_excel("data_pull_EU_metabolomics.xlsx",sheet="Annotation",
                          col_names = TRUE )
as.character(dat.met.names[3,2])
dat.met.names[1:5,1]=="id5"
? corrplot

seq(2,434,by=2)
unlist(strsplit(as.character(dat.met.names[3,2]),"1-hr "))

dat.met.names.1<-apply(dat.met.names[2:nrow(dat.met.names),2],1,function(x){unlist(strsplit(as.character(x),"Fasting "))[2]})

dat.met.names11_fasting<-cbind(dat.met.names[2:nrow(dat.met.names),1],dat.met.names.1)

dat.met.names12_fasting<-dat.met.names11_fasting[!is.na(dat.met.names.1),]


dim(dat.met.names12_fasting)

dat.met_names_fasting_final<-dat.met.names12_fasting[,1]

names(dat.met_names_fasting_final)<-dat.met.names12_fasting[,2]


dat.met_names_fasting_final[pathway_met_fasting[,1]]

pathway_met_fasting[1:4,]
pred_met_fasting<-cbind(id=1:nrow(pathway_met_fasting),pathway_met_fasting[,c(1,8,9)],original_name=dat.met_names_fasting_final[pathway_met_fasting[,1]]
)

pred_met_fasting[1:4,]
dat.met.names.11

adj_idx_met_fasting<-sapply(1:nrow(pred_met_fasting),function(x){if(pred_met_fasting[x,3]!=""
) {c1<-as.numeric(unlist(strsplit(as.character(pred_met_fasting[x,3]),",")));
return(cbind(rep(pred_met_fasting[x,1] ,length(c1)),c1))}
else return(NA)

})

lapply(adj_idx_met_fasting,function(x){if(!is.na(x)) return(x)})

adj_idx_met_fasting1<-adj_idx_met_fasting[!sapply(adj_idx_met_fasting, function(x) all(is.na(x)))]


    
x<-adj_idx_met_fasting1[[1]]
for (i in 2:length(adj_idx_met_fasting1)){
x<-rbind(x,adj_idx_met_fasting1[[i]])       
}     

colnames(metobolites.final.fasting)

dim(pred_met_fasting)
pred_met_fasting[1:5,]

pred_met_fasting[,5]


met.names.other<-setdiff(colnames(metobolites.final.fasting),pred_met_fasting[,5])

colnames()
#list of pair of ids in adjacent matrix from fasting predictors with KEGG ids in fasting metabolites 
adj_pair_lst_met_fasting<-x       
             
adj_mat_kegg<-matrix(0,ncol=128,nrow=128)




#adjusted met id, firstly 119 from pathway analysis of KEGG database,
#other 9 no KEGG ids
colnames(adj_mat_kegg)<-c(pred_met_fasting[,5],met.names.other)
rownames(adj_mat_kegg)<-c(pred_met_fasting[,5],met.names.other)
adj_mat_kegg[adj_pair_lst_met_fasting]<-1



#######################################################################################
colnames()



dat.met.names.1hour<-dat.met.names[apply(dat.met.names[,2],1,function(x){grepl("1-hr",x)}),]

dat.met.names.fasting<-dat.met.names[apply(dat.met.names[,2],1,function(x){grepl("Fasting",x)}),2]


dat.met.names.fasting



library(visNetwork)
nodes <- data.frame(id = 1:3,label=c("Antler","Bear","Canoe"))
edges <- data.frame(from = c(1,2), to = c(2,3))
visNetwork(nodes, edges, width = "100%")

devtools::install_github("SimonLarsen/grandforest")


library(MetaboAnalystR)

library(devtools)
devtools::install_github('OSS-Lab/MetQy',subdir = 'MetQy_1.1.0',dependencies = TRUE) 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MAIT")

BiocManager::install("fgsea")
library(BiocManager)
BiocManager::install("KEGGgraph")

library(fgsea)
data(examplePathways)
data(exampleRanks)
fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=10000, maxSize=500)
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
                                      examplePathways, exampleRanks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]


library(fgsea)
library(data.table)
library(ggplot2)

## ----echo=FALSE---------------------------------------------------------------
library(BiocParallel)
register(SerialParam())

## -----------------------------------------------------------------------------
data(examplePathways)

head(exampleRanks)

data(exampleRanks)
set.seed(42)

## -----------------------------------------------------------------------------
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)
#####################################################################
library(readxl)
setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")
dat.met<-read_excel("data_pull_EU_metabolomics.xlsx")
summary(dat.met)

#metabolitic enrichment set analysis
dat_met_mesa<-read.csv("msea_ora_result_network.csv")

summary(dat_met_mesa)

dat_met_mesa[1,c(10,11,9,12)]

mat_met_mesa<-matrix(as.matrix(dat_met_mesa[,c(9:12)]),ncol=4,nrow=nrow(dat_met_mesa))

fisher.test(matrix(as.matrix(dat_met_mesa[1,c(10,11,9,12)]),nrow=2))$p.value

p.value_mesa<-sapply(1:nrow(dat_met_mesa), function(x){fisher.test(matrix(as.matrix(dat_met_mesa[x,c(10,11,9,12)]),nrow=2))$p.value
})
fdr_mesa<-round(p.adjust(p.value_mesa, "BH"), 3)

dat_mesa_tot<-cbind(dat_met_mesa,p.value_mesa,fdr_mesa)

dat_met_mesa_tot<-write.csv(dat_mesa_tot,file="msea_ora_result_network_tot.csv")


#pathway analysis on whole network
dat_met_pathway_tot<-read.csv("pathway_results_network.csv")

dat_met_pathway_subnetwork<-read.csv("pathway_results_subnetwork_e.csv")

summary(dat_met_pathway_tot)

dat_met_pathway_tot1<-dat_met_pathway_tot[,c(1,4)]

dat_met_pathway_tot1$un_hit_tot<-56-dat_met_pathway_tot1[,2]

colnames(dat_met_pathway_tot1)<-c("Pathway","Hit_tot","Un_hit_tot")

dat_met_pathway_subnetwork1<-dat_met_pathway_subnetwork[,c(1,4)]

dat_met_pathway_subnetwork1$un_hit_subnetwork<-18-dat_met_pathway_subnetwork1[,2]

colnames(dat_met_pathway_subnetwork1)<-c("Pathway","Hit_subnetwork","Un_hit_subnetwork")

dat_met_pathway_tot2<-merge(dat_met_pathway_tot1, dat_met_pathway_subnetwork1, by.x="Pathway", by.y="Pathway",all=T)


dat_met_pathway_tot2$Hit_subnetwork[is.na(dat_met_pathway_tot2$Hit_subnetwork)]<-0

dat_met_pathway_tot2$Un_hit_subnetwork[is.na(dat_met_pathway_tot2$Un_hit_subnetwork)]<-18

dat_met_pathway_tot2$Hit_other<-  dat_met_pathway_tot2$Hit_tot-dat_met_pathway_tot2$Hit_subnetwork

dat_met_pathway_tot2$Un_hit_other<-  dat_met_pathway_tot2$Un_hit_tot-dat_met_pathway_tot2$Un_hit_subnetwork

summary(dat_met_pathway_tot2)


fisher.test(matrix(as.matrix(dat_met_pathway_tot2[1,c(4,6,5,7)]),nrow=2))$p.value

p.value_pathway<-sapply(1:nrow(dat_met_pathway_tot2), function(x){fisher.test(matrix(as.matrix(dat_met_pathway_tot2[x,c(4,6,5,7)]),nrow=2))$p.value
})
fdr_pathway<-round(p.adjust(p.value_pathway, "BH"), 3)

dat_pathway_tot<-cbind(dat_met_pathway_tot2,p.value_pathway,fdr_pathway)

met_pathway_result_network_tot<-write.csv(dat_pathway_tot,file="met_pathway_result_network_tot.csv")

min(p.value_pathway)

? merge




p.value_mesa[order(p.value_mesa)]

fdr_mesa<-round(p.adjust(p.value_mesa, "BH"), 3)

p.adjust(p.value_mesa, "fdr")

Convictions <- matrix(c(2, 10, 15, 3), nrow = 2,
                      dimnames =
                        list(c("Dizygotic", "Monozygotic"),
                             c("Convicted", "Not convicted")))
Convictions


sapply()


# Exactly the same p-value, as Cochran's conditions are never met:
fisher.test(MP6, hybrid=TRUE)$p.value




dim(dat.met)
head(dat.met)
ncol(dat.met)
colnames(dat.met)


dat.met1<-dat.met[,-1]
cout.mis<-apply(dat.met1,2,function(x){sum(is.na(x))})

#cout.mis<-sapply(1:ncol(dat.met1),function(x){sum(is.na(as.matrix(dat.met1)[,x]))})

names(cout.mis)<-colnames(dat.met1)
table(cout.mis)

#the number of missing in metabolite variables
cout.mis1<-cout.mis[-c(length(cout.mis):(length(cout.mis)-12))]
table(cout.mis1)

sum(cout.mis1>=200)/length(cout.mis1)
summary(cout.mis1)
#the percentage of missing in metabolite variables
perc.count.mis<-cout.mis1/nrow(dat.met)
hist(cout.mis1,breaks=100,col="skyblue",xlab="Count of Missing")
summary(perc.count.mis)

hist(perc.count.mis,breaks=100,col="skyblue",xlab="Percentage of Missing",main="Histogram")


sum(perc.count.mis>0.1)/length(perc.count.mis)

sum(perc.count.mis>0.5)

sum(perc.count.mis>0.6)

sum(perc.count.mis>0.8)

sum(perc.count.mis>0.9)

sum(perc.count.mis==1)/length(cout.mis1)

#set a cutooff value
perc.mis.cutoff<-0.1
metobolites<-dat.met1[,-c(length(cout.mis):(length(cout.mis)-12))]
dim(metobolites)
metobolites.nonmis<-metobolites[,perc.count.mis<=0.1]
dim(metobolites.nonmis)

sapply(1:5,function(x){abs(x)})
min(1:6,na.rm=T)
#impute missing values with minimum value/2
impute.min.5<-function(x){
  min.x<-min(x,na.rm=T)
  sapply(x,function(y){ifelse(is.na(y),min.x/2.0,y)})
}
metobolites.nonmis.impute<-apply(metobolites.nonmis,2,impute.min.5)

sum(apply(metobolites.nonmis.impute,2,function(z){any(is.na(z))}))
#metabolite with log transformed
metobolites.nonmis.impute.log<-metobolites.nonmis.impute[,grepl("log",colnames(metobolites.nonmis.impute)
)]
dim(metobolites.nonmis.impute.log)



#metabolite without log transformed
metobolites.nonmis.impute.nolog<-metobolites.nonmis.impute[,!grepl("log",colnames(metobolites.nonmis.impute)
)]
dim(metobolites.nonmis.impute.nolog)
hist(metobolites.nonmis.impute.nolog[,102],breaks=100,col="skyblue",xlab="Percentage of Missing",main="Histogram")


#standardize the metabolite covariates
metobolites.nonmis.impute.nolog.std<-apply(metobolites.nonmis.impute.nolog,2,function(x){(log(x)-mean(log(x)))/sd(log(x))})

metobolites.nonmis.impute.log.std<-apply(metobolites.nonmis.impute.log,2,function(x){(x-mean(x))/sd(x)})

metobolites.final<-cbind(metobolites.nonmis.impute.nolog.std,metobolites.nonmis.impute.log.std)
apply(metobolites.final,2,mean)

apply(metobolites.final.fasting,2,sd)

dim()
# 1 hour
metobolites.final.1hour<-metobolites.final[,grepl("_12",colnames(metobolites.final))]
#fasting
metobolites.final.fasting<-metobolites.final[,grepl("_01",colnames(metobolites.final))]

colnames(metobolites.final.fasting)
dim(metobolites.final.fasting)
#generate P values for the nodes
#regression adjusting for field center(field_center), gestational age(ga_ogtt_wks), maternal age(age_ogtt),
#maternal BMI(bmi_ogtt), and mean arterial pressure at OGTT(mapm_ogtt)
#, newborn sex(nn_gender), and sample storage time(storetime_nt_mom_yrs)
#outcome Y-new birth weights (nn_bthwt)
cov.other<-dat.met1[,c("field_center","bmi_ogtt","ga_ogtt_wks","age_ogtt",
                       "mapm_ogtt","nn_gender","storetime_nt_mom_yrs")]

dat.analy.fasting<-data.frame(metobolites.final.fasting,cov.other,dat.met1[,"nn_bthwt"])
Y<-as.vector(unlist(dat.met1[,"nn_bthwt"]))


#export to text file
write.table(metobolites.final.fasting, file = "meta_fasting_dat.txt", sep = "\t",
            row.names = F,col.names = F)



library(igraph)
library(BioNet)
getwd()
setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")
source("edge_score_func.R")

source("pval_perm_corr_edge.R")

names.met.fasting1<-gsub("_", ".", colnames(metobolites.final.fasting))

metobolites.final.fasting1<-metobolites.final.fasting
colnames(metobolites.final.fasting1)<-names.met.fasting1

#pval.edge.met.fasting1<-pval.perm.corr(dat=metobolites.final.fasting1,nsim=100000,MatchId=NA,do.parallel=T,no_cores=4)
names(pval.met.fasting)<-colnames(metobolites.final.fasting1)

#pval.edge.met.1hour<-pval.perm.corr(dat=metobolites.final.1hour1,nsim=100000,MatchId=NA,do.parallel=T,no_cores=4)



range(edge.score.fasting)
edge.score.fasting<-uniform.beta.edge.score(pval=pval.edge.met.fasting1,fdr=0.01)



load("/Users/yubingyao/Google Drive/Network analysis/R code/pval_edge_met_1hour.RData")
load("/Users/yubingyao/Google Drive/Network analysis/R code/pval_edge_met_fasting.RData")

edge.score.1hour<-uniform.beta.edge.score(pval=pval.edge.met.1hour,fdr=0.01)
names(edge.score.1hour)

edge.score.fasting<-uniform.beta.edge.score(pval=pval.edge.met.fasting1,fdr=0.01)

#contruct the test network with the format-igraph based on the dataset of the covariates and edge scores
da.igraph.test<-induced.graph.data.frame(dat=metobolites.final.fasting1,node.score=NA,edge.score=edge.score.fasting,node.weight=NA,edge.weight=NA)


#get the node score based on the network and p.values of the nodes in the network
node.scores.fasting<-node.score(da.igraph=da.igraph.test,pval=pval.met.fasting,fdr=0.2)
node.scores.fasting[node.scores.fasting>0]
network.fasting<-induced.graph.data.frame(dat=metobolites.final.fasting1,node.score=node.scores.fasting,edge.score=edge.score.fasting,node.weight=NA,edge.weight=NA)



names.met.1hour1<-gsub("_", ".", colnames(metobolites.final.1hour))
metobolites.final.1hour1<-metobolites.final.1hour
colnames(metobolites.final.1hour1) <-names.met.1hour1
#contruct the test network with the format-igraph based on the dataset of the covariates and edge scores
da.igraph.test2<-induced.graph.data.frame(dat=metobolites.final.1hour1,node.score=NA,edge.score=edge.score.1hour,node.weight=NA,edge.weight=NA)


##############################################################
###without age adjusted
pvalue.met32<-function(i){
  
  dat.analy.fasting_temp2<-data.frame(met=metobolites.final.fasting[,i],cov.other)
  lm1<-lm(met~bmi_ogtt,data=dat.analy.fasting_temp2)
  summary(lm1)$coefficients[2,4]
  
}
pval.met.fasting32<-sapply(1:ncol(metobolites.final.fasting),pvalue.met32)

names(pval.met.fasting32)<-colnames(metobolites.final.fasting1)
node.scores.fasting32<-node.score(da.igraph=da.igraph.test,pval=pval.met.fasting32,fdr=0.05)
node.scores.fasting32[node.scores.fasting32>0]
length(node.scores.fasting32)
length(node.scores.fasting32[node.scores.fasting32>0])


##############################################################
###without age adjusted
pvalue.met32<-function(i){
  
  
  
  dat.analy.fasting_temp2<-data.frame(met=metobolites.final.fasting[,i],cov.other)
  lm1<-lm(bmi_ogtt~met+age_ogtt,data=dat.analy.fasting_temp2)
  summary(lm1)$coefficients[2,4]
  
  
  
}

coeff.met32<-function(i){
  
  
  
  dat.analy.fasting_temp2<-data.frame(met=metobolites.final.fasting[,i],cov.other)
  lm1<-lm(bmi_ogtt~met+age_ogtt,data=dat.analy.fasting_temp2)
  summary(lm1)$coefficients[2,1]
  
  
  
}


pval.met.fasting32<-sapply(1:ncol(metobolites.final.fasting),pvalue.met32)

coeff.met.fasting32<-sapply(1:ncol(metobolites.final.fasting),coeff.met32)


#######################################################
##1 hour metabolites-metabolite~maternal BMI(+age)
#######################################################
pvalue.met42<-function(i){
  
  dat.analy.1hour_temp2<-data.frame(met=metobolites.final.1hour[,i],cov.other)
  lm1<-lm(met~bmi_ogtt,data=dat.analy.1hour_temp2)
  summary(lm1)$coefficients[2,4]
  
}
pval.met.fasting42<-sapply(1:ncol(metobolites.final.1hour1),pvalue.met42)

names(pval.met.fasting42)<-colnames(metobolites.final.1hour1)
node.scores.fasting42<-node.score(da.igraph=da.igraph.test2,pval=pval.met.fasting42,fdr=0.05)
node.scores.fasting42[node.scores.fasting42>0]


pvalue.met31<-function(i){
  
  dat.analy.fasting_temp2<-data.frame(met=metobolites.final.fasting[,i],cov.other)
  lm1<-lm(met~bmi_ogtt+age_ogtt,data=dat.analy.fasting_temp2)
  summary(lm1)$coefficients[2,4]
  
}
pval.met.fasting31<-sapply(1:ncol(metobolites.final.fasting),pvalue.met31)

names(pval.met.fasting31)<-colnames(metobolites.final.fasting1)
node.scores.fasting31<-node.score(da.igraph=da.igraph.test,pval=pval.met.fasting31,fdr=0.05)
node.scores.fasting31[node.scores.fasting31>0]

median(metobolites.final.fasting1$bmi_ogtt)
#contruct the test network with the format-igraph based on the dataset of the covariates and edge scores
da.igraph.test31<-induced.graph.data.frame(dat=metobolites.final.fasting1,node.score=node.scores.fasting31,edge.score=edge.score.fasting,node.weight=NA,edge.weight=NA)

da.igraph.test32<-induced.graph.data.frame(dat=metobolites.final.fasting1,node.score=node.scores.fasting32,edge.score=edge.score.fasting,node.weight=NA,edge.weight=NA)

##############################################################

#######################################################
##1 hour metabolites-metabolite~maternal BMI(+age)
#######################################################
pvalue.met41<-function(i){
  
  dat.analy.1hour_temp2<-data.frame(met=metobolites.final.1hour[,i],cov.other)
  lm1<-lm(met~bmi_ogtt+age_ogtt,data=dat.analy.1hour_temp2)
  summary(lm1)$coefficients[2,4]
  
}
pval.met.fasting41<-sapply(1:ncol(metobolites.final.1hour),pvalue.met41)

names(pval.met.fasting41)<-colnames(metobolites.final.1hour1)
node.scores.fasting41<-node.score(da.igraph=da.igraph.test2,pval=pval.met.fasting41,fdr=0.05)
node.scores.fasting41[node.scores.fasting41>0]



dat.met.names<-read_excel("data_pull_EU_metabolomics.xlsx",sheet="Annotation",
                          col_names = TRUE )
as.character(dat.met.names[3,2])
dat.met.names[1:5,1]=="id5"
? corrplot

seq(2,434,by=2)
unlist(strsplit(as.character(dat.met.names[3,2]),"1-hr "))

dat.met.names.1<-apply(dat.met.names[2:nrow(dat.met.names),2],1,function(x){unlist(strsplit(as.character(x)," "))[2]})

dat.met.names.11

dat.met.names.1hour<-dat.met.names[apply(dat.met.names[,2],1,function(x){grepl("1-hr",x)}),]

dat.met.names.fasting<-dat.met.names[apply(dat.met.names[,2],1,function(x){grepl("Fasting",x)}),]


met.fasting.shortname<-apply(dat.met.names.fasting[,1],1,function(x){gsub("_",".",x)})
met.fasting.fullname<-apply(dat.met.names.fasting[,2],1,function(x){unlist(strsplit(as.character(x),"Fasting "))[2]})
names(met.fasting.fullname)<-met.fasting.shortname
met.fasting.shortname.final<-node.scores.fasting31[met.fasting.shortname][na.omit(node.scores.fasting31[met.fasting.shortname])]

met.1hr.shortname<-apply(dat.met.names.1hour[,1],1,function(x){gsub("_",".",x)})
met.1hr.fullname<-apply(dat.met.names.1hour[,2],1,function(x){unlist(strsplit(as.character(x),"1-hr "))[2]})
names(met.1hr.fullname)<-met.1hr.shortname

format(x, digits=2, nsmall=2)

met.fasting.fullname.final<-met.fasting.fullname[names(na.omit(node.scores.fasting31[met.fasting.shortname]))]
met.fasting.name.nodescore<-paste(met.fasting.fullname.final,as.character(format(na.omit(node.scores.fasting31[met.fasting.shortname]),digits=1, nsmall=2)),sep=" ")
names(met.fasting.name.nodescore)<-names(na.omit(node.scores.fasting31[met.fasting.shortname]))

met.1hr.fullname.final<-met.1hr.fullname[names(na.omit(node.scores.fasting41[met.1hr.shortname]))]
met.1hr.name.nodescore<-paste(met.1hr.fullname.final,as.character(format(na.omit(node.scores.fasting41[met.1hr.shortname]),digits=1, nsmall=2)),sep=" ")
names(met.1hr.name.nodescore)<-names(na.omit(node.scores.fasting41[met.1hr.shortname]))

#subnetwork of selected nodes and node scores, edge scores
as.matrix(met.fasting.fullname.final,ncol=1)

write.csv(met.fasting.fullname.final,"met_fasting_names.csv", row.names = FALSE)


names(met.fasting.fullname)

names.cluster1_fasting<-names(met.fasting.fullname)[sapply(c("alpha-Tocopherol","Methyl palmitate","Methyl stearate","Methyl linoleate","Erythronic acid/Threonic acid",
  "Glycerol 1-phosphate","Glyceric acid","Myoinositol","Triglycerides"),function(i){which(met.fasting.fullname==i)})]


library(wActNet)
v.id.cluster1<-names.cluster1_fasting
E(da.igraph.test31)$name<-names(edge.score.fasting)
V(da.igraph.test31)$score

V(da.igraph.test31)$name<-met.fasting.fullname.final[V(da.igraph.test31)$name]

node.scores.fasting31

network.cluster1.fasting <- subnetwork.e(da.igraph.test31, vid = v.id.cluster1, eid = names(edge.score.fasting), remove.vertex = F)

V(network.cluster1.fasting)$score<-node.scores.fasting31[V(network.cluster1.fasting)$name]


V(network.cluster1.fasting)$name<-met.fasting.fullname.final[V(network.cluster1.fasting)$name]

V(network.cluster1.fasting)$weight

plotModule(network.cluster1.fasting)

node.scores.fasting42[names.cluster1_fasting]

V(network.cluster1.fasting)$name<-c("alpha-Tocopherol","Methyl palmitate","Methyl stearate","Methyl linoleate","Erythronic acid/Threonic acid",
                                     "Glycerol 1-phosphate","Glyceric acid","Myoinositol","Triglycerides")
library(rgexf)


#correlation threshold-0.2,0.4,0.6 to construct metabolic network
#apply Dittrich's subnetwork detection method-node score only

##########################################################
###correlation threshold to decide edges of network
##########################################################
#with cutoff 0.6,0.7,0.8
##########################################################
cutoff.thre<-0.4
data.corr<-metobolites.final.fasting1
colnames(metobolites.final.fasting1)

corr.thre<-function(data.corr,cutoff.thre){
  
  
  #cor.metobolites.fasting<-cor(metobolites.final.fasting1)
  cor.metobolites.fasting1<-sapply(1:(ncol(data.corr)-1),
                                   function(i){sapply((i+1):ncol(data.corr),function(j){cor(data.corr[,i],data.corr[,j])})})
  
  cor.metobolites.fasting.1<-unlist(cor.metobolites.fasting1)
  
  names.cor.met.fasting1<-sapply(1:(ncol(data.corr)-1),
                                 function(i){sapply((i+1):ncol(data.corr),function(j){paste(colnames(data.corr)[i],colnames(data.corr)[j],sep="_")})})
  
  names(cor.metobolites.fasting.1)<-unlist(names.cor.met.fasting1)
  cor.metobolites.fasting.11<-cor.metobolites.fasting.1
  cor.metobolites.fasting.11[abs(cor.metobolites.fasting.1)<cutoff.thre]<-NA
  
  cor.metobolites.fasting.2<-cor.metobolites.fasting.11[!is.na(cor.metobolites.fasting.11)]
  
  cor.metobolites.fasting.2
}

cor.met.fasting.thre1<-corr.thre(metobolites.final.fasting1,0.25)

cor.met.1hour.thre1<-corr.thre(metobolites.final.1hour1,0.25)


cor.met.fasting.thre2<-corr.thre(metobolites.final.fasting1,0.5)

cor.met.1hour.thre2<-corr.thre(metobolites.final.1hour1,0.5)
length(cor.met.fasting.thre2)

#cutoff 0.4
cor.met.fasting.thre3<-corr.thre(metobolites.final.fasting1,0.4)

cor.met.1hour.thre3<-corr.thre(metobolites.final.1hour1,0.4)
length(cor.met.fasting.thre3)

#cutoff 0.6
cor.met.fasting.thre4<-corr.thre(metobolites.final.fasting1,0.6)

cor.met.1hour.thre4<-corr.thre(metobolites.final.1hour1,0.6)
length(cor.met.fasting.thre4)


#cutoff 0.7
cor.met.fasting.thre5<-corr.thre(metobolites.final.fasting1,0.7)

cor.met.1hour.thre5<-corr.thre(metobolites.final.1hour1,0.7)
length(cor.met.fasting.thre5)
library(wActNet)

network.fasting.cut.corr4<-induced.graph.data.frame(dat=metobolites.final.fasting1,node.score=node.scores.fasting32,edge.score=cor.met.fasting.thre4,node.weight=NA,edge.weight=NA)

network.fasting.cut.corr3<-induced.graph.data.frame(dat=metobolites.final.fasting1,node.score=node.scores.fasting32,edge.score=cor.met.fasting.thre3,node.weight=NA,edge.weight=NA)

lst.mod.fasting.n32_c06<-MultiModuleFind(network.fasting.cut.corr4,node.scores.fasting32,cor.met.fasting.thre4,ncluster=4,method="NodeOnly")



lst.mod.fasting.n32_c04<-MultiModuleFind(network.fasting.cut.corr3,node.scores.fasting32,cor.met.fasting.thre3,ncluster=4,method="NodeOnly")

V(lst.mod.fasting.n32_c04[[1]])$name<-met.fasting.fullname.final[V(lst.mod.fasting.n32_c04[[1]])$name]

V(lst.mod.fasting.n32_c04[[2]])$name<-met.fasting.fullname.final[V(lst.mod.fasting.n32_c04[[2]])$name]


V(lst.mod.fasting.n32_c06[[1]])$name<-met.fasting.fullname.final[V(lst.mod.fasting.n32_c06[[1]])$name]

V(lst.mod.fasting.n32_c06[[2]])$name<-met.fasting.fullname.final[V(lst.mod.fasting.n32_c06[[2]])$name]


####################################################################################################
subnetwork.e<-function(graph,vid,eid,remove.vertex=F){
  
  
  if(is.null(vid)){
    warning("No nodes for subnetwork")  
    break
  }
  
  
  
  
  if(!is.null(V(graph)$score)){
    node.scores<-V(graph)$score
    
    
    
    
    names(node.scores)<-V(graph)$name
    node.score.sub<-node.scores[vid]
    
  }
  if(!is.null(V(graph)$weight)){
    node.weight<-V(graph)$weight
    
    names(node.weight)<-V(graph)$name
    node.weight.sub<-node.weight[vid]
    
  }
  
  
  if(!is.null(E(graph)$weight)){
    edge.weight<-E(graph)$weight
    names(edge.weight)<-E(graph)$name}
  
  
  
  
  if(is.null(eid)){
    warning("No edges for subnetwork")  
    if(is.null(V(graph)$score) && is.null(V(graph)$weight) ){
      pos.nodes.dat<-data.frame(name=vid)}
    
    if(is.null(V(graph)$score) && !is.null(V(graph)$weight) ){
      pos.nodes.dat<-data.frame(name=vid,weight= node.weight.sub)} 
    
    if(!is.null(V(graph)$score) && is.null(V(graph)$weight) ){
      pos.nodes.dat<-data.frame(name=vid,score= node.score.sub)}
    g <- graph_from_data_frame(edge.pos.nodes.dat, directed=F, vertices=pos.nodes.dat)
    return(g)
    break
  }
  
  
  if(!is.null(E(graph)$score)){
    edge.scores<-E(graph)$score
    names(edge.scores)<-E(graph)$name
    edge.score.sub<-edge.scores[eid]
    
  }
  from.name.sub<-unlist(strsplit(E(graph)$name,"_"))[seq(1,2*length(E(graph)$name),by=2)]
  to.name.sub<-unlist(strsplit(E(graph)$name,"_"))[seq(2,2*length(E(graph)$name),by=2)]
  edge.name.sub<-as.matrix(cbind(from.name.sub,to.name.sub))
  if(dim( edge.name.sub)[1]>=1){
    edge.name.sub.node<-matrix(NA,nrow=dim(edge.name.sub)[1],ncol=dim(edge.name.sub)[2])
    
    
    
    for(i in 1:dim( edge.name.sub)[1]){
      if(length(intersect(vid,edge.name.sub[i,]))==2){
        edge.name.sub.node[i,] <-edge.name.sub[i,]
      }
      else
        edge.name.sub.node[i,] <-c(NA,NA)
      
    }
    edge.name.sub.node<-na.omit(edge.name.sub.node)
    from.name.edge.sub.nodes<-edge.name.sub.node[,1]
    to.name.edge.sub.nodes<-edge.name.sub.node[,2]
  }
  
  if( remove.vertex){
    if(length(unique(c( from.name.edge.sub.nodes, to.name.edge.sub.nodes)))<length(unique(vid))){
      vid<-unique(c( from.name.edge.sub.nodes, to.name.edge.sub.nodes))
      if ( dim(edge.name.sub.node)[1] >= 1){
        edge.name.sub.node2<-matrix(NA,nrow=dim(edge.name.sub.node)[1],ncol=dim(edge.name.sub.node)[2])
        
        for(i in 1:dim( edge.name.sub.node)[1]){
          if(length(intersect(vid,edge.name.sub.node[i,]))==2){
            edge.name.sub.node2[i,] <-edge.name.sub[i,]
          }
          else
            edge.name.sub.node2[i,] <-c(NA,NA)
          
        }
        edge.name.sub.node2<-na.omit(edge.name.sub.node2)
        from.name.edge.sub.nodes<-edge.name.sub.node2[,1]
        to.name.edge.sub.nodes<-edge.name.sub.node2[,2]
      }
    }
  }
  
  
  
  if(!is.null(V(graph)$score)){
    
    node.score.sub<-node.scores[vid]
    names(node.score.sub)<-vid
    
  }
  if(!is.null(V(graph)$weight)){
    
    node.weight.sub<-node.weight[vid]
    names(node.weight.sub)<-vid
  }
  
  
  
  
  
  
  
  
  names.edge.sub.nodes<-paste( from.name.edge.sub.nodes,to.name.edge.sub.nodes,sep="_")
  
  if(!is.null(E(graph)$score)){
    edge.score.sub.nodes<-edge.scores[names.edge.sub.nodes]
  }
  
  
  if(!is.null(E(graph)$weight)){
    edge.weight.sub.nodes<-edge.weight[names.edge.sub.nodes]}
  
  
  if(is.null(V(graph)$score) && is.null(V(graph)$weight) ){
    pos.nodes.dat<-data.frame(name=vid)}
  
  if(is.null(V(graph)$score) && !is.null(V(graph)$weight) ){
    pos.nodes.dat<-data.frame(name=vid,weight= node.weight.sub)} 
  
  if(!is.null(V(graph)$score) && is.null(V(graph)$weight) ){
    pos.nodes.dat<-data.frame(name=vid,score= node.score.sub)}
  
  if(!is.null(V(graph)$score) && !is.null(V(graph)$weight) ){
    pos.nodes.dat<-data.frame(name=vid,score= node.score.sub,weight=node.weight.sub)}
  
  if(!is.null(E(graph)$weight) && !is.null(E(graph)$score)){
    edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                   score=edge.score.sub.nodes,name= names.edge.sub.nodes,weight=edge.weight.sub.nodes)}
  
  if(is.null(E(graph)$weight) && !is.null(E(graph)$score) ){
    edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                   score=edge.score.sub.nodes,name= names.edge.sub.nodes)}
  
  
  if(is.null(E(graph)$weight) && is.null(E(graph)$score) ){
    edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                   name= names.edge.sub.nodes)}
  
  
  if(!is.null(E(graph)$weight) && is.null(E(graph)$score) ){
    edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                   weight=edge.weight.sub.nodes,name= names.edge.sub.nodes)}
  
  
  
  g <- graph_from_data_frame(edge.pos.nodes.dat, directed=F, vertices=pos.nodes.dat)
  return(g)}


#
#check any existing isolated nodes in some cluster
#if exisiting some existing isolated nodes, generate one edge between isolated nodes to some random  
#connected nodes within the cluster
#input 1: g-the igraph class object with names of nodes defined as "Varid" index of nodes
#with edges names defined by "Varid1-Varid2" where id1,id2 is the indices of two connecting nodes
#ionput 2: clusterid-the cluter index list in vector form
#output:the new igraph class object with no isolated nodes within the cluster
make_full_connected_cluster<-function(g,clusterid){
  library(igraph)
  
  g.edge.id<-cbind(as.numeric(unlist(strsplit(get.edgelist(g,names=T)[,1],"Var"))[seq(2,2*dim(get.edgelist(g,names=T))[1],2)]),
                   
                   as.numeric(unlist(strsplit(get.edgelist(g,names=T)[,2],"Var"))[seq(2,2*dim(get.edgelist(g,names=T))[1],2)]))
  
  
  ind.cluster1<-unlist(apply(g.edge.id,1,function(x) {length(intersect(x,clusterid))==2}))
  
  
  
  connected.nodes.cluster1<-unique(as.vector(g.edge.id[ind.cluster1==1,]))
  
  if(length(connected.nodes.cluster1)==0){
    
    
    #if no connected edges within all the nodes in the cluster  #make miminal number connections-cluster.size-1 to make all nodes  #connected in the cluster
    clusterid.new<-sample(clusterid,length(clusterid))
    new.cluster.edge<- t(sapply(1:(n.size-1),function(x){c(clusterid.new[x],clusterid.new[x+1])}))
    new.edgelist.cluster11<-t(apply(new.cluster.edge,1,function(x){paste("Var",x,sep="")}))
    #new.cluster.edge1.name<-paste(new.cluster.edge1[,1],new.cluster.edge1[,2],sep="_")
    
  } 
  
  if (length(connected.nodes.cluster1)==length(clusterid)){
    
    return(g)}
  
  
  else if(0 < length(connected.nodes.cluster1) & length(connected.nodes.cluster1) <length(clusterid)){
    
    #identity the subgraph with the nodes with clusterid
    g2<-subnetwork.e(g,vid=V(g)$name[clusterid],eid=E(g)$name)
    
    #identify the separated connected components within the subgraph
    conn.comp.graph <- decompose.graph(g2)
    
    
    #new.nodelist.cluster1<-setdiff(clusterid,connected.nodes.cluster1)
    new.edgelist.cluster11<-t(sapply(1:(length(conn.comp.graph)-1),function(x){c(sample(V(conn.comp.graph[[x]])$name,1)
                                                                                 , sample(V(conn.comp.graph[[x+1]])$name,1))}))
    
  }
  
  
  
  new.edgelist.not.cluster1<-get.edgelist(g,names=T)
  new.edgelist.cluster1<-rbind(new.edgelist.not.cluster1,new.edgelist.cluster11)
  
  from.new.edge.list<-new.edgelist.cluster1[,1]
  to.new.edge.list<-new.edgelist.cluster1[,2]
  names.new.edge.list<-paste(from.new.edge.list,to.new.edge.list,sep="_")
  
  
  
  edge.dat<-data.frame(from= from.new.edge.list,to=to.new.edge.list,name= names.new.edge.list)
  
  nodes.dat<-data.frame(name=V(g)$name)
  
  interactome2 <- graph_from_data_frame(edge.dat, directed=F, vertices=nodes.dat)
  
  
  # interactome2 <- graph.edgelist(new.edgelist.cluster1, F)
  # E(g)$name<-paste(get.edgelist(interactome2,name=T)[,1],get.edgelist(interactome2,name=T)[,2],sep="_")
  
  return(interactome2)
  
}


MultiModuleFind_cos_dat<-function(network,node.scores,edge.scores,ncluster,clustersize){
  
  #subnetwork to induced subgraph from graph object considering both the vertices and edges,
  #with remove.vertex=T,F; if True, then remove the vertices which is not in the edges, otherwise keep them
  #vid-the names of the vertices in the graph
  #eid-the names of the edges in the graph, with the format-Name(from_vertex)_Name(to_vertex)
  #output:the subgraph containing only the vertices of vid and the edges
  library(COSINE)
  if(is.null(E(network)$name)){
    E(network)$name<-names(edge.scores)}
  
  num.cluster<-1
  node.scores0<-node.scores
  edge.score0<-edge.scores
  lst.modules<-list()
  mat_e_idx<-matrix(unlist(strsplit(names(edge.scores),"_")),ncol=2,nrow=length(edge.scores),byrow=T)
  
  GA_result10<-GA_search_PPI(lambda=0.5,node.scores,edge.scores,mat_e_idx,
                             num_iter=20000, muCh=0.1, zToR=10, minsize=clustersize)
  module.10<-GA_result10$Subnet
  
  module.1<-GA_result10$Subnet
  
  while( !(is.null(module.1))  && num.cluster<ncluster ){
    
    
    lst.modules[[num.cluster]]<-module.10
    
    vid.after.cluster1<-setdiff(names(node.scores)
                                ,names(node.scores)[module.1])
    
    network.after.cluster1<-subnetwork.e(network,vid=vid.after.cluster1,eid=names(edge.scores),remove.vertex=F)
    
    node.scores2<-node.scores[vid.after.cluster1]
    eid.after.cluster1<-E(network.after.cluster1)$name
    
    
    
    edge.scores2<-edge.scores[eid.after.cluster1]
    
    mat_e_idx2<-matrix(unlist(strsplit(names(edge.scores2),"_")),ncol=2,nrow=length(edge.scores2),byrow=T)
    GA_result20<-GA_search_PPI(lambda=0.5,node.scores2,edge.scores2,mat_e_idx2,
                               num_iter=20000, muCh=0.1, zToR=10, minsize=clustersize)
    
    
    module.20<-sapply(names(node.scores2)[GA_result20$Subnet],function(i){which(names(node.scores0)==i)})
    
    module.2<-GA_result20$Subnet
    
    num.cluster<-num.cluster+1
    
    if(num.cluster==ncluster && (!is.null(module.2)) ){
      lst.modules[[num.cluster]]<-module.20
    }
    network<-network.after.cluster1
    node.scores<-node.scores2
    
    edge.scores<-edge.scores2
    module.1<-module.2
    module.10<-module.20
    
  }
  
  lst.modules
  
}
out_dir<-"/Users/yubingyao/Google Drive/Network analysis/R code/"
module.fasting1.1_cos<-MultiModuleFind_cos_dat(network=da.igraph.test32, node.scores=node.scores.fasting32,edge.scores=edge.score.fasting,ncluster=4,clustersize=10)
save(module.fasting1.1_cos,file=paste(out_dir,"fasting32_",as.character(10),"n_iter",as.character(4),".Rdata",sep=""))

library(igraph)
load("fasting32_10n_iter4.Rdata")

V(module.fasting1.1_cos)$name

node.scores.fasting32[module.fasting1.1_cos[[1]]]

node.scores.fasting32[module.fasting1.1_cos[[2]]]

node.scores.fasting32[module.fasting1.1_cos[[3]]]

node.scores.fasting32[module.fasting1.1_cos[[4]]]

vid_cos_fasting<-c(module.fasting1.1_cos[[1]],module.fasting1.1_cos[[2]],module.fasting1.1_cos[[3]])
unlist(module.fasting1.1_cos)
node.scores.fasting32[vid_cos_fasting]
V(da.igraph.test32)$name[vid_cos_fasting]

E(da.igraph.test32)$name<-names(edge.score.fasting)

module.fasting_cos<-subnetwork.e(da.igraph.test32,vid=V(da.igraph.test32)$name[vid_cos_fasting],eid=names(edge.score.fasting),remove.vertex=F)



V(module.fasting_cos)$name<-met.fasting.fullname.final[V(module.fasting_cos)$name]  

V(module.fasting_cos)$score

V(da.igraph.test32)$name
E(da.igraph.test32)$score
? igraph.to.gexf
# Converts the given igraph object to GEXF format and saves it at the given filepath location
#     g: input igraph object to be converted to gexf format
#     filepath: file location where the output gexf file should be saved
#
saveAsGEXF = function(g, filepath="converted_graph.gexf")
{
  require(igraph)
  require(rgexf)
  require(XML)
  
  # gexf nodes require two column data frame (id, label)
  # check if the input vertices has label already present
  # if not, just have the ids themselves as the label
  if(is.null(V(g)$name))
    V(g)$label <- as.character(V(g))
  if(!is.null(V(g)$name))
    V(g)$label <- V(g)$name
  
  # similarily if edges does not have weight, add default 1 weight
  if(is.null(E(g)$score))
    E(g)$weight <- rep.int(1, ecount(g))
  if(!is.null(E(g)$score))
    E(g)$weight <-E(g)$score
  
  nodes <- data.frame(cbind(V(g), V(g)$label))
  edges <- t(Vectorize(get.edge, vectorize.args='id')(g, 1:ecount(g)))
  
  # combine all node attributes into a matrix (and take care of & for xml)
  vAttrNames <- setdiff(list.vertex.attributes(g), "label") 
  nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&",get.vertex.attribute(g, attr))))
  
  # combine all edge attributes into a matrix (and take care of & for xml)
  eAttrNames <- setdiff(list.edge.attributes(g), "weight") 
  edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&",get.edge.attribute(g, attr))))
  
  # combine all graph attributes into a meta-data
  graphAtt <- sapply(list.graph.attributes(g), function(attr) sub("&", "&",get.graph.attribute(g, attr)))
  
  # generate the gexf object
  output <- write.gexf(nodes, edges, 
                       edgesWeight=E(g)$weight,
                       edgesAtt = edgesAtt,
                       nodesAtt = nodesAtt,
                       meta=c(list(creator="Yubing Yao", description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt),
                       output=filepath)
  
  print(output, filepath, replace=T)
}
saveAsGEXF(network.cluster1.fasting, filepath=paste0("/Users/yubingyao/Google Drive/Network analysis/R code/","fasting_cluster",as.character(1),".gexf",sep=""))


V(dat.network[[1]])$name
E(dat.network[[1]])$score

? order
#construct new optimized subnetwork from our wActNet algorithm 
#removing those edges with edge score less than or equal to 0
dat.network.sub<-dat.network[[1]]


dat.network.sub1<-delete_edges(dat.network.sub,which(E(dat.network.sub)$score<=0))

V(dat.network.sub1)
library(visNetwork)
nodes <- data.frame(id = V(dat.network.sub1),label=V(dat.network.sub1)$name)
edges <- data.frame(from = c(1,2), to = c(2,3))
visNetwork(nodes, edges, width = "100%")


# get data and plot :
data <- toVisNetworkData(dat.network.sub1)
visNetwork(nodes = data$nodes, edges = data$edges, height = "500px")
#4 nodes without any connecting edges with positive edge scores
#select connecting edges with minimum edge scores

Isolated = which(degree(dat.network.sub1)==0)

not.Isolated = which(degree(dat.network.sub1)>0)

v.niso_name<-V(dat.network.sub1)$name[not.Isolated]


v.iso_name<-V(dat.network.sub1)$name[Isolated]

V(dat.network.sub1)$score[Isolated]



which(E(dat.network.sub)$name<=0)

E(dat.network.sub)$name[1:2]

edge.name.full<-matrix(unlist(strsplit(E(dat.network.sub)$name,"_")),ncol=2,nrow=length(E(dat.network.sub)$name),byrow = T)



library(readxl)
setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")
dat.met<-read_excel("data_pull_EU_metabolomics.xlsx")
dim(dat.met)
head(dat.met)
ncol(dat.met)

dat.met.names<-read_excel("data_pull_EU_metabolomics.xlsx",sheet="Annotation",
                          col_names = TRUE )


edge.name.full[1,]

dat.met.names.fasting<-dat.met.names[apply(dat.met.names[,2],1,function(x){grepl("Fasting",x)}),]


met.fasting.shortname<-apply(dat.met.names.fasting[,1],1,function(x){gsub("_",".",x)})
met.fasting.fullname<-apply(dat.met.names.fasting[,2],1,function(x){unlist(strsplit(as.character(x),"Fasting "))[2]})
names(met.fasting.fullname)<-met.fasting.shortname

met.fasting.fullname[edge.name.full[1,]]

edge.name.full_mat<-apply(edge.name.full,2,function(x){met.fasting.fullname[x]})

edge.name.full.new<-paste(edge.name.full_mat[,1], edge.name.full_mat[,2], sep = "_")

idx.2iso<-cbind(edge.name.full_mat[,1] %in% v.iso_name,edge.name.full_mat[,2] %in% v.iso_name)

idx_iso_niso_max_edge<-sapply(v.iso_name,function(x){which.max(E(dat.network.sub)$score[(apply(idx.2iso,1,sum)!=2) & grepl(x, edge.name.full.new)])})

idx_iso_niso_max_edge2<-sapply(v.iso_name,function(x){which.max(E(dat.network.sub)$score[grepl(x, edge.name.full.new)])})


dat.network.sub.score1<-E(dat.network.sub)$score[(apply(idx.2iso,1,sum)!=2) & grepl(v.iso_name[1], edge.name.full.new)]

dat.network.sub.score2<-E(dat.network.sub)$score[(apply(idx.2iso,1,sum)!=2) & grepl(v.iso_name[2], edge.name.full.new)]

dat.network.sub.score3<-E(dat.network.sub)$score[(apply(idx.2iso,1,sum)!=2) & grepl(v.iso_name[3], edge.name.full.new)]

dat.network.sub.score4<-E(dat.network.sub)$score[(apply(idx.2iso,1,sum)!=2) & grepl(v.iso_name[4], edge.name.full.new)]



dat.network.sub.name1<-E(dat.network.sub)$name[(apply(idx.2iso,1,sum)!=2) & grepl(v.iso_name[1], edge.name.full.new)]

dat.network.sub.name2<-E(dat.network.sub)$name[(apply(idx.2iso,1,sum)!=2) & grepl(v.iso_name[2], edge.name.full.new)]

dat.network.sub.name3<-E(dat.network.sub)$name[(apply(idx.2iso,1,sum)!=2) & grepl(v.iso_name[3], edge.name.full.new)]

dat.network.sub.name4<-E(dat.network.sub)$name[(apply(idx.2iso,1,sum)!=2) & grepl(v.iso_name[4], edge.name.full.new)]


#scores of 4 added edges
e.score.iso.niso<-c(dat.network.sub.score1[idx_iso_niso_max_edge[1]],dat.network.sub.score2[idx_iso_niso_max_edge[2]],
  dat.network.sub.score3[idx_iso_niso_max_edge[3]],dat.network.sub.score4[idx_iso_niso_max_edge[4]])
#names of 4 added edges
e.name.iso.niso<-c(dat.network.sub.name1[idx_iso_niso_max_edge[1]],dat.network.sub.name2[idx_iso_niso_max_edge[2]],
  dat.network.sub.name3[idx_iso_niso_max_edge[3]],dat.network.sub.name4[idx_iso_niso_max_edge[4]])

unlist(sapply(e.name.iso.niso,function(x){strsplit(x,"_")}))


dat.network.sub11<-delete_edges(dat.network.sub,which(E(dat.network.sub)$score<=0))


idx.vnode<-sapply(unlist(sapply(e.name.iso.niso,function(x){strsplit(x,"_")})),function(x){which(V(dat.network.sub11)$name==x)})

dat.network.sub12<-add_edges(dat.network.sub11,idx.vnode)

E(dat.network.sub12)$name<-c(E(dat.network.sub11)$name,e.name.iso.niso)

E(dat.network.sub12)$score<-c(E(dat.network.sub11)$score,e.score.iso.niso)
names(E(dat.network.sub12)$score)<-c(E(dat.network.sub11)$name,e.name.iso.niso)


save(dat.network.sub12, file = "network_fasting_w_plot.RData")

visNetwork(nodes, edges, height = "500px", width = "100%")

getwd()
setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")
load("network_fasting_w_plot.RData")
library(visNetwork)
library(igraph)
# get data and plot :
E(dat.network.sub12)$dashes<-ifelse(E(dat.network.sub12)$score>0,FALSE,TRUE)
E(dat.network.sub12)$width<-(as.integer(abs(E(dat.network.sub12)$score))+1)
V(dat.network.sub12)$color<-ifelse(V(dat.network.sub12)$score>0,"#E7298A","purple")
V(dat.network.sub12)$size<-(as.integer(abs(V(dat.network.sub12)$score)))+20

data <- toVisNetworkData(dat.network.sub12)
visNetwork(nodes = data$nodes, edges = data$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))

dash<-ifelse(E(dat.network.sub12)$score>0,T,F)

width = as.integer(abs(E(dat.network.sub12)$score))+1
"#33F00"

V(dat.network.sub12)$color<-ifelse(V(dat.network.sub12)$score>0,"#33F00","purple")

library("RColorBrewer")
display.brewer.all()
?? visNetwork

#illustrative example of 4 clusters in simulation setting,20 nodes,1 to 20-node name
#cluster 1-large node score, large edge score
g11 <- barabasi.game(20,m=5,power=1.5,directed=F)
V(g11)$name<-as.character(1:20)
V(g11)$score<-runif(20,min=4,max=6)
E(g11)$score<-runif(length(E(g11)),min=4,max=6)
#cluster 2-moderate node score, large edge score
g12 <- barabasi.game(20,m=5,power=1.5,directed=F)
V(g12)$name<-as.character(1:20)
V(g12)$score<-runif(20,min=0.1,max=2)
E(g12)$score<-runif(length(E(g12)),min=4,max=6)
#cluster 3-large node score, moderate edge score
g13 <- barabasi.game(20,m=5,power=1.5,directed=F)
V(g13)$name<-as.character(1:20)
V(g13)$score<-runif(20,min=4,max=6)
E(g13)$score<-runif(length(E(g13)),min=0.1,max=2)
#cluster 4-moderate node score, moderate edge score
g14 <- barabasi.game(20,m=5,power=1.5,directed=F)
V(g14)$name<-as.character(1:20)
V(g14)$score<-runif(20,min=0.1,max=2)
E(g14)$score<-runif(length(E(g14)),min=0.1,max=2)

plot(g14)
#E(g)$name<-paste(as.character(1:20),get.edgelist(g,name=T)[,2],sep="_")
# get data and plot :
E(g14)$dashes<-ifelse(E(g14)$score>0,FALSE,TRUE)
E(g14)$width<-E(g14)$score-1
V(g14)$color<-ifelse(V(g14)$score>0,"#E7298A","purple")
V(g14)$size<-V(g14)$score+4

data14 <- toVisNetworkData(g14)
visNetwork(nodes = data14$nodes, edges = data14$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))


E(g12)$dashes<-ifelse(E(g12)$score>0,FALSE,TRUE)
E(g12)$width<-E(g12)$score-1
V(g12)$color<-ifelse(V(g12)$score>0,"#E7298A","purple")
V(g12)$size<-V(g12)$score+4

data12 <- toVisNetworkData(g12)
visNetwork(nodes = data12$nodes, edges = data12$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))



E(g13)$dashes<-ifelse(E(g13)$score>0,FALSE,TRUE)
E(g13)$width<-E(g13)$score-1
V(g13)$color<-ifelse(V(g13)$score>0,"#E7298A","purple")
V(g13)$size<-V(g13)$score+12

data13 <- toVisNetworkData(g13)
visNetwork(nodes = data13$nodes, edges = data13$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))

E(g11)$dashes<-ifelse(E(g11)$score>0,FALSE,TRUE)
E(g11)$width<-E(g11)$score-1
V(g11)$color<-ifelse(V(g11)$score>0,"#E7298A","purple")
V(g11)$size<-V(g11)$score+12

data11 <- toVisNetworkData(g11)
visNetwork(nodes = data11$nodes, edges = data11$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))

title1="Cluster 1"
title2="Cluster 2"
title3="Cluster 3"
title4="Cluster 4"

network.test<-network
dev.off()
par(mfrow=c(2,2), mar=c(1,0.5,1,0.5))

visNetwork(nodes = data11$nodes, edges = data11$edges, main = "Cluster 1", height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))

visNetwork(nodes = data12$nodes, edges = data12$edges, main = "Cluster 2", height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))

visNetwork(nodes = data13$nodes, edges = data13$edges, main = "Cluster 3", height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))

visNetwork(nodes = data14$nodes, edges = data14$edges, main = "Cluster 4",height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))

library(visNetwork)

visDocumentation()




meta_name_cluster1_corr0p4<-c("AC C2",
"AC C4-OH",
"AC C8:1",
"AC C12:1",
"AC C14:2",
"AC C16:1",
"AC C16:1-OH/C14:1-DC",
"AC C16-OH/C14-DC",
"AC C18:1",
"AC C20-OH/C18-DC",
"3-Hydroxybutyric acid",
"beta-hydroxybutyrate",
"Palmitoleic acid",
"Glycerol")

V(lst.mod.fasting.n32_c06[[1]])$score

E(lst.mod.fasting.n32_c04[[2]])$score


g_d_corr04_c1<-lst.mod.fasting.n32_c04[[1]]

g_d_corr04_c2<-lst.mod.fasting.n32_c04[[2]]


g_d_corr06_c1<-lst.mod.fasting.n32_c06[[1]]

g_d_corr06_c2<-lst.mod.fasting.n32_c06[[2]]


E(g_d_corr04_c1)$dashes<-ifelse(E(g_d_corr04_c1)$score>0,FALSE,TRUE)
E(g_d_corr04_c1)$width<-E(g_d_corr04_c1)$score*5
V(g_d_corr04_c1)$color<-ifelse(V(g_d_corr04_c1)$score>0,"#E7298A","purple")
V(g_d_corr04_c1)$size<-V(g_d_corr04_c1)$score+12

d_corr04_c1 <- toVisNetworkData(g_d_corr04_c1)
visNetwork(nodes = d_corr04_c1$nodes, edges = d_corr04_c1$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))


E(g_d_corr04_c2)$dashes<-ifelse(E(g_d_corr04_c2)$score>0,FALSE,TRUE)
E(g_d_corr04_c2)$width<-E(g_d_corr04_c2)$score*5
V(g_d_corr04_c2)$color<-ifelse(V(g_d_corr04_c2)$score>0,"#E7298A","purple")
V(g_d_corr04_c2)$size<-V(g_d_corr04_c2)$score+12

d_corr04_c2 <- toVisNetworkData(g_d_corr04_c2)
visNetwork(nodes = d_corr04_c2$nodes, edges = d_corr04_c2$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))




E(g_d_corr06_c1)$dashes<-ifelse(E(g_d_corr06_c1)$score>0,FALSE,TRUE)
E(g_d_corr06_c1)$width<-E(g_d_corr06_c1)$score*5
V(g_d_corr06_c1)$color<-ifelse(V(g_d_corr06_c1)$score>0,"#E7298A","purple")
V(g_d_corr06_c1)$size<-V(g_d_corr06_c1)$score+12

d_corr06_c1 <- toVisNetworkData(g_d_corr06_c1)
visNetwork(nodes = d_corr06_c1$nodes, edges = d_corr06_c1$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))


E(g_d_corr06_c2)$dashes<-ifelse(E(g_d_corr06_c2)$score>0,FALSE,TRUE)
E(g_d_corr06_c2)$width<-E(g_d_corr06_c2)$score*5
V(g_d_corr06_c2)$color<-ifelse(V(g_d_corr06_c2)$score>0,"#E7298A","purple")
V(g_d_corr06_c2)$size<-V(g_d_corr06_c2)$score+12

d_corr06_c2 <- toVisNetworkData(g_d_corr06_c2)
visNetwork(nodes = d_corr06_c2$nodes, edges = d_corr04_c2$edges, height = "600px")%>% visEdges(color = list(color = "#66A61E", highlight = "red"))





meta_name_cluster1_corr0p4<-c(
"Aspartic acid",
"Asparagine/Aspartic acid",
"Arginine",
"Erythronic acid/Threonic acid",
"Glutamic acid",
"Threonine",
"Myoinositol",
"Glutamine/Glutamic acid",
"Leucine/Isoleucine",
"Phenylalanine",
"Serine",
"Valine")


network.cluster.plot(network.sim.ba.m24.power15.clustersize20,1,20,title1)
network.cluster.plot(network.sim.ba.m24.power15.clustersize20,floor(vcount(network.test)/2),20,title2)



V(dat.network.sub11)$name

g <- make_empty_graph(n = 5) %>%
  add_edges(c(1,2, 2,3, 3,4, 4,5)) %>%
  set_edge_attr("color", value = "red") %>%
  add_edges(c(5,1), color = "green")
E(g)[[]]
plot(g)


dat.network.sub.score2[idx_iso_niso_max_edge[2]]

dat.network.sub.name2[idx_iso_niso_max_edge[2]]

dat.network.sub.name1[idx_iso_niso_max_edge[1]]



E(dat.network.sub)$score[grepl(v.iso_name[1], edge.name.full.new)]

E(dat.network.sub)$name[grepl(v.iso_name[1], edge.name.full.new)]


E(dat.network.sub)$name[idx_iso_niso_max_edge2]

E(dat.network.sub)$score[idx_iso_niso_max_edge2]




edge.name.full_mat[which(apply(idx.2iso,1,sum)==2),]


!(apply(idx.2iso,1,sum)==2) & grepl(v.iso_name[1], edge.name.full.new)

v.iso_name

E(dat.network.sub)$name[which.max(E(dat.network.sub)$score[grepl(v.iso_name[4], edge.name.full.new)])]
  
E(dat.network.sub)$score[which.max(E(dat.network.sub)$score[grepl(v.iso_name[4], edge.name.full.new)])]

names(E(dat.network.sub)$score)<-edge.name.full.new

E(dat.network.sub)$name<-edge.name.full.new


edge.name.full_mat



edge.name.full.new[grepl(v.iso_name[1], edge.name.full.new, fixed = TRUE)]

V()


v.iso_name

gsub()


E(dat.network.sub)$name

E(dat.network.sub)$name

E(dat.network.sub)$name



? matrix

as.vector(E(dat.network[[1]])$name[E(dat.network[[1]])$score>0])

sort(as.vector(E(dat.network[[1]])$name[E(dat.network[[1]])$score>0]))

sort(E(dat.network.sub1)$name)




index(E(dat.network[[1]])$score,E(dat.network[[1]])$score<=0)
E(dat.network.sub)$score<-E(dat.network[[1]])$score[E(dat.network[[1]])$score>0]
E(dat.network.sub)$name<-E(dat.network[[1]])$name[E(dat.network[[1]])$score>0]


#gene enrichment analysis 
metabolite_set<-c("Hartnup Disease",
                  "Alzheimer's Disease",
"Sarcosinemia",
"Dementia")

FDR_adj<-c(2.13E-04,
           0.00913,
           0.0802,
           0.0802)
trans_FDR_set=-log(FDR_adj,10)

#pathway analysis
pathway_kegg<-c("Aminoacyl-tRNA biosynthesis",
"Arginine biosynthesis",
"Glyoxylate and dicarboxylate metabolism",
"Glycine, serine and threonine metabolism")
FDR_path<-c(2.05E-08,
            0.00039037,
            0.0082749,
            0.10252
)

trans_FDR_path=-log(FDR_path,10)

library(ggplot2) 
df_enrich<-data.frame(metabolite_set,trans_FDR_set)

p1 <- ggplot(df_enrich, aes(x = metabolite_set, y = trans_FDR_set))+ geom_col() + xlab("Metabolite Set in SMPDB") + ylab("-log10(FDR) in Metolite Enrichment Analysis")+geom_bar(stat="identity",  fill="#56B4E9")


p1 + coord_flip()


df_path<-data.frame(pathway_kegg,trans_FDR_path)

p2 <- ggplot(df_path, aes(x = pathway_kegg, y = trans_FDR_path))+ geom_col() + xlab("Pathway in KEGG") + ylab("-log10(FDR) in Pathway Analysis")+geom_bar(stat="identity",  fill="#56B4E9")+scale_x_discrete(guide = guide_axis(n.dodge = 2)) 
p2 + geom_hline(yintercept=2, linetype="dashed", color = "purple", size=1.5)

+ geom_hline(yintercept=1.30103, linetype="dashed", color = "purple", size=1.5)



+ coord_flip()


+opts(axis.text.y = theme_text(family = "sans", face = "bold", size = 12)) +opts(axis.text.x = theme_text(family = "sans", face = "bold", size = 12))



