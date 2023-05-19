#####################################################################################################
#please comment the below code out if you already install the package-wActNet from github
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")
install_github("yuy113/wActNet")
###############################################################################################
library(readxl)
library(wActNet)
library(igraph)
library(BioNet)
#setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")
#read and load HAPO real dataset
dat.met<-read_excel("data_pull_EU_metabolomics.xlsx")
summary(dat.met)
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


#sum(perc.count.mis>0.1)/length(perc.count.mis)

#sum(perc.count.mis>0.5)

#sum(perc.count.mis>0.6)

#sum(perc.count.mis>0.8)

#sum(perc.count.mis>0.9)

#sum(perc.count.mis==1)/length(cout.mis1)

#set a cutooff value
perc.mis.cutoff<-0.1
metobolites<-dat.met1[,-c(length(cout.mis):(length(cout.mis)-12))]
metobolites.nonmis<-metobolites[,perc.count.mis<=0.1]

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



#metabolite without log transformed
metobolites.nonmis.impute.nolog<-metobolites.nonmis.impute[,!grepl("log",colnames(metobolites.nonmis.impute)
)]
dim(metobolites.nonmis.impute.nolog)
hist(metobolites.nonmis.impute.nolog[,102],breaks=100,col="skyblue",xlab="Percentage of Missing",main="Histogram")


#standardize the metabolite covariates
metobolites.nonmis.impute.nolog.std<-apply(metobolites.nonmis.impute.nolog,2,function(x){(log(x)-mean(log(x)))/sd(log(x))})

metobolites.nonmis.impute.log.std<-apply(metobolites.nonmis.impute.log,2,function(x){(x-mean(x))/sd(x)})

metobolites.final<-cbind(metobolites.nonmis.impute.nolog.std,metobolites.nonmis.impute.log.std)
#apply(metobolites.final,2,mean)

#apply(metobolites.final.fasting,2,sd)


# 1 hour
metobolites.final.1hour<-metobolites.final[,grepl("_12",colnames(metobolites.final))]
#fasting
metobolites.final.fasting<-metobolites.final[,grepl("_01",colnames(metobolites.final))]

#colnames(metobolites.final.fasting)
#dim(metobolites.final.fasting)
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


#setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")
#source("edge_score_func.R")

#source("pval_perm_corr_edge.R")

names.met.fasting1<-gsub("_", ".", colnames(metobolites.final.fasting))

metobolites.final.fasting1<-metobolites.final.fasting
colnames(metobolites.final.fasting1)<-names.met.fasting1

pval.edge.met.fasting1<-pval.perm.corr(dat=metobolites.final.fasting1,nsim=100000,MatchId=NA,do.parallel=T,no_cores=4)
names(pval.met.fasting)<-colnames(metobolites.final.fasting1)

pval.edge.met.1hour<-pval.perm.corr(dat=metobolites.final.1hour1,nsim=100000,MatchId=NA,do.parallel=T,no_cores=4)



edge.score.fasting<-uniform.beta.edge.score(pval=pval.edge.met.fasting1,fdr=0.01)



#load("/Users/yubingyao/Google Drive/Network analysis/R code/pval_edge_met_1hour.RData")
#load("/Users/yubingyao/Google Drive/Network analysis/R code/pval_edge_met_fasting.RData")

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
#contruct the  network with the format-igraph based on the dataset of the covariates and edge scores
da.igraph.test2<-induced.graph.data.frame(dat=metobolites.final.1hour1,node.score=NA,edge.score=edge.score.1hour,node.weight=NA,edge.weight=NA)


##############################################################

##############################################################
###with age adjusted
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


dat.met.names<-read_excel("data_pull_EU_metabolomics.xlsx",sheet="Annotation",
                          col_names = TRUE )

dat.met.names.1<-apply(dat.met.names[2:nrow(dat.met.names),2],1,function(x){unlist(strsplit(as.character(x)," "))[2]})


#dat.met.names.1hour<-dat.met.names[apply(dat.met.names[,2],1,function(x){grepl("1-hr",x)}),]

dat.met.names.fasting<-dat.met.names[apply(dat.met.names[,2],1,function(x){grepl("Fasting",x)}),]


met.fasting.shortname<-apply(dat.met.names.fasting[,1],1,function(x){gsub("_",".",x)})
met.fasting.fullname<-apply(dat.met.names.fasting[,2],1,function(x){unlist(strsplit(as.character(x),"Fasting "))[2]})
names(met.fasting.fullname)<-met.fasting.shortname
met.fasting.shortname.final<-node.scores.fasting31[met.fasting.shortname][na.omit(node.scores.fasting31[met.fasting.shortname])]

#met.1hr.shortname<-apply(dat.met.names.1hour[,1],1,function(x){gsub("_",".",x)})
#met.1hr.fullname<-apply(dat.met.names.1hour[,2],1,function(x){unlist(strsplit(as.character(x),"1-hr "))[2]})
#names(met.1hr.fullname)<-met.1hr.shortname

format(x, digits=2, nsmall=2)

met.fasting.fullname.final<-met.fasting.fullname[names(na.omit(node.scores.fasting31[met.fasting.shortname]))]
met.fasting.name.nodescore<-paste(met.fasting.fullname.final,as.character(format(na.omit(node.scores.fasting31[met.fasting.shortname]),digits=1, nsmall=2)),sep=" ")
names(met.fasting.name.nodescore)<-names(na.omit(node.scores.fasting31[met.fasting.shortname]))



network.fasting2.n<-induced.graph.data.frame(dat=metobolites.final.fasting1,node.score=node.scores.fasting32,edge.score=edge.score.fasting.pos,node.weight=NA,edge.weight=NA)

module.fasting2.1n<-runFastHeinz(network.fasting2.n, node.scores.fasting32)
#plotModule(module.fasting2.1n)

module.fasting2<-runFastHeinz.e(network.fasting, node.scores.fasting32,edge.score.fasting)

lst.mod.fasting.e<-MultiModuleFind(network.fasting1,node.scores.fasting32,edge.score.fasting,ncluster=4,method="NodeEdge")


#met.1hr.fullname.final<-met.1hr.fullname[names(na.omit(node.scores.fasting41[met.1hr.shortname]))]
#met.1hr.name.nodescore<-paste(met.1hr.fullname.final,as.character(format(na.omit(node.scores.fasting41[met.1hr.shortname]),digits=1, nsmall=2)),sep=" ")
#names(met.1hr.name.nodescore)<-names(na.omit(node.scores.fasting41[met.1hr.shortname]))

#subnetwork of selected nodes and node scores, edge scores
#as.matrix(met.fasting.fullname.final,ncol=1)

#write.csv(met.fasting.fullname.final,"met_fasting_names.csv", row.names = FALSE)


#names(met.fasting.fullname)

names.cluster1_fasting<-names(met.fasting.fullname)[sapply(c("alpha-Tocopherol","Methyl palmitate","Methyl stearate","Methyl linoleate","Erythronic acid/Threonic acid",
  "Glycerol 1-phosphate","Glyceric acid","Myoinositol","Triglycerides"),function(i){which(met.fasting.fullname==i)})]


library(wActNet)
v.id.cluster1<-names.cluster1_fasting
E(da.igraph.test31)$name<-names(edge.score.fasting)

V(da.igraph.test31)$name<-met.fasting.fullname.final[V(da.igraph.test31)$name]


network.cluster1.fasting <- subnetwork.e(da.igraph.test31, vid = v.id.cluster1, eid = names(edge.score.fasting), remove.vertex = F)

V(network.cluster1.fasting)$score<-node.scores.fasting31[V(network.cluster1.fasting)$name]


V(network.cluster1.fasting)$name<-met.fasting.fullname.final[V(network.cluster1.fasting)$name]

#V(network.cluster1.fasting)$weight

#plotModule(network.cluster1.fasting)

#node.scores.fasting42[names.cluster1_fasting]

V(network.cluster1.fasting)$name<-c("alpha-Tocopherol","Methyl palmitate","Methyl stearate","Methyl linoleate","Erythronic acid/Threonic acid",
                                     "Glycerol 1-phosphate","Glyceric acid","Myoinositol","Triglycerides")



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

#cor.met.fasting.thre1<-corr.thre(metobolites.final.fasting1,0.25)

#cor.met.1hour.thre1<-corr.thre(metobolites.final.1hour1,0.25)


#cor.met.fasting.thre2<-corr.thre(metobolites.final.fasting1,0.5)

#cor.met.1hour.thre2<-corr.thre(metobolites.final.1hour1,0.5)
#length(cor.met.fasting.thre2)

#cutoff 0.4
cor.met.fasting.thre3<-corr.thre(metobolites.final.fasting1,0.4)

#cor.met.1hour.thre3<-corr.thre(metobolites.final.1hour1,0.4)
length(cor.met.fasting.thre3)

#cutoff 0.6
cor.met.fasting.thre4<-corr.thre(metobolites.final.fasting1,0.6)

#cor.met.1hour.thre4<-corr.thre(metobolites.final.1hour1,0.6)
#length(cor.met.fasting.thre4)


#cutoff 0.7
#cor.met.fasting.thre5<-corr.thre(metobolites.final.fasting1,0.7)

#cor.met.1hour.thre5<-corr.thre(metobolites.final.1hour1,0.7)
#length(cor.met.fasting.thre5)
#library(wActNet)

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
#out_dir<-"/Users/yubingyao/Google Drive/Network analysis/R code/"
module.fasting1.1_cos<-MultiModuleFind_cos_dat(network=da.igraph.test32, node.scores=node.scores.fasting32,edge.scores=edge.score.fasting,ncluster=4,clustersize=10)
#save(module.fasting1.1_cos,file=paste(out_dir,"fasting32_",as.character(10),"n_iter",as.character(4),".Rdata",sep=""))

library(igraph)
#load("fasting32_10n_iter4.Rdata")

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


#construct new optimized subnetwork from our wActNet algorithm 
#removing those edges with edge score less than or equal to 0
dat.network.sub<-dat.network[[1]]


dat.network.sub1<-delete_edges(dat.network.sub,which(E(dat.network.sub)$score<=0))


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

#met.fasting.fullname[edge.name.full[1,]]

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

#unlist(sapply(e.name.iso.niso,function(x){strsplit(x,"_")}))


dat.network.sub11<-delete_edges(dat.network.sub,which(E(dat.network.sub)$score<=0))


idx.vnode<-sapply(unlist(sapply(e.name.iso.niso,function(x){strsplit(x,"_")})),function(x){which(V(dat.network.sub11)$name==x)})

dat.network.sub12<-add_edges(dat.network.sub11,idx.vnode)

E(dat.network.sub12)$name<-c(E(dat.network.sub11)$name,e.name.iso.niso)

E(dat.network.sub12)$score<-c(E(dat.network.sub11)$score,e.score.iso.niso)
names(E(dat.network.sub12)$score)<-c(E(dat.network.sub11)$name,e.name.iso.niso)


save(dat.network.sub12, file = "network_fasting_w_plot.RData")

visNetwork(nodes, edges, height = "500px", width = "100%")

#getwd()
#setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")
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






library(visNetwork)

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

#V(lst.mod.fasting.n32_c06[[1]])$score

#E(lst.mod.fasting.n32_c04[[2]])$score


g_d_corr04_c1<-lst.mod.fasting.n32_c04[[1]]

g_d_corr04_c2<-lst.mod.fasting.n32_c04[[2]]


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





library(visNetwork)

g_d_corr06_c1<-lst.mod.fasting.n32_c06[[1]]

g_d_corr06_c2<-lst.mod.fasting.n32_c06[[2]]



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


dat.network.sub.score2[idx_iso_niso_max_edge[2]]

dat.network.sub.name2[idx_iso_niso_max_edge[2]]

dat.network.sub.name1[idx_iso_niso_max_edge[1]]



E(dat.network.sub)$score[grepl(v.iso_name[1], edge.name.full.new)]

E(dat.network.sub)$name[grepl(v.iso_name[1], edge.name.full.new)]


E(dat.network.sub)$name[idx_iso_niso_max_edge2]

E(dat.network.sub)$score[idx_iso_niso_max_edge2]




edge.name.full_mat[which(apply(idx.2iso,1,sum)==2),]



E(dat.network.sub)$name[which.max(E(dat.network.sub)$score[grepl(v.iso_name[4], edge.name.full.new)])]
  
E(dat.network.sub)$score[which.max(E(dat.network.sub)$score[grepl(v.iso_name[4], edge.name.full.new)])]

names(E(dat.network.sub)$score)<-edge.name.full.new

E(dat.network.sub)$name<-edge.name.full.new


edge.name.full_mat



edge.name.full.new[grepl(v.iso_name[1], edge.name.full.new, fixed = TRUE)]

V()



