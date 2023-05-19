 
  library(visNetwork)
  library(igraph)
  #2 illustrative clusters at simulation scenario B
  #cluster 3-large node score, moderate edge score
  g11 <- barabasi.game(20,m=5,power=1.5,directed=F)
  V(g11)$name<-as.character(1:20)
  V(g11)$score<-runif(20,min=4,max=6)
  E(g11)$score<-runif(length(E(g11)),min=4,max=6)
  #cluster 2-moderate node score, large edge score
  g12 <- barabasi.game(20,m=5,power=1.5,directed=F)
  V(g12)$name<-as.character(31:50)
  V(g12)$score<-runif(20,min=0.1,max=2)
  E(g12)$score<-runif(length(E(g12)),min=4,max=6)
 
  
  EL  =get.edgelist(g11)
  EL1 = get.edgelist(g12)
  ELU = rbind(EL, EL1)
  ELU = ELU[!duplicated(ELU),]
  GU = graph_from_edgelist(ELU, directed=FALSE)
  V(GU)$score<-c(runif(20,min=4,max=6),)
  
  g13 <- barabasi.game(10,m=2,power=1.5,directed=F)
  V(g13)$name<-as.character(21:30)
  plot(g13)
  
  EL2 = get.edgelist(g13)
  
  ELU2 = rbind(ELU, EL2)
  ELU2 = ELU2[!duplicated(ELU2),]
  
  
  ELU2<- rbind(ELU2,as.character(c(12,23)))
  
  ELU2<- rbind(ELU2,as.character(c(8,23)))
  
  ELU2<-rbind(ELU2,as.character(c(3,23)))
  
  ELU2<-rbind(ELU2,as.character(c(18,24)))
  
  ELU2<-rbind(ELU2,as.character(c(17,27)))
  
  
  ELU2<-rbind(ELU2,as.character(c(22,43)))
  
  ELU2<-rbind(ELU2,as.character(c(28,43)))
  
  ELU2<-rbind(ELU2,as.character(c(23,43)))
  
  ELU2<-rbind(ELU2,as.character(c(28,34)))
  
  ELU2<-rbind(ELU2,as.character(c(27,34)))
  
  GU2 = graph_from_edgelist(ELU2, directed=FALSE)
  
  V(GU2)$name
  V(GU2)$score<-c(runif(20,min=4,max=6),runif(20,min=0.1,max=1),runif(10,min=-6,max=-4))
  names( V(GU2)$score)<-V(GU2)$name

  E(GU2)$score<-c(runif(85,min=0.1,max=1),runif(85,min=4,max=6),runif(27,min=-4,max=2))
  
  
  
  
   E(GU2)$dashes<-ifelse(E(GU2)$score>0,FALSE,FALSE)
  E(GU2)$width<-(E(GU2)$score+3)*0.9
  V(GU2)$color<-ifelse(V(GU2)$score>0,"#E7298A","grey")
  V(GU2)$size<-c((V(GU2)$score[1:20]+12)*1.8,(V(GU2)$score[21:50]+12)*1.2)
  
  E(GU2)$color<-ifelse(E(GU2)$score>0,"#66A61E","grey")
 #########################################################################################
  # setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")
  
 # save(GU2,file="network2_2clusters_plot.RData")
  
  
  
  #load("network_2clusters_plot.RData")
  
  data_gu2 <- toVisNetworkData(GU2)
  
  visNetwork(nodes = data_gu2$nodes, edges = data_gu2$edges, height = "600px")%>% visNodes(font=list(size=30))
