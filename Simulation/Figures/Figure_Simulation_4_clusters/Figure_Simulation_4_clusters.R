

#illustrative example of 4 clusters in simulation setting,20 nodes,1 to 20-node name
#cluster 1-large node score, large edge score
library(visNetwork)
library(igraph)
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



