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
