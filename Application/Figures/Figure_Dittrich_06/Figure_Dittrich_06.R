
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
