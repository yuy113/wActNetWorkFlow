#please comment the below code out if you already install the package-wActNet from github
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")
install_github("yuy113/wActNet")
############################################################3
library(wActNet)
library(igraph)
#R package BioNet from BioConductor 
library(BioNet)
#Construct the graph from some real dataset
#set up the directory for the code
#change to your corresponding directory

#setwd("/Users/yubingyao/Google Drive/Network analysis/R code/")

#out_dir<-"/Users/yubingyao/Google Drive/Network analysis/R code/"
node_size<-500
network <- erdos.renyi.game(node_size, 0.3)
################################################################################################
MultiModuleFind_cos<-function(network,node.scores,edge.scores,ncluster,clustersize){
  
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
                             num_iter=200, muCh=0.1, zToR=10, minsize=clustersize)
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
                               num_iter=200, muCh=0.1, zToR=10, minsize=clustersize)
    
    
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

####################################################################################################
#functions used for simulated clustered network
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

make_noconnection_cluster1_cluster2<-function(g,cluster1.id,cluster2.id){
  
  
  #make sure no connecting edges between two clusters
  #input 1: g-the igraph class object with names of nodes defined as "Varid" index of nodes
  #with edges names defined by "Varid1-Varid2" where id1,id2 is the indices of two connecting nodes
  #ionput 2: cluster1.id-the 1st cluter index list in vector form
  #ionput 3: cluster2.id-the 2nd cluter index list in vector form
  
  #output:the new igraph class object with no connection edges across the two clusters
  
  condition.cluster1.cluster2<-function(x,cluster1.id,cluster2.id){
    
    
    ( (length(unique(intersect(x,cluster1.id)))==1) & (length(unique(intersect(x,cluster2.id)))==1) ) 
    
    
  }
  #ind.cluster2<-apply(get.edgelist(g,names=F),1,function(x){ ((length(unique(intersect(x,cluster1.id)))==1) & (length(unique(intersect(x,cluster2.id)))==1)) })
  
  
  
  ind.cluster1<-apply(get.edgelist(g,names=F),1,condition.cluster1.cluster2,cluster1.id=cluster1.id,cluster2.id=cluster2.id)
  
  
  if (sum(ind.cluster1) == 0 ){
    return(g)
  }
  if (sum(ind.cluster1) > 0 ){
    
    new.edgelist.name<-get.edgelist(g,names=F)
    
    new.edgelist.name[ind.cluster1==1,]<-c(NA,NA)
    
    
    new.edgelist.name2<-na.omit(new.edgelist.name)
    
    new.edgelist.varname2<-t(apply(new.edgelist.name2,1,function(x){paste("Var",as.character(x),sep="")}))
    
    
    from.new.edge.list<-new.edgelist.varname2[,1]
    to.new.edge.list<-new.edgelist.varname2[,2]
    names.new.edge.list<-paste(from.new.edge.list,to.new.edge.list,sep="_")
    
    
    
    edge.dat<-data.frame(from= from.new.edge.list,to=to.new.edge.list,name= names.new.edge.list)
    
    nodes.dat<-data.frame(name=V(g)$name)
    
    interactome2 <- graph_from_data_frame(edge.dat, directed=F, vertices=nodes.dat)
    
    
    #interactome2 <- graph.edgelist(new.edgelist.varname2, T)
    
    
    #E(interactome2)$name<-paste(get.edgelist(interactome2,name=T)[,1],get.edgelist(interactome2,name=T)[,2],sep="_")
    
    
    return(interactome2)}
  
}






#subnetwork to induced subgraph from graph object considering both the vertices and edges,
#with remove.vertex=T,F; if True, then remove the vertices which is not in the edges, otherwise keep them
#vid-the names of the vertices in the graph
#eid-the names of the edges in the graph, with the format-Name(from_vertex)_Name(to_vertex)
#output:the subgraph containing only the vertices of vid and the edges


subnetwork.e<-function(graph,vid,eid,remove.vertex=F){
  
  vid<-unique(vid)
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
###################################################################################################################
id.node.module<-function(g){
  
  if(!is.null(V(g)$name)){
    
    as.numeric(unlist(strsplit(V(g)$name,"r"))[seq(2,2*length(V(g)$name),2)])
    
  }
}



#simulation scenario A2
#under ER random graph model
##############################################################################################################
#generate the random graph using the number of nodes in the network from real dataset
#p:the probability of connecting edges
#n.size-the size of positive scoring cluster,node scores in positive scoring clusters-uniform(0.1,1) or uniform(2,4)
#total 4 positive scoring clusters-within first one edge scores all positive-uniform(2,4) with high node score
#node scores in other nodes in the output random graph-uniform(-4,-2)
#edge scores in other edges in the output random graph-uniform(-4,1)
rand.posneg.network.new12<-function(network.size,p,n.size){
  library(igraph)
  library(BioNet)
  g <- erdos.renyi.game(network.size, p)
  V(g)$name<-paste("Var",as.character(1:network.size),sep="")
  
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  #make sure each node of two signal clusters is not isolated within the clusters 
  n.start.cluster2<-floor(vcount(g)/2)
  cluster2.id<-seq(n.start.cluster2,n.start.cluster2+n.size-1,1)
  cluster3.id<-seq((vcount(g)-3*n.size+1),(vcount(g)-2*n.size),1)
  
  # g1<-subnetwork.e(g,vid=cluster2.id,eid=E(g)$name)
  #  plotModule(g1)
  
  
  cluster1.id<-1:n.size
 
  cluster4.id<-(vcount(g)-n.size+1):vcount(g)

  
  
  
    #reorder all cluster ids
  
  id.clusters.tot<-c(cluster1.id,cluster2.id,cluster3.id,cluster4.id)
  
  id.clusters<-sample(id.clusters.tot,length(id.clusters.tot))
 
  cluster1.id<-id.clusters[1:n.size]
  cluster2.id<-id.clusters[(n.size+1):(2*n.size)]
  cluster3.id<-id.clusters[(2*n.size+1):(3*n.size)]
  cluster4.id<-id.clusters[(3*n.size+1):(4*n.size)]
  

   
  
  
  g<-make_full_connected_cluster(g,cluster1.id)
  g<-make_full_connected_cluster(g,cluster2.id)
  g<-make_full_connected_cluster(g,cluster3.id)
  
  g<-make_full_connected_cluster(g,cluster4.id)
  
  
  #make sure the two clusters don't have direct connections
  g<-make_noconnection_cluster1_cluster2(g,cluster1.id,cluster2.id)
  
  g<-make_noconnection_cluster1_cluster2(g,cluster2.id,cluster4.id)
  
  g<-make_noconnection_cluster1_cluster2(g,cluster3.id,cluster4.id)
  
  g<-make_noconnection_cluster1_cluster2(g,cluster2.id,cluster3.id)
  
  g<-make_noconnection_cluster1_cluster2(g,cluster1.id,cluster3.id)
  
  g<-make_noconnection_cluster1_cluster2(g,cluster1.id,cluster4.id)
  
  
  
  
  
  V(g)$score<-runif(vcount(g),min=-4,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  
  
  g.edge.id<-cbind(as.numeric(unlist(strsplit(get.edgelist(g,names=T)[,1],"Var"))[seq(2,2*dim(get.edgelist(g,names=T))[1],2)]),
                   
                   as.numeric(unlist(strsplit(get.edgelist(g,names=T)[,2],"Var"))[seq(2,2*dim(get.edgelist(g,names=T))[1],2)]))
  
  
  ind.cluster1<-unlist(apply(g.edge.id,1,function(x) return(sum(length(intersect(x,cluster1.id)))==2)))
  new.edgelist.cluster1<-get.edgelist(g,names=T)[ind.cluster1==1,]
  edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
  
  
  
  n.edge1=length(edge.name.cluster1)
  score.edge1<-runif(n.edge1,min=2,max=4)
  names(score.edge1)<-edge.name.cluster1
  score.node.cluster1<-runif(n.size,min=2,max=4)
  names(score.node.cluster1)<-V(g)$name[cluster1.id]
  V(g)$score[cluster1.id]<-score.node.cluster1
  E(g)$score[ind.cluster1==1]<-score.edge1
  ##################################################################################
  
  
  #moderate node score,moderate edge score cluster
  ind.cluster4<-unlist(apply(g.edge.id,1,function(x) return(sum(length(intersect(x,cluster4.id)))==2)))
  new.edgelist.cluster4<-get.edgelist(g,names=T)[ind.cluster4==1,]
  edge.name.cluster4<-paste(new.edgelist.cluster4[,1],new.edgelist.cluster4[,2],sep="_")
  
  n.edge4=length(edge.name.cluster4)
  score.edge4<-runif(n.edge4,min=0.1,max=2)
  names(score.edge4)<-edge.name.cluster4
  score.node.cluster4<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster4)<-V(g)$name[cluster4.id]
  V(g)$score[cluster4.id]<-score.node.cluster4
  E(g)$score[ind.cluster4==1]<-score.edge4
  
  
  
  
  #high node score,moderate edge score cluster
  ind.cluster3<-unlist(apply(g.edge.id,1,function(x) return(sum(length(intersect(x,cluster3.id)))==2)))
  new.edgelist.cluster3<-get.edgelist(g,names=T)[ind.cluster3==1,]
  edge.name.cluster3<-paste(new.edgelist.cluster3[,1],new.edgelist.cluster3[,2],sep="_")
  n.edge3=length(edge.name.cluster3)
  score.edge3<-runif(n.edge3,min=0.1,max=2)
  names(score.edge3)<-edge.name.cluster3
  score.node.cluster3<-runif(n.size,min=2,max=4)
  names(score.node.cluster3)<-V(g)$name[cluster3.id]
  V(g)$score[cluster3.id]<-score.node.cluster3
  E(g)$score[ind.cluster3==1]<-score.edge3
  
  
  
  
  ##moderate node score,high edge score cluster
  
  ind.cluster2<-unlist(apply(g.edge.id,1,function(x) return(sum(length(intersect(x,cluster2.id)))==2)))
  new.edgelist.cluster2<-get.edgelist(g,names=T)[ind.cluster2==1,]
  edge.name.cluster2<-paste(new.edgelist.cluster2[,1],new.edgelist.cluster2[,2],sep="_")
  n.edge2=length(edge.name.cluster2)
  score.edge2<-runif(n.edge2,min=2,max=4)
  names(score.edge2)<-edge.name.cluster2
  score.node.cluster2<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster2)<-V(g)$name[cluster2.id]
  V(g)$score[cluster2.id]<-score.node.cluster2
  E(g)$score[ind.cluster2==1]<-score.edge2
  
  list(g,id.clusters)
  
}
###########################################################################


#simulation scenario A2
#under ER random graph model
p=0.3
n.size=50
n.cluster=4
lst.nodenum.c<-list()
lst.id.clusters<-list()


for (i in 1:n.sim){
  
  lst_network_er<-rand.posneg.network.new12(network.size=500,n.size=n.size,p=p)
  network.test<-lst_network_er[[1]]
  id.cluster.true<-lst_network_er[[2]]
  
  node.scores<-V(network.test)$score
  names(node.scores)<-V(network.test)$name
  
  edge.scores<-E(network.test)$score
  names(edge.scores)<-E(network.test)$name
  
  
  
  
  
  
  #find all posible optimized subnetworks by Dittrich node only method and our method
  modules.cos<-MultiModuleFind_cos(network.test,node.scores,edge.scores,ncluster=n.cluster,clustersize=n.size)
  
  node_size<-length(unique(V(network.test)$score))
  id.cluster1<-id.cluster.true[1:n.size]
  id.cluster2<-id.cluster.true[(n.size+1):(2*n.size)]
  id.cluster3<-id.cluster.true[(2*n.size+1):(3*n.size)]
  id.cluster4<-id.cluster.true[(3*n.size+1):(4*n.size)]
  
  id.noise<-setdiff(1:node_size,id.cluster.true)
  
  
  lst.modules.c.nodelist<-modules.cos
  
  
  
  lst.nodenum.c[[i]]<-lapply(lst.modules.c.nodelist,function(x){c(length(x),length(intersect(x,id.cluster1)),length(intersect(x,id.cluster2)),
                                                                  length(intersect(x,id.cluster3)),length(intersect(x,id.cluster4)),
                                                                  length(intersect(x,id.noise)))})
  
  lst.id.clusters[[i]]<-id.cluster.true
  
}
result2.sim11<-list(lst.nodenum.c,lst.id.clusters)
save(result2.sim11,file=paste(out_dir,"ERnew_clustersize",as.character(n.size),"_p",as.character(p*10),".Rdata",sep=""))

#simulation scenario A2
#under ER random graph model
p=0.1
n.size=50
n.cluster=4
lst.nodenum.c<-list()
lst.id.clusters<-list()

g<-network

for (i in 1:n.sim){
  
  
  lst_network_er<-rand.posneg.network.new12(network.size=500,n.size=n.size,p=p)
  network.test<-lst_network_er[[1]]
  id.cluster.true<-lst_network_er[[2]]
  
  node.scores<-V(network.test)$score
  names(node.scores)<-V(network.test)$name
  
  edge.scores<-E(network.test)$score
  names(edge.scores)<-E(network.test)$name
  
  
  
  
  
  
  #find all posible optimized subnetworks by Dittrich node only method and our method
  modules.cos<-MultiModuleFind_cos(network.test,node.scores,edge.scores,ncluster=n.cluster,clustersize=n.size)
  
  node_size<-length(unique(V(network.test)$score))
  id.cluster1<-id.cluster.true[1:n.size]
  id.cluster2<-id.cluster.true[(n.size+1):(2*n.size)]
  id.cluster3<-id.cluster.true[(2*n.size+1):(3*n.size)]
  id.cluster4<-id.cluster.true[(3*n.size+1):(4*n.size)]
  
  id.noise<-setdiff(1:node_size,id.cluster.true)
  
  
  lst.modules.c.nodelist<-modules.cos
  
  
  
  lst.nodenum.c[[i]]<-lapply(lst.modules.c.nodelist,function(x){c(length(x),length(intersect(x,id.cluster1)),length(intersect(x,id.cluster2)),
                                                                  length(intersect(x,id.cluster3)),length(intersect(x,id.cluster4)),
                                                                  length(intersect(x,id.noise)))})
  
  lst.id.clusters[[i]]<-id.cluster.true
  
}
result2.sim12<-list(lst.nodenum.c,lst.id.clusters)
save(result2.sim12,file=paste(out_dir,"ERnew_clustersize",as.character(n.size),"_p",as.character(p*10),".Rdata",sep=""))



#simulation scenario A2
#under ER random graph model

p=0.3
n.size=20
n.cluster=4
lst.nodenum.c<-list()
lst.id.clusters<-list()

g<-network


n.sim=40
for (i in 1:n.sim){
  
  
  lst_network_er<-rand.posneg.network.new12(network.size=500,n.size=n.size,p=p)
  network.test<-lst_network_er[[1]]
  id.cluster.true<-lst_network_er[[2]]
  
  
  node.scores<-V(network.test)$score
  names(node.scores)<-V(network.test)$name
  
  edge.scores<-E(network.test)$score
  names(edge.scores)<-E(network.test)$name
  
  
  
  
  
  
  
  #find all posible optimized subnetworks by Dittrich node only method and our method
  modules.cos<-MultiModuleFind_cos(network.test,node.scores,edge.scores,ncluster=n.cluster,clustersize=n.size)
  
  node_size<-length(unique(V(network.test)$score))
  id.cluster1<-id.cluster.true[1:n.size]
  id.cluster2<-id.cluster.true[(n.size+1):(2*n.size)]
  id.cluster3<-id.cluster.true[(2*n.size+1):(3*n.size)]
  id.cluster4<-id.cluster.true[(3*n.size+1):(4*n.size)]
  
  id.noise<-setdiff(1:node_size,id.cluster.true)
  
  
  lst.modules.c.nodelist<-modules.cos
  
  

  lst.nodenum.c[[i]]<-lapply(lst.modules.c.nodelist,function(x){c(length(x),length(intersect(x,id.cluster1)),length(intersect(x,id.cluster2)),
                                                                  length(intersect(x,id.cluster3)),length(intersect(x,id.cluster4)),
                                                                  length(intersect(x,id.noise)))})
  
  lst.id.clusters[[i]]<-id.cluster.true
  
}
result2.sim21<-list(lst.nodenum.c,lst.id.clusters)
save(result2.sim21,file=paste(out_dir,"ERnew_clustersize",as.character(n.size),"_p",as.character(p*10),".Rdata",sep=""))


#simulation scenario A2
#under ER random graph model

p=0.1
n.size=20
n.cluster=4
lst.nodenum.c<-list()
lst.id.clusters<-list()

g<-network

for (i in 1:n.sim){
  
  
  lst_network_er<-rand.posneg.network.new12(network.size=500,n.size=n.size,p=p)
  network.test<-lst_network_er[[1]]
  id.cluster.true<-lst_network_er[[2]]
  
  node.scores<-V(network.test)$score
  names(node.scores)<-V(network.test)$name
  
  edge.scores<-E(network.test)$score
  names(edge.scores)<-E(network.test)$name
  
  
  
  
  
  #find all posible optimized subnetworks by Dittrich node only method and our method
  modules.cos<-MultiModuleFind_cos(network.test,node.scores,edge.scores,ncluster=n.cluster,clustersize=n.size)
  
  node_size<-length(unique(V(network.test)$score))
  id.cluster1<-id.cluster.true[1:n.size]
  id.cluster2<-id.cluster.true[(n.size+1):(2*n.size)]
  id.cluster3<-id.cluster.true[(2*n.size+1):(3*n.size)]
  id.cluster4<-id.cluster.true[(3*n.size+1):(4*n.size)]
  
  id.noise<-setdiff(1:node_size,id.cluster.true)
  
  
  lst.modules.c.nodelist<-modules.cos
  
  
  
  lst.nodenum.c[[i]]<-lapply(lst.modules.c.nodelist,function(x){c(length(x),length(intersect(x,id.cluster1)),length(intersect(x,id.cluster2)),
                                                                  length(intersect(x,id.cluster3)),length(intersect(x,id.cluster4)),
                                                                  length(intersect(x,id.noise)))})
  
  lst.id.clusters[[i]]<-id.cluster.true
  
}
result2.sim22<-list(lst.nodenum.c,lst.id.clusters)
save(result2.sim22,file=paste(out_dir,"ERnew_clustersize",as.character(n.size),"_p",as.character(p*10),".Rdata",sep=""))

