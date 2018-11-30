subS2B_version3=function(seed_graph,index1,index2,meandist){
  betweencount=rep(0,igraph::gorder(seed_graph))
  seedmat1=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index1))
  seedmat2=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index2))
  sp1_in=igraph::distances(seed_graph,v=igraph::V(seed_graph),to=index1, mode = "out") # element sp1_in[i,j] indicates the shortest path length going out of node i of seed_graph to seed j of index1
  sp2_in=igraph::distances(seed_graph,v=igraph::V(seed_graph),to=index2, mode = "out") # matrix similar to sp1_in but reffering to index 2 seeds
  sp1_in[sp1_in==Inf]=igraph::vcount(seed_graph) # if there is no shortest path, the value is set to inf
  sp2_in[sp2_in==Inf]=igraph::vcount(seed_graph)
  maxbc=0
  for (i in 1:length(index1)){
    for (j in 1:length(index2)){
      sumsp=sp1_in[,i]+sp2_in[,j] # sum of the shortests distances going in to seed i and in to seed j, from each node in seed_graph
      m1=min(sp1_in[,i][sp1_in[,i]>0]) # "shortest path" between candidates and seed i
      m2=min(sp2_in[,j][sp2_in[,j]>0]) # "shortest path" between candidates and seed j
      if (m1>0 & m1<meandist & m2>0 & m2<meandist){
        maxbc=maxbc+1
        nodelist=which(sumsp == (m1+m2)) # the seed_graph nodes that are present in a bidirected shortest path linking seed i to seed j as the diverging node are added to the list
        nodelist=nodelist[!nodelist %in% c(index1[i],index2[j])] # removing seeds from nodelist
        betweencount[nodelist]=betweencount[nodelist]+1
        seedmat1[nodelist,i]=1
        seedmat2[nodelist,j]=1
      }
    }
  }
  betweencount=betweencount/maxbc
  list(allcount=betweencount,smat1=seedmat1,smat2=seedmat2,maxS2B=maxbc)
}