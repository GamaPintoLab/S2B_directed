subS2B_version2=function(seed_graph,index1,index2,meandist){
  betweencount=rep(0,igraph::gorder(seed_graph))
  seedmat1=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index1))
  seedmat2=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index2))
  sp1_out=igraph::distances(seed_graph,v=index1,to=igraph::V(seed_graph), mode = "out") # element sp1_out[i,j] indicates the shortest path length going out of seed i of index1 to node j of seed_graph
  sp2_out=igraph::distances(seed_graph,v=index2,to=igraph::V(seed_graph), mode = "out") # matrix similar to sp1_out but reffering to index 2 seeds
  sp1_out[sp1_out==Inf]=igraph::vcount(seed_graph) # if there is no shortest path, the value is set to inf
  sp2_out[sp2_out==Inf]=igraph::vcount(seed_graph)
  maxbc=0
  for (i in 1:length(index1)){
    for (j in 1:length(index2)){
      sumsp=sp1_out[i,]+sp2_out[j,] # sum of the shortests distances going out of seed i and out of seed j, to each node in seed_graph
      m1=min(sp1_out[i,][sp1_out[i,]>0]) # "shortest path" between seed i and candidates
      m2=min(sp2_out[j,][sp2_out[j,]>0]) # "shortest path" between seed j and candidates
      if (m1>0 & m1<meandist & m2>0 & m2<meandist){
        maxbc=maxbc+1
        nodelist=which(sumsp == (m1+m2)) # the seed_graph nodes that are present in a bidirected shortest path linking seed i to seed j as the converging node are added to the list
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