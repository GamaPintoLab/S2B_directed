subS2B_version1=function(seed_graph,index1,index2,meandist){
  betweencount=rep(0,igraph::gorder(seed_graph))
  seedmat1=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index1))
  seedmat2=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index2))
  sp1_out=igraph::distances(seed_graph,v=index1,to=igraph::V(seed_graph), mode = "out") # element sp1_out[i,j] indicates the shortest path length going out of seed i of index1 to node j of seed_graph
  sp1_in=igraph::distances(seed_graph,v=igraph::V(seed_graph),to=index1, mode = "out") # element sp1_in[i,j] indicates the shortest path length going out of node i of seed_graph to seed j of index1
  sp2_out=igraph::distances(seed_graph,v=index2,to=igraph::V(seed_graph), mode = "out") # matrix similar to sp1_out but reffering to index 2 seeds
  sp2_in=igraph::distances(seed_graph,v=igraph::V(seed_graph),to=index2, mode = "out") # matrix similar to sp1_in but reffering to index 2 seeds
  sp1_out[sp1_out==Inf]=igraph::vcount(seed_graph) # if there is no shortest path, the value is set to inf
  sp1_in[sp1_in==Inf]=igraph::vcount(seed_graph)
  sp2_out[sp2_out==Inf]=igraph::vcount(seed_graph)
  sp2_in[sp2_in==Inf]=igraph::vcount(seed_graph)
  sp_1to2=sp1_out[,index2] # matrix with shortests paths going out of index 1 seeds to index 2 seeds
  sp_2to1=sp2_out[,index1] # matrix with shortests paths going out of index 2 seeds to index 1 seeds
  maxbc=sum(sp_1to2>0 & sp_1to2<meandist) + sum(sp_2to1>0 & sp_2to1<meandist)
  for (i in 1:length(index1)){
    for (j in 1:length(index2)){
      m1=sp_1to2[i,j] # shortest path goinf out of seed i of index 1 to seed j of index 2
      m2=sp_2to1[j,i] # shortest path going out of seed j of index 2 to seed i of index 1
      if (m1<meandist | m2<meandist){
        if (m1<meandist){
          sumsp=sp1_out[i,]+sp2_in[,j] # sum of the shortests distances going out of seed i and in to seed j, to and from each node in seed_graph
          nodelist=which(sumsp==m1) # the seed_graph nodes that are present in a shortest path linking seed i to seed j are added to the list
          nodelist=nodelist[!nodelist %in% c(index1[i],index2[j])] # removing seeds from nodelist
          betweencount[nodelist]=betweencount[nodelist]+1
          seedmat1[nodelist,i]=1
          seedmat2[nodelist,j]=1
        }
        if (m2<meandist){
          sumsp=sp1_in[,i]+sp2_out[j,] # sum of the shortest distances going out of seed j and in to seed i, to and from each node in seed_graph
          nodelist=which(sumsp==m2) # the seed_graph nodes that are present in a shortest path linking seed j to seed i are added to the list
          nodelist=nodelist[c(-index1[i],-index2[j])] # removing seeds from nodelist
          betweencount[nodelist]=betweencount[nodelist]+1
          seedmat1[nodelist,i]=1
          seedmat2[nodelist,j]=1
        }
      }
    }
  }
  betweencount=betweencount/maxbc
  list(allcount=betweencount,smat1=seedmat1,smat2=seedmat2,maxS2B=maxbc)
}