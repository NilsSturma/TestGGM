library(igraph)

g = graph_from_literal(1--6--7--8--4, 2--6--7--8--3, 7--5)

E(g) # this gives all edges
length(E(g))
E(g)$weight <- c(1,1,1,1,1,1,1) 
E(g)$weight

V(g)
length(V(g))
V(g)$type <- c(1,2,2,2,1,1,1,1) #1=observed
V(g)$type

colors <- c("tomato", "gray50")
V(g)$color <- colors[V(g)$type]
V(g)$color
plot(g)
g


g[]  # matrix with weights (could use this to calculate a covariance matrix given a network)

degree(g, mode="all")  # get degree og all nodes (important to distinguish leaves)


# Compute intersection of two paths
# !!! look at the function shortestPath from gdistance !!!
path1 = shortest_paths(g, 
               from=V(g)[name==1], 
               to =V(g)[name==5],
               output="epath",
               weight=NA)
path2 = shortest_paths(g, 
                      from=V(g)[name==2], 
                      to=V(g)[name==3],
                      output="epath",       # "epath", "vpath" or "both"
                      weight=NA)              # no weights are used

intersection = intersect(unlist(path1$epath), unlist(path2$epath))
intersection  #=2
E(g)[2]




# # get weights
# E(g)[unlist(path_1_5$epath)]$weight


##############################################################################################
##############################################################################################
library(CombMSC) # subsets

m = sum(V(g)$type==1)
sub_sets = subsets(5,4,1:m)

Q = list()
not_Q = list()
res = NONE


for (j in 1:nrow(sub_sets)){
  p = sub_sets[j,1]
  others = sub_sets[j,2:4]
  count = 0
  # Find the pair that gives an empty set when the intersection of its component paths is taken
  for (i in 1:3){
    path1 = shortest_paths(g, from=V(g)[name==p], to=V(g)[name==others[i]], output="epath", weight=NA)
    path2 = shortest_paths(g, from=V(g)[name==others[-i][1]], to=V(g)[name==others[-i][2]], output="epath", weight=NA)  
    intersection = intersect(unlist(path1$epath), unlist(path2$epath))
    if (length(intersection)==0){
      res = c(p, others[i], others[-i])
      count = count + 1
    }
    if (count>1){
      not_Q = c(not_Q, list( sub_sets[j,1:4]))
      break
    }
  }
  Q = c(Q, list(res))
}

# convert to matrix
Q = matrix(unlist(Q), nrow = length(Q), byrow = TRUE)
# not_Q = matrix(unlist(not_Q), nrow = length(Q), byrow = TRUE)




