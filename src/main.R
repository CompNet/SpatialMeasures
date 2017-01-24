############################################################################
# This script shows how to process the discrete and continuous averages 
# of the Straightness measure.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/main.R")
############################################################################
library("igraph")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")




############################################################################
# Preparation
############################################################################
# first you need to provide a graph taking the form of an igraph object
# for instance, by loading it from a file. the graph must have two nodal
# attributes x and y describing the spatial position of the nodes
#g <- read.graph(my.graph,format="graphml")

# here we just generate it for the example
g <- erdos.renyi.game(n=100,p.or.m=0.4,type="gnp",directed=FALSE,loops=FALSE)
V(g)$x <- runif(n=gorder(g),0,1)
V(g)$y <- runif(n=gorder(g),0,1)

# plot the graph
#plot(g)




############################################################################
# Node-to-node Straightness
############################################################################
# process the vertex-to-vertex straightness for each pair of nodes
res <- straightness.nodes(graph=g)
cat("Node-to-node straightness for each pair of nodes in the graph:")
print(res)
# plot the histogram, just by curiosity...
#hist(res)

# discrete (traditional) average straightness between each node and the rest of the graph
res <- mean.straightness.nodes(graph=g, v=V(g))
cat("\n\nAverage node-to-node straightness between each node and the rest of the graph, and associated standard deviation:")
print(res)
hist(res)

# discrete (traditional) average straightness over the whole graph
res <- mean.straightness.nodes(graph=g, v=NA)
cat("\n\nAverage node-to-node straightness over the whole graph: ",res[1]," standard deviation: ",res[2],"\n",sep="")




############################################################################
# Point-to-point Straightness
############################################################################
mean.straightness.graph(graph=g)
mean.straightness.links.graph(graph=g, e=1:ecount(g))
mean.straightness.link.link(graph=g, e1, e2)
mean.straightness.nodes.graph(graph=g, u=1:vcount(g))
mean.straightness.nodes.link(graph=g, u=1:vcount(g), e)
