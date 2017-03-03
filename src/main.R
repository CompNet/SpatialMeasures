############################################################################
# This script shows how to process the discrete and continuous averages 
# of the Straightness measure.
#
# Some plotting instructions are commented, you need to un-comment them to
# see the corresponding graphical representations of the graph.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/main.R")
############################################################################
library("igraph")

source("src/misc/plot.R")
source("src/misc/transformations.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")




############################################################################
# Preparation
############################################################################
# first you need to provide a graph taking the form of an igraph object
# for instance, by loading it from a file. the graph must have two nodal
# attributes x and y describing the spatial position of the nodes
#g <- read.graph("mygraph.graphml",format="graphml")


# here we just generate it for the example
n <- 10
g <- graph.empty(n, directed=FALSE)
V(g)$x <- runif(vcount(g),min=-1,max=1)
V(g)$y <- runif(vcount(g),min=-1,max=1)
g <- connect.triangulation(g)
g <- distances.as.weights(g)
g <- delete_edges(graph=g,edges=sample(E(g),n/2,FALSE))				


# plot the graph
#plot(g)




############################################################################
# Node-to-node Straightness
############################################################################
# process the vertex-to-vertex straightness for each pair of nodes
res <- straightness.nodes(graph=g)
#cat("Node-to-node straightness for each pair of nodes in the graph:");print(res)
#hist(res) # plot the histogram, just by curiosity...


# discrete (traditional) average straightness between each node and the rest of the graph
res <- mean.straightness.nodes(graph=g, v=V(g))
#cat("\n\nAverage node-to-node straightness between each node and the rest of the graph, and associated standard deviation:");print(res)
#hist(res)	# to plot their histogram
#myplot.graph(g, node.str=res[,1], link.str=NA, large=TRUE, filename=NA, out.folder=NA, export=FALSE, formats=NA) # represent straightness with colors


# discrete (traditional) average straightness over the whole graph
res <- mean.straightness.nodes(graph=g, v=NA)
cat("\n\nAverage node-to-node straightness over the whole graph: ",res[1]," standard deviation: ",res[2],"\n",sep="")




############################################################################
# Point-to-point Straightness
############################################################################
# continuous average straightness between some nodes and a link of interest
edge <- 1
res <- mean.straightness.nodes.link(graph=g, u=1:vcount(g), e=edge)
E(g)$marked <- 1:ecount(g)==edge
#myplot.graph(g, node.str=res, link.str=NA, large=TRUE, filename=NA, out.folder=NA, export=FALSE, formats=NA) # represent straightness with colors
E(g)$marked <- FALSE


# continuous average straightness between some nodes and the whole graph
res <- mean.straightness.nodes.graph(graph=g, u=1:vcount(g))
#myplot.graph(g, node.str=res, link.str=NA, large=TRUE, filename=NA, out.folder=NA, export=FALSE, formats=NA) # represent straightness with colors


# continuous average straightness between two links
res <- mean.straightness.link.link(graph=g, e1=1, e2=2)
cat("\n\nContinuous average Straightness between edges 1 and 2: ",res,"\n",sep="")


# continuous average straightness between some links and the graph
res <- mean.straightness.links.graph(graph=g, e=1:ecount(g))
#myplot.graph(g, node.str=NA, link.str=res, large=TRUE, filename=NA, out.folder=NA, export=FALSE, formats=NA) # represent straightness with colors


# continuous average straightness between all pairs of points in the graph
res <- mean.straightness.graph(graph=g)
cat("\n\nContinuous average Straightness over the whole graph: ",res,"\n",sep="")
