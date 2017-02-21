############################################################################
# Functions related to distance processing.
#
# Vincent Labatut 02/2017
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/misc/distances.R")
############################################################################
library("igraph")




############################################################################
# Processes the Euclidean distance between the specified nodes and the rest
# of the graph. If v is NA, then all the nodes are considered for v.
#
# Note that, paradoxically, processing all distances can be faster than
# focusing only on certain nodes. This is due to the fact that these cases
# are handled using different functions. Similarly, if v contains more than
# hald the nodes of the considered graph, the dist object returned when 
# processing all nodes will be smaller than the rectangular matrix returned 
# when focusing on the nodes from v.
#
# g: the concerned graphs.
# v: the source nodes (the function processes the distance between them
#	 and all the nodes in the graph, including themselves).
#
# returns: either a rectangular matrix (if v is neither NA nor V(g)) or a dist
#		   object if v is NA or V(g) (i.e. we want the distance between all 
#		   pairs of nodes in the graph).
############################################################################
process.euclidean.distance <- function(g, v=NA)
{	if(all(is.na(v)))
		v <- V(g)
	
	pos <- cbind(vertex_attr(g, name="x"),vertex_attr(g, name="y"))
		
	# all nodes at once
	if(length(v)==gorder(g))
		result <- dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2)
	
	# only certain nodes
	else
		result <- t(sapply(v, function(n) sqrt(apply(cbind(pos[,1]-pos[n,1],pos[,2]-pos[n,2])^2, 1, sum))))
	
	return(result)
}




############################################################################
# Returns the distance extracted from a dist object, for the specified nodes.
#
# u,v: the concerned nodes.
# d: the dist object.
#
# returns: the distance between u and v as represented by d.
############################################################################
get.dist <- function(u, v, d)
{	if(u==v)
		res <- 0
	else
	{	n <- attr(d, "Size")
		if(u>v)
			res <- d[n*(v-1) - v*(v-1)/2 + u-v]
		else
			res <- d[n*(u-1) - u*(u-1)/2 + v-u]
	}
	
	return(res)
}
