############################################################################
# Processing the classic discrete average straightness, based on summations.
# There are several variants, described in the article (see the readme file).
#
# Vincent Labatut 05/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/straightness/discrete.R")
############################################################################
library("igraph")


source("src/common/misc.R")
source("src/common/transformations.R")


############################################################################
# Processes the straightness for each pair of nodes:
#  1) Get the spatial distance between each pair of nodes.
#  2) Get the spatial length of the shortest path between each pair of nodes.
#  3) Divide the first by the second.
# If the node are not connected by some path, their straightness is zero.
#
# graph: the graph to process.
# v: ids of the nodes whose straightness must be returned. The values are
#    processed by considering all the nodes in the graph, though.
# 
# returns: a matrix containing length(v) rows and vcount(g) columns, whose 
#		   (i,j) element represents the straightness between nodes i and j.
############################################################################
straightness.nodes <- function(graph, v=V(graph))
{	# process spatial distances
	pos <- cbind(V(graph)$x,V(graph)$y)
	m <- as.matrix(dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2))
	m <- m[v,]
	
	# process geodesic distances
#	graph <- distances.as.weights(graph)
	msp <- shortest.paths(graph=graph, v=v, to=V(graph), weights=E(graph)$dist)
	
	# process the ratio
	res <- m / msp
	
	# straightness between a node and itself must be 1
	for(i in 1:length(v))
		res[i,v[i]] <- 1
	
	# unconnected nodes get a 0 straightness
	res[is.na(res)] <- 0
	
	return(res)
}


############################################################################
# Processes the average straightness between the specified nodes and the rest
# of the graph, or between all pairs of nodes in the graph. Also returns the 
# corresponding standard deviations. If self=TRUE, the straigthness between one 
# nodeand itself is taken into account (otherwise, it is ignored).
# 
# g: graph to consider.
# v: if NA, then the average is processed over all pairs of nodes. If a numerical
#    vector, then the average straightness is processed for each node, over the rest
#	 of the graph nodes.
# self: whether or not to consider the straightness between a node and itself, 
#       when processing the average.
# 
# returns: if v is NA, a vector of two values (average straightness and standard deviation),
#		   otherwise a matrix of length(v) rows and 2 columns (average straightnes and
#		   standard deviation for each node specified in v).
############################################################################
mean.straightness.nodes <- function(graph, v=NA, self=TRUE)
{	# global average
	if(all(is.na(v)))
	{	# process the straightness values
		strn <- straightness.nodes(graph)
		strn <- strn[upper.tri(strn,diag=self)]
		# average them
		res1 <- mean(strn)
		res2 <- sd(strn)
		# setup result variable
		result <- c(res1,res2)
		names(result) <- c("AVG","STDEV")
	}
	
	# individual averages
	else
	{	# process the straightness values
		strn <- straightness.nodes(graph,v)
		# average them
		result <- sapply(1:length(v),function(i)
				{	vals <- strn[i,]
					if(!self)
						vals <- vals[-v[i]]
					res1 <- mean(vals)
					res2 <- sd(vals)
					res <- c(res1,res2)
					return(res)
				})
		result <- t(result)
		colnames(result) <- c("AVG","STDEV")
	}
	
	return(result)
}
