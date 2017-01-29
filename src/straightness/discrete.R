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

source("src/straightness/continuous.R")




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
# e.dist: Euclidean distance between all pairs of nodes. Must be a dist object.
#		  If missing, it processed by the function (processing once before 
#		  calling this function can speed up the processing, if there are
#		  several calls).
# g.dist: graph distance between all pairs of nodes. Same remark than for
#		  parameter e.dist.
# 
# returns: a matrix containing length(v) rows and vcount(graph) columns, whose 
#		   (i,j) element represents the straightness between nodes i and j.
############################################################################
straightness.nodes <- function(graph, v=V(graph), e.dist, g.dist)
{	# process spatial distances
	if(missing(e.dist))
	{	pos <- cbind(V(graph)$x,V(graph)$y)
		e.dist <- dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2)
	}
	
	# process geodesic distances
#	graph <- distances.as.weights(graph)
	denominators <- shortest.paths(graph=graph, v=v, to=V(graph), weights=E(graph)$dist)
	
	# process the ratio
	res <- matrix(NA,nrow=nrow(denominators),ncol=ncol(denominators))
	for(i in 1:(nrow(res)-1))
	{	for(j in (i+1):ncol(res))
		{	res[i,j] <- get.dist(i,j,e.dist) / denominators[i,j]
			res[j,i] <- res[i,j]
		}
	}
	
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
# e.dist: Euclidean distance between all pairs of nodes. Must be a dist object.
#		  If missing, it processed by the function (processing once before 
#		  calling this function can speed up the processing, if there are
#		  several calls).
# g.dist: graph distance between all pairs of nodes. Same remark than for
#		  parameter e.dist.
# 
# returns: if v is NA, a vector of two values (average straightness and standard deviation),
#		   otherwise a matrix of length(v) rows and 2 columns (average straightnes and
#		   standard deviation for each node specified in v).
############################################################################
mean.straightness.nodes <- function(graph, v=NA, self=TRUE, e.dist, g.dist)
{	# global average
	if(all(is.na(v)))
	{	# process the straightness values
		strn <- straightness.nodes(graph, v=1:V(g), e.dist, g.dist)
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
		strn <- straightness.nodes(graph, v, e.dist, g.dist)
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
