############################################################################
# Various functions used by various scripts.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/common/misc.R")
############################################################################
library("igraph")


############################################################################
# Processes the total (spatial) length of all edges in the graph.
#
# Note : if the graph already has a "dist" edge attribute, it is used (not
# processed again, so it must be up to date).
#
# g: the graph to process.
# 
# returns: the total length summed over all edges.
############################################################################
total.length <- function(g)
{	# possibly processes the spatial distances
	ea <- list.edge.attributes(graph=g)
	if(!("dist" %in% ea))
		g <- distances.as.weights(g)
	
	# sum the distances associated to all links
	res <- sum(E(g)$dist)
	
	return(res)
}


############################################################################
# Checks if the three specified points are aligned.
#
# x1,y1: coordinates of the first point.
# x2,y2: coordinates of the second point.
# x3,y3: coordinates of the third point.
#
# returns: TRUE iff the points are aligned (colinear).
############################################################################
check.alignment <- function(x1, y1, x2, y2, x3, y3)
{	result <- FALSE
	
#cat("x1==x2: ",abs(x1-x2)<1e-10," x1==x3: ",abs(x1-x3)<1e-10,"\n",sep="")	

	# possibly vertical
	if(abs(x1-x2)<1e-10)
		result <- (abs(x1-x3)<1e-10)
	
	# not vertical
	else
	{	if(abs(x1-x3)>=1e-10)
		{	slope1 <- (y1-y2)/(x1-x2)
			slope2 <- (y1-y3)/(x1-x3)
			result <- abs(slope1-slope2)<1e-10
#cat("slope1: ",(y1-y2)/(x1-x2)," slope2: ",(y1-y3)/(x1-x3)," slope1==slope2: ",abs(slope1-slope2)<1e-10,"\n",sep="")
		}
	}
	
	return(result)
}


############################################################################
# Processes the Euclidean distance between the two specified nodes (not
# necessarily neighbors).
# 
# g: the considered graph.
# i: index of the first node.
# j: index of the second node.
#
# returns: Euclidean distance between the specified nodes.
############################################################################
dist.nodes <- function(g, i, j)
{	sqrt((V(g)$x[i]-V(g)$x[j])^2 + (V(g)$y[i]-V(g)$y[j])^2)	
}


############################################################################
# Process the mean and standard deviation of the internode distance for
# graph "g". Depending on the "connected" parameter, all pairs of nodes are 
# considered, even if not connected (FALSE), or only nodes connected by
# links (TRUE). In the latter case, the distances are supposed to be the
# weights of the links.
#
# g: the graph to process.
# connected: TRUE to consider only connected pairs of nodes, FALSE to consider
#			 all possible pairs of nodes.
#
# returns: a vector containing the mean distance and the associated standard deviation.
############################################################################
mean.distance <- function(g, connected=FALSE)
{	# process spatial distances
	if(connected)
	{	dists <- E(g)$dist	
	}	
	else
	{	pos <- cbind(V(g)$x,V(g)$y)
		m <- as.matrix(dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p = 2))
		dists <- m[upper.tri(m)]
	}
	
	# mean and stdev values
	result <- c(mean(dists),sd(dists))
	return(result)
}
