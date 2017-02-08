############################################################################
# Functions related to distance processing.
#
# Vincent Labatut 02/2017
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/misc/distancess.R")
############################################################################
library("igraph")




############################################################################
# Processes the Euclidean distance between various sets of nodes.
#
# u,v: the concerned nodes.
# d: the dist object.
#
# returns: the distance between u and v as represented by d.
############################################################################
process.euclidean.distance <- function()
{
	
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
