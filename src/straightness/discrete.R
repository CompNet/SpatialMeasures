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
# If the nodes are not connected by some path, their straightness is zero.
#
# graph: the graph to process. Can be omitted if both distances are specified.
# v: ids of the nodes whose straightness must be returned. The values are
#    processed by considering all the nodes in the graph, though. If this 
# 	 parameter is NA, then we consider all the nodes (v=V(g)).
# e.dist: Euclidean distance between all pairs of nodes. Must be a dist object.
#		  If missing, the function will process only the necessary distances
#		  with respect to v. If the available memory allows it, processing all
#		  the distances once and for all can significantly speeds things up.
# g.dist: graph distance between all pairs of nodes, as a square matrix. If 
#		  this parameter is missing, it is processed on the fly as explained 
#		  for parameter e.dist.
# slow: which mode to use. TRUE to select the method that uses much less memory,
# 		but is much slower, or FALSE (the default) to use the method that needs
# 		much more memory, but is much faster.
# 
# returns: a matrix containing length(v) rows and vcount(graph) columns, whose 
#		   (i,j)th element represents the straightness between nodes i and j.
############################################################################
straightness.nodes <- function(graph, v=NA, e.dist, g.dist, slow=FALSE)
{	disp <- FALSE
	# for some reason, displaying text helps measuring the memory usage!
#	if(disp) for(i in 1:1000) cat("x")
	
	# check if the link lengths are present in the graph
	if(missing(g.dist))
	{	eatt <- list.edge.attributes(graph)
		if(!("dist" %in% eatt))
			stop("You need to either provide the graph distances (using parameter g.dist) or a graph possessing an edge attribute called \"dist\"")
	}
	
	# set up the number of nodes
	if(missing(graph))
		n <- ncol(g.dist)
	else
		n <- gorder(graph)
	
	# possibly set up v
	if(all(is.na(v)))
	{	if(missing(graph))
			v <- 1:attr(e.dist, "Size")
		else
			v <- V(graph)
	}
	
	# set up the distances for when we need to process all the nodes
	if(length(v)==n)
	{	# possibly process the Euclidean distances
		if(missing(e.dist))
		{	if(disp) tlog(2,"Processing the Euclidean distances")
			# since we process all nodes, we necessarily want a dist object
			e.dist <- process.euclidean.distance(graph, v)
		}
		# possibly convert dist to matrix (requires more memory)
		if(!slow)
			e.dist <- as.matrix(e.dist)
		
		# possibly process the geodesic distances
		if(missing(g.dist))
		{	if(disp) tlog(2,"Processing the graph distances")
			g.dist <- shortest.paths(graph=graph, v=v, to=V(graph), weights=E(graph)$dist)
		}
	}
	# otherwise, set up the distances for when we focus only on certain nodes
	else
	{	# possibly process the Euclidean distances
		if(missing(e.dist))
		{	if(disp) tlog(2,"Processing the Euclidean distances for ",length(v)," node(s) / ",n)
			if(slow)
			{	# if v contains more than half the nodes, it's worth using dist 
				if(length(v)>(n/2))
					e.dist <- process.euclidean.distance(graph, V(graph))
				# otherwise, the slower approach
				else
					e.dist <- process.euclidean.distance(graph, v)
			}
			else
				e.dist <- process.euclidean.distance(graph, V(graph))
		}
		# possibly convert dist to matrix
		if(class(e.dist)=="dist" && !slow)
			e.dist <- as.matrix(e.dist)[v,,drop=FALSE]
		
		# possibly process the geodesic distances
		if(missing(g.dist))
		{	if(disp) tlog(2,"Processing the graph distances")
			g.dist <- shortest.paths(graph=graph, v=v, to=V(graph), weights=E(graph)$dist)
		}
		else
			g.dist <- g.dist[v,,drop=FALSE]
	}
	
	# if the Euclidean distance is still a dist object (slower)
	if(class(e.dist)=="dist")
	{	if(disp) tlog(2,"Working on the dist object")
		res <- matrix(NA,nrow=nrow(g.dist),ncol=ncol(g.dist))
		for(i in 1:nrow(res))
		{	for(j in 1:ncol(res))
			{	if(j==v[i])
					res[i,j] <- 1
				else
				{	val <- get.dist(v[i],j,e.dist) / g.dist[i,j]
					if(is.na(val))
						res[i,j] <- 0
					else
						res[i,j] <- val
				}
			}
		}
	}
	# otherwise, if the Euclidean distance is now a matrix (faster) 
	else
	{	if(disp) tlog(2,"Working on the matrices")
		
		# process the ratio
		if(disp) tlog(2,"Processing the ratio")
		res <- e.dist / g.dist
		
		# straightness between a node and itself must be 1
		if(disp) tlog(2,"Correcting self straightness")
		for(i in 1:length(v))
		{	#if(disp) tlog("i=",i," v[i]=",v[i]," nrow=",nrow(res)," ncol=",ncol(res))
			res[i,v[i]] <- 1
		}
		
		# unconnected nodes get a 0 straightness
		if(disp) tlog(2,"Correcting NA values")
		res[is.na(res)] <- 0
	}
	
	# should not happen... but still happens (rounding problem)
	res[res<0] <- 0
	res[res>1] <- 1
	
	return(res)
}



############################################################################
# Processes the average straightness between the specified nodes and the rest
# of the graph, or between all pairs of nodes in the graph. Also returns the 
# corresponding standard deviations. If self=TRUE, the straigthness between one 
# node and itself is taken into account (otherwise, it is ignored).
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
# slow: which mode to use. TRUE to select the method that uses much less memory,
# 		but is much slower, or FALSE (the default) to use the method that needs
# 		much more memory, but is much faster.
# 
# returns: if v is NA, a vector of two values (average straightness and standard deviation),
#		   otherwise a matrix of length(v) rows and 2 columns (average straightnes and
#		   standard deviation for each node specified in v).
############################################################################
mean.straightness.nodes <- function(graph, v=NA, self=TRUE, e.dist, g.dist, slow=FALSE)
{	# global average
	if(all(is.na(v)))
	{	# process the straightness values
		strn <- straightness.nodes(graph, v, e.dist, g.dist, slow)
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
		strn <- straightness.nodes(graph, v, e.dist, g.dist, slow)
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
