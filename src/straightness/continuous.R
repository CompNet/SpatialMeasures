############################################################################
# Processing the continuous average straightness, based on integrations.
# There are several variants, described in the article (see the readme file).
#
# Vincent Labatut 05/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/straightness/continuous.R")
############################################################################
library("igraph")


source("src/common/misc.R")
source("src/common/transformations.R")



############################################################################
# Auxiliary functions processing the lambda_u and lambda_v break-even distances 
# for the two specified links.
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# u2,v2: end nodes of the second link.
#
# returns: a vector containing lambda_u and lambda_v.
############################################################################
process.lambdauv <- function(e.dist, g.dist, u1, v1, u2, v2)
{	lambdau <- (g.dist[u2,v1] - g.dist[u2,u1] + e.dist[u1,v1])/2
	lambdav <- (g.dist[v2,v1] - g.dist[v2,u1] + e.dist[u1,v1])/2
	result <- c(lambdau,lambdav)
	return(result)
}


############################################################################
# Auxiliary functions processing the lambda_2 threshold for specified links
# and point.
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# u2,v2: end nodes of the second link.
# lambda_u, lambda_v: pre-processed break-even distances.
#
# returns: the lambda_2 break-even distance.
############################################################################
process.lambda2 <- function(e.dist, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
{	result <- NA
	if(ellp1<=lambdau && ell1<=lambdav)
		result <- (g.dist[u1,v2] - g.dist[u1,u2] + e.dist[u2,v2])/2
	else if(ellp1<=lambdau && ell1>lambdav)
		result <- (g.dist[v1,v2] - g.dist[u1,u2] + e.dist[u1,v1] + e.dist[u2,v2] - 2*ellp1)/2
	else if(ellp1>lambdau && ell1<=lambdav)
		result <- (g.dist[u1,v2] - g.dist[v1,u2] - e.dist[u1,v1] + e.dist[u2,v2] + 2*ellp1)/2
	else if(ellp1>lambdau && ell1>lambdav)
		result <- (g.dist[v1,v2] - g.dist[v1,u2] + e.dist[u2,v2])/2
	return(result)
}


############################################################################
# Processes the straightness between two points (not necessarily nodes).
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# u2,v2: end nodes of the second link.
# ellp2: position of p_2 on the second link.
# lambda_u, lambda_v: pre-processed break-even distances.
# lambda_2: pre-processed break-even distance.
# 
# returns: straightness value between p_1 and p_2. 
############################################################################
straightness.point.point <- function(e.dist, g.dist, u1, v1, ellp1, u2, v2, ellp2, lambdau, lambdav, lambda2)
{	# get the node coordinates
	xu1 <- V(g)$x[u1]
	yu1 <- V(g)$y[u1]
	xu2 <- V(g)$x[u2]
	yu2 <- V(g)$y[u2]
	xv1 <- V(g)$x[v1]
	yv1 <- V(g)$y[v1]
	xv2 <- V(g)$x[v2]
	yv2 <- V(g)$y[v2]
	
	# euclidean distance between the points
	edp1p2 <- sqrt((xu2 + ellp2/e.dist[u2,v2]*(xv2-xu2) - xu1 - ellp1/E.dist[u1,v1]*(xv1-xu1))^2 
				+ (yu2 + ellp2/e.dist[u2,v2]*(yv2-yu2) - yu1 - ellp1/E.dist[u1,v1]*(yv1-yu1))^2)
	
	# process the straightness
	result <- NA
	# if the points are disconnected
	if(is.infinite(g.dist[u1,u2]))
		result <- 0
	# if p1 and p2 are lying on the same link
	else if(setequal(c(u1,v1),c(u2,v2)))
		result <- 1
	# general case
	else
	{	# u1u2
		if(ellp1<=lambdau && ell1<=lambdav && ellp2<=lambda2
			|| ellp1<=lambdau && ell1>lambdav && ellp2<=lambda2)
			result <- edp1p2 / (ellp1 + g.dist[u1,u2] + ellp2)
		# u1v2
		else if(ellp1<=lambdau && ell1<=lambdav && ellp2>lambda2
			|| ellp1>lambdau && ell1<=lambdav && ellp2>lambda2)
			result <- edp1p2 / (ellp1 + g.dist[u1,v2] + e.dist[u2,v2] - ellp2)
		# v1u2
		else if(ellp1>lambdau && ell1<=lambdav && ellp2<=lambda2
			|| ellp1>lambdau && ell1>lambdav && ellp2<=lambda2)
			result <- edp1p2 / (e.dist[u1,v1] - ellp1 + g.dist[v1,u2] + ellp2)
		# v1v2
		else if(ellp1<=lambdau && ell1>lambdav && ellp2>lambda2
			|| ellp1>lambdau && ell1>lambdav && ellp2>lambda2)
			result <- edp1p2 / (e.dist[u1,v1] - ellp1 + g.dist[v1,v2] + e.dist[u2,v2] - ellp2)
	}
		
	return(result)
}


############################################################################
# Processes the total straightness between a point and a link (i.e. all the 
# points constituting this link).
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# u2,v2: end nodes of the second link.
# lambda_u, lambda_v: pre-processed break-even distances.
# 
# returns: total straightness between p_1 and (u_2,v_2). 
############################################################################
total.straightness.point.link <- function(e.dist, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
{	# process the appropriate lambda2
	lambda2 <- process.lambda2(e.dist, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
			
	# process the straightness
	result <- NA
	# if the points are disconnected
	if(is.infinite(g.dist[u1,u2]))
		result <- 0
	# if p1 and p2 are lying on the same link
	else if(setequal(c(u1,v1),c(u2,v2)))
		result <- e.dist[u2,v2]
	# general case
	else
	{	# define the function to integrate
		fun <- Vectorize(function(x)
		{	res <- straightness.point.point(e.dist, g.dist, u1, v1, ellp1, u2, v2, ellp2=x, lambdau, lambdav, lambda2)
			return(res)
		})
		# perform the approximate integration
		result <- integrate(f=fun,lower=0,upper=e.dist[u2,v2],abs.tol=1e-15)
	}
	
	return(result)
}


############################################################################
# Processes the average straightness between a point and a link (i.e. all the 
# points constituting this link).
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# u2,v2: end nodes of the second link.
# lambda_u, lambda_v: pre-processed break-even distances.
# 
# returns: average straightness between p_1 and (u_2,v_2). 
############################################################################
mean.straightness.point.link <- function(e.dist, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
{	# get the total straightness between the point and the link
	total <- total.straightness.point.link(e.dist, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
	# normalize to get the mean
	result <- total / e.dist[u2,v2]
	
	return(result)
}


############################################################################
# Processes the average straightness between each specified node and a link 
# (i.e. all the points constituting this link).
# 
# graph: considered graph.
# u: vector of nodes of interest.
# e: link of interest.
# 
# returns: vector containing the average straightness between each node of u
#          and the link e=(u_2,v_2). 
############################################################################
mean.straightness.nodes.link <- function(graph, u=V(graph), e)
{	# init the result vector
	result <- c()
	
	# init the distances
	pos <- cbind(V(graph)$x,V(graph)$y)
	e.dist <- as.matrix(dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2))
	g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
	
	# get the second link
	el <- get.edgelist(graph)
	u2 <- el[e,1]
	v2 <- el[e,2]
	
	# process each specified node
	for(i in u)
	{	# if the node is an isolate, no need to go further
		if(degree(graph,u)==0)
			str <- 0
		
		# otherwise, we process the mean straightness for the node
		else
		{	# get the relative position
			tmp <- neighbors(g,1)[1]
			if(tmp<u)
			{	u1 <- tmp
				v1 <- u
				ellp1 <- e.dist[u1,v1]
			}
			else
			{	u1 <- u
				v1 <- tmp
				ellp1 <- 0
			}
			
			# process the lambdas
			temp <- process.lambdauv(e.dist, g.dist, u1, v1, u2, v2)
			lambdau <- temp[1]
			lambdav <- temp[2]
			
			# get the total straightness between the point and the link
			str <- mean.straightness.point.link(e.dist, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
		}
		
		# add to the result vector
		result <- c(result,str)
	}
	
	return(result)
}


############################################################################
# Processes the average straightness between a point and the graph (i.e. all 
# the points constituting the graph).
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# 
# returns: average straightness between p_1 and all the points constituting the graph. 
############################################################################
mean.straightness.point.graph <- function(graph, e.dist, g.dist, u1, v1, ellp1)
{	total.str <- 0
	total.lgt <- 0
	
	# process each link
	el <- get.edgelist(graph)
	for(i in 1:nrow(el))
	{	# get the second link
		u2 <- el[i,1]
		v2 <- el[i,2]

		# process the lambdas
		temp <- process.lambdauv(e.dist, g.dist, u1, v1, u2, v2)
		lambdau <- temp[1]
		lambdav <- temp[2]
		
		# get the mean straightness between the point and the link
		total.str <- total.str + total.straightness.point.link(e.dist, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
		total.lgt <- total.lgt + e.dist[u2,v2]
	}
	
	# normalize to get the mean
	result <- total.str / total.lgt
	
	return(result)
}


############################################################################
# Processes the average straightness between each one of the specified nodes 
# and the graph (i.e. all the points constituting the graph). 
# 
# This can be interpreted as a centrality measure corresponding to the 
# average accessibility of the node from (or to) the rest of the graph.
# 
# graph: considered graph.
# u: vector of nodes of interest.
# 
# returns: vector containing the average straightness between each node of u
#          and the points constituting the graph.
############################################################################
mean.straightness.nodes.graph <- function(graph, u)
{	# init the result vector
	result <- c()
	
	# init the distances
	pos <- cbind(V(graph)$x,V(graph)$y)
	e.dist <- as.matrix(dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2))
	g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
	
	# process each specified node
	for(i in 1:length(u))
	{	# if the node is an isolate, no need to go further
		if(degree(graph,u[i])==0)
			str <- 0
		
		# otherwise, we process the mean straightness for the node
		else
		{	# get the relative position
			tmp <- neighbors(g,u[i])[1]
			if(tmp<u[i])
			{	u1 <- tmp
				v1 <- u[i]
				ellp1 <- e.dist[u1,v1]
			}
			else
			{	u1 <- u[i]
				v1 <- tmp
				ellp1 <- 0
			}
			
			# get the mean straightness between the point and the graph
			str <- mean.straightness.point.graph(e.dist, g.dist, u1, v1, ellp1)
		}
		
		# add to the result vector
		result <- c(result,str)
	}
	
	return(result)
}


############################################################################
# Processes the total straightness between two links (i.e. all the pairs of points 
# such that each point of the pair is located on a different link).
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# u2,v2: end nodes of the second link.
# lambda_u, lambda_v: pre-processed break-even distances.
# 
# returns: total straightness between (u_1,v_1) and (u_2,v_2). 
############################################################################
total.straightness.link.link <- function(e.dist, g.dist, u1, v1, u2, v2, lambdau, lambdav)
{	# process the straightness
	result <- NA
	# if the links are disconnected
	if(is.infinite(g.dist[u1,u2]))
		result <- 0
	# if the links are the same
	else if(setequal(c(u1,v1),c(u2,v2)))
		result <- e.dist[u1,v1]^2/2
	# general case
	else
	{	# define the function to integrate
		fun <- Vectorize(function(x)
				{	res <- total.straightness.point.link(e.dist, g.dist, u1, v1, ellp1=x, u2, v2, lambdau, lambdav)
					return(res)
				})
		
		# define the integral bounds
		b1 <- min(lambda_u,lmabda_v)
		b2 <- max(lambda_u,lmabda_v)
		
		# perform the approximate integration
		res1 <- integrate(f=fun,lower=0,upper=b1,abs.tol=1e-15)
		res2 <- integrate(f=fun,lower=b1,upper=b2,abs.tol=1e-15)
		res3 <- integrate(f=fun,lower=b2,upper=e.dist[u1,v1],abs.tol=1e-15)
		result <- res1 + res2 + res3
	}
	
	return(result)
}
	

############################################################################
# Processes the average straightness between two links (i.e. all the pairs of 
# points such that each point of the pair is located on a different link).
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# u2,v2: end nodes of the second link.
# lambda_u, lambda_v: pre-processed break-even distances.
# 
# returns: average straightness between (u_1,v_1) and (u_2,v_2). 
############################################################################
mean.straightness.link.link <- function(e.dist, g.dist, u1, v1, u2, v2, lambdau, lambdav)
{	result <- NA
	
	# if the links are the same
	if(setequal(c(u1,v1),c(u2,v2)))
		result <- 1
	# general case
	else
	{	# get the total straightness between the links
		total <- total.straightness.link.link(e.dist, g.dist, u1, v1, u2, v2, lambdau, lambdav)
		# normalize to get the mean
		result <- total / (e.dist[u2,v2]*e.dist[u1,v1])
	}
	
	return(result)
}


############################################################################
# Processes the average straightness between a link and the whole graph (i.e. 
# all the pairs of points such that one point is on the link and the other is
# on the graph).
# 
# e.dist: pre-processed Euclidean distances for all pairs of nodes in the graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the link of interest.
# u2,v2: end nodes of the second link.
# exclude.self: whether or not the straightness between the link and itself 
#				should be considered in the average.
# 
# returns: average straightness between (u_1,v_1) and the rest of the graph.
############################################################################
mean.straightness.link.graph <- function(e.dist, g.dist, u1, v1, exclude.self=FALSE)
{	el <- get.edgelist(graph)
	
	# process all links (possibly excluding (u_1,v_1)
	total.str <- 0
	total.length <- 0
	for(i in 1:nrow(el))
	{	# get the second link end-nodes
		u2 <- el[i,1]
		v2 <- el[i,2]
		
		# process the lambdas
		temp <- process.lambdauv(e.dist, g.dist, u1, v1, u2, v2)
		lambdau <- temp[1]
		lambdav <- temp[2]
		
		# check if the same than (u_1,v_1)
		if(!setequal(c(u1,v1),c(u2,v2)) || !exclude.self)
		{	# process the total straightness
			str <- total.straightness.link.link(e.dist, g.dist, u1, v1, u2, v2, lambdau, lambdav)
			# add to result
			total.str <- total.str + str
			total.length <- total.length + e.dist[u1,v1]*e.dist[u2,v2]
		}
	}
	
	# process the normalizing value
	if(exclude.self)
		denom <- total.length
	else
		denom <- total.length - e.dist[u1,v1]^2/2
	
	# process the average
	result <- total.str / denom
	return(result)
}
	
	
############################################################################
# Processes the average straightness between each one of the specified links 
# and the graph (i.e. all the points constituting the graph). 
# 
# This can be interpreted as a centrality measure corresponding to the 
# average accessibility of the link from (or to) the rest of the graph.
# 
# graph: considered graph.
# e: vector of nodes of interest.
# exclude.self: whether or not the straightness between one link and itself 
#				should be considered in the average.
# 
# returns: vector containing the average straightness between each node of u
#          and the points constituting the graph.
############################################################################
mean.straightness.links.graph <- function(graph, e, exclude.self=FALSE)
{	# init the result vector
	result <- c()
	
	# init the distances
	pos <- cbind(V(graph)$x,V(graph)$y)
	e.dist <- as.matrix(dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2))
	g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
	
	# process averages based on the previous totals
	for(i in e)
	{	# get the link end-nodes
		u1 <- el[i,1]
		v1 <- el[i,2]
		
		# get the total straightness between the point and the link
		str <- mean.straightness.link.graph(e.dist, g.dist, u1, v1, exclude.self)
		
		# add to the result vector
		result <- c(result,str)
	}
	
	return(result)
}


mean.straightness.graph <- function(graph, exclude.self=FALSE)
{	result <- NA
	el <- get.edgelist(graph)
	
	# init the distances
	pos <- cbind(V(graph)$x,V(graph)$y)
	e.dist <- as.matrix(dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2))
	g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
	
	# process all pairs of links
	total.str <- 0
	total.length <- 0
	for(i in 1:nrow(el))
	{	# get the first link end-nodes
		u1 <- el[i,1]
		v1 <- el[i,2]
		
		for(j in i:nrow(el))
		{	# check if the links are the same one
			if(i!=j || !exclude.self)
			{	# get the second link end-nodes
				u2 <- el[j,1]
				v2 <- el[j,2]
				
				# process the lambdas
				temp <- process.lambdauv(e.dist, g.dist, u1, v1, u2, v2)
				lambdau <- temp[1]
				lambdav <- temp[2]
				
				# process straightness
				str <- total.straightness.link.link(e.dist, g.dist, u1, v1, u2, v2, lambdau, lambdav)
				# add to result
				total.str <- total.str + str
				total.length <- total.length + e.dist[u1,v1]*e.dist[u2,v2]
			}
		}
	}
	
	# possibly substract superfluous values
	if(exclude.self)
		denom <- total.length
	else
	{	denom <- total.length
		temp <- sapply(1:nrow(el), function(i) e.dist[el[i,1],el[i,2]]^2/2)
		denom <- denom - sum(temp)
	}
	
	# process the average
	result <- total.str / denom
	return(result)
}
