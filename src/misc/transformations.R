############################################################################
# Functions applying various transformations to existing graphs.
#
# Vincent Labatut 09/2015
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/misc/transformations.R")
############################################################################
library("igraph")
library("geometry")	# used for the triangulation when generating random planar graphs

source("src/misc/distances.R")
source("src/misc/log.R")




############################################################################
# Processes the spatial distances between each pair of connected nodes of the 
# graph, and stores them as an edge numerical attribute named "dist", which 
# can be used as an edge weight by igraph for various calculations.
#
# g : the graph to process.
# slow: which mode to use. TRUE to select the method that uses much less
#		memory, but is much slower, or FALSE to use the method that uses
#		much more memory, but is much faster.
#
# returns: the modified graph.
############################################################################
distances.as.weights <- function(g, slow=FALSE)
{	# slow but memory-cheap (for large graphs)
	if(slow)
	{	# retrieve the list of links
		el <- get.edgelist(graph=g, names=FALSE)
		# process each one individually
		weights <- apply(el, 1, function(nodes)
		# old way: much slower
		#			sqrt((V(g)$x[nodes[1]] - V(g)$x[nodes[2]])^2
		#				+ (V(g)$y[nodes[1]] - V(g)$y[nodes[2]])^2))
					sqrt((vertex_attr(g, name="x", index=nodes[1]) - vertex_attr(g, name="x", index=nodes[2]))^2
						+ (vertex_attr(g, name="y", index=nodes[1]) - vertex_attr(g, name="y", index=nodes[2]))^2))
	}
	
	# fast but memory-expansive (for not so large graphs)
	else
	{	# process all the Euclidean distances at once
		pos <- cbind(vertex_attr(g, name="x"),vertex_attr(g, name="y"))
		e.dist <- dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2)
		
		# process the link weights
		el <- get.edgelist(graph=g, names=FALSE)
		weights <- apply(el, 1, function(nodes) get.dist(nodes[1],nodes[2],e.dist))
	}
	
	# update the graph
	g <- set.edge.attribute(graph=g, name="dist", value=weights)
	
	#print(E(g)$dist)
	#plot(g,edge.label=E(g)$dist) 
	
	return(g)
}



############################################################################################
# From the existing graph, builds a new one by adding extra nodes on existing edges. These
# nodes are added so that the resulting edges are not longer than a given granularity value.
# 
# Note that if the original graph does not have the edge attribute dist indicating the 
# spatial length position of each link, it is automatically added.
# 
# g: the original graph.
# granularity: approximate length of the edges in the returned graph (the smaller, the better).
#			   Zero means no change at all.
# slow: which mode to use. TRUE to select the method that uses much less memory,
# 		but is much slower, or FALSE (the default) to use the method that needs
# 		much more memory, but is much faster.
#
# returns: the same graph, with additional nodes.
############################################################################################
add.intermediate.nodes <- function(g, granularity, slow=FALSE)
{	disp <- TRUE
	
	if(granularity==0)
		g2 <- g
	else
	{	# possibly process the spatial distances in the original graph
		eatt <- list.edge.attributes(g)
		if("dist" %in% eatt)
			g <- distances.as.weights(g)
		
		# create an empty graph (same nodes, no link)
		g2 <- delete.edges(graph=g,edges=E(g))
		V(g2)$type <- "original"
		
		# create empty data structures
		new.pos.x <- c()
		new.pos.y <- c()
		new.links <- c()
		new.lengths <- c()
		next.node <- vcount(g2) + 1
		
		# consider each original link independently
		edgelist <- get.edgelist(g)
		dcst <- nrow(edgelist)%/%10
		for(r in 1:nrow(edgelist))
		{	if(disp && r%%dcst==0) tlog(2,sprintf("%.2f",r/nrow(edgelist)*100),"% done")
			
			# get the nodes and their spatial distance
			n.from <- edgelist[r,1]
			n.to <- edgelist[r,2]
			d <- E(g)[n.from %--% n.to]$dist
			#cat(n.from,"-(",d,")-",n.to,"\n",sep="")
			
			# approximate the desired granularity
			k <- trunc(d / granularity  + 0.000001) # workaround for some weird rounding problem
			# check if the link should be split or not
			if(k==0)
			{	# just add the existing link
				g2 <- add.edges(graph=g2, edges=c(n.from,n.to), attr=list(dist=d))
			}
			# split the link
			else
			{	delta <- d / k
				delta.x <- (vertex_attr(g, name="x", index=n.to) - vertex_attr(g, name="x", index=n.from)) / k
				delta.y <- (vertex_attr(g, name="y", index=n.to) - vertex_attr(g, name="y", index=n.from)) / k
				#cat(sprintf("%.10f", d/granularity),"\n")		
				#cat("d/g=",d/granularity," rnd=",round(d/granularity)," int=",as.integer(d/granularity)," trc=",trunc(d/granularity)," k=",k," delta=",delta," delta.x=",delta.x," delta.y=",delta.y,"\n",sep="")
				
				# add the corresponding nodes and links
				prev.node <- n.from
				if(slow)
				{	new.links <- c()
					new.lengths <- c()
				}
				if(k>1)
				{	pos.x <- vertex_attr(g, name="x", index=n.from)
					pos.y <- vertex_attr(g, name="y", index=n.from)
					if(slow)
					{	new.pos.x <- c()
						new.pos.y <- c()
					}
					for(i in 1:(k-1))
					{	# process the spatial position of the new node
						pos.x <- pos.x + delta.x
						new.pos.x <- c(new.pos.x, pos.x)
						pos.y <- pos.y + delta.y
						new.pos.y <- c(new.pos.y, pos.y)
						# set up the new link
						new.links <- c(new.links, prev.node, next.node)
						new.lengths <- c(new.lengths,delta)
						prev.node <- next.node
						next.node <- next.node + 1
					}
					# add the new nodes
					if(slow)
					{	new.node.nbr <- length(new.pos.x)
						g2 <- add.vertices(graph=g2, nv=new.node.nbr, attr=list(x=new.pos.x, y=new.pos.y, type=rep("extra",new.node.nbr)))
					}
				}
				# add the new links
				new.links <- c(new.links,prev.node,n.to)
				new.lengths <- c(new.lengths,delta)
				if(slow)
					g2 <- add.edges(graph=g2, edges=new.links, attr=list(dist=new.lengths))
			}
		}
		
		# add all nodes and links at once
		if(!slow)
		{	new.node.nbr <- length(new.pos.x)
			g2 <- add.vertices(graph=g2, nv=new.node.nbr, attr=list(x=new.pos.x, y=new.pos.y, type=rep("extra",new.node.nbr)))
			g2 <- add.edges(graph=g2, edges=new.links, attr=list(dist=new.lengths))
		}
	}
	
	return(g2)
}



############################################################################################
# Previous version of add.intermediate.nodes, kept for testing purposes.
############################################################################################
add.intermediate.nodes.old <- function(g, granularity)
{	disp <- TRUE
	
	if(granularity==0)
		g2 <- g
	else
	{	# possibly process the spatial distances in the original graph
		eatt <- list.edge.attributes(g)
		if("dist" %in% eatt)
			g <- distances.as.weights(g)
		
		# create an empty graph (same nodes, no link)
		g2 <- delete.edges(graph=g,edges=E(g))
		V(g2)$type <- "original"
		
		# consider each original link independently
		edgelist <- get.edgelist(g)
		dcst <- nrow(edgelist)%/%10
		for(r in 1:nrow(edgelist))
		{	if(disp && r%%dcst==0) tlog(2,sprintf("%.2f",r/nrow(edgelist)*100),"% done")
			
			# get the nodes and their spatial distance
			n.from <- edgelist[r,1]
			n.to <- edgelist[r,2]
			d <- E(g)[n.from %--% n.to]$dist
			#cat(n.from,"-(",d,")-",n.to,"\n",sep="")
			
			# approximate the desired granularity
			k <- trunc(d / granularity  + 0.000001) # workaround for some weird rounding problem
			# check if the link should be split or not
			if(k==0)
			{	# just add the existing link
				g2 <- add.edges(graph=g2, edges=c(n.from,n.to), attr=list(dist=d))
			}
			# split the link
			else
			{	delta <- d / k
				delta.x <- (vertex_attr(g, name="x", index=n.to) - vertex_attr(g, name="x", index=n.from)) / k
				delta.y <- (vertex_attr(g, name="y", index=n.to) - vertex_attr(g, name="y", index=n.from)) / k
				#cat(sprintf("%.10f", d/granularity),"\n")		
				#cat("d/g=",d/granularity," rnd=",round(d/granularity)," int=",as.integer(d/granularity)," trc=",trunc(d/granularity)," k=",k," delta=",delta," delta.x=",delta.x," delta.y=",delta.y,"\n",sep="")
				
				# add the corresponding nodes and links
				prev.node <- n.from
				if(k>1)
				{	pos.x <- vertex_attr(g, name="x", index=n.from)
					pos.y <- vertex_attr(g, name="y", index=n.from)
					for(i in 1:(k-1))
					{	# process the spatial position of the new node
						pos.x <- pos.x + delta.x
						pos.y <- pos.y + delta.y
						# create the new node
						g2 <- add.vertices(graph=g2, nv=1, attr=list(x=pos.x, y=pos.y, type="extra"))
						current.node <- vcount(g2)
						# connect it to the previous node
						g2 <- add.edges(graph=g2, edges=c(prev.node,current.node), attr=list(dist=delta))
						prev.node <- current.node
					}
				}
				# connect to the very last node
				g2 <- add.edges(graph=g2, edges=c(prev.node,n.to), attr=list(dist=delta))
			}
		}
	}
	
	return(g2)
}



############################################################################################
# Adds links in the graph using Delonay's triangulation, and possibly removing
# extra links in order to meet the specified average degree k.
#
# g: existing graph (only nodes, no links).
# k: approximate average degree of the returned graph,
#    or NA to keep the full triangulation result.
#
# returns: the modified graph.
############################################################################################
connect.triangulation <- function(g, k=NA)
{	# triangulation
	pos <- cbind(vertex_attr(g, name="x"),vertex_attr(g, name="y"))
	triangles <- delaunayn(p=pos)
	links <- rbind(triangles[,1:2],triangles[,2:3],triangles[,c(3,1)])
	
	# add links to graph
	links <- c(t(links))
	g <- add.edges(graph=g, edges=links)
	# remove multiple links
	g <- simplify(graph=g,remove.multiple=TRUE)
	
	# remove surnumerary links
	if(!is.na(k))
	{	m <- vcount(g)*k/2
		if(ecount(g)>m)
			g <- delete.edges(graph=g,edges=sample(x=1:ecount(g),size=(ecount(g)-m)))
	}
	
	return(g)
}

# test of the discretization function
#n <- 10
#g <- graph.empty(n=n, directed=FALSE)
#V(g)$x <- runif(vcount(g),min=-1,max=1)
#V(g)$y <- runif(vcount(g),min=-1,max=1)
#g <- connect.triangulation(g)
#g <- distances.as.weights(g)# add inter-node distances as link attributes
#V(g)$label <- 1:vcount(g)
#g$duration <- 0
#g$size <- vcount(g)
#g$granularity <- NA
#g$avgseg <- 0
#start.time <- Sys.time();g1 <- add.intermediate.nodes.old(g,0.01);end.time <- Sys.time();duration1 <- difftime(end.time,start.time,units="s");print(duration1)
#start.time <- Sys.time();g2 <- add.intermediate.nodes(g,0.01,slow=TRUE);end.time <- Sys.time();duration2 <- difftime(end.time,start.time,units="s");print(duration2)
#start.time <- Sys.time();g3 <- add.intermediate.nodes(g,0.01,slow=FALSE);end.time <- Sys.time();duration3 <- difftime(end.time,start.time,units="s");print(duration3)
