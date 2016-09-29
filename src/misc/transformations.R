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



############################################################################
# Processes the spatial distances between each pair of connected nodes of the 
# graph, and stores them as an edge numerical attribute named "dist", which 
# can be used as an edge weight by igraph for various calculations.
#
# g : the graph to process.
#
# returns: the modified graph.
############################################################################
distances.as.weights <- function(g)
{	# process spatial distances
	pos <- cbind(V(g)$x,V(g)$y)
	m <- as.matrix(dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p = 2))
	
	# select the values corresponding to the existing links
	el <- get.edgelist(graph=g, names=FALSE)
	weights <- apply(el, 1, function(nodes) m[nodes[1],nodes[2]])
	
	# set them in the graph
	g <- set.edge.attribute(graph=g, name="dist", value=weights)
	#print(E(g)$dist)
	#plot(g,edge.label=E(g)$dist)
	
	return(g)
}



############################################################################################
# From the existing graph, builds a new one by adding extra nodes on existing edges. These
# nodes are added so that the resulting edges are not longer than a given granularity value.
# Note that if the original graph does not have the edge attribute dist indicating the 
# spatial length position of each link, it is automatically added.
# 
# g: the original graph.
# granularity: approximate size of the edges in the returned graph (the smaller, the better).
#			   Zero means no change at all.
#
# returns: the same graph, with additional nodes.
############################################################################################
add.intermediate.nodes <- function(g, granularity)
{	if(granularity==0)
		g2 <- g
	else
	{	# process the spatial distances in the original graph
		g <- distances.as.weights(g)
		# create an empty graph (same nodes, no link)
		g2 <- delete.edges(graph=g,edges=E(g))
		V(g2)$type <- "original"
		
		# consider each original link independently
		edgelist <- get.edgelist(g)
		for(r in 1:nrow(edgelist))
		{	# get the nodes and their spatial distance
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
				delta.x <- (V(g)[n.to]$x - V(g)[n.from]$x) / k
				delta.y <- (V(g)[n.to]$y - V(g)[n.from]$y) / k
				#cat(sprintf("%.10f", d/granularity),"\n")		
				#cat("d/g=",d/granularity," rnd=",round(d/granularity)," int=",as.integer(d/granularity)," trc=",trunc(d/granularity)," k=",k," delta=",delta," delta.x=",delta.x," delta.y=",delta.y,"\n",sep="")
				
				# add the corresponding nodes and links
				prev.node <- n.from
				if(k>1)
				{	pos.x <- V(g)[n.from]$x
					pos.y <- V(g)[n.from]$y
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
	pos <- cbind(V(g)$x,V(g)$y)
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