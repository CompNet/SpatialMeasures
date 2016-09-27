############################################################################
# Generates various types of spatial networks.
#
# Vincent Labatut 09/2015
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/common/generation.R")
############################################################################
library("splancs") # used to process polygon areas

#source("src/raw-generation.R")
source("src/common/transformations.R")


# constant used as a workaround to rounding problems
EPSILON <- 0.0001



############################################################################################
# Connects each node to its k closest neighbors.
# Each must have two attributes x and y representing
# its spatial position. It is possible to define a maximal
# limit on the distance of the considered neighbors.
#
# g: graph to process.
# k: number of closest neighbors to link to.
# dist.fun: spatial distance function (Euclidean distance by default).
# min.pow: power of Minkowski's distance (default: 2).
# max.dist: maximal distance for the neighborhood (default: NA=no-limit)
#
# returns: the modified graph.
############################################################################################
connect.closest <- function(g, k, dist.fun="euclidean", min.pow=2, max.dist=NA)
{ # process spatial distances
	m <- as.matrix(dist(x=cbind(V(g)$x,V(g)$y), method=dist.fun, p=min.pow))
#	print(m)
	
	# get closest neighbors as a k*n matrix
	tgt.mat <- sapply(1:vcount(g),function(i)
			{   sorted.row <- sort(m[i,])
				if(!is.na(max.dist))
					sorted.row <- sorted.row[sorted.row<=max.dist]
				sorted.neighs <- names(sorted.row)
				sorted.neighs <- as.integer(sorted.neighs[sorted.neighs!=i])
				if(length(sorted.neighs)<k)
					sorted.neighs <- c(sorted.neighs,rep(i,k-length(sorted.neighs)))
				else
					sorted.neighs <- sorted.neighs[1:k]
				return(sorted.neighs)
			}
	)
	
	# build link matrix
	src.mat <- t(matrix(ncol=k,nrow=vcount(g),1:vcount(g)))
	link.mat <- cbind(src.mat,tgt.mat)
	link.mat <- link.mat[,c(matrix(1:ncol(link.mat), nrow=2, byrow = TRUE))]
	link.vect <- as.vector(t(link.mat))
	
	# add links
	g <- add.edges(graph=g,edges=link.vect)
	g <- simplify(graph=g, remove.multiple=TRUE, remove.loops=TRUE)
	
	return(g)
}



############################################################################################
# Generates a spatial graph whose nodes have the specified positions.
#
# x: list of vectors of x positions. Each vector corresponds to a distinct row. 
# y: vector of y positions, each one associated to a specific row.
# reps: number of times the rows from x must be generated.
# gap: vertical space between the last and first rows of two consecutive repetitions.
# directed: whether the graph should be directed or not.
#
# returns: the created graph.
############################################################################################
generate.regular.graph <- function(x, y, reps, gap, directed=FALSE)
{	# compute node positions
	x.coords <- c()
	y.coords <- c()
	cum.gap <- 0
	for(r in 1:reps)
	{	for(i in 1:length(x))
		{	x.coords <- c(x.coords,x[[i]])
			y.coords <- c(y.coords, cum.gap + rep(y[i], length(x[[i]])))
		}
		cum.gap <- cum.gap + y[length(x)] + gap
	}
	
	if(length(x)>1)
	{	x.coords <- c(x.coords,x[[1]])
		y.coords <- c(y.coords, cum.gap + rep(y[1], length(x[[1]])))
	}
	
#	print(x.coords)
#	print(y.coords)
	
	# create the empty graph
	g <- graph.empty(n=length(x.coords),directed=directed)
	g$name <- "regular"
	V(g)$x <- x.coords
	V(g)$y <- y.coords
	
	
	return(g)
}



############################################################################################
# Generates a graph of squares.
#
# n: number of squares in one row/column.
# area: area covered by the network (by default: 1).
#
# returns: a graph made of n*n squares, with approximately the specified area.
############################################################################################
produce.square.graph <- function(n, area=1)
{	# init scale
	ov.side <- sqrt(area) 	# total width of the graph
	side <- ov.side/n		# side of the squares
	
	# generate graph
	g <- generate.regular.graph(
		x=list(seq(from=0,to=ov.side,by=side)),
		y=0, gap=side, reps=(n+1)
	)
	g <- connect.closest(g=g,k=4,max.dist=side+EPSILON)
	
	# set graph attributes
	g$area <- area
	g$side <- side
	g$name <- "squares"
	g$title <- paste("Squares n=",n,sep="") 
	
	return(g)
}


	
############################################################################################
# Generates a radioconcentric graph.
#
# r: number of rays.
# s: number of spires.
# area: area covered by the network (by default: 1).
#
# returns: a radioconcentric graph with approximately the specified area.
############################################################################################
produce.radiocentric.graph <- function(r, s, area=1)
{	# generate graph
	g <- generate.radiocentric.graph(ray.nbr=r, spire.nbr=s)
	
	# translate to (0,0)
	dx <- min(V(g)$x)
	V(g)$x <- V(g)$x - dx
	dy <- min(V(g)$y)
	V(g)$y <- V(g)$y - dy
	
	# normalize area
	env <- chull(V(g)$x,V(g)$y) #convex hull of the network nodes
	true.area <- areapl(cbind(V(g)$x[env],V(g)$y[env]))
	coef <- sqrt(area/true.area)
	coord <- cbind(V(g)$x,V(g)$y)
	coord <- coord %*% (coef*diag(2))
	V(g)$x <- coord[,1]
	V(g)$y <- coord[,2]
	
	# set graph attributes
	g$area <- area
	g$rays <- r
	g$spires <- s
	g$name <- "radioconcentric"
	g$title <- paste("Radioconcentric rays=",r," spires=",s,sep="") 
	
	return(g)
}

