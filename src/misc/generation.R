############################################################################
# Generates various types of spatial networks.
#
# Vincent Labatut 09/2015
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/misc/generation.R")
############################################################################
library("splancs") # used to process polygon areas (needed when generating spiderweb graphs)


source("src/misc/transformations.R")




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
	g <- distances.as.weights(g)
	
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
	g <- distances.as.weights(g)
	
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
	g <- distances.as.weights(g)
	
	return(g)
}


	
############################################################################################
# Generates a graph of triangles.
#
# n: number of triangles in one row.
# area: area covered by the network (by default: 1).
#
# returns: a graph made of triangles, with approximately the specified area.
############################################################################################
produce.triangle.graph <- function(n, area=1)
{	# init scale
	m <- round(n*sqrt(3)/2+EPSILON)			# number of triangles on a column
	side <- sqrt((4*area)/(sqrt(3)*m*n))	# side of the equilateral triangles
	h <- sqrt(3)*side/2						# height of the equilateral triangle
	ov.side <- n*h	 						# total width of the graph
	
	# generate graph
	g <- generate.regular.graph(
			x=list(
					seq(from=0, to=ov.side, by=2*h),
					seq(from=h, to=ov.side, by=2*h)
			),
			y=c(0,side/2), gap=side/2, reps=m
	)
	g <- connect.closest(g=g,k=6,max.dist=side+EPSILON)
	
	# set graph attributes
	g$area <- area
	g$side <- side
	g$name <- "triangles"
	g$title <- paste("Triangles n=",n,sep="") 
	g <- distances.as.weights(g)
	
	return(g)
}



############################################################################################
# Generates a graph of hexagons.
#
# m: number of hexagons in one column.
# area: area covered by the network (by default: 1).
#
# returns: a graph made of hexagons, with approximately the specified area.
############################################################################################
produce.hexagon.graph <- function(m, area=1)
{	# init scale
	n <- trunc(2*(m*sqrt(3)-1/2)/3+EPSILON)					# number of hexagons on one row
	if((n %% 2) == 0)										# this number must be odd, for construction reasons
		n <- n + 1 											
	side <- sqrt(2*area / ((n*(m-1)+n%/%2+1)*3*sqrt(3)))	# side of the hexagons
	h <- sqrt(3)*side/2										# height of the "corner *right* triangle" of the hexagons
	ov.side <- (n*3/2+1/2)*side								# total width of the graph
	
	# generate graph
	g <- generate.regular.graph(
			x=list(
					sort(c(seq(from=0.5*side,to=ov.side-1.5*side,by=3*side),
									seq(from=1.5*side,to=ov.side-0.5*side,by=3*side)
							)),
					sort(c(seq(from=0,to=ov.side-2*side,by=3*side),
									seq(from=2*side,to=ov.side,by=3*side)
							))
			),
			y=c(0,h), gap=h, reps=m
	)
	g <- connect.closest(g=g,k=3,max.dist=side+EPSILON)
	
	# set graph attributes
	g$area <- area
	g$side <- side
	g$name <- "hexagons"
	g$title <- paste("Hexagons m=",m,sep="")
	g <- distances.as.weights(g)
	
	return(g)
}	



############################################################################################
# Generates a graph of octagons.
#
# m: number of octagons in one row/column.
# area: area covered by the network (by default: 1).
#
# returns: a graph made of octagons, with approximately the specified area.
############################################################################################
produce.octagon.graph <- function(n, area=1)
{	# init scale
	side <- sqrt(area / (n^2*2*(1+sqrt(2))+(n-1)^2))
	h <- side / sqrt(2)
	ov.side <- n*(side+2*h)
	
	# generate graph
	g <- generate.regular.graph(
			x=list(
					sort(c(seq(from=h,to=ov.side-h-side,by=2*h+side),
									seq(from=h+side,to=ov.side-h,by=2*h+side)
							)),
					seq(from=0,to=ov.side,by=2*h+side),
					seq(from=0,to=ov.side,by=2*h+side)
			),
			y=c(0,h,h+side), gap=h, reps=n)
	g <- connect.closest(g=g,k=4,max.dist=side+EPSILON)
	
	# set graph attributes
	g$area <- area
	g$side <- side
	g$name <- "octagons"
	g$title <- paste("Octagons n=",n,sep="") 
	g <- distances.as.weights(g)
	
	return(g)
}	



############################################################################################
# Generates a radio-concentric graph.
#
# ray.nbr: total number of rays in the graph.
# spire.nbr: number of spires connected to one ray.
# directed: whether the graph should be directed or not.
# empty: produces a graph without any link.
#
# returns: the created graph.
############################################################################################
generate.radioconcentric.graph <- function(ray.nbr, spire.nbr, directed=FALSE, empty=FALSE)
{	# init
	n <- spire.nbr * ray.nbr + 1
	x.coords <- rep(NA,n)
	y.coords <- rep(NA,n)
	links <- c()
	
	# set coordinates
	x.coords[1] <- 0
	y.coords[1] <- 0
	i <- 2
	for(r in 1:ray.nbr)
	{	angle <- pi/2 + (r-1) * (2*pi/ray.nbr)
		prev <- 1
		for(s in 1:spire.nbr)
		{	x.coords[i] <- cos(angle)*s/spire.nbr
			y.coords[i] <- sin(angle)*s/spire.nbr
			links <- c(links, c(prev,i))
			if(r>1)
				links <- c(links, c(i-spire.nbr,i))
			else
				links <- c(links, c(i+spire.nbr*(ray.nbr-1),i))
			prev <- i
			i <- i + 1
		}
	}
	
	# build graph
	g <- graph.empty(n=length(x.coords),directed=directed)
	g$name <- "radioconcentric"
	if(!empty)
		g <- add.edges(graph=g,edges=links)
	V(g)$x <- x.coords
	V(g)$y <- y.coords
	g <- distances.as.weights(g)
	
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
produce.radioconcentric.graph <- function(r, s, area=1)
{	# generate graph
	g <- generate.radioconcentric.graph(ray.nbr=r, spire.nbr=s)
	
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
	g <- distances.as.weights(g)
	
	return(g)
}



############################################################################################
# Generates an orb-web graph.
#
# ray.nbr: total number of rays in the graph.
# spire.nbr: approximate number of spires connected to one ray.
# directed: whether the graph should be directed or not.
# empty: produces a graph without any link.
#
# returns: the created graph.
############################################################################################
generate.orbweb.graph <- function(ray.nbr, spire.nbr, directed=FALSE, empty=FALSE)
{	# init
#	n <- spire.nbr * ray.nbr + 1
	x.coords <- c()
	y.coords <- c()
	links <- c()
	
	# set coordinates
	x.coords[1] <- 0
	y.coords[1] <- 0
	delta <- 1 / (ray.nbr*spire.nbr)
	#cat("delta:",delta,"\n")      
	d <- 0
	i <- 2
	for(s in 1:spire.nbr)
	{	for(r in 1:ray.nbr)
		{	angle <- pi/2 + (r-1) * (2*pi/ray.nbr)
			d <- d + delta
			#cat(i,":",d,"\n")      
			x.coords[i] <- cos(angle)*d
			y.coords[i] <- sin(angle)*d
			# add spiral links
			if(i>2)
				links <- c(links, c(i-1,i))
			# add ray links
			if(s>1)
				links <- c(links, c(i-ray.nbr,i))
			else
				links <- c(links, c(1,i))
			
			i <- i + 1
		}
	}
	
	# build graph
	g <- graph.empty(n=length(x.coords),directed=directed)
	g$name <- "orbweb"
	if(!empty)
		g <- add.edges(graph=g,edges=links)
	V(g)$x <- x.coords
	V(g)$y <- y.coords
	g <- distances.as.weights(g)
	
	return(g)
}



############################################################################################
# Generates an orbweb graph.
#
# r: number of rays.
# s: number of spires.
# area: area covered by the network (by default: 1).
#
# returns: an orbweb graph with approximately the specified area.
############################################################################################
produce.orbweb.graph <- function(r, s, area=1)
{	# generate graph
	g <- generate.orbweb.graph(ray.nbr=r, spire.nbr=s)
	
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
	g$name <- "orbweb"
	g$title <- paste("Orb-web rays=",r," spires=",s,sep="") 
	g <- distances.as.weights(g)
	
	return(g)
}

