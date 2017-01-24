############################################################################
# Generate the networks for performance evaluation.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/evaluation/common.R")
############################################################################
source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")




############################################################################
# Creates a new graph corresponding to the specified parameters, or loads
# it from the file, if it exists. When created, the graph is also recorded
# as a Graphml file.
#
# n: number of nodes.
# type: randplanar (planar random graph) or erdosrenyi (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# folder: folder in which to record the graph.
#
# returns: the created or loaded graph.
############################################################################
init.graph <- function(n=5, type="randplanar", iteration=1, folder)
{	tlog(2,"Initializing the graph")
	
	# init file name
	graph.folder <- file.path(folder,type,paste0("n=",n),paste0("it=",iteration))
	graph.filename <- "disc=0"
	gf <- file.path(graph.folder,paste0(graph.filename,".graphml"))
	
	# graph already exists as a file
	if(file.exists(gf))
	{	tlog(4,"The graph already exists: loading file \"",gf,"\"")
		g <- read.graph(gf,format="graphml")
	}
	
	# graph must be created
	else
	{	tlog(4,"The graph is generated and recorded in file \"",gf,"\"")
		
		# create a planar graph
		if(type=="randplanar")
		{	g <- graph.empty(n=n, directed=FALSE)									# create empty graph
			V(g)$x <- runif(vcount(g),min=-1,max=1)									# setup the node spatial positions
			V(g)$y <- runif(vcount(g),min=-1,max=1)
			g <- connect.triangulation(g)											# use triangulation to init the links
		}
		
		# create an Erdös-Rényi graph
		else if(type=="erdosrenyi")
		{	g <- erdos.renyi.game(n=n,p.or.m=0.1,directd=FALSE)
			V(g)$x <- runif(vcount(g),min=-1,max=1)									# setup the node spatial positions
			V(g)$y <- runif(vcount(g),min=-1,max=1)
		}
		
		# common steps
		g <- distances.as.weights(g)											# add inter-node distances as link attributes
		V(g)$label <- 1:vcount(g)
		g$duration <- 0
		g$size <- vcount(g)
		g$granularity <- NA
		g$avgseg <- 0
		myplot.graph(g, node.str=NA, link.str=NA, large=TRUE, filename=graph.filename, out.folder=graph.folder, export=TRUE, formats=c("pdf",NA))
	}
	
	return(g)
}
