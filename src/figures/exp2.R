############################################################################
# Script used to produce some of the paper figures from the experiments section.
# This one focuses on a selection of lattices (regular graphs).
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/figures/exp2.R")
############################################################################
source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")
source("src/misc/generation.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")


graph.types <- c(
		"hexagons",
		"octogons"
#		"orbitele",
#		"radiocentric",
#		"squares",
#		"triangles"
)


########################################
# Generate the graphs
# Note: this part of the script is separated from the rest so that it can be easily disabled, 
# if the graphs have already been generated. 
########################################
for(gtype in graph.types)
{	# possibly create the folder
	out.folder <- file.path("data","figures",gtype)
	dir.create(path=out.folder, showWarnings=FALSE, recursive=TRUE)
	
	if(gtype=="hexagons")
	{	tlog("Generate a graph of hexagons")
		g <- produce.hexagon.graph(m=11,area=20)
	}
	else if(gtype=="octogons")
	{	tlog("Generate a graph of octogons")
		g <- produce.octogon.graph(n=9,area=30)
	}
	else if(gtype=="orbitele")
	{	tlog("Generate a spider graph")
		g <- produce.orbitele.graph(r=8,s=7,area=20)
	}
	else if(gtype=="radiocentric")
	{	tlog("Generate a radio-concentric graph")
		g <- produce.radiocentric.graph(r=8,s=10,area=20)
	}
	else if(gtype=="squares")
	{	tlog("Generate a graph of squares")
		g <- produce.square.graph(n=14,area=20)
	}
	else if(gtype=="triangles")
	{	tlog("Generate a graph of triangles")
		g <- produce.triangle.graph(n=14,area=20)
	}
	
	# record the graph as a graphml file
	write.graph(g,file.path(out.folder,"graph.graphml"),format="graphml")
}



########################################
# Process all the average straightness variants
########################################
for(gtype in graph.types)
{	########################################
	# Load the (previously created) graph
	########################################
	tlog("Load the ",gtype," graph")
	out.folder <- file.path("data","figures",gtype)
	g <- read.graph(file.path(out.folder,"graph.graphml"),format="graphml")
	
	########################################
	# Process the link-graph and node-graph averages
	########################################
	# process several variants of the average straightness
	tlog(2,"Process the average straightness between each link and the rest of the graph")
	link.str <- mean.straightness.links.graph(graph=g)
	tlog(2,"Process the average straightness between each node and the rest of the graph")
	node.str <- mean.straightness.nodes.graph(graph=g)
	
	# plot them
	tlog(2,"Generate the corresponding plots")
	myplot.graph(g, node.str=NA, link.str=link.str, large=TRUE, filename="link-graph", out.folder=out.folder, export=FALSE, formats="pdf")
	myplot.graph(g, node.str=node.str, link.str=NA, large=TRUE, filename="node-graph", out.folder=out.folder, export=FALSE, formats="pdf")
	
	# record them as text files
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=link.str,file=file.path(out.folder,"link-graph.txt"),row.names=FALSE,col.names=FALSE)
	write.table(x=node.str,file=file.path(out.folder,"node-graph.txt"),row.names=FALSE,col.names=FALSE)
	
	
	
	########################################
	# Process the node-link and link-node averages
	########################################
	# process the average node-link straightness
	tlog(2,"Process average straightness between each node and each link")
	nl.str <- matrix(NA,nrow=ecount(g),ncol=vcount(g)) 
	for(e in 1:ecount(g))
	{	tlog(4,"Process link ",e,"/",ecount(g))
		str <- mean.straightness.nodes.link(graph=g, u=1:vcount(g), e=e)
		#print(str)
		nl.str[e,] <- str
	}
	
	# plot the values for each node
	tlog(2,"Plot them relatively to the nodes")
	for(u in 1:vcount(g))
	{	tlog(4,"Process node ",u,"/",vcount(g))
		V(g)$marked <- 1:vcount(g)==u
		myplot.graph(g, node.str=NA, link.str=nl.str[,u], large=TRUE, filename=paste0("node=",u,"-link"), out.folder=out.folder, export=FALSE, formats="pdf")
	}
	V(g)$marked <- FALSE
	
	# plot the values for each link 
	tlog(2,"Plot them relatively to the links")
	el <- get.edgelist(g)
	for(e in 1:ecount(g))
	{	tlog(4,"Process link ",e,"/",ecount(g))
		E(g)$marked <- 1:ecount(g)==e
		myplot.graph(g, node.str=nl.str[e,], link.str=NA, large=TRUE, filename=paste0("link=",e,"-node"), out.folder=out.folder, export=FALSE, formats="pdf")
	}
	E(g)$marked <- FALSE
	
	# record them as a text file
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=nl.str,file=file.path(out.folder,"node-link.txt"),row.names=FALSE,col.names=FALSE)
	
	
	
	########################################
	# Process the link-link averages
	########################################
	# process the average link-link straightness
	tlog(2,"Process average straightness between each link and each other link")
	ll.str <- matrix(NA,nrow=ecount(g),ncol=ecount(g))
	for(e1 in 1:ecount(g))
	{	tlog(4,"Process link ",e1,"/",ecount(g))
		for(e2 in 1:ecount(g))
		{	tlog(6,"Process link ",e2,"/",ecount(g))
			str <- mean.straightness.link.link(graph=g, e1=e1, e2=e2)
			#print(str)
			ll.str[e1,e2] <- str
			ll.str[e2,e1] <- ll.str[e1,e2]
		}
	}
	
	# plot the values for each link 
	tlog(2,"Plot them relatively to the links")
	el <- get.edgelist(g)
	for(e in 1:ecount(g))
	{	tlog(4,"Process link ",e,"/",ecount(g))
		E(g)$marked <- 1:ecount(g)==e
		myplot.graph(g, node.str=NA, link.str=ll.str[e,], large=TRUE, filename=paste0("link=",e,"-link"), out.folder=out.folder, export=FALSE, formats="pdf")
	}
	E(g)$marked <- FALSE
	
	# record them as a text file
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=ll.str,file=file.path(out.folder,"link-link.txt"),row.names=FALSE,col.names=FALSE)

	# check the values
	if(length(which(is.infinite(link.str)))>0)
		error("Straightness pb in link.str")
	if(length(which(is.infinite(node.str)))>0)
		error("Straightness pb in node.str")
	if(length(which(is.infinite(nl.str)))>0)
		error("Straightness pb in nl.str")
	if(length(which(is.infinite(ll.str)))>0)
		error("Straightness pb in ll.str")
}
	


# myplot.graph(g, node.str=NA, link.str=NA, large=TRUE, filename=NA, out.folder=NA, export=FALSE, formats=NA)
