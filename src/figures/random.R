############################################################################
# Script used to produce some of the paper figures from the experiments section.
# This one focuses on a small random planar graph (Fig.7-left in the paper).
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/figures/exp1.R")
############################################################################
source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")
source("src/misc/generation.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")




########################################
# Load/generate the graphs
########################################
gl <- list()
# load the original graph
tlog("Load the original graph")
gl[[1]] <- read.graph("data/figures/randplan-original/graph.graphml",format="graphml")
# remove a few nodes from the graph (to highlight the consequence on the average straightness)
tlog("Remove nodes 16-18 from the graph")
gl[[2]] <- delete.vertices(gl[[1]],16:18)



########################################
# Process all the average straightness variants
########################################
for(i in 1:length(gl))
{	tlog("Process graph #",i)
	g <- gl[[i]]
	out.folder <- file.path("data","figures","randplan-original",paste0("small_",i))
	dir.create(path=out.folder, showWarnings=FALSE, recursive=TRUE)
	
	########################################
	# Process the link-graph and node-graph averages
	########################################
	# process several variants of the average straightness
	tlog(2,"Process the average straightness between each node and the rest of the graph")
	node.str <- mean.straightness.nodes.graph(graph=g)
	tlog(2,"Process the average straightness between each link and the rest of the graph")
	link.str <- mean.straightness.links.graph(graph=g)
	
	# plot them
	tlog(2,"Generate the corresponding plots")
	myplot.graph(g, node.str=node.str, link.str=NA, large=TRUE, filename="node-graph", out.folder=out.folder, export=FALSE, formats="pdf")
	myplot.graph(g, node.str=NA, link.str=link.str, large=TRUE, filename="link-graph", out.folder=out.folder, export=FALSE, formats="pdf")
	
	# record them as text files
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=node.str,file=file.path(out.folder,"node-graph.txt"),row.names=FALSE,col.names=FALSE)
	write.table(x=link.str,file=file.path(out.folder,"link-graph.txt"),row.names=FALSE,col.names=FALSE)
	
	
	
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
		{	#tlog(6,"Process link ",e2,"/",ecount(g))
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



########################################
# Check if the values are within the theoretical ranges (used when debugging)
########################################
for(i in 1:length(gl))
{	tlog("Processing graph #",i)
	out.folder <- file.path("data","figures","randplan-original",paste0("small_",i))
	for(d in c("node-graph","link-graph","node-link","link-link"))
	{	tlog(2,"Processing ",d)
		data <- as.matrix(read.table(file=file.path(out.folder,paste0(d,".txt")),header=FALSE))
		tlog(4,"           Infinite values: ",length(which(is.infinite(data))))
		tlog(4,"                NaN values: ",length(which(is.nan(data))))
		tlog(4,"                 NA values: ",length(which(is.na(data))))
		tlog(4,"     Other negative values: ",length(which(!is.infinite(data) & !is.nan(data) & !is.na(data) & data<0)))
		tlog(4,"Other values larger than 1: ",length(which(!is.infinite(data) & !is.nan(data) & !is.na(data) & data>1)))
	}
}



# this command allows to plot the network on screen (useful for visual verifications)
# myplot.graph(g, node.str=NA, link.str=NA, large=TRUE, filename=NA, out.folder=NA, export=FALSE, formats=NA)
