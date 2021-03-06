############################################################################
# Script used to produce some of the paper figures from the experiments section.
# This one focuses on a selection of regular graphs, cf. Fig.5-7 in the paper.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/figures/regular.R")
############################################################################
source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")
source("src/misc/generation.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")




graph.types <- c(
	"hexagons",
	"octagons",
	"orbweb",
	"radioconcentric",
	"squares",
	"triangles"
#	"randplan-original",
#	"randplan-yifanhu",
#	"randplan-fruchtergold"
)

data.folder <- "data"
figures.folder <- file.path(data.folder,"figures")




########################################
# This function generates a random planar graph.
#
# gtype: type of layout to apply to the graph.
# returns: the produced graph, as an igraph object.
########################################
produce.random.graph <- function(gtype)
{	out.folder <- file.path(figures.folder,gtype)
	net.file <- file.path(out.folder,"graph.graphml")
	
	if(gtype=="randplan-original")
	{	g <- graph.empty(n=50, directed=FALSE)
		V(g)$x <- runif(vcount(g),min=-1,max=1)
		V(g)$y <- runif(vcount(g),min=-1,max=1)
		g <- connect.triangulation(g)
		g <- distances.as.weights(g)
		V(g)$label <- 1:vcount(g)
	}
	else
	{	net.file0 <- file.path(figures.folder,"randplan-original","graph.graphml")
		if(file.exists(net.file0))
			g <- read.graph(net.file0,format="graphml")
		else
		{	g <- produce.random.graph(gtype="randplan-original")
			write.graph(g,net.file0,format="graphml")
		}
		if(gtype=="randplan-yifanhu")
		{	lay <- layout_with_kk(g)
			V(g)$x <- lay[,1]
			V(g)$y <- lay[,2]
		}
		else if(gtype=="randplan-fruchtergold")
		{	lay <- layout_with_fr(g)
			V(g)$x <- lay[,1]
			V(g)$y <- lay[,2]
		}
	}
	
	return(g)
}



########################################
# Generates the graphs
# Note: this part of the script is separated from the rest so that it can be easily disabled, 
# if the graphs have already been generated. 
########################################
for(gtype in graph.types)
{	# possibly create the folder
	out.folder <- file.path(figures.folder,gtype)
	dir.create(path=out.folder, showWarnings=FALSE, recursive=TRUE)
	net.file <- file.path(out.folder,"graph.graphml")
	
	if(file.exists(net.file))
		g <- read.graph(net.file,format="graphml")
	else
	{	if(gtype=="hexagons")
		{	tlog("Generate a graph of hexagons")
			g <- produce.hexagon.graph(m=7,area=7)
			write.graph(g,net.file,format="graphml")
		}
		else if(gtype=="octagons")
		{	tlog("Generate a graph of octagons")
			g <- produce.octagon.graph(n=7,area=15)
			write.graph(g,net.file,format="graphml")
		}
		else if(gtype=="orbweb")
		{	tlog("Generate an orb-web graph")
			g <- produce.orbweb.graph(r=8,s=7,area=20)
			write.graph(g,net.file,format="graphml")
		}
		else if(gtype=="radioconcentric")
		{	tlog("Generate a radio-concentric graph")
			g <- produce.radioconcentric.graph(r=8,s=10,area=20)
			write.graph(g,net.file,format="graphml")
		}
		else if(gtype=="squares")
		{	tlog("Generate a graph of squares")
			g <- produce.square.graph(n=14,area=20)
			write.graph(g,net.file,format="graphml")
		}
		else if(gtype=="triangles")
		{	tlog("Generate a graph of triangles")
			g <- produce.triangle.graph(n=7,area=4)
			write.graph(g,net.file,format="graphml")
		}
		#
		else if(gtype=="randplan-original")
		{	tlog("Generate a random planar graph")
			g <- produce.random.graph(gtype)
			write.graph(g,net.file,format="graphml")
		}
		else if(gtype=="randplan-fruchtergold")
		{	tlog("Generate a random planar graph with a Fruchterman-Reingold layout")
			g <- produce.random.graph(gtype)
			write.graph(g,net.file,format="graphml")
		}
		else if(gtype=="randplan-yifanhu")
		{	tlog("Generate a random planar graph with a Yifan-Hu layout")
			g <- produce.random.graph(gtype)
			write.graph(g,net.file,format="graphml")
		}
	}
}



########################################
# Process all the average straightness variants
########################################
for(gtype in graph.types)
{	########################################
	# Load the (previously created) graph
	########################################
	tlog("Load the ",gtype," graph")
	out.folder <- file.path(figures.folder,gtype)
	g <- read.graph(file.path(out.folder,"graph.graphml"),format="graphml")
	
	
	
	########################################
	# Process the node-node discrete straightness (not an average)
	########################################
	# process the discrete straightness
	tlog(2,"Process the discrete straightness between each pair of nodes")
	nn.str <- straightness.nodes(graph=g, v=V(g))
	
	# plot the values for each node
	tlog(2,"Plot them relatively to the nodes")
	for(u in 1:vcount(g))
	{	tlog(4,"Process node ",u,"/",vcount(g))
		V(g)$marked <- 1:vcount(g)==u
		myplot.graph(g, node.str=nn.str[u,]	, link.str=NA, large=TRUE, filename=paste0("disc-node=",u,"-node"), out.folder=out.folder, export=FALSE, formats="pdf")
	}
	V(g)$marked <- FALSE
	
	# record them as text files
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=nn.str,file=file.path(out.folder,"disc-node-node.txt"),row.names=FALSE,col.names=FALSE)
	
	
	
	########################################
	# Process the node-graph discrete averages
	########################################
	# process the average straightness
	tlog(2,"Process the discrete average straightness between each node and the rest of the graph")
	disc.node.str <- mean.straightness.nodes(graph=g, v=V(g))[,1]
	
	# plot them
	tlog(2,"Generate the corresponding plots")
	myplot.graph(g, node.str=disc.node.str, link.str=NA, large=TRUE, filename="disc-node-graph", out.folder=out.folder, export=FALSE, formats="pdf")
	
	# record them as text files
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=disc.node.str,file=file.path(out.folder,"disc-node-graph.txt"),row.names=FALSE,col.names=FALSE)
	
	
	
	########################################
	# Process the link-graph and node-graph continuous averages
	########################################
	# process several variants of the average straightness
	tlog(2,"Process the continuous average straightness between each node and the rest of the graph")
	node.str <- mean.straightness.nodes.graph(graph=g)
	tlog(2,"Process the continuous average straightness between each link and the rest of the graph")
	link.str <- mean.straightness.links.graph(graph=g)
	
	# plot them
	tlog(2,"Generate the corresponding plots")
	myplot.graph(g, node.str=node.str, link.str=NA, large=TRUE, filename="cont-node-graph", out.folder=out.folder, export=FALSE, formats="pdf")
	myplot.graph(g, node.str=NA, link.str=link.str, large=TRUE, filename="cont-link-graph", out.folder=out.folder, export=FALSE, formats="pdf")
	rank.diff.barplot(disc.vals=disc.node.str, cont.vals=node.str, out.folder=out.folder, formats="pdf")
	
	# record them as text files
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=node.str,file=file.path(out.folder,"cont-node-graph.txt"),row.names=FALSE,col.names=FALSE)
	write.table(x=link.str,file=file.path(out.folder,"cont-link-graph.txt"),row.names=FALSE,col.names=FALSE)
	
	
	
	########################################
	# Process the node-link and link-node continuous averages
	########################################
	# process the average node-link straightness
	tlog(2,"Process continuous average straightness between each node and each link")
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
		myplot.graph(g, node.str=NA, link.str=nl.str[,u], large=TRUE, filename=paste0("cont-node=",u,"-link"), out.folder=out.folder, export=FALSE, formats="pdf")
	}
	V(g)$marked <- FALSE
	
	# plot the values for each link 
	tlog(2,"Plot them relatively to the links")
	el <- get.edgelist(g)
	for(e in 1:ecount(g))
	{	tlog(4,"Process link ",e,"/",ecount(g))
		E(g)$marked <- 1:ecount(g)==e
		myplot.graph(g, node.str=nl.str[e,], link.str=NA, large=TRUE, filename=paste0("cont-link=",e,"-node"), out.folder=out.folder, export=FALSE, formats="pdf")
	}
	E(g)$marked <- FALSE
	
	# record them as a text file
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=nl.str,file=file.path(out.folder,"cont-node-link.txt"),row.names=FALSE,col.names=FALSE)
	
	
	
	########################################
	# Process the link-link continuous averages
	########################################
	# process the average link-link straightness
	tlog(2,"Process continuous average straightness between each link and each other link")
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
		myplot.graph(g, node.str=NA, link.str=ll.str[e,], large=TRUE, filename=paste0("cont-link=",e,"-link"), out.folder=out.folder, export=FALSE, formats="pdf")
	}
	E(g)$marked <- FALSE
	
	# record them as a text file
	tlog(2,"Record the numerical results for later consultation")
	write.table(x=ll.str,file=file.path(out.folder,"cont-link-link.txt"),row.names=FALSE,col.names=FALSE)
	
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
for(gtype in graph.types)
{	tlog("Processing graph type ",gtype)
	out.folder <- file.path(figures.folder,gtype)
	for(d in c("cont-node-graph","cont-link-graph","cont-node-link","cont-link-link"))
	{	tlog(2,"Processing ",d)
		data <- as.matrix(read.table(file=file.path(out.folder,paste0(d,".txt")),header=FALSE))
		tlog(4,"           Infinite values: ",length(which(is.infinite(data))))
		tlog(4,"                NaN values: ",length(which(is.nan(data))))
		tlog(4,"                 NA values: ",length(which(is.na(data))))
		tlog(4,"     Other negative values: ",length(which(!is.infinite(data) & !is.nan(data) & !is.na(data) & data<0)))
		tlog(4,"Other values larger than 1: ",length(which(!is.infinite(data) & !is.nan(data) & !is.na(data) & data>1)))
	}
}



########################################
# Additional computation for my HDR
########################################
for(gtype in graph.types)
{	cat("Processing ",gtype,"\n",sep="")
	
	# load the (previously created) graph
	tlog("Load the ",gtype," graph")
	out.folder <- file.path(figures.folder, gtype)
	g <- read.graph(file.path(out.folder,"graph.graphml"),format="graphml")
	tlog(2,"Number of nodes: ",vcount(g))
	myplot.graph(g, large=FALSE, filename=paste0("neutral"), out.folder=out.folder, export=FALSE, formats="pdf", autoscale=FALSE)
	
	# discretize the edges
	gran <- max(E(g)$dist)/10
	g2 <- add.intermediate.nodes(g, granularity=gran, slow=FALSE)
	
	#########################
	# compute all node-to-node straightness values
	vals1 <- straightness.nodes(graph=g, v=V(g))
	vals2 <- straightness.nodes(graph=g2, v=V(g2))
	
	# plot these values as a histogram
	h1 <- hist(vals1, freq=FALSE, plot=FALSE, breaks=seq(0,1,1/40))
	h2 <- hist(vals2, freq=FALSE, plot=FALSE, breaks=seq(0,1,1/40))
	pdf(file.path(out.folder,"histo_straightness.pdf"))
		if(max(h1$density)>max(h2$density))
		{	plot(h1, col=rgb(0,0,1,1/4), xlab="Straightness", main=NA, xlim=c(0.4,1), freq=FALSE)
			plot(h2, col=rgb(1,0,0,1/4), xlim=c(0.4,1), add=TRUE, freq=FALSE)
		}else
		{	plot(h2, col=rgb(1,0,0,1/4), xlab="Straightness", main=NA, xlim=c(0.4,1), freq=FALSE)
			plot(h1, col=rgb(0,0,1,1/4), xlim=c(0.4,1), add=TRUE, freq=FALSE)
		}
		legend(x="topright", legend=c("Vertex-to-vertex","Point-to-point"), fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)))
	dev.off()
	
	# display average straightness
	tlog(2,"Average v2v straightness: ",mean(vals1))
	tlog(2,"Average p2p straightness: ",mean(vals2))

	#########################
	# compute all node-to-node betweenness values
	vals1 <- betweenness(graph=g, v=V(g), directed=FALSE, weights=E(g)$dist, normalized=TRUE)
	vals2 <- betweenness(graph=g2, v=V(g2), directed=FALSE, weights=E(g2)$dist, normalized=TRUE)
	tlog(2,"Betweenness: min=",min(c(vals1,vals2)), " max=",max(c(vals1,vals2)))	
	
	# plot these values as a histogram
	h1 <- hist(vals1, freq=FALSE, plot=FALSE, breaks=seq(0,0.4,0.4/40))
	h2 <- hist(vals2, freq=FALSE, plot=FALSE, breaks=seq(0,0.4,0.4/40))
	pdf(file.path(out.folder,"histo_betweenness.pdf"))
	if(max(h1$density)>max(h2$density))
	{	plot(h1, col=rgb(0,0,1,1/4), xlab="Straightness", main=NA, xlim=c(0,0.4), freq=FALSE)
		plot(h2, col=rgb(1,0,0,1/4), xlim=c(0,0.4), add=TRUE, freq=FALSE)
	}else
	{	plot(h2, col=rgb(1,0,0,1/4), xlab="Straightness", main=NA, xlim=c(0,0.4), freq=FALSE)
		plot(h1, col=rgb(0,0,1,1/4), xlim=c(0,0.4), add=TRUE, freq=FALSE)
	}
	legend(x="topleft", legend=c("Vertex-to-vertex","Point-to-point"), fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)))
	dev.off()
	
	# display average straightness
	tlog(2,"Average v2v betweenness: ",mean(vals1))
	tlog(2,"Average p2p betweenness: ",mean(vals2))
	
#	readline(prompt="Press [enter] to continue")
}



# this command allows to plot the network on screen (useful for visual verifications)
# myplot.graph(g, node.str=NA, link.str=NA, large=TRUE, filename=NA, out.folder=NA, export=FALSE, formats=NA)
