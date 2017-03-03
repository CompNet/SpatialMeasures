############################################################################
# Script used to produce some of the paper figures from the experiments section.
# This one focuses on the road networks.
#
# Vincent Labatut 01/2017
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/figures/urban.R")
############################################################################
source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")
source("src/misc/generation.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")




cities <- list(
	abidjan=2103,
	karlskrona=2091,
	soustons=3092,
	maastricht=3296,
	troisrivieres=7032,
	alicesprings=483,
	sfax=1565,
	avignon=18756,
	liverpool=7612,
	ljubljana=30766,
	lisbon=283,
	dakar=1006,
	hongkong=26602,
	beijin=49588,
	manhattan=10098,
##	"newyork",
	istanbul=62899,
	roma=52599,
	bordeaux=63565,
	stpetersburg=105809,
	tokyo=76391
)

data.folder <- "data"
urban.folder <- file.path(data.folder,"urban")




########################################
# Process all the average straightness variants
########################################
# we suppose the graphs were already extracted before (see src/misc/extraction.R)
for(c in 1:length(cities))
{	city <- names(cities)[c]
	node <- cities[[c]]
	tlog("Process city ",city)
	
	########################################
	# Load the (previously created) graph
	########################################
	tlog(2,"Load the ",city," graph")
	out.folder <- file.path(urban.folder,city)
	g <- read.graph(file.path(out.folder,"graph.graphml"),format="graphml")
	
	
	
#	########################################
#	# process the graph distance values, or load them if previously processed
#	########################################
#	dist.file <- file.path(out.folder,"g-dist-original.data")
#	if(file.exists(dist.file))
#	{	tlog(2,"Loading cached graph distances")
#		load(dist.file)
#		tlog(4,"Graph distances loaded")
#	}
#	else
#	{	tlog(2,"Processing graph distances")
#		g.dist <- shortest.paths(graph=g, weights=E(g)$dist)
#		tlog(4,"Recording graph distances")
#		save(g.dist,file=dist.file)
#	}
#	
#	
#	
#	########################################
#	# Process the node-graph discrete averages
#	########################################
#	out.file <- file.path(out.folder,"disc-node-graph.txt")
#	# load the cached average straightness values
#	if(file.exists(out.file))
#	{	tlog(2,"Load the discrete average straightness between each node and the rest of the graph")
#		disc.node.str <- as.matrix(read.table(file=out.file,header=FALSE))[,1]
#	}
#	# process the average straightness values
#	else
#	{	tlog(2,"Process the discrete average straightness between each node and the rest of the graph")
#		disc.node.str <- mean.straightness.nodes(graph=g, v=V(g))
#		# record them as text files
#		tlog(2,"Record the numerical results for later consultation")
#		write.table(x=disc.node.str,file=out.file,row.names=FALSE,col.names=FALSE)
#	}
#	
#	# plot them
#	tlog(2,"Generate the corresponding plots")
#	myplot.graph(g, node.str=disc.node.str, link.str=NA, large=FALSE, filename="disc-node-graph", out.folder=out.folder, export=FALSE, formats="pdf", autoscale=TRUE)
#	
#	
#	
#	########################################
#	# Process the node-graph continuous averages
#	########################################
#	out.file <- file.path(out.folder,"cont-node-graph.txt")
#	# load the cached average straightness values
#	if(file.exists(out.file))
#	{	tlog(2,"Load the continuous average straightness between each node and the rest of the graph")
#		node.str <- as.matrix(read.table(file=out.file,header=FALSE))
#	}
#	# process the average straightness values
#	else
#	{	# process the values
#		tlog(2,"Process the continuous average straightness between each node and the rest of the graph")
#		node.str <- mean.straightness.nodes.graph(graph=g, g.dist=g.dist)
#		# record them as text files
#		tlog(4,"Record the numerical results for later consultation")
#		write.table(x=node.str,file=out.file,row.names=FALSE,col.names=FALSE)
#	}
#	
#	# plot them
#	tlog(2,"Generate the corresponding plots")
#	myplot.graph(g, node.str=node.str, link.str=NA, large=FALSE, filename="cont-node-graph", out.folder=out.folder, export=FALSE, formats="pdf", autoscale=TRUE)
#	
#	
#	
#	########################################
#	# Process the link-graph continuous averages
#	########################################
#	out.file <- file.path(out.folder,"cont-link-graph.txt")
#	# load the cached average straightness values
#	if(file.exists(out.file))
#	{	tlog(2,"Load the continuous average straightness between each link and the rest of the graph")
#		link.str <- as.matrix(read.table(file=out.file,header=FALSE))
#	}
#	# process the average straightness values
#	else
#	{	# process the values
#		tlog(2,"Process the continuous average straightness between each link and the rest of the graph")
#		link.str <- mean.straightness.links.graph(graph=g, g.dist=g.dist)
#		# record them as text files
#		tlog(4,"Record the numerical results for later consultation")
#		write.table(x=link.str,file=out.file,row.names=FALSE,col.names=FALSE)
#	}
#	
#	# plot them
#	tlog(2,"Generate the corresponding plots")
#	myplot.graph(g, node.str=NA, link.str=link.str, large=FALSE, filename="cont-link-graph", out.folder=out.folder, export=FALSE, formats="pdf", autoscale=FALSE)
	
	
	
	########################################
	# Process the node-link and link-node continuous averages
	########################################
#	out.file <- file.path(out.folder,"cont-node-link.txt")
#	# load the cached average node-link straightness
#	if(file.exists(out.file))
#	{	tlog(2,"Load the continuous average straightness between each node and each link")
#		nl.str <- as.matrix(read.table(file=out.file,header=FALSE))
#	}
#	# process the average node-link straightness
#	else
#	{	# process the values
#		tlog(2,"Process continuous average straightness between each node and each link")
#		nl.str <- matrix(NA,nrow=ecount(g),ncol=vcount(g)) 
#		for(e in 1:ecount(g))
#		{	tlog(4,"Process link ",e,"/",ecount(g))
#			str <- mean.straightness.nodes.link(graph=g, u=1:vcount(g), e=e, g.dist=g.dist)
#			#print(str)
#			nl.str[e,] <- str
#		}
#		# record them as a text file
#		tlog(4,"Record the numerical results for later consultation")
#		write.table(x=nl.str,file=out.file,row.names=FALSE,col.names=FALSE)
#	}
#	
#	# plot the values for each node
#	tlog(2,"Plot them relatively to the nodes")
#	for(u in 1:vcount(g))
#	{	tlog(4,"Process node ",u,"/",vcount(g))
#		V(g)$marked <- 1:vcount(g)==u
#		myplot.graph(g, node.str=NA, link.str=nl.str[,u], large=FALSE, filename=paste0("cont-node=",u,"-link"), out.folder=out.folder, export=FALSE, formats="pdf", autoscale=FALSE)
#	}
#	V(g)$marked <- FALSE
#	
#	# plot the values for each link 
#	tlog(2,"Plot them relatively to the links")
#	el <- get.edgelist(g)
#	for(e in 1:ecount(g))
#	{	tlog(4,"Process link ",e,"/",ecount(g))
#		E(g)$marked <- 1:ecount(g)==e
#		myplot.graph(g, node.str=nl.str[e,], link.str=NA, large=FALSE, filename=paste0("cont-link=",e,"-node"), out.folder=out.folder, export=FALSE, formats="pdf", autoscale=TRUE)
#	}
#	E(g)$marked <- FALSE
	
	
	
	########################################
	# Process the link-link continuous averages
	########################################
#	out.file <- file.path(out.folder,"cont-link-link.txt")
#	# load the cached average link-link straightness
#	if(file.exists(out.file))
#	{	tlog(2,"Load the continuous average straightness between each link and the rest of the graph")
#		ll.str <- as.matrix(read.table(file=out.file,header=FALSE))
#	}
#	# process the average link-link straightness
#	else
#	{	# process the values
#		tlog(2,"Process continuous average straightness between each link and each other link")
#		ll.str <- matrix(NA,nrow=ecount(g),ncol=ecount(g))
#		for(e1 in 1:ecount(g))
#		{	tlog(4,"Process link ",e1,"/",ecount(g))
#			for(e2 in 1:ecount(g))
#			{	#tlog(6,"Process link ",e2,"/",ecount(g))
#				str <- mean.straightness.link.link(graph=g, e1=e1, e2=e2, g.dist=g.dist)
#				#print(str)
#				ll.str[e1,e2] <- str
#				ll.str[e2,e1] <- ll.str[e1,e2]
#			}
#		}
#		# record them as a text file
#		tlog(4,"Record the numerical results for later consultation")
#		write.table(x=ll.str,file=out.file,row.names=FALSE,col.names=FALSE)
#	}
#	
#	# plot the values for each link 
#	tlog(2,"Plot them relatively to the links")
#	el <- get.edgelist(g)
#	for(e in 1:ecount(g))
#	{	tlog(4,"Process link ",e,"/",ecount(g))
#		E(g)$marked <- 1:ecount(g)==e
#		myplot.graph(g, node.str=NA, link.str=ll.str[e,], large=FALSE, filename=paste0("cont-link=",e,"-link"), out.folder=out.folder, export=FALSE, formats="pdf", autoscale=FALSE)
#	}
#	E(g)$marked <- FALSE
#	
#	# check the values
#	if(length(which(is.infinite(link.str)))>0)
#		error("Straightness pb in link.str")
#	if(length(which(is.infinite(node.str)))>0)
#		error("Straightness pb in node.str")
#	if(length(which(is.infinite(nl.str)))>0)
#		error("Straightness pb in nl.str")
#	if(length(which(is.infinite(ll.str)))>0)
#		error("Straightness pb in ll.str")
	
	
	
	########################################
	# generate node centric plot
	########################################
	nl.str <- rep(NA,ecount(g)) 
	for(e in 1:ecount(g))
	{	tlog("Process link ",e,"/",ecount(g))
		nl.str[e] <- mean.straightness.nodes.link(graph=g, u=node, e=e)
	}
	V(g)$marked <- 1:vcount(g)==node
	myplot.graph(g, node.str=NA, link.str=nl.str, large=FALSE, filename=paste0("cont-node=",node,"-links"), out.folder=out.folder, export=FALSE, formats="pdf", autoscale=FALSE)
}



########################################
# Check if the values are within the theoretical ranges (used when debugging)
########################################
#for(city in names(cities))
#{	tlog("Processing graph type ",city)
#	out.folder <- file.path(urban.folder,city)
#	for(d in c("cont-node-graph","cont-link-graph","cont-node-link","cont-link-link"))
#	{	tlog(2,"Processing ",d)
#		data <- as.matrix(read.table(file=file.path(out.folder,paste0(d,".txt")),header=FALSE))
#		tlog(4,"           Infinite values: ",length(which(is.infinite(data))))
#		tlog(4,"                NaN values: ",length(which(is.nan(data))))
#		tlog(4,"                 NA values: ",length(which(is.na(data))))
#		tlog(4,"     Other negative values: ",length(which(!is.infinite(data) & !is.nan(data) & !is.na(data) & data<0)))
#		tlog(4,"Other values larger than 1: ",length(which(!is.infinite(data) & !is.nan(data) & !is.na(data) & data>1)))
#	}
#}




# this command allows to plot the network on screen (useful for visual verifications)
# myplot.graph(g, node.str=NA, link.str=NA, large=FALSE, filename=NA, out.folder=NA, export=FALSE, formats=NA)

##vertex.sizes=2
