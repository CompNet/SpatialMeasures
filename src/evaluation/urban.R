############################################################################
# Generate the time- and memory-related plots from the article cited in the 
# project readme file. This script focuses on the real-world road networks.
#
# Vincent Labatut 01/2017
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/evaluation/urban.R")
############################################################################
source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")




data.folder <- "data"
urban.folder <- file.path(data.folder,"urban")
res.graph.filename <- "table-graph.txt"
res.nodes.filename <- "table-nodes.txt"
sample.size <- 25




############################################################################
# Column headers
############################################################################
ORIG_GDIST_TIME <- "Orig-Gdist-Time"
DISCRETIZATION_TIME <- "Discretization-Time"
DISC_EDIST_TIME <- "Disc-Edist-Time"
DISC_GDIST_TIME <- "Disc-Gdist-Time"
CONT_STRAIGHT <- "Cont-Straight"
CONT_PROC_TIME <- "Cont-Proc-Time"
CONT_MEM_USE <- "Cont-Mem-Use"
DISC_STRAIGHT <- "Disc-Straight"
DISC_PROC_TIME <- "Disc-Proc-Time"
DISC_MEM_USE <- "Disc-Mem-Use"
STRAIGHT_DIFFERENCE <- "Straight-Diff"
PROP_NODE_NBR <- "Node-Nbr"
PROP_LINK_NBR <- "Link-Nbr"
PROP_DENSITY <- "Density"
NODE_ID <- "Node-ID"
NETWORK_ID <- "Network-ID"




############################################################################
# Initializes the tables that will contain the results. If they already exist
# as files, they are loaded. Otherwise, they are initialized.
#
# g: the concerned graph.
# city.folder: folder of the currently processed city.
# light.process: whether the whole graph (FALSE) or only a sample of the nodes
#				 (TRUE) should be processed.
#
# returns: a list containing two tables: one for the whole graph, and one for
#		   the individual nodes.
############################################################################
init.result.tables <- function(g, city.folder, light.process)
{	tlog(2,"Initializing the stats table for city ",basename(city.folder))
	res.tables <- list()
	
	# whole graph table
	res.graph.file <- file.path(city.folder,res.graph.filename)
	if(file.exists(res.graph.file))
	{	tlog(4,"The graph table already exists >> we load it")
		res.tables$graph <- as.matrix(read.table(res.graph.file,header=TRUE,check.names=FALSE))
	}
	else
	{	tlog(4,"The graph table does no exist yet >> we create it")
		cnames <- c(ORIG_GDIST_TIME,DISCRETIZATION_TIME,DISC_EDIST_TIME,DISC_GDIST_TIME,
				CONT_STRAIGHT,CONT_PROC_TIME,CONT_MEM_USE,DISC_STRAIGHT,DISC_PROC_TIME,DISC_MEM_USE,STRAIGHT_DIFFERENCE,
				PROP_NODE_NBR,PROP_LINK_NBR,PROP_DENSITY)
		res.tables$graph <- matrix(NA,nrow=1,ncol=length(cnames))
		colnames(res.tables$graph) <- cnames
		res.tables$graph[1,PROP_NODE_NBR] <- gorder(g)
		res.tables$graph[1,PROP_LINK_NBR] <- gsize(g)
		res.tables$graph[1,PROP_DENSITY] <- graph.density(g)
		write.table(x=res.tables$graph,file=res.graph.file,row.names=FALSE,col.names=TRUE)
	}
	
	# individual nodes table
	res.nodes.file <- file.path(city.folder,res.nodes.filename)
	if(file.exists(res.nodes.file))
	{	tlog(4,"The nodes table already exists >> we load it")
		res.tables$nodes <- as.matrix(read.table(res.nodes.file,header=TRUE,check.names=FALSE))
	}
	else
	{	tlog(4,"The nodes table does no exist yet >> we create it")
		cnames <- c(CONT_STRAIGHT,CONT_PROC_TIME,CONT_MEM_USE,DISC_STRAIGHT,DISC_PROC_TIME,DISC_MEM_USE,STRAIGHT_DIFFERENCE)
		if(light.process)
		{	cnames <- c(NODE_ID,cnames)
			res.tables$nodes <- matrix(NA,nrow=sample.size,ncol=length(cnames))
			colnames(res.tables$nodes) <- cnames
			res.tables$nodes[,NODE_ID] <- sample(gorder(g),sample.size)
		}
		else
		{	res.tables$nodes <- matrix(NA,nrow=gorder(g),ncol=length(cnames))
			colnames(res.tables$nodes) <- cnames
		}
		write.table(x=res.tables$nodes,file=res.nodes.file,row.names=FALSE,col.names=TRUE)
	}
	
	return(res.tables)
}




############################################################################
# Process both the Euclidean and graph distance between all pairs of nodes
# in the specified graph, or load them if they have already been processed
# before (and cached).
#
# g: concerned graph.
# city.folder: folder of the currently processed city.
# res.tables: result tables (possibly contain the distance processing times).
# light.process: whether the whole graph (FALSE) or only a sample of the nodes
#				 (TRUE) should be processed.
# 
# returns: an updated graph table.
############################################################################
init.distances <- function(g, city.folder, res.tables, light.process)
{	if(!light.process)
	{	tlog(2,"Initializing both distances for city ",basename(city.folder))
		
		# process the graph distances if needed
		g.file <- file.path(city.folder,"g-dist-original.data")
		if(!file.exists(g.file))
		{	tlog(4,"Processing graph distances")
			start.time <- Sys.time()
				g.dist <- shortest.paths(graph=g, weights=E(g)$dist)
			end.time <- Sys.time()
			g.duration <- difftime(end.time,start.time,units="s")
			tlog(6,"Recording graph distances (",g.duration,"s)")
			save(g.dist,file=g.file)
			# update table
			res.tables$graph[1,ORIG_GDIST_TIME] <- g.duration
			table.file <- file.path(city.folder,res.graph.filename)
			write.table(x=res.tables$graph,file=table.file,row.names=FALSE,col.names=TRUE)
		}
	}
	
	return(res.tables)
}




############################################################################
# Discretizes the graph and processes the corresponding Euclidean and graph
# distances between all pairs of nodes. Also caches the resulting data.
#
# g: concerned graph.
# city.folder: folder of the currently processed city.
# res.tables: result tables (possibly contain the distance processing times).
# light.process: whether the whole graph (FALSE) or only a sample of the nodes
#				 (TRUE) should be processed.
# avgseg: granularity of the discretization threshold (only one is used here, by 
#		opposition to the random networks, for which we sequentially used several 
#		values).
#
# returns: an updated graph table.
############################################################################
init.discretization <- function(g, city.folder, res.tables, light.process, avgseg=50)
{	tlog(2,"Discretizing the graph links for city ",basename(city.folder))
	
	# process the discretized graph
	graph.file <- file.path(city.folder,"graph_discretized.graphml")
	if(file.exists(graph.file))
	{	tlog(4,"Loading the discretized graph")
		g2 <- read.graph(graph.file, format="graphml")
	}
	else
	{	# discretize the graph
		gran <- mean(E(g)$dist)/avgseg
		tlog(4,"Discretizing the graph once and for all (gran=",gran,")")
		start.time <- Sys.time()
			g2 <- add.intermediate.nodes(g, granularity=gran)
		end.time <- Sys.time()
		duration <- difftime(end.time,start.time,units="s")
		write.graph(graph=g2,file=graph.file,format="graphml")
		# record the result table
		res.tables$graph[1,DISCRETIZATION_TIME] <- duration
		table.file <- file.path(city.folder,res.graph.filename)
		write.table(x=res.tables$graph,file=table.file,row.names=FALSE,col.names=TRUE)
		# display for debug
		tlog(6,"Number of nodes: ",gorder(g2)," - Duration: ",duration," s - Average segmentation: ",avgseg," - Recorded as \"",graph.file,"\"")
	}
	
	tlog(2,"Initializing both distances for city ",basename(city.folder))
	
	# process the Euclidean distances on the discretized graph
	if(!light.process)
	{	e.file <- file.path(city.folder,"e-dist-discr.data")
		if(!file.exists(e.file))
		{	tlog(4,"Processing Euclidean distances")
			start.time <- Sys.time()
				pos <- cbind(vertex_attr(g2, name="x"),vertex_attr(g2, name="y"))
				e.dist <- dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2)
			end.time <- Sys.time()
			e.duration <- difftime(end.time,start.time,units="s")
			tlog(6,"Recording Euclidean distances (",e.duration,"s)")
			save(e.dist,file=e.file)
			# update table
			res.tables$graph[1,DISC_EDIST_TIME] <- e.duration
			table.file <- file.path(city.folder,res.graph.filename)
			write.table(x=res.tables$graph,file=table.file,row.names=FALSE,col.names=TRUE)
		}
	}
	
	# process the graph distances on the discretized graph
	if(!light.process)
	{	g.file <- file.path(city.folder,"g-dist-discr.data")
		if(!file.exists(g.file))
		{	tlog(4,"Processing graph distances")
			start.time <- Sys.time()
				g.dist <- shortest.paths(graph=g2, weights=E(g2)$dist)
			end.time <- Sys.time()
			g.duration <- difftime(end.time,start.time,units="s")
			tlog(6,"Recording graph distances (",g.duration,"s)")
			save(g.dist,file=g.file)
			# update table
			res.tables$graph[1,DISC_GDIST_TIME] <- g.duration
			table.file <- file.path(city.folder,res.graph.filename)
			write.table(x=res.tables$graph,file=table.file,row.names=FALSE,col.names=TRUE)
		}
	}
	
	return(res.tables)
}




############################################################################
# Processes the continuous version of the average straightness, for the
# specified graph.
#
# g: the original graph.
# city.folder: folder of the currently processed city.
# res.tables: list containing both result tables.
# light.process: whether the whole graph (FALSE) or only a sample of the nodes
#				 (TRUE) should be processed.
#
# returns: the list of tables containing the average straightness and corresponding
#		   durations values for the graph ($graph) and each node ($nodes).
############################################################################
process.continuous.straightness <- function(g, city.folder, res.tables, light.process)
{	tlog(2,"Processing the continuous average straightness for city ",basename(city.folder))
	gc()
	if(light.process)
		additional.duration <- 0
	else
		additional.duration <- res.tables$graph[1,ORIG_GDIST_TIME]
	
	# possibly load the previously processed distances
	if(!light.process)
	{	g.file <- file.path(city.folder,"g-dist-original.data")
		tlog(4,"Loading cached graph distances")
		load(g.file)
		tlog(6,"Graph distances processing time: ",res.tables$graph[1,ORIG_GDIST_TIME]," s")
	}
	
	# possibly process the whole graph
	if(light.process)
		tlog(4,"Light process, so not dealing with the whole graph")
	else
	{	# check if results already exist, in which case we don't process them again
		if(is.na(res.tables$graph[1,CONT_STRAIGHT]))
		{	# process straightness
			tlog(4,"Processing average over the whole graph")
			start.time <- Sys.time()
				value <- mean.straightness.graph(graph=g, g.dist=g.dist)
			end.time <- Sys.time()
			duration <- difftime(end.time,start.time,units="s") + additional.duration
			tlog(6,"Continuous average straightness: ",value," - Duration: ",duration," s")
			gc()
			
			# record the result table
			res.tables$graph[1,CONT_STRAIGHT] <- value
			res.tables$graph[1,CONT_PROC_TIME] <- duration
			table.file <- file.path(city.folder,res.graph.filename)
			write.table(x=res.tables$graph,file=table.file,row.names=FALSE,col.names=TRUE)
		}
		else
		{	value <- res.tables$graph[1,CONT_STRAIGHT]
			duration <- res.tables$graph[1,CONT_PROC_TIME]
			tlog(4,"Continuous average straightness: already processed (",value," - ",duration," s)")
		}
	}
	
	# process the nodes individually
	if(is.na(res.tables$nodes[nrow(res.tables$nodes),CONT_STRAIGHT]))
	{	tlog(2,"Processing continuous average for the individual node")
		imax <- nrow(res.tables$nodes)
		for(i in 1:imax)
		{	if(is.na(res.tables$nodes[i,CONT_STRAIGHT]))
			{	# process straightness
				start.time <- Sys.time()
					if(light.process)
						value <- mean.straightness.nodes.graph(graph=g, u=res.tables$nodes[i,NODE_ID], use.primitive=TRUE)
					else
						value <- mean.straightness.nodes.graph(graph=g, u=i, g.dist=g.dist, use.primitive=TRUE)
				end.time <- Sys.time()
				duration <- difftime(end.time,start.time,units="s") + additional.duration
				tlog(4,"Continuous average straightness for node ",i,"/",imax,": ",value," - Duration: ",duration," s")
				gc()
				
				# record the result table
				res.tables$nodes[i,CONT_STRAIGHT] <- value
				res.tables$nodes[i,CONT_PROC_TIME] <- duration
				table.file <- file.path(city.folder,res.nodes.filename)
				write.table(x=res.tables$nodes,file=table.file,row.names=FALSE,col.names=TRUE)
			}
			else
			{	value <- res.tables$nodes[i,CONT_STRAIGHT] 
				duration <- res.tables$nodes[i,CONT_PROC_TIME] 
				tlog(4,"Continuous average straightness for node ",i,"/",imax,": already processed (",value," - ",duration," s)")
			}
		}
	}
	
	return(res.tables)
}




############################################################################
# Processes the discrete version of the average straightness, for the
# specified graph.
#
# g: the original graph.
# city.folder: folder of the currently processed city.
# res.tables: list containing both result tables.
# light.process: whether the whole graph (FALSE) or only a sample of the nodes
#				 (TRUE) should be processed.
#
# returns: the list of tables containing the average straightness and corresponding
#		   durations values for the graph ($graph) and each node ($nodes).
############################################################################
process.discrete.straightness <- function(g, city.folder, res.tables, light.process)
{	tlog(2,"Processing the discrete approximation of the average straightness for city ",basename(city.folder))
	gc()
	if(light.process)
		additional.duration <- 0
	else
		additional.duration <- res.tables$graph[1,DISC_EDIST_TIME] + res.tables$graph[1,DISC_GDIST_TIME]
	
	# retrieve the discretized graph
	graph.file <- file.path(city.folder,"graph_discretized.graphml")
	tlog(4,"Loading the discretized graph")
	g2 <- read.graph(graph.file, format="graphml")
	additional.duration <- additional.duration + res.tables$graph[1,DISCRETIZATION_TIME]
	
	# possibly load the previously processed distances
	if(!light.process)
	{	e.file <- file.path(city.folder,"e-dist-discr.data")
		tlog(4,"Loading cached Euclidean distances")
		load(e.file)
		tlog(6,"Euclidean distances processing time: ",res.tables$graph[1,DISC_EDIST_TIME]," s")
		g.file <- file.path(city.folder,"g-dist-discr.data")
		tlog(4,"Loading cached graph distances")
		load(g.file)
		tlog(6,"Graph distances processing time: ",res.tables$graph[1,DISC_GDIST_TIME]," s")
	}
	
	# possibly process the whole graph
	if(light.process)
		tlog(4,"Light process, so not dealing with the whole graph")
	else
	{	# check if results already exist, in which case we don't process them again
		if(is.na(res.tables$graph[1,DISC_STRAIGHT]))
		{	# process straightness
			tlog(4,"Processing average over the whole graph")
			start.time <- Sys.time()
				value <- mean.straightness.nodes(graph=g2, v=NA, e.dist=e.dist, g.dist=g.dist, slow=TRUE)[1]
			end.time <- Sys.time()
			duration <- difftime(end.time,start.time,units="s") + additional.duration
			diff <- value - res.tables$graph[1,CONT_STRAIGHT]
			tlog(6,"Discrete average straightness: ",value," (Difference: ",diff,") - Duration: ",duration," s")
			gc()
			
			# record the result table
			res.tables$graph[1,DISC_STRAIGHT] <- value
			res.tables$graph[1,DISC_PROC_TIME] <- duration
			res.tables$graph[1,STRAIGHT_DIFFERENCE] <- diff
			table.file <- file.path(city.folder,res.graph.filename)
			write.table(x=res.tables$graph,file=table.file,row.names=FALSE,col.names=TRUE)
		}
		else
		{	value <- res.tables$graph[1,DISC_STRAIGHT]
			duration <- res.tables$graph[1,DISC_PROC_TIME]
			diff <- res.tables$graph[1,STRAIGHT_DIFFERENCE]
			tlog(4,"Continuous average straightness: already processed (",value," - ",duration," s)")
		}
	}
	
	# process the nodes individually
	if(is.na(res.tables$nodes[nrow(res.tables$nodes),DISC_STRAIGHT]))
	{	tlog(4,"Processing discrete average for each individual node")
		imax <- nrow(res.tables$nodes)
		for(i in 1:imax)
		{	if(is.na(res.tables$nodes[i,DISC_STRAIGHT]))
			{	# process straightness
				start.time <- Sys.time()
					if(light.process)
						value <- mean.straightness.nodes(graph=g2, v=res.tables$nodes[i,NODE_ID], slow=TRUE)[1,1]
					else
						value <- mean.straightness.nodes(graph=g2, v=i, e.dist=e.dist, g.dist=g.dist, slow=TRUE)[1,1]
				end.time <- Sys.time()
				duration <- difftime(end.time,start.time,units="s") + additional.duration
				diff <- value - res.tables$nodes[i,CONT_STRAIGHT]
				tlog(6,"Discrete average straightness for node ",i,"/",imax,": ",value," (Difference: ",diff,") - Duration: ",duration," s")
				gc()
				
				# record the result table
				res.tables$nodes[i,DISC_STRAIGHT] <- value
				res.tables$nodes[i,DISC_PROC_TIME] <- duration
				res.tables$nodes[i,STRAIGHT_DIFFERENCE] <- diff
				table.file <- file.path(city.folder,res.nodes.filename)
				write.table(x=res.tables$nodes,file=table.file,row.names=FALSE,col.names=TRUE)
			}
			else
			{	value <- res.tables$nodes[i,DISC_STRAIGHT] 
				duration <- res.tables$nodes[i,DISC_PROC_TIME] 
				diff <- res.tables$nodes[i,STRAIGHT_DIFFERENCE] 
				tlog(6,"Discrete average straightness for node ",i,"/",imax,": already processed (",value," - ",duration," s)")
			}
		}
	}
	
	g2 <- NULL
	gc()
	return(res.tables)
}




############################################################################
# Generates the plots for the current city only.
#
# city.folder: the folder of the current city.
# res.tables: previously processed stats.
############################################################################
generate.city.plots <- function(city.folder, res.tables)
{	tlog(2,"Generating plots for city ",basename(city.folder))
	
	# cont duration vs. disc duration
	plot.file <- file.path(city.folder,paste0("comparison-durations.pdf"))
	pdf(file=plot.file)
		plot(x=res.tables$nodes[,DISC_PROC_TIME],
			y=res.tables$nodes[,CONT_PROC_TIME],
			xlab="Discrete average processing time (s)",ylab="Continuous average processing time (s)"
		)
	dev.off()
	
	# cont straightness vs. disc straightness
	plot.file <- file.path(city.folder,paste0("comparison-straightness.pdf"))
	pdf(file=plot.file)
		plot(x=res.tables$nodes[,DISC_STRAIGHT],
			y=res.tables$nodes[,CONT_STRAIGHT],
			xlab="Discrete average Straightness",ylab="Continuous average Straightness"
		)
	dev.off()
	
#	# cont memory vs. disc memory
#	plot.file <- file.path(city.folder,paste0("comparison-memory.pdf"))
#	pdf(file=plot.file)
#		plot(x=res.tables$nodes[,RES_DISC_MEM],
#			y=res.tables$nodes[,RES_CONT_MEM],
#			xlab="Discrete average memory usage (MB)",ylab="Continuous average memory usage (MB)"
#		)
#	dev.off()
	
	# histogram of the processing times for cont and disc straightness processing
	plot.file <- file.path(city.folder,paste0("histo-durations.pdf"))
	pdf(file=plot.file)
		multi.hist(x1=res.tables$nodes[,CONT_PROC_TIME], 
				x2=res.tables$nodes[,DISC_PROC_TIME], 
				breaks=20,
				x.label="Processing time (s)", 
				series.names=c("Discrete average","Continuous average"), 
				leg.pos="topleft"
			)	
	dev.off()
	
	# histogram of the straightness differences
	plot.file <- file.path(city.folder,paste0("histo-differences.pdf"))
	pdf(file=plot.file)
		hist(x=res.tables$nodes[,STRAIGHT_DIFFERENCE],
			xlab="Straightness difference",ylab="Frequency",
			main=""
		)
	dev.off()
	
	# histogram of the memory usages for cont and disc straightness processing
#	plot.file <- file.path(city.folder,paste0("histo-memory.pdf"))
#	pdf(file=plot.file)
#	multi.hist(x1=res.tables$nodes[,CONT_MEM_USE], 
#			x2=res.tables$nodes[,DISC_MEM_USE], 
#			breaks=20,
#			x.label="Memory usage (MB)", 
#			series.names=c("Discrete average","Continuous average"), 
#			leg.pos="topleft"
#	)	
#	dev.off()
}




############################################################################
# Generates the plots comparing the results obtained for the collection of
# cities.
#
# cities: list of considered cities.
############################################################################
generate.overall.plots <- function(cities)
{	tlog(2,"Generating overall plots, considering all cities")
	
	# get the city names
	city.names <- c()
	for(city in cities)
		city.names <- c(city.names, city[[2]])
	# colors taken from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
	city.colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")[1:length(cities)]
	
	# retrieve all the necessary data
	tlog(4,"Retrieve the data for each city")
	for(c in 1:length(cities))
	{	# get city folder
		city <- names(cities)[c]
		tlog(6,"Dealing with city ",city)
		city.folder <- file.path(urban.folder,city)
		
		# retrieve the graph data 
		tlog(8,"Loading the graph-related data")
		res.graph.file <- file.path(city.folder,res.graph.filename)
		temp <- as.matrix(read.table(res.graph.file,header=TRUE,check.names=FALSE))
		if(c==1)
			graph.table <- temp
		else
			graph.table <- rbind(graph.table, temp)

		# retrieve the nodes data
		tlog(8,"Loading the node-related data")
		res.nodes.file <- file.path(city.folder,res.nodes.filename)
		temp <- as.matrix(read.table(res.nodes.file,header=TRUE,check.names=FALSE))
		if(nrow(temp)>sample.size)
		{	idx <- sample(nrow(temp),sample.size)
			temp <- temp[idx,]
			temp <- cbind(idx,temp)
			colnames(temp)[1] <- NODE_ID
		}
		temp <- cbind(temp,rep(c,nrow(temp)))
		colnames(temp)[length(colnames(temp))] <- NETWORK_ID
		if(c==1)
			nodes.table <- temp
		else
			nodes.table <- rbind(nodes.table, temp)
	}
	
	# generate the comparison plots
	tlog(4,"Generate discr. perf vs. cont perf plots")
	for(yaxis in c("duration","straightness"))#,"memory"))
	{	if(yaxis=="duration")
		{	d.col <- DISC_PROC_TIME
			c.col <- CONT_PROC_TIME
			x.lab <- "Discrete processing time (s)"
			y.lab <- "Continuous processing time (s)"
			leg.pos <- "topleft"
			inset <- 0.3
		}
		else if(yaxis=="straightness")
		{	d.col <- DISC_STRAIGHT
			c.col <- CONT_STRAIGHT
			x.lab <- "Discrete average Straightness"
			y.lab <- "Continuous average Straightness"
			leg.pos <- "topleft"
			inset <- 0.3
		}
		else if(yaxis=="memory")
		{	d.col <- DISC_PROC_TIME
			c.col <- CONT_PROC_TIME
			x.lab <- "Discrete memory usage (MB)"
			y.lab <- "Continuous memory usage (MB)"
			leg.pos <- "topleft"
			inset <- 0.3
		}
		
		plot.file <- file.path(urban.folder,paste0("comparison-",yaxis,".pdf"))
		pdf(file=plot.file)
			plot(x=nodes.table[which(nodes.table[,NETWORK_ID]==1),d.col],
				y=nodes.table[which(nodes.table[,NETWORK_ID]==1),c.col],
				xlab=x.lab,ylab=y.lab,
				col=city.colors[1]
			)
			if(nrow(graph.table)>1)
			{	for(net in 2:nrow(graph.table))
				{	points(x=nodes.table[which(nodes.table[,NETWORK_ID]==net),d.col],
						y=nodes.table[which(nodes.table[,NETWORK_ID]==net),c.col],
						col=city.colors[net]
					)
				}
			}
			legend(x=leg.pos,legend=city.names,
				inset=inset,
				fill=city.colors
			)
		dev.off()
	}
	
	# generate the perf vs. prop plots
	tlog(4,"Generate perf vs. prop plots")
	for(prop in c(PROP_NODE_NBR,PROP_DENSITY,PROP_LINK_NBR))
	{	xvals <- rep(graph.table[,prop],each=sample.size)
		x.lim <- c(min(graph.table[,prop]),max(graph.table[,prop]))
		
		if(prop==PROP_NODE_NBR)
		{	xaxis <- "nodes"
			xlab <- "Number of nodes"
			str.inset <- 0.03
			str.leg.pos <- "topleft"
			dur.inset <- 0.03
			dur.leg.pos <- "topleft"
			mem.inset <- 0.03
			mem.leg.pos <- "topleft"
			diff.inset <- 0.03
			diff.leg.pos <- "topleft"
		}
		else if(prop==PROP_DENSITY)
		{	xaxis <- "density"
			xlab <- "Density"
			str.inset <- 0.03
			str.leg.pos <- "topleft"
			dur.inset <- 0.03
			dur.leg.pos <- "topleft"
			mem.inset <- 0.03
			mem.leg.pos <- "topleft"
			diff.inset <- 0.03
			diff.leg.pos <- "topleft"
		} 
		else if(prop==PROP_LINK_NBR)
		{	xaxis <- "links"
			xlab <- "Number of edges"
			str.inset <- 0.03
			str.leg.pos <- "topleft"
			dur.inset <- 0.03
			dur.leg.pos <- "topleft"
			mem.inset <- 0.03
			mem.leg.pos <- "topleft"
			diff.inset <- 0.03
			diff.leg.pos <- "topleft"
		} 
		
		for(yaxis in c("straightness","duration"))#,"memory"))
		{	if(yaxis=="straightness")
			{	ylab <- "Straightness"
				disc.yvals <- nodes.table[,DISC_STRAIGHT]
				cont.yvals <- nodes.table[,CONT_STRAIGHT]
				inset <- str.inset
			}
			else if(yaxis=="duration")
			{	ylab <- "Processing time (s)"
				disc.yvals <- nodes.table[,DISC_PROC_TIME]
				cont.yvals <- nodes.table[,CONT_PROC_TIME]
				inset <- dur.inset
			}
			else if(yaxis=="memory")
			{	ylab <- "Memory usage (MB)"
				disc.yvals <- nodes.table[,DISC_MEM_USE]
				cont.yvals <- nodes.table[,CONT_MEM_USE]
				inset <- mem.inset
			}
			y.lim <- c(min(c(disc.yvals,cont.yvals)),max(c(disc.yvals,cont.yvals)))
			
			plot.file <- file.path(urban.folder,paste0(yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(x=xvals, 
					y=disc.yvals,
					xlab=xlab, ylab=ylab,
					col="BLUE",
					xlim=x.lim,
					ylim=y.lim
				)
				points(x=xvals,
					y=cont.yvals,
					pch=4,	#cross instead of circle (1)
					col="RED"
				)
				legend(x=str.leg.pos,legend=c("Discrete average","Continuous average"),
					inset=inset,
					fill=c("BLUE","RED")
				)
			dev.off()
		}
		
		{	yaxis <- "difference"
			ylab <- "Straightness difference"
			yvals <- nodes.table[,STRAIGHT_DIFFERENCE]
			inset <- diff.inset
			y.lim <- c(min(yvals),max(yvals))
			
			plot.file <- file.path(urban.folder,paste0(yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(x=xvals[which(nodes.table[,NETWORK_ID]==1)], 
					y=yvals[which(nodes.table[,NETWORK_ID]==1)],
					xlab=xlab, ylab=ylab,
					col=city.colors[1],
					xlim=x.lim,
					ylim=y.lim
				)
				if(nrow(graph.table)>1)
				{	for(net in 2:nrow(graph.table))
					{	points(x=xvals[which(nodes.table[,NETWORK_ID]==net)],
							y=yvals[which(nodes.table[,NETWORK_ID]==net)],
							col=city.colors[net]
						)
					}
				}
				legend(x=str.leg.pos,legend=city.names,
					inset=inset,
					fill=city.colors
				)
			dev.off()
		}
	}
}




############################################################################
# Process the time usage stats on a collection of real-world road networks,
# and generates the corresponding tables and plots.
#
# city: name of the graph.
############################################################################
monitor.time <- function(cities)
{	for(c in 1:length(cities))
	{	city <- names(cities)[c]
		light.process <- cities[[c]][[1]]
		tlog("Processing city ",city,if(light.process) " (light process)" else " (complete process)")
		gc()
		
		# retrieve the graph
		tlog(2,"Retrieve the graph for city ",city)
		city.folder <- file.path(urban.folder,city)
		net.file <- file.path(city.folder,"graph.graphml")
		g <- read.graph(file=net.file,format="graphml")
		
		# init/retrieve the result tables
		res.tables <- init.result.tables(g, city.folder, light.process)		
		
		# process both distances once and for all
		res.tables <- init.distances(g, city.folder, res.tables, light.process)
		# process the discretized graph once and for all / init the corresponding distances
		res.tables <- init.discretization(g, city.folder, res.tables, light.process)
		
		# deal with the continuous version
		res.tables <- process.continuous.straightness(g, city.folder, res.tables, light.process) 
		gc()
		
		# deal with the discrete version
		res.tables <- process.discrete.straightness(g, city.folder, res.tables, light.process)
		gc()
		
		# generate plots for the current city only
		generate.city.plots(city.folder, res.tables)
	}
	
	# generate plots considering all cities
#	generate.overall.plots(cities)
}




# the boolean value controls the processing: complete (FALSE) or light (TRUE)
# "light" means: only certain nodes, and not the whole graph (appropriate for large graphs)
cities <- list(
	abidjan=list(TRUE,"Abidjan")
#	alicesprings=list(TRUE,"Alice Springs"),
#	avignon=list(TRUE,"Avignon"),
#	beijin=list(TRUE,"Beijin"),
#	bordeaux=list(TRUE,"Bordeaux"),
#	dakar=list(TRUE,"Dakar"),
#	hongkong=list(TRUE,"Hong Kong"),
#	istanbul=list(TRUE,"Istanbul"),
#	karlskrona=list(TRUE,"Karlskrona"),
#	lisbon=list(TRUE,"Lisbon"),
#	liverpool=list(TRUE,"Liverpool"),
#	ljubljana=list(TRUE,"Ljubljana"),
#	maastricht=list(TRUE,"Maastricht"),
#	manhattan=list(TRUE,"Manhattan"),
###	newyork=list(TRUE,"New York"),
#	stpetersburg=list(TRUE,"St-Petersburg"),
#	roma=list(TRUE,"Roma"),
#	sfax=list(TRUE,"Sfax"),
#	soustons=list(TRUE,"Soustons")
#	tokyo=list(TRUE,"Tokyo"),
#	troisrivieres=list(FALSE,"Trois-Rivières")
)

monitor.time(cities)
