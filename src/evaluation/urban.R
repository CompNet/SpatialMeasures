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
			res.tables$nodes[,NODE_ID] <- sample(sample.size)
		}
		else
			res.tables$nodes <- matrix(NA,nrow=gorder(g),ncol=length(cnames))
		colnames(res.tables$nodes) <- cnames
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
# 
# returns: an updated graph table.
############################################################################
init.distances <- function(g, city.folder, res.tables)
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
	
	return(res.tables)
}




############################################################################
# Discretizes the graph and processes the corresponding Euclidean and graph
# distances between all pairs of nodes. Also caches the resulting data.
#
# g: concerned graph.
# city.folder: folder of the currently processed city.
# res.tables: result tables (possibly contain the distance processing times).
# avgseg: granularity of the discretization threshold (only one is used here, by 
#		opposition to the random networks, for which we sequentially used several 
#		values).
#
# returns: an updated graph table.
############################################################################
init.discretization <- function(g, city.folder, res.tables, avgseg=50)
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
	e.file <- file.path(city.folder,"e-dist-discr.data")
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
	
	# process the graph distances on the discretized graph
	g.file <- file.path(city.folder,"g-dist-discr.data")
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
	
	return(res.tables)
}




############################################################################
# Processes the continuous version of the average straightness, for the
# specified graph.
#
# g: the original graph.
# city.folder: folder of the currently processed city.
# dists: previously retrieved distances (both of them: Euclidean and graph).
# res.tables: list containing both result tables.
# light.process: whether the whole graph (FALSE) or only a sample of the nodes
#				 (TRUE) should be processed.
#
# returns: the list of tables containing the average straightness and corresponding
#		   durations values for the graph ($graph) and each node ($nodes).
############################################################################
process.continuous.straightness <- function(g, city.folder, dists, res.tables, light.process)
{	tlog(2,"Processing the continuous average straightness for city ",basename(city.folder))
	gc()
	if(light.process)
		additional.duration <- res.tables$graph[1,ORIG_GDIST_TIME]
	else
		additional.duration <- 0
	
	# possibly load the previously processed distances
	if(light.process)
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
				value <- mean.straightness.graph(graph=g, g.dist=dists$g.dist)
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
						value <- mean.straightness.nodes.graph(graph=g, u=i, e.dist=dists$e.dist, g.dist=dists$g.dist, use.primitive=TRUE)
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
# dists: previously retrieved distances (both of them: Euclidean and graph).
# res.tables: list containing both result tables.
# light.process: whether the whole graph (FALSE) or only a sample of the nodes
#				 (TRUE) should be processed.
#
# returns: the list of tables containing the average straightness and corresponding
#		   durations values for the graph ($graph) and each node ($nodes).
############################################################################
process.discrete.straightness <- function(g, city.folder, dists, res.tables, light.process)
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
	if(light.process)
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
				value <- mean.straightness.nodes(graph=g2, v=NA, e.dist=dists$e.dist, g.dist=dists$g.dist)[1]
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
			duration <- res.tables$graph[1,RES_DISCT_TIME]
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
						value <- mean.straightness.nodes(graph=g2, v=res.tables$nodes[i,NODE_ID])[1,1]
					else
						value <- mean.straightness.nodes(graph=g2, v=i, e.dist=dists$e.dist, g.dist=dists$g.dist)[1,1]
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
	
	# cont memory vs. disc memory
#	plot.file <- file.path(city.folder,paste0("comparison-memory.pdf"))
#	pdf(file=plot.file)
#		plot(x=res.tables$nodes[,RES_DISC_MEM],
#			y=res.tables$nodes[,RES_CONT_MEM],
#			xlab="Discrete average memory usage (MB)",ylab="Continuous average memory usage (MB)"
#		)
#	dev.off()
	
	# histogram of the durations for cont and disc straightness processing
	plot.file <- file.path(city.folder,paste0("histo-durations.pdf"))
	pdf(file=plot.file)
		hist(x=res.tables$nodes[,STRAIGHT_DIFFERENCE],
			xlab="Processing time (s)",ylab="Frequency"
		)
	dev.off()
	
	# histogram of the straightness differences
	plot.file <- file.path(city.folder,paste0("histo-differences.pdf"))
	pdf(file=plot.file)
		hist(x=res.tables$nodes[,STRAIGHT_DIFFERENCE],
			xlab="Straightness difference",ylab="Frequency"
		)
	dev.off()
}




############################################################################
# Generates the plots for each iteration.
#
# n: number of nodes.
# type: randplanar (planar random graph) or erdosrenyi (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# disc.table: discretization table.
# cont.tables: list of tables for the continuous average straightness.
# disc.tables: list of tables for the discrete average straightness.
# light.process: whether the whole graph (FALSE) or only a sample of the nodes
#				 (TRUE) should be processed.
#
# returns: list of tables corresponding to a more compact representation of
#		   the previously processed values.
############################################################################
generate.rep.plots <- function(n=5, type="randplanar", iteration=1, disc.table, cont.tables, disc.tables, light.process)
{	tlog("Generating plots and tables for the iteration ",iteration)
	it.folder <- file.path(urban.folder,type,paste0("n=",n),paste0("it=",iteration))
	nm <- paste0("d=",0:(nrow(disc.table)-1))
	
	# build the graph table
	tlog(2,"Building the graph table")
	graph.all <- matrix(NA,nrow=nrow(disc.table),ncol=3)
	colnames(graph.all) <- c("Straightness","Difference","Duration")
	for(d in 0:(nrow(disc.table)-1))
		graph.all[d+1,] <- c(disc.tables$graph[[as.character(d)]])
	rownames(graph.all) <- nm
	
	# record the graph table
	table.file <- file.path(it.folder,paste0("discrete-graph.txt"))
	write.table(graph.all,file=table.file,col.names=TRUE,row.names=TRUE)
	
	# build the nodes tables
	tlog(2,"Building the nodes tables")
	nodes.values <- matrix(NA,ncol=nrow(disc.table),nrow=disc.table[1,"Nodes"])
	colnames(nodes.values) <- nm
	nodes.differences <- matrix(NA,ncol=nrow(disc.table),nrow=disc.table[1,"Nodes"])	
	colnames(nodes.differences) <- nm
	nodes.durations <- matrix(NA,ncol=nrow(disc.table),nrow=disc.table[1,"Nodes"])	
	colnames(nodes.durations) <- nm
	for(d in 0:(nrow(disc.table)-1))
	{	nodes.values[,d+1] <- disc.tables$nodes[[as.character(d)]][,"Straightness"]
		nodes.differences[,d+1] <- disc.tables$nodes[[as.character(d)]][,"Difference"]
		nodes.durations[,d+1] <- disc.tables$nodes[[as.character(d)]][,"Duration"]
	}
	
	# record the nodes tables
	table.file <- file.path(it.folder,paste0("discrete-nodes-straightness.txt"))
	write.table(nodes.values,file=table.file,col.names=TRUE,row.names=FALSE)
	table.file <- file.path(it.folder,paste0("discrete-nodes-differences.txt"))
	write.table(nodes.differences,file=table.file,col.names=TRUE,row.names=FALSE)
	table.file <- file.path(it.folder,paste0("discrete-nodes-durations.txt"))
	write.table(nodes.durations,file=table.file,col.names=TRUE,row.names=FALSE)
	
	# generate the plots
	tlog(2,"Generating the plots for iteration ",iteration)
	for(xaxis in c("nodes","avgseg","granularity"))
	{	if(xaxis=="nodes")
		{	xlab <- "Total number of nodes"
			log.axes <- ""
			xvals <- disc.table[,"Nodes"]
		}
		else if(xaxis=="avgseg")
		{	xlab <- "Average segmentation"
			log.axes <- "x"
			xvals <- disc.table[,"AverageSegmentation"] + 1
		} 
		else if(xaxis=="granularity")
		{	xlab <- "Granularity"
			log.axes <- ""
			xvals <- disc.table[,"Granularity"]
		} 
		
		# straightness
		{	yaxis <- "straightness"
			ylab <- "Straightness"
			graph.yvals <- graph.all[,"Straightness"]
			graph.cont.val <- cont.tables$graph[1,"Straightness"]
			nodes.yvals <- nodes.values
			nodes.cont.vals <- cont.tables$nodes[,"Straightness"]
			
			# graph plot
			plot.file <- file.path(it.folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(x=xvals, y=graph.yvals,
					xlab=xlab, ylab=ylab,
					col="BLUE",
					log=log.axes,
					ylim=c(min(c(graph.yvals,graph.cont.val)),max(c(graph.yvals,graph.cont.val)))
			)
			lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
					y=rep(graph.cont.val,2),
					col="RED"
			)
			legend(x="bottomright",legend=c("Discrete average","Continuous average"),
					inset=0.03,
					fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			for(j in 1:length(nodes.cont.vals))
			{	plot.file <- file.path(it.folder,paste0("node=",j,"-",yaxis,"-vs-",xaxis,".pdf"))
				pdf(file=plot.file)
				plot(x=xvals, y=nodes.yvals[j,],
						xlab=xlab, ylab=ylab,
						col="BLUE",
						log=log.axes,
						ylim=c(min(c(nodes.yvals[j,],nodes.cont.vals[j])),max(c(nodes.yvals[j,],nodes.cont.vals[j])))
				)
				lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(nodes.cont.vals[j],2),
						col="RED"
				)
				legend(x="bottomright",legend=c("Discrete average","Continuous average"),
						inset=0.03,
						fill=c("BLUE","RED"))
				dev.off()
			}
		}
		
		# durations
		{	yaxis <- "durations"
			ylab <- "Time (s)"
			graph.yvals <- graph.all[,"Duration"]
			graph.cont.val <- cont.tables$graph[1,"Duration"]
			nodes.yvals <- nodes.durations
			nodes.cont.vals <- cont.tables$nodes[,"Duration"]
			
			# graph plots
			plot.file <- file.path(it.folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(x=xvals, y=graph.yvals,
					xlab=xlab, ylab=ylab,
					col="BLUE",
					log=log.axes,
					ylim=c(min(c(graph.yvals,graph.cont.val)),max(c(graph.yvals,graph.cont.val)))
			)
			lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
					y=rep(graph.cont.val,2),
					col="RED"
			)
			legend(x="bottomright",legend=c("Discrete average","Continuous average"),
					inset=0.03,
					fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			plot.file <- file.path(it.folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(x=rep(xvals,nrow(nodes.yvals)), y=c(t(nodes.yvals)),
					col="BLUE",#add.alpha("BLUE", 0.25),pch=20,
					xlab=xlab, ylab=ylab,
					log=log.axes,
					ylim=c(min(c(nodes.yvals,nodes.cont.vals)),max(c(nodes.yvals,nodes.cont.vals)))
			)
			for(j in 1:length(nodes.cont.vals))
			{	lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(nodes.cont.vals[j],2),
						col="RED"#add.alpha("RED", 0.25)
				)
			}
			legend(x="bottomright",legend=c("Discrete average","Continuous average"),
					inset=0.03,
					fill=c("BLUE","RED"))
			dev.off()
		}
		
		# straightness differences
		{	yaxis <- "differences"
			ylab <- "Straightness difference"
			graph.yvals <- graph.all[,"Difference"]
			nodes.yvals <- nodes.differences
			
			# graph plots
			plot.file <- file.path(it.folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(NULL,ylim=c(min(graph.yvals),max(graph.yvals)),
					xlim=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)+10),
					log=log.axes,
					xlab=xlab, ylab=ylab
			)
			lines(x=c(min(xvals,na.rm=TRUE), max(xvals,na.rm=TRUE)+10), 
					y=c(0,0), 
					col="BLACK", lty=2
			)
			points(x=xvals, y=graph.yvals,
					col="BLUE"
			)
			dev.off()
			
			# node plots
			plot.file <- file.path(it.folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(NULL,ylim=c(min(nodes.yvals),max(nodes.yvals)),
					xlim=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)+10),
					log=log.axes,
					xlab=xlab, ylab=ylab
			)
			lines(x=c(min(xvals,na.rm=TRUE), max(xvals,na.rm=TRUE)+10), 
					y=c(0,0), 
					col="BLACK", lty=2
			)
			points(x=rep(xvals,nrow(nodes.yvals)), y=c(t(nodes.yvals)),
					col="BLUE"#add.alpha("BLUE", 0.25),pch=20					
			)
			dev.off()
		}
	}
	
	result <- 	list(graph=graph.all, 
			nodes=list(straightness=nodes.values,
					differences=nodes.differences,
					durations=nodes.durations))
	return(result)
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
		light.process <- cities[[c]][1]
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
		if(!light.process)
			res.tables <- init.distances(g, city.folder, res.tables)
		# process the discretized graph once and for all / init the corresponding distances
		res.tables <- init.discretization(g, city.folder, res.tables)
		
		# deal with the continuous version
		res.tables <- process.continuous.straightness(g, city.folder, dists, res.tables, light.process) 
		gc()
		
		# deal with the discrete version
		res.tables <- process.discrete.straightness(g, city.folder, dists, res.tables, light.process)
		gc()
		
		# generate plots for the current city only
		generate.city.plots(city.folder, res.tables)
	}
	
	# generate plots considering all cities
#	generate.overall.plots(cities, light.process)
}


# the boolean value controls the processing: complete (FALSE) or light (TRUE)
# "light" means: only certain nodes, and not the whole graph (appropriate for large graphs)
cities <- list(
#	abidjan=c(FALSE),
#	alicesprings=c(TRUE),
#	avignon=c(TRUE),
#	beijin=c(TRUE),
#	bordeaux=c(TRUE),
#	dakar=c(TRUE),
#	hongkong=c(TRUE),
#	istanbul=c(TRUE),
#	karlskrona=c(TRUE),
#	lisbon=c(TRUE),
#	liverpool=c(TRUE),
#	ljubljana=c(TRUE),
#	maastricht=c(TRUE),
#	manhattan=c(TRUE),
###	newyork=c(FALSE),
#	stpetersburg=c(TRUE),
#	roma=c(TRUE),
#	sfax=c(TRUE),
	soustons=c(FALSE)
#	tokyo=c(TRUE),
#	troisrivieres=c(FALSE)
)

monitor.time(cities)



# for each city
#   load the graph
#   process the continuous variants while measuring the time
#   record these times
#   process the discrete variants while also measuring the time
#   record these values
# data structure:
#   Graphs: matrix with 
#		- rows=cities
#		- cols=cont strgtns, cont time, cont mem, disc strgtns, disc time,  disc mem, node number, link number, density
#   Nodes: matrix with
#		- rows=nodes (cities are placed sequentially)
#		- cols=city name, cont strgtns, cont time, cont mem, disc strgtns, disc time, disc mem


# TODO : plot discrete vs. continuous values?
# like: for each one, time disc vs. time cont
# same for memory
# and same for graphs (?)

#TODO only process the values that have not been processed yet
#(so need to record them)
# maybe just focus on certain nodes ? dealing with all distances seems too much


#TODO must process only the necessary distances when specifying nodes
#TODO possility to restrict the processing:
#	- do not preprocess the distances
#	- not the whole graph
#	- only a sample of the nodes (not all of them)


# setwd("~/eclipse/workspaces/Networks/SpatialMeasures");source("src/evaluation/urban.R")