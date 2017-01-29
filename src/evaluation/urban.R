############################################################################
# Generate the time- and memory-related plots from the article cited in the 
# project readme file. This script focuses on the real-world road networks.
#
# Vincent Labatut 01/2017
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/evaluation/time.R")
############################################################################
source("src/straightness/continuous.R")
source("src/straightness/discrete.R")




data.folder <- "data"
urban.folder <- file.path(data.folder,"urban")
res.graph.filename <- "table-graph.txt"
res.nodes.filename <- "table-nodes.txt"




RES_EDIST_TIME <- "Edist-Time"
RES_GDIST_TIME <- "Gdist-Time"
RES_DISCRETIZATION_TIME <- "Dscrtz-Time"
RES_CONT_STR <- "Cont-Str"
RES_CONT_TIME <- "Cont-Time"
RES_CONT_MEM <- "Cont-Mem"
RES_DISC_STR <- "Disc-Str"
RES_DISC_TIME <- "Disc-Time"
RES_DISC_MEM <- "Disc-Mem"
RES_NODE_NBR <- "Node-Nbr"
RES_LINK_NBR <- "Link-Nbr"
RES_DENSITY <- "Density"




############################################################################
# Initializes the tables that will contain the results. If they already exist
# as files, they are loaded. Otherwise, they are initialized.
#
# g: the concerned graph.
# city.folder: folder of the currently processed city.
#
# returns: a list containing two tables: one for the whole graph, and one for
#		   the individual nodes.
############################################################################
init.result.tables <- function(g, city.folder)
{	res.tables <- list()
	
	# whole graph table
	res.graph.file <- file.path(city.folder,res.graph.filename)
	if(file.exists(res.graph.file))
		res.tables$graph <- as.matrix(read.table(res.graph.file,header=TRUE))
	else
	{	cnames <- c(EDIST_TIME,GDIST_TIME,RES_CONT_STR,RES_CONT_TIME,RES_CONT_MEM,RES_DISC_STR,RES_DISC_TIME,RES_DISC_MEM,RES_NODE_NBR,RES_LINK_NBR,RES_DENSITY)
		res.tables$graph <- matrix(NA,nrow=1,ncol=length(cnames))
		res.tables$graph[1,RES_NODE_NBR] <- gorder(g)
		res.tables$graph[1,RES_LINK_NBR] <- gsize(g)
		res.tables$graph[1,RES_DENSITY] <- graph.density(g)
	}
	
	# individual nodes table
	res.nodes.file <- file.path(city.folder,res.nodes.filename)
	if(file.exists(res.nodes.file))
		res.tables$nodes <- as.matrix(read.table(res.nodes.file,header=TRUE))
	else
	{	cnames <- c(RES_CONT_STR,RES_CONT_TIME,RES_CONT_MEM,RES_DISC_STR,RES_DISC_TIME,RES_DISC_MEM)
		res.tables$nodes <- matrix(NA,nrow=gorder(g),ncol=length(cnames))
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
# returns: a list containing the distances and the corresponding processing times.
############################################################################
init.distances <- function(g, city.folder, res.tables)
{	tlog(2,"Initializing both distances")
	
	# Euclidean distance
	e.file <- file.path(city.folder,"e-dist.txt")
	if(file.exists(e.file))
	{	tlog(2,"Loading Euclidean distances")
		e.dist <- load(e.file)
		e.duration <- res.tables$graph[1,RES_EDIST_TIME]
	}
	else
	{	tlog(4,"Processing Euclidean distances")
		start.time <- Sys.time()
			pos <- cbind(V(g)$x,V(g)$y)
			e.dist <- dist(x=pos, method="euclidean", diag=FALSE, upper=TRUE, p=2)
		end.time <- Sys.time()
		e.duration <- difftime(end.time,start.time,units="s")
		tlog(4,"Recording Euclidean distances")
		save(e.dist,file=e.file)
	}
	
	# graph distance
	g.file <- file.path(city.folder,"g-dist.txt")
	if(file.exists(g.file))
	{	tlog(2,"Loading graph distances")
		g.dist <- as.matrix(read.table(g.file,header=FALSE))
		g.duration <- res.tables$graph[1,RES_GDIST_TIME]
	}
	else
	{	tlog(4,"Processing graph distances")
		start.time <- Sys.time()
			g.dist <- shortest.paths(graph=g, weights=E(g)$dist)
		end.time <- Sys.time()
		g.duration <- difftime(end.time,start.time,units="s")
		tlog(4,"Recording graph distances")
		write.table(g.dist,g.file,row.names=FALSE,col.names=FALSE)
	}
	
	result <- list(g.dist=g.dist, e.dist=e.dist, g.time=g.duration, e.time=e.duration)
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
#
# returns: the list of tables containing the average straightness and corresponding
#		   durations values for the graph ($graph) and each node ($nodes).
############################################################################
process.continuous.straightness <- function(g, city.folder, dists, res.tables)
{	tlog("Processing the continuous average straightness")
	additional.duration <- res.tables$graph[1,RES_EDIST_TIME] + res.tables$graph[1,RES_GDIST_TIME]
	
	# check if results already exist, in which case we don't process them again
	if(is.na(res.tables$graph[1,RES_CONT_STR]))
	{	# process straightness
		tlog(2,"Processing average over the whole graph")
		start.time <- Sys.time()
			value <- mean.straightness.graph(graph=g, e.dist=dists$e.dist, g.dist=dists$g.dist)
		end.time <- Sys.time()
		duration <- difftime(end.time,start.time,units="s") + additional.duration
		tlog(4,"Continuous average straightness: ",value," - Duration: ",duration," s")
		gc()
		
		# record the result table
		res.tables$graph[1,RES_CONT_STR] <- value
		res.tables$graph[1,RES_CONT_TIME] <- duration
		table.file <- file.path(city.folder,res.graphs.filename)
		write.table(x=res.tables$graph,file=table.file,row.names=FALSE,col.names=TRUE)
	}
	
	# process each node in the graph
	if(is.na(res.tables$nodes[nrow(res.tables$graph),RES_CONT_STR]))
	{	tlog(2,"Processing continuous average for each individual node")
		for(i in 1:vcount(g))
		{	if(is.na(res.tables$nodes[i,RES_CONT_STR]))
			{	# process straightness
				start.time <- Sys.time()
					value <- mean.straightness.nodes.graph(graph=g, u=i, e.dist=dists$e.dist, g.dist=dists$g.dist) #TODO if we decide to switch to mean only: would be faster to process all nodes at once
				end.time <- Sys.time()
				duration <- difftime(end.time,start.time,units="s") + additional.duration
				tlog(4,"Continuous average straightness for node ",i,": ",value," - Duration: ",duration," s")
				gc()
				
				# record the result table
				res.tables$nodes[i,RES_CONT_STR] <- value
				res.tables$nodes[i,RES_CONT_TIME] <- duration
				table.file <- file.path(city.folder,res.nodes.filename)
				write.table(x=res.tables$nodes,file=table.file,row.names=FALSE,col.names=TRUE)
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
# avgseg: granularity of the discretization threshold (only one is used here, by 
#		opposition to the random networks, for which we sequentially used several 
#		values). 
#
# returns: the list of tables containing the average straightness and corresponding
#		   durations values for the graph ($graph) and each node ($nodes).
############################################################################
process.discrete.straightness <- function(g, city.folder, dists, res.tables, avgseg=50)
{	tlog("Processing the discrete approximation of the average straightness")
	additional.duration <- res.tables$graph[1,RES_EDIST_TIME] + res.tables$graph[1,RES_GDIST_TIME]
	
	# retrieve the discretized graph
	graph.file <- file.path(city.folder,"discretized.graphml")
	if(file.exists(graph.file))
	{	tlog(2,"Loading the discretized graph")
		g2 <- read.graph(graph.file, format="graphml")
		additional.duration <- additional.duration + res.tables$graph[1,RES_DISCRETIZATION_TIME]
	}
	else
	{	# discretize the graph
		tlog(2,"Discretizing the graph once and for all")
		start.time <- Sys.time()
			g2 <- add.intermediate.nodes(g, granularity=gran)
# TODO should actually use the average segmentation instead.
# maybe normalize first depending on the graph spread?
avgseg <- mean(E(g)$dist/grans[d])
			
		end.time <- Sys.time()
		duration <- difftime(end.time,start.time,units="s")
		additional.duration <- additional.duration + duration
		# record the graph
		write.graph(graph=g2,file=graph.file,format="graphml")
		# record the result table
		res.tables$graph[1,RES_DISCRETIZATION_TIME] <- duration
		# display for debug
		tlog(6,"Number of nodes: ",gorder(g2)," - Duration: ",duration," s - Average segmentation: ",avgseg," - Recorded as \"",graph.file,"\"")
	}
	
	
	
	
	
	tlog(4,"Loading the previously discretized graph from file \"",graph.file,"\"")
	g2 <- read.graph(file=graph.file,format="graphml")
	seg.duration <- disc.table[d+1,"SegmentationDuration"]
	
	# for the whole graph
	table.file <- file.path(it.folder,paste0("disc=",d,"-graph.txt"))
	# check if the table already exists, in which case it is simply loaded
	if(file.exists(table.file))
	{	tlog(4,"The results for the whole graph are already available: we load them from file \"",table.file,"\"")
		graph.table <- as.matrix(read.table(file=table.file,header=TRUE))
	}
	# otherwise, we process and record it 
	else
	{	tlog(4,"Processing average over the whole graph")
		start.time <- Sys.time()
			value <- mean.straightness.nodes(graph=g2, v=NA)[1]
		end.time <- Sys.time()
		duration <- difftime(end.time,start.time,units="s") + seg.duration
		diff <- value - cont.tables$graph[1,"Straightness"]
		tlog(6,"Discrete average straightness: ",value," (Difference: ",diff,") - Duration: ",duration," s")
		gc()
		
		# record the data as a text file
		tlog(4,"Record results as file \"",table.file,"\"")
		graph.table <- matrix(c(value,diff,duration),nrow=1)
		colnames(graph.table) <- c("Straightness","Difference","Duration")
		write.table(x=graph.table,file=table.file,row.names=FALSE,col.names=TRUE)
	}
	graph.tables[[as.character(d)]] <- graph.table
	
	# process each node in the graph
	table.file <- file.path(it.folder,paste0("disc=",d,"-nodes.txt"))
	# check if the table already exists, in which case it is simply loaded
	if(file.exists(table.file))
	{	tlog(4,"The results for the nodes are already available: we load them from file \"",table.file,"\"")
		nodes.table <- as.matrix(read.table(file=table.file,header=TRUE))
	}
	# otherwise, we process and record it 
	else
	{	tlog(4,"Processing discrete average for each individual node")
		values <- c()
		durations <- c()
		diffs <- c()
		for(i in 1:vcount(g))
		{	# process the straightness
			start.time <- Sys.time()
				value <- mean.straightness.nodes(graph=g2,v=i)[1,1]
			end.time <- Sys.time()
			values <- c(values,value)
			duration <- difftime(end.time,start.time,units="s") + seg.duration
			durations <- c(durations,duration)
			diff <- value - cont.tables$nodes[i,"Straightness"]
			diffs <- c(diffs,diff)
			tlog(6,"Discrete average straightness for node ",i,": ",value," (Difference: ",diff,") - Duration: ",duration," s")
			gc()
		}
		
		# record the data as a text file
		tlog(4,"Record results as file \"",table.file,"\"")
		nodes.table <- cbind(values,diffs,durations)
		colnames(nodes.table) <- c("Straightness","Difference","Duration")
		write.table(x=nodes.table,file=table.file,row.names=FALSE,col.names=TRUE)
	}
	nodes.tables[[as.character(d)]] <- nodes.table
	g2 <- NULL
	gc()
		
	return(result)
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
#
# returns: list of tables corresponding to a more compact representation of
#		   the previously processed values.
############################################################################
generate.rep.plots <- function(n=5, type="randplanar", iteration=1, disc.table, cont.tables, disc.tables)
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
			xvals <- disc.table[,"Nodes"]
		}
		else if(xaxis=="avgseg")
		{	xlab <- "Average segmentation"
			xvals <- disc.table[,"AverageSegmentation"]
		} 
		else if(xaxis=="granularity")
		{	xlab <- "Granularity"
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
					ylim=c(min(c(graph.yvals,graph.cont.val)),max(c(graph.yvals,graph.cont.val)))
			)
			lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),y=rep(graph.cont.val,2),col="RED")
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			for(j in 1:length(nodes.cont.vals))
			{	plot.file <- file.path(it.folder,paste0("node=",j,"-",yaxis,"-vs-",xaxis,".pdf"))
				pdf(file=plot.file)
				plot(x=xvals, y=nodes.yvals[j,],
						xlab=xlab, ylab=ylab,
						col="BLUE",
						ylim=c(min(c(nodes.yvals[j,],nodes.cont.vals[j])),max(c(nodes.yvals[j,],nodes.cont.vals[j])))
				)
				lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(nodes.cont.vals[j],2),
						col="RED"
				)
				legend(x="bottomright",legend=c("Approximation","Exact value"),
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
					ylim=c(min(c(graph.yvals,graph.cont.val)),max(c(graph.yvals,graph.cont.val)))
			)
			lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),y=rep(graph.cont.val,2),col="RED")
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			plot.file <- file.path(it.folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(x=rep(xvals,nrow(nodes.yvals)), y=c(t(nodes.yvals)),
					col="BLUE",#add.alpha("BLUE", 0.25),pch=20,
					xlab=xlab, ylab=ylab,
					ylim=c(min(c(nodes.yvals,nodes.cont.vals)),max(c(nodes.yvals,nodes.cont.vals)))
			)
			for(j in 1:length(nodes.cont.vals))
			{	lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(nodes.cont.vals[j],2),
						col="RED"#add.alpha("RED", 0.25)
				)
			}
			legend(x="bottomright",legend=c("Approximation","Exact value"),
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
					xlim=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
					xlab=xlab, ylab=ylab
			)
			lines(x=c(min(xvals,na.rm=TRUE), max(xvals,na.rm=TRUE)), y=c(0,0), col="BLACK", lty=2)
			points(x=xvals, y=graph.yvals,
					col="BLUE"
			)
			dev.off()
			
			# node plots
			plot.file <- file.path(it.folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(NULL,ylim=c(min(nodes.yvals),max(nodes.yvals)),
					xlim=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
					xlab=xlab, ylab=ylab
			)
			lines(x=c(min(xvals,na.rm=TRUE), max(xvals,na.rm=TRUE)), y=c(0,0), col="BLACK", lty=2)
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
{	for(city in cities)
	{	gc()
	
		# retrieve the graph
		city.folder <- file.path(urban.folder,city)
		net.file <- file.path(city.folder,"graph.graphml")
		g <- read.graph(file=net.file,format="graphml")
		
		# init/retrieve the result tables
		init.result.tables(g, city.folder)
		
		# init/retrieve both distances once and for all
		dists <- init.distances(g, city.folder)
		res.tables$graph[city,RES_EDIST_TIME] <- dists$e.time
		res.tables$graph[city,RES_GDIST_TIME] <- dists$g.time
		
		# deal with the continuous version
		res.tables <- process.continuous.straightness(g, city.folder, dists, res.tables) 
			
		# deal with the discrete version
		res.tables <- process.discrete.straightness(g, city.folder, dists, res.tables)
			
		# generate plots
		generate.plots(city.folder, res.tables)
	}
}


cities <- c(
#	"abidjan",
	"alicesprings",
	"avignon",
	"beijin",
	"bordeaux",
	"dakar",
	"hongkong",
	"istanbul",
	"karlskrona",
	"lisbon",
	"liverpool",
	"ljubljana",
	"maastricht",
	"manhattan",
	"stpetersburg",
	"roma",
	"sfax",
	"soustons",
	"tokyo",
	"troisrivieres"
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
