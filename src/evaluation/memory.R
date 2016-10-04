############################################################################
# Generate the memory-related plots from the article cited in the project
# readme file. This script must be executed *after* the one dealing with
# processing times, because it takes advantages of certain files it produces. 
# 
# Note: the memory profiler of R is kind of buggy. I had to launch this script
# a bunch of time to gather enough data.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/evaluation/memory.R")
############################################################################
source("src/evaluation/common.R")

MEM_FILE <- "data/profiling.txt"
MEM_INTER <- 0.00002
MEM_SLEEP <- 0.5 	# in seconds



############################################################################
# Discretizes the original graph and record all the resulting graphs as 
# Graphml files for later use. A table summarizing the process is recorded
# too, also for later use. If the table already exists as a file, we assume 
# the graph files also exist, and load the table then return it.
#
# n: number of nodes.
# type: RAND_PLANAR (planar random graph) or ERDOS_RENYI (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# g: the original graph.
#
# returns: the table summarizing the discretization process.
############################################################################
load.disc.table <- function(n=5, type="RAND_PLANAR", iteration=1, g)
{	tlog(2,"Loading the discretization table")
	
	# init file names
	it.folder <- file.path("data",type,paste0("n=",n),paste0("it=",iteration))
	disc.file <- file.path(it.folder,"discretizations.txt")

	# get the table
	if(file.exists(disc.file))
		disc.table <- as.matrix(read.table(file=disc.file,header=TRUE))
	else
		stop("Discretization table not found. You must first execute the time.R script.")
	nbr.disc <- nrow(disc.table) - 1
	
	# process each granularity value
	memuses <- c()
	for(d in 0:nbr.disc)
	{	tlog(4,"Granularity: ",disc.table[d+1,"Granularity"]," (",d,"/",nbr.disc,")")
		memuse <- 0
		if(d>0)
		{	again <- TRUE
			while(again)
			{	again <- FALSE
				gc()
				Rprof(MEM_FILE, memory.profiling=TRUE, interval=MEM_INTER)
#					Sys.sleep(MEM_SLEEP)
					g2 <- add.intermediate.nodes(g, granularity=disc.table[d+1,"Granularity"])
#					Sys.sleep(MEM_SLEEP)
				Rprof(NULL)
#				Sys.sleep(MEM_SLEEP)
				# update memory usage
				mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[["\"source\""]],
					error=function(e){again <<- TRUE})
				if(again)
					tlog(4,"Error while trying to process memory usage. Trying again.")
				else
					memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
			}
		}
		memuses <- c(memuses,memuse)
		tlog(6,"Memory: ",memuse," MB")
		
		# reset to free memory
		g2 <- NULL
	}
	
	# record the table
	if("SegmentationMemory" %in% colnames(disc.table))
		disc.table[,"SegmentationMemory"] <- memuses
	else
	{	disc.table <- cbind(disc.table,memuses)
		colnames(disc.table)[ncol(disc.table)] <- "SegmentationMemory"
	}
	write.table(disc.table,file=disc.file,col.names=TRUE,row.names=FALSE)
	
	return(disc.table)
}



############################################################################
# Processes the continuous version of the average straightness, for the
# specified graphs.
#
# n: number of nodes.
# type: RAND_PLANAR (planar random graph) or ERDOS_RENYI (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# g: the original graph.
#
# returns: the list of tables containing the memory usage.
############################################################################
process.continuous.straightness <- function(n=5, type="RAND_PLANAR", iteration=1, g)
{	tlog("Processing the continuous average straightness")
	it.folder <- file.path("data",type,paste0("n=",n),paste0("it=",iteration))
	
	# for the whole graph
	table.file <- file.path(it.folder,"continuous-graph.txt")
	if(file.exists(table.file))
		graph.table <- as.matrix(read.table(file=table.file,header=TRUE))
	else
		stop("Whole graph table not found (\"",table.file,"\"). You must first execute the time.R script.")
	tlog(2,"Processing average over the whole graph")
	again <- TRUE
	while(again)
	{	again <- FALSE
		gc()
		Rprof(MEM_FILE, memory.profiling=TRUE, interval=MEM_INTER)
			value <- mean.straightness.graph(graph=g)
		Rprof(NULL)
		mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[["\"source\""]],
			error=function(e){again <<- TRUE})
		if(again)
			tlog(4,"Error while trying to process memory usage. Trying again.")
		else
			memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
	}
	tlog(4,"Continuous average straightness: ",value," - Memory: ",memuse," MB")
	
	# record the data as a text file
	tlog(4,"Record results as file \"",table.file,"\"")
	if("Memory" %in% colnames(graph.table))
		graph.table[,"Memory"] <- memuse
	else
	{	graph.table <- cbind(graph.table,memuse)
		colnames(graph.table)[ncol(graph.table)] <- "Memory"
	}
	write.table(x=graph.table,file=table.file,row.names=FALSE,col.names=TRUE)
	
	# process each node in the graph
	table.file <- file.path(it.folder,"continuous-nodes.txt")
	if(file.exists(table.file))
		nodes.table <- as.matrix(read.table(file=table.file,header=TRUE))
	else
		stop("Individual nodes table not found (\"",table.file,"\"). You must first execute the time.R script.")
	tlog(2,"Processing continuous average for each individual node")
	memuses <- c()
	for(i in 1:vcount(g))
	{	again <- TRUE
		while(again)
		{	again <- FALSE
			gc()
			Rprof(MEM_FILE, memory.profiling=TRUE, interval=MEM_INTER)
				value <- mean.straightness.nodes.graph(graph=g, u=i)
			Rprof(NULL)
			mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[["\"source\""]],
				error=function(e){again <<- TRUE})
			if(again)
				tlog(4,"Error while trying to process memory usage. Trying again.")
			else
				memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
		}
		memuses <- c(memuses,memuse)
		tlog(4,"Continuous average straightness for node ",i,": ",value," - Memory: ",memuse," MB")
	}
	
	# record the data as a text file
	tlog(2,"Record results as file \"",table.file,"\"")
	if("Memory" %in% colnames(nodes.table))
		nodes.table[,"Memory"] <- memuses
	else
	{	nodes.table <- cbind(nodes.table,memuses)
		colnames(nodes.table)[ncol(nodes.table)] <- "Memory"
	}
	write.table(x=nodes.table,file=table.file,row.names=FALSE,col.names=TRUE)
	
	result <- list(graph=graph.table, nodes=nodes.table)
	return(result)
}



############################################################################
# Processes the discrete version of the average straightness, for the
# specified graphs.
#
# n: number of nodes.
# type: RAND_PLANAR (planar random graph) or ERDOS_RENYI (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# g: the original graph.
#
# returns: the list of lists of tables containing the memory usage.
############################################################################
process.discrete.straightness <- function(n=5, type="RAND_PLANAR", iteration=1, g, cont.tables)
{	tlog("Processing the discrete approximation of the average straightness")
	it.folder <- file.path("data",type,paste0("n=",n),paste0("it=",iteration))
	
	# load the discretization table
	disc.file <- file.path(it.folder,"discretizations.txt")
	tlog(2,"Loading the discretization table from file \"",disc.file,"\"")
	disc.table <- as.matrix(read.table(file=disc.file,header=TRUE))
	nbr.disc <- nrow(disc.table) - 1
	
	# process each discretization step
	graph.tables <- list()
	nodes.tables <- list()
	for(d in 0:nbr.disc)
	{	tlog(2,"Processing discretization number ",d,"/",nbr.disc)
		
		# load the discretized graph
		graph.file <- file.path(it.folder,paste0("disc=",d,".graphml"))
		tlog(4,"Loading the previously discretized graph from file \"",graph.file,"\"")
		g2 <- read.graph(file=graph.file,format="graphml")
		seg.mem <- disc.table[d+1,"SegmentationMemory"]
		
		# for the whole graph
		table.file <- file.path(it.folder,paste0("disc=",d,"-graph.txt"))
		if(file.exists(table.file))
			graph.table <- as.matrix(read.table(file=table.file,header=TRUE))
		else
			stop("Whole graph table not found (\"",table.file,"\"). You must first execute the time.R script.")
		tlog(4,"Processing average over the whole graph")
		again <- TRUE
		while(again)
		{	again <- FALSE
			gc()
			Rprof(MEM_FILE, memory.profiling=TRUE, interval=MEM_INTER)
				value <- mean.straightness.nodes(graph=g2, v=NA)[1]
			Rprof(NULL)
			mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[["\"source\""]],
				error=function(e){again <<- TRUE})
			if(again)
				tlog(4,"Error while trying to process memory usage. Trying again.")
			else
				memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
		}
		memuse <- max(seg.mem,memuse)
		tlog(6,"Discrete average straightness: ",value," - Memory: ",memuse," MB")
		
		# record the data as a text file
		tlog(4,"Record results as file \"",table.file,"\"")
		if("Memory" %in% colnames(graph.table))
			graph.table[1,"Memory"] <- memuse
		else
		{	graph.table <- cbind(graph.table,memuse)
			colnames(graph.table)[ncol(graph.table)] <- "Memory"
		}
		write.table(x=graph.table,file=table.file,row.names=FALSE,col.names=TRUE)

		graph.tables[[as.character(d)]] <- graph.table
		
		# process each node in the graph
		table.file <- file.path(it.folder,paste0("disc=",d,"-nodes.txt"))
		if(file.exists(table.file))
			nodes.table <- as.matrix(read.table(file=table.file,header=TRUE))
		# otherwise, we process and record it 
		else
			stop("Individual nodes table not found (\"",table.file,"\"). You must first execute the time.R script.")
		tlog(4,"Processing discrete average for each individual node")
		memuses <- c()
		for(i in 1:vcount(g))
		{	again <- TRUE
			while(again)
			{	again <- FALSE
				gc()
				Rprof(MEM_FILE, memory.profiling=TRUE, interval=MEM_INTER)
					value <- mean.straightness.nodes(graph=g2,v=i)[1,1]
				Rprof(NULL)
				mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[["\"source\""]],
					error=function(e){again <<- TRUE})
				if(again)
					tlog(4,"Error while trying to process memory usage. Trying again.")
				else
					memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
			}
			memuse <- max(seg.mem,memuse)
			memuses <- c(memuses,memuse)
			tlog(6,"Discrete average straightness for node ",i,": ",value," - Memory: ",memuse," MB")
		}
		
		# record the data as a text file
		tlog(4,"Record results as file \"",table.file,"\"")
		if("Memory" %in% colnames(nodes.table))
			nodes.table[,"Memory"] <- memuses
		else
		{	nodes.table <- cbind(nodes.table,memuses)
			colnames(nodes.table)[ncol(nodes.table)] <- "Memory"
		}
		write.table(x=nodes.table,file=table.file,row.names=FALSE,col.names=TRUE)
		
		nodes.tables[[as.character(d)]] <- nodes.table
		g2 <- NULL
		gc()
	}
	
	result <- list(graph=graph.tables, nodes=nodes.tables)
	return(result)
}



############################################################################
# Generates the plots for each iteration.
#
# n: number of nodes.
# type: RAND_PLANAR (planar random graph) or ERDOS_RENYI (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# disc.table: discretization table.
# cont.tables: list of tables for the continuous average straightness.
# disc.tables: list of tables for the discrete average straightness.
#
# returns: list of tables corresponding to a more compact representation of
#		   the previously processed memory usages.
############################################################################
generate.rep.plots <- function(n=5, type="RAND_PLANAR", iteration=1, disc.table, cont.tables, disc.tables)
{	tlog("Generating plots and tables for the iteration ",iteration)
	it.folder <- file.path("data",type,paste0("n=",n),paste0("it=",iteration))
	nm <- paste0("d=",0:(nrow(disc.table)-1))
	
	# build the graph table
	tlog(2,"Building the graph table")
	graph.all <- matrix(NA,nrow=nrow(disc.table),ncol=4)
	colnames(graph.all) <- c("Straightness","Difference","Duration","Memory")
	for(d in 0:(nrow(disc.table)-1))
		graph.all[d+1,] <- c(disc.tables$graph[[as.character(d)]])
	rownames(graph.all) <- nm
	
	# record the graph table
	table.file <- file.path(it.folder,paste0("discrete-graph.txt"))
	write.table(graph.all,file=table.file,col.names=TRUE,row.names=TRUE)
		
	# build the nodes tables
	tlog(2,"Building the nodes tables")
	nodes.memory <- matrix(NA,ncol=nrow(disc.table),nrow=disc.table[1,"Nodes"])
	colnames(nodes.memory) <- nm
	for(d in 0:(nrow(disc.table)-1))
		nodes.memory[,d+1] <- disc.tables$nodes[[as.character(d)]][,"Memory"]
	
	# record the nodes tables
	table.file <- file.path(it.folder,paste0("discrete-nodes-memory.txt"))
	write.table(nodes.memory,file=table.file,col.names=TRUE,row.names=FALSE)
	
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
		
		# memory
		{	yaxis <- "memory"
			ylab <- "Memory (MB)"
			graph.yvals <- graph.all[,"Memory"]
			graph.cont.val <- cont.tables$graph[1,"Memory"]
			nodes.yvals <- nodes.memory
			nodes.cont.vals <- cont.tables$nodes[,"Memory"]
			
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
					pch=20,
					xlab=xlab, ylab=ylab,
					col=add.alpha("BLUE", 0.25),
					ylim=c(min(c(nodes.yvals,nodes.cont.vals)),max(c(nodes.yvals,nodes.cont.vals)))
			)
			for(j in 1:length(nodes.cont.vals))
			{	lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(nodes.cont.vals[j],2),
						col=add.alpha("RED", 0.25)
				)
			}
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
		}
	}
	
	result <- 	list(graph=graph.all, 
			nodes=list(memuses=nodes.memory))
	return(result)
}



############################################################################
# Generates the plots for all iteration.
#
# n: number of nodes.
# type: RAND_PLANAR (planar random graph) or ERDOS_RENYI (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# discretizations: list of discretization tables.
# data.cont: list of lists of tables for the continuous average straightness.
# data.disc: list of list of tables for the discrete average straightness.
############################################################################
generate.overall.plots <- function(n=10, type="RAND_PLANAR", discretizations, data.cont, data.disc)
{	tlog("Generating plots and tables for all the repetitions")
	folder <- file.path("data",type,paste0("n=",n))
	
	# collecting the data
	ngran <- nrow(data.disc[[1]]$graph)
	graph.cont.memory <- data.cont[[1]]$graph[1,"Memory"]
	nodes.cont.memory <- data.cont[[1]]$nodes[,"Memory"]
	graph.disc.memory <- data.disc[[1]]$graph[,"Memory"]
	nodes.disc.memory <- data.disc[[1]]$nodes$memuses
	granularities <- discretizations[[1]][,"Granularity"]
	nodes <- discretizations[[1]][,"Nodes"]
	avgseg <- discretizations[[1]][,"AverageSegmentation"]
	for(r in 2:length(data.disc))
	{	# continuous results
		graph.cont.memory <- c(graph.cont.memory,data.cont[[r]]$graph[1,"Memory"])
		nodes.cont.memory <- c(nodes.cont.memory,data.cont[[r]]$nodes[,"Memory"])
		# discrete results
		graph.disc.memory <- c(graph.disc.memory,data.disc[[r]]$graph[,"Memory"])
		nodes.disc.memory <- cbind(nodes.disc.memory,data.disc[[r]]$nodes$memuses)
		# x values
		granularities <- c(granularities,discretizations[[r]][,"Granularity"])
		nodes <- c(nodes,discretizations[[r]][,"Nodes"])
		avgseg <- c(avgseg,discretizations[[r]][,"AverageSegmentation"])
	}
	
	# generating the plots
	tlog(2,"Generating the plots for all iterations")
	for(xaxis in c("nodes","avgseg","granularity"))
	{	if(xaxis=="nodes")
		{	xlab <- "Total number of nodes"
			xvals <- nodes
		}
		else if(xaxis=="avgseg")
		{	xlab <- "Average segmentation"
			xvals <- avgseg
		} 
		else if(xaxis=="granularity")
		{	xlab <- "Granularity"
			xvals <- granularities
		} 
		
		# memory
		{	yaxis <- "memory"
			ylab <- "Memory (MB)"
			
			# graph plots
			plot.file <- file.path(folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(x=xvals, y=graph.disc.memory,
					pch=20,
					xlab=xlab, ylab=ylab,
					col=add.alpha("BLUE", 0.25),
					ylim=c(min(c(graph.disc.memory,graph.cont.memory)),max(c(graph.disc.memory,graph.cont.memory)))
			)
			for(j in 1:length(graph.cont.memory))
			{	lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(graph.cont.memory[j],2),
						col=add.alpha("RED", 0.25)
				)
			}
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			plot.file <- file.path(folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(x=rep(xvals,nrow(nodes.disc.memory)), y=c(t(nodes.disc.memory)),
					pch=20,
					xlab=xlab, ylab=ylab,
					col=add.alpha("BLUE", 0.25),
					ylim=c(min(c(nodes.disc.memory,nodes.cont.memory)),max(c(nodes.disc.memory,nodes.cont.memory)))
			)
			for(j in 1:length(nodes.cont.memory))
			{	lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(nodes.cont.memory[j],2),
						col=add.alpha("RED", 0.25)
				)
			}
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
		}
	}
}



############################################################################
# Process the time usage stats on a collection of randomly generated networks,
# and generates the corresponding tables and plots.
#
# n: number of nodes in the graph.
# type: type of graph, RAND_PLANAR (planar random graph) or ERDOS_RENYI 
#		(Erdös-Rényi random graph).
# repetitions: number of instances of the graph to generate and process.
############################################################################
monitor.memory <- function(n=5, type="RAND_PLANAR", repetitions=10)
{	gc()
	
	# process the specified number of repetitions
	data.disc <- list()
	data.cont <- list()
	discretizations <- list()
	for(r in 1:repetitions)
	{	# retrieve or create the graph
		g <- init.graph(n, type, iteration=r)
		
		# retrieve the discretization table
		disc.table <- load.disc.table(n, type, iteration=r, g)
		
		# deal with the continuous version
		cont.tables <- process.continuous.straightness(n,type,iteration=r,g)
		
		# deal with the discrete version
		disc.tables <- process.discrete.straightness(n,type,iteration=r,g,cont.tables)
		
		# generate plots
		data.disc[[r]] <- generate.rep.plots(n, type, iteration=r, disc.table, cont.tables, disc.tables)

		data.cont[[r]] <- cont.tables
		discretizations[[r]] <- disc.table
		gc()
	}
	
	# generate plots using the repetitions
	generate.overall.plots(n, type, discretizations, data.cont, data.disc)
}

monitor.memory(n=10, type="RAND_PLANAR", repetitions=10)
#monitor.memory(n=25, type="RAND_PLANAR", repetitions=10)
#monitor.memory(n=50, type="RAND_PLANAR", repetitions=10)
#monitor.memory(n=100, type="RAND_PLANAR", repetitions=10)
