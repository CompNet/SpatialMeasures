############################################################################
# Generate the time-related plots from the article cited in the project
# readme file.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/time.R")
############################################################################
source("src/common/misc.R")
source("src/common/plot.R")
source("src/common/transformations.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")



############################################################################
# Creates a new graph corresponding to the specified parameters, or loads
# it from the file, if it exists. When created, the graph is also recorded
# as a Graphml file.
#
# n: number of nodes.
# type: RAND_PLANAR (planar random graph) or ERDOS_RENYI (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
#
# returns: the created or loaded graph.
############################################################################
init.graph <- function(n=5, type="RAND_PLANAR", iteration=1)
{	tlog(2,"Initializing the graph")
	
	# init file name
	graph.folder <- file.path("data",type,paste0("n=",n),paste0("it=",iteration))
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
		if(type=="RAND_PLANAR")
		{	g <- graph.empty(n=n, directed=FALSE)									# create empty graph
			V(g)$x <- runif(vcount(g),min=-1,max=1)									# setup the node spatial positions
			V(g)$y <- runif(vcount(g),min=-1,max=1)
			g <- connect.triangulation(g)											# use triangulation to init the links
		}
		
		# create an Erdös-Rényi graph
		else if(type=="ERDOS_RENYI")
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



############################################################################
# Discretizes the original graph and record all the resulting graphs as 
# Graphml files for later use. A table summarizing the process is also
# recorded, also for later use. If the table already exists as a file, we 
# assime the graph files also exist, and load the table then return it.
#
# n: number of nodes.
# type: RAND_PLANAR (planar random graph) or ERDOS_RENYI (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# g: the original graph.
#
# returns: the table summarizing the discretization process.
############################################################################
init.disc.table <- function(n=5, type="RAND_PLANAR", iteration=1, g)
{	tlog(2,"Initializing the discretization table")
	
	# init file names
	it.folder <- file.path("data",type,paste0("n=",n),paste0("it=",iteration))
	disc.file <- file.path(it.folder,"discretizations.txt")
	
	# if the table already exists, we load it
	if(file.exists(disc.file))
	{	tlog(2,"The discretization table is already available: we load it from file \"",disc.file,"\"")
		disc.table <- as.matrix(read.table(file=disc.file,header=TRUE))
	}
	
	# other wise, we init it (which may take a while)
	else
	{	# init discretization table
		disc.table <- matrix(nrow=1,data=c(NA,n,0,0))
		colnames(disc.table) <- c("Granularity","Nodes","AverageSegmentation","SegmentationDuration")
		
		# init granularity values
		grans <- seq(from=max(E(g)$dist)/2,to=0.004,by=-0.001)#seq(from=0.10,to=0.004,by=-0.0001))
		
		# process each granularity value
		prev.size <- 0
		for(d in 1:length(grans))
		{	tlog(4,"Iteration "	,d,"/",length(grans)," granularity: ",grans[d])
			start.time <- Sys.time()
			g2 <- add.intermediate.nodes(g, granularity=grans[d])
			if(vcount(g2)!=prev.size)
			{	# update duration
				end.time <- Sys.time()
				duration <- difftime(end.time,start.time,units="s")
				g2$duration <- duration
				# update size
				size <- vcount(g2)
				g2$size <- size
				prev.size <- size
				# update granularity
				g2$granularity <- grans[d]
				# update average segmentation
				avgseg <- mean(E(g)$dist/grans[d])
				g2$avgseg <- avgseg
				
				# record graph for later
				graph.file <- file.path(it.folder,paste0("disc=",nrow(disc.table),".graphml"))
				write.graph(graph=g2,file=graph.file,format="graphml")
				
				# update discretization table
				row <- c(grans[d],size,avgseg,duration)
				disc.table <- rbind(disc.table,row)
				
				# display for debug
				tlog(6,"Number of nodes: ",size," - Duration: ",duration," s - Average segmentation: ",avgseg," - Recorded as \"",graph.file,"\"")
			}
			
			# reset to free memory
			g2 <- NULL
			gc()
			
			# record the current table
			write.table(disc.table,file=disc.file,col.names=TRUE,row.names=FALSE)
		}
	}
	
	return(disc.table)
}



############################################################################
############################################################################
process.continuous.straightness <- function(n=5, type="RAND_PLANAR", iteration=1, g)
{	tlog("Processing the continuous average straightness")
	it.folder <- file.path("data",type,paste0("n=",n),paste0("it=",iteration))
	
	# for the whole graph
	table.file <- file.path(it.folder,"continuous-graph.txt")
	# check if the table already exists, in which case it is simply loaded
	if(file.exists(table.file))
	{	tlog(2,"The results for the whole graph are already available: we load them from file \"",table.file,"\"")
		graph.table <- as.matrix(read.table(file=table.file,header=TRUE))
	}
	# otherwise, we process and record it 
	else
	{	tlog(2,"Processing average over the whole graph")
		start.time <- Sys.time()
			value <- mean.straightness.graph(graph=g)
		end.time <- Sys.time()
		duration <- difftime(end.time,start.time,units="s")
		tlog(4,"Continuous average straightness: ",value," - Duration: ",duration," s")
		gc()
		
		# record the data as a text file
		tlog(4,"Record results as file \"",table.file,"\"")
		graph.table <- matrix(c(value,duration),nrow=1)
		colnames(graph.table) <- c("Straightness","Duration")
		write.table(x=graph.table,file=table.file,row.names=FALSE,col.names=TRUE)
	}
	
	# process each node in the graph
	table.file <- file.path(it.folder,"continuous-nodes.txt")
	# check if the table already exists, in which case it is simply loaded
	if(file.exists(table.file))
	{	tlog(2,"The results for the nodes are already available: we load them from file \"",table.file,"\"")
		nodes.table <- as.matrix(read.table(file=table.file,header=TRUE))
	}
	# otherwise, we process and record it 
	else
	{	tlog(2,"Processing continuous average for each individual node")
		values <- c()
		durations <- c()
		for(i in 1:vcount(g))
		{	# process the straightness
			start.time <- Sys.time()
				value <- mean.straightness.nodes.graph(graph=g, u=i)
			end.time <- Sys.time()
			values <- c(values,value)
			duration <- difftime(end.time,start.time,units="s")
			durations <- c(durations,duration)
			tlog(4,"Continuous average straightness for node ",i,": ",value," - Duration: ",duration," s")
			gc()
		}
	
		# record the data as a text file
		tlog(2,"Record results as file \"",table.file,"\"")
		nodes.table <- cbind(values,durations)
		colnames(nodes.table) <- c("Straightness","Duration")
		write.table(x=nodes.table,file=table.file,row.names=FALSE,col.names=TRUE)
	}
	
	result <- list(graph=graph.table, nodes=nodes.table)
	return(result)
}



############################################################################
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
	}
	
	result <- list(graph=graph.tables, nodes=nodes.tables)
	return(result)
}



############################################################################
############################################################################
generate.rep.plots <- function(n=5, type="RAND_PLANAR", iteration=1, disc.table, cont.tables, disc.tables)
{	tlog("Generating plots and tables for the iteration ",iteration)
	it.folder <- file.path("data",type,paste0("n=",n),paste0("it=",iteration))
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
					pch=20,
					col=add.alpha("BLUE", 0.25)
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
############################################################################
generate.overall.plots <- function(n=10, type="RAND_PLANAR", discretizations, data.cont, data.disc)
{	tlog("Generating plots and tables for all the repetitions")
	folder <- file.path("data",type,paste0("n=",n))
	
	# collecting the data
	ngran <- nrow(data.disc[[1]]$graph)
	graph.cont.durations <- data.cont[[1]]$graph[1,"Duration"]
	nodes.cont.durations <- data.cont[[1]]$nodes[,"Duration"]
	graph.disc.durations <- data.disc[[1]]$graph[,"Duration"]
	graph.disc.differences <- data.disc[[1]]$graph[,"Difference"]
	nodes.disc.durations <- data.disc[[1]]$nodes$durations	#TODO pas le même nbre de discretizations
	nodes.disc.differences <- data.disc[[1]]$nodes$differences
	granularities <- discretizations[[1]][,"Granularity"]
	nodes <- discretizations[[1]][,"Nodes"]
	avgseg <- discretizations[[1]][,"AverageSegmentation"]
	for(r in 2:length(data.disc))
	{	# continuous results
		graph.cont.durations <- c(graph.cont.durations,data.cont[[r]]$graph[1,"Duration"])
		nodes.cont.durations <- c(nodes.cont.durations,data.cont[[r]]$nodes[,"Duration"])
		# discrete results
		graph.disc.durations <- c(graph.disc.durations,data.disc[[r]]$graph[,"Duration"])
		graph.disc.differences <- c(graph.disc.differences,data.disc[[r]]$graph[,"Difference"])
		nodes.disc.durations <- cbind(nodes.disc.durations,data.disc[[r]]$nodes$durations)
		nodes.disc.differences <- cbind(nodes.disc.differences,data.disc[[r]]$nodes$differences)
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
		
		# durations
		{	yaxis <- "durations"
			ylab <- "Time (s)"
			
			# graph plots
			plot.file <- file.path(folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(x=xvals, y=graph.disc.durations,
					pch=20,
					xlab=xlab, ylab=ylab,
					col=add.alpha("BLUE", 0.25),
					ylim=c(min(c(graph.disc.durations,graph.cont.durations)),max(c(graph.disc.durations,graph.cont.durations)))
			)
			for(j in 1:length(graph.cont.durations))
			{	lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(graph.cont.durations[j],2),
						col=add.alpha("RED", 0.25)
				)
			}
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			plot.file <- file.path(folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(x=rep(xvals,nrow(nodes.disc.durations)), y=c(t(nodes.disc.durations)),
					pch=20,
					xlab=xlab, ylab=ylab,
					col=add.alpha("BLUE", 0.25),
					ylim=c(min(c(nodes.disc.durations,nodes.cont.durations)),max(c(nodes.disc.durations,nodes.cont.durations)))
			)
			for(j in 1:length(nodes.cont.durations))
			{	lines(x=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
						y=rep(nodes.cont.durations[j],2),
						col=add.alpha("RED", 0.25)
				)
			}
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
		}
		
		# straightness differences
		{	yaxis <- "differences"
			ylab <- "Straightness difference"
			
			# graph plots
			plot.file <- file.path(folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(NULL,ylim=c(min(graph.disc.differences),max(graph.disc.differences)),
				xlim=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
				xlab=xlab, ylab=ylab
			)
			lines(x=c(min(xvals,na.rm=TRUE), max(xvals,na.rm=TRUE)), y=c(0,0), col="BLACK", lty=2)
			points(x=xvals, y=graph.disc.differences,
					pch=20,
					col=add.alpha("BLUE", 0.25)
			)
			dev.off()
			
			# node plots
			plot.file <- file.path(folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
			plot(NULL,ylim=c(min(nodes.disc.differences),max(nodes.disc.differences)),
					xlim=c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE)),
					xlab=xlab, ylab=ylab
			)
			lines(x=c(min(xvals,na.rm=TRUE), max(xvals,na.rm=TRUE)), y=c(0,0), col="BLACK", lty=2)
			points(x=rep(xvals,nrow(nodes.disc.differences)), y=c(t(nodes.disc.differences)),
					pch=20,
					col=add.alpha("BLUE", 0.25)
			)
			dev.off()
		}
	}
}


############################################################################
############################################################################
mymain <- function(n=5, type="RAND_PLANAR", repetitions=10)
{	gc()
	
	# process the specified number of repetitions
	data.disc <- list()
	data.cont <- list()
	discretizations <- list()
	for(r in 1:repetitions)
	{	# retrieve or create the graph
		g <- init.graph(n, type, iteration=r)
		
		# retrieve the discretization table exists
		disc.table <- init.disc.table(n, type, iteration=r, g)
		
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

mymain(n=10, type="RAND_PLANAR", repetitions=10)
