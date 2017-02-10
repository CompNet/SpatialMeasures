############################################################################
# Generate the time-related plots from the article cited in the project
# readme file. This script must be executed before memory.R, which takes
# advantage of certain outputs of this script. It focuses on the random
# graphs, the road networks being processed with the script urban.R.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/evaluation/time.R")
############################################################################
source("src/evaluation/common.R")




data.folder <- "data"
urban.folder <- file.path(data.folder,"eval")




############################################################################
# Discretizes the original graph and record all the resulting graphs as 
# Graphml files for later use. A table summarizing the process is recorded
# too, also for later use. If the table already exists as a file, we assume 
# the graph files also exist, and load the table then return it.
#
# n: number of nodes.
# type: randplanar (planar random graph) or erdosrenyi (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# g: the original graph.
#
# returns: the table summarizing the discretization process.
############################################################################
init.disc.table <- function(n=5, type="randplanar", iteration=1, g)
{	tlog(2,"Initializing the discretization table")
	
	# init file names
	it.folder <- file.path(urban.folder,type,paste0("n=",n),paste0("it=",iteration))
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
# Processes the continuous version of the average straightness, for the
# specified graphs.
#
# n: number of nodes.
# type: randplanar (planar random graph) or erdosrenyi (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# g: the original graph.
#
# returns: the list of tables containing the average straightness and corresponding
#		   durations values for the graph ($graph) and each node ($nodes).
############################################################################
process.continuous.straightness <- function(n=5, type="randplanar", iteration=1, g)
{	tlog("Processing the continuous average straightness")
	it.folder <- file.path(urban.folder,type,paste0("n=",n),paste0("it=",iteration))
	
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
# Processes the discrete version of the average straightness, for the
# specified graphs.
#
# n: number of nodes.
# type: randplanar (planar random graph) or erdosrenyi (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# g: the original graph.
#
# returns: the list of lists of tables containing the average straightness 
#		   and corresponding durations for the graph ($graph) and each node 
#		   ($nodes). Each element in the sublists correspond to a given 
#		   granularity in the discretization process.
############################################################################
process.discrete.straightness <- function(n=5, type="randplanar", iteration=1, g, cont.tables)
{	tlog("Processing the discrete approximation of the average straightness")
	it.folder <- file.path(urban.folder,type,paste0("n=",n),paste0("it=",iteration))
	
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
	cnames <- c("Straightness","Difference","Duration")
	colnames(graph.all) <- cnames
	for(d in 0:(nrow(disc.table)-1))
	{	#print(c(disc.tables$graph[[as.character(d)]])) # debug
		graph.all[d+1,] <- c(disc.tables$graph[[as.character(d)]])[1:3]
	}
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
			x.lim <- c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE))
			log.axes <- "" 
			str.inset <- 0.03
			dur.inset <- 0.03
			str.leg.pos <- "topright"
			dur.leg.pos <- "topleft"
		}
		else if(xaxis=="avgseg")
		{	xlab <- "Average segmentation"
			xvals <- disc.table[,"AverageSegmentation"]
			xvals.zero <- which(xvals==0)
			xvals.positive <- which(xvals>0)
			#x.lim <- c(min(xvals[xvals.positive],na.rm=TRUE),max(xvals[xvals.positive],na.rm=TRUE))
			x.lim <- c(0.5,max(xvals[xvals.positive],na.rm=TRUE))
			log.axes <- "x" 
			str.inset <- 0.03
			dur.inset <- c(0.10, 0.03)
			str.leg.pos <- "topright"
			dur.leg.pos <- "topleft"
		} 
		else if(xaxis=="granularity")
		{	xlab <- "Granularity"
			xvals <- disc.table[,"Granularity"]
			x.lim <- c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE))
			log.axes <- "" 
			str.inset <- 0.03
			dur.inset <- 0.03
			str.leg.pos <- "topleft"
			dur.leg.pos <- "topright"
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
				plot(x=if(log.axes=="x") xvals[xvals.positive] else xvals, 
						y=if(log.axes=="x") graph.yvals[xvals.positive] else graph.yvals,
						xlab=xlab, ylab=ylab,
						col="BLUE",
						log=log.axes,
						xlim=x.lim,
						ylim=c(min(c(graph.yvals,graph.cont.val)),max(c(graph.yvals,graph.cont.val)))
				)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)),
						y=graph.yvals[xvals.zero],
						col="BLUE"
				)
				lines(x=x.lim,
						y=rep(graph.cont.val,2),
						col="RED"
				)
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
				legend(x=str.leg.pos,legend=c("Discrete average","Continuous average"),
						inset=str.inset,
						fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			for(j in 1:length(nodes.cont.vals))
			{	plot.file <- file.path(it.folder,paste0("node=",j,"-",yaxis,"-vs-",xaxis,".pdf"))
				pdf(file=plot.file)
					plot(x=if(log.axes=="x") xvals[xvals.positive] else xvals, 
							y=if(log.axes=="x") nodes.yvals[j,xvals.positive] else nodes.yvals[j,],
							xlab=xlab, ylab=ylab,
							col="BLUE",
							log=log.axes,
							xlim=x.lim,
							ylim=c(min(c(nodes.yvals[j,],nodes.cont.vals[j])),max(c(nodes.yvals[j,],nodes.cont.vals[j])))
					)
					if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)),
								y=nodes.yvals[j,xvals.zero],
								col="BLUE"
						)
					lines(x=x.lim,
							y=rep(nodes.cont.vals[j],2),
							col="RED"
					)
					if(log.axes=="x") axis.break2(1,0.6,style="gap")
					legend(x=str.leg.pos,legend=c("Discrete average","Continuous average"),
							inset=str.inset,
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
				plot(x=if(log.axes=="x") xvals[xvals.positive] else xvals, 
						y=if(log.axes=="x") graph.yvals[xvals.positive] else graph.yvals,
						xlab=xlab, ylab=ylab,
						col="BLUE",
						log=log.axes,
						xlim=x.lim,
						ylim=c(min(c(graph.yvals,graph.cont.val)),max(c(graph.yvals,graph.cont.val)))
				)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)),
							y=graph.yvals[xvals.zero],
							col="BLUE"
					)
				lines(x=x.lim,
						y=rep(graph.cont.val,2),
						col="RED"
				)
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
				legend(x=dur.leg.pos,legend=c("Discrete average","Continuous average"),
						inset=dur.inset,
						fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			plot.file <- file.path(it.folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(x=if(log.axes=="x") rep(xvals[xvals.positive],nrow(nodes.yvals)) else rep(xvals,nrow(nodes.yvals)), 
						y=if(log.axes=="x") c(t(nodes.yvals[,xvals.positive])) else c(t(nodes.yvals)),
						col="BLUE",#add.alpha("BLUE", 0.25),pch=20,
						xlab=xlab, ylab=ylab,
						log=log.axes,
						xlim=x.lim,
						ylim=c(min(c(nodes.yvals,nodes.cont.vals)),max(c(nodes.yvals,nodes.cont.vals)))
				)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)*nrow(nodes.yvals)),
							y=c(t(nodes.yvals[,xvals.zero])),
							col="BLUE"
					)
				for(j in 1:length(nodes.cont.vals))
				{	lines(x=x.lim,
							y=rep(nodes.cont.vals[j],2),
							col="RED"#add.alpha("RED", 0.25)
					)
				}
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
				legend(x=dur.leg.pos,legend=c("Discrete average","Continuous average"),
						inset=dur.inset,
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
				plot(NULL,
						xlim=c(x.lim[1],x.lim[2]*1.1),
						ylim=c(min(graph.yvals),max(graph.yvals)),
						log=log.axes,
						xlab=xlab, ylab=ylab
				)
				lines(x=c(x.lim[1],x.lim[2]*1.1), 
						y=c(0,0), 
						col="BLACK", lty=2)
				points(x=if(log.axes=="x") xvals[xvals.positive] else xvals, 
						y=if(log.axes=="x") graph.yvals[xvals.positive] else graph.yvals,
						col="BLUE"
				)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)),
							y=graph.yvals[xvals.zero],
							col="BLUE"
					)
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
			dev.off()
			
			# node plots
			plot.file <- file.path(it.folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(NULL,
						xlim=c(x.lim[1],x.lim[2]*1.1),
						ylim=c(min(nodes.yvals),max(nodes.yvals)),
						log=log.axes,
						xlab=xlab, ylab=ylab
				)
				lines(x=c(x.lim[1],x.lim[2]*1.1), 
						y=c(0,0), 
						col="BLACK", lty=2)
				points(x=if(log.axes=="x") rep(xvals[xvals.positive],nrow(nodes.yvals)) else rep(xvals,nrow(nodes.yvals)), 
						y=if(log.axes=="x") c(t(nodes.yvals[,xvals.positive])) else c(t(nodes.yvals)),
						col="BLUE"#add.alpha("BLUE", 0.25),pch=20					
				)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)*nrow(nodes.yvals)),
							y=c(t(nodes.yvals[,xvals.zero])),
							col="BLUE"
					)
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
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
# Generates the plots for all iteration.
#
# n: number of nodes.
# type: randplanar (planar random graph) or erdosrenyi (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# discretizations: list of discretization tables.
# data.cont: list of lists of tables for the continuous average straightness.
# data.disc: list of list of tables for the discrete average straightness.
############################################################################
generate.overall.plots <- function(n=10, type="randplanar", discretizations, data.cont, data.disc)
{	tlog("Generating plots and tables for all the repetitions")
	folder <- file.path(urban.folder,type,paste0("n=",n))
	
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
			x.lim <- c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE))
			log.axes <- ""
			dur.inset <- 0.03
			dur.leg.pos <- "topleft"
		}
		else if(xaxis=="avgseg")
		{	xlab <- "Average segmentation"
			xvals <- avgseg
			xvals.zero <- which(xvals==0)
			xvals.positive <- which(xvals>0)
			#x.lim <- c(min(xvals[xvals.positive],na.rm=TRUE),max(xvals[xvals.positive],na.rm=TRUE))
			x.lim <- c(0.5,max(xvals[xvals.positive],na.rm=TRUE))
			log.axes <- "x"
			dur.inset <- c(0.10, 0.03)
			dur.leg.pos <- "topleft"
		} 
		else if(xaxis=="granularity")
		{	xlab <- "Granularity"
			xvals <- granularities
			x.lim <- c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE))
			log.axes <- "" 
			dur.inset <- 0.03
			dur.leg.pos <- "topright"
		} 
		
		# durations
		{	yaxis <- "durations"
			ylab <- "Time (s)"
			
			# graph plots
			plot.file <- file.path(folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(x=if(log.axes=="x") xvals[xvals.positive] else xvals, 
						y=if(log.axes=="x") graph.disc.durations[xvals.positive] else graph.disc.durations,
						col="BLUE",#add.alpha("BLUE", 0.25),pch=20,
						xlab=xlab, ylab=ylab,
						log=log.axes,
						xlim=x.lim,
						ylim=c(min(c(graph.disc.durations,graph.cont.durations)),max(c(graph.disc.durations,graph.cont.durations)))
				)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)),
							y=graph.disc.durations[xvals.zero],
							col="BLUE"
					)
				for(j in 1:length(graph.cont.durations))
				{	lines(x=x.lim,
							y=rep(graph.cont.durations[j],2),
							col="RED"#add.alpha("RED", 0.25)
					)
				}
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
				legend(x=dur.leg.pos,legend=c("Discrete average","Continuous average"),
						inset=dur.inset,
						fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			plot.file <- file.path(folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(x=if(log.axes=="x") rep(xvals[xvals.positive],nrow(nodes.disc.durations)) else rep(xvals,nrow(nodes.disc.durations)), 
						y=if(log.axes=="x") c(t(nodes.disc.durations[,xvals.positive])) else c(t(nodes.disc.durations)),
						col="BLUE",#add.alpha("BLUE", 0.25),pch=20,
						xlab=xlab, ylab=ylab,
						log=log.axes,
						xlim=x.lim,
						ylim=c(min(c(nodes.disc.durations,nodes.cont.durations)),max(c(nodes.disc.durations,nodes.cont.durations)))
				)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)*nrow(nodes.disc.durations)),
							y=c(t(nodes.disc.durations[,xvals.zero])),
							col="BLUE"
					)
				for(j in 1:length(nodes.cont.durations))
				{	lines(x=x.lim,
							y=rep(nodes.cont.durations[j],2),
							col="RED"#add.alpha("RED", 0.25)
					)
				}
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
				legend(x=dur.leg.pos,legend=c("Discrete average","Continuous average"),
						inset=dur.inset,
						fill=c("BLUE","RED"))
				dev.off()
		}
		
		# straightness differences
		{	yaxis <- "differences"
			ylab <- "Straightness difference"
			
			# graph plots
			plot.file <- file.path(folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(NULL,
					xlim=c(x.lim[1],x.lim[2]*1.1),
					ylim=c(min(graph.disc.differences),max(graph.disc.differences)),
					log=log.axes,
					xlab=xlab, ylab=ylab
				)
				lines(x=c(x.lim[1],x.lim[2]*1.1), 
						y=c(0,0), 
						col="BLACK", lty=2
				)
				points(x=if(log.axes=="x") xvals[xvals.positive] else xvals, 
						y=if(log.axes=="x") graph.disc.differences[xvals.positive] else graph.disc.differences,
						col="BLUE"#add.alpha("BLUE", 0.25),pch=20
					)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)),
							y=graph.disc.differences[xvals.zero],
							col="BLUE"
					)
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
			dev.off()
			
			# node plots
			plot.file <- file.path(folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(NULL,
						xlim=c(x.lim[1],x.lim[2]*1.1),
						ylim=c(min(nodes.disc.differences),max(nodes.disc.differences)),
						log=log.axes,
						xlab=xlab, ylab=ylab
				)
				lines(x=c(x.lim[1],x.lim[2]*1.1), 
						y=c(0,0), 
						col="BLACK", lty=2
				)
				points(x=if(log.axes=="x") rep(xvals[xvals.positive],nrow(nodes.disc.differences)) else rep(xvals,nrow(nodes.disc.differences)), 
						y=if(log.axes=="x") c(t(nodes.disc.differences[,xvals.positive])) else c(t(nodes.disc.differences)),
						col="BLUE"#add.alpha("BLUE", 0.25),pch=20
					)
				if(log.axes=="x") points(x=rep(0.5,length(xvals.zero)*nrow(nodes.disc.differences)),
							y=c(t(nodes.disc.differences[,xvals.zero])),
							col="BLUE"
					)
				if(log.axes=="x") axis.break2(1,0.6,style="gap")
			dev.off()
		}
	}
}



############################################################################
# Process the time usage stats on a collection of randomly generated networks,
# and generates the corresponding tables and plots.
#
# n: number of nodes in the graph.
# type: type of graph, randplanar (planar random graph) or erdosrenyi 
#		(Erdös-Rényi random graph).
# repetitions: number of instances of the graph to generate and process.
############################################################################
monitor.time <- function(n=5, type="randplanar", repetitions=10)
{	gc()
	
	# process the specified number of repetitions
	data.disc <- list()
	data.cont <- list()
	discretizations <- list()
	for(r in 1:repetitions)
#	r <- 1
	{	# retrieve or create the graph
		g <- init.graph(n, type, iteration=r, folder=urban.folder)
		
		# retrieve the discretization table, or init it if it doesn't exist
		disc.table <- init.disc.table(n, type, iteration=r, g)
		
		# deal with the continuous version
		cont.tables <- process.continuous.straightness(n, type, iteration=r, g)
		
		# deal with the discrete version
		disc.tables <- process.discrete.straightness(n, type, iteration=r, g, cont.tables)
		
		# generate plots
		data.disc[[r]] <- generate.rep.plots(n, type, iteration=r, disc.table, cont.tables, disc.tables)

		data.cont[[r]] <- cont.tables
		discretizations[[r]] <- disc.table
		gc()
	}
	
	# generate plots using the repetitions
	generate.overall.plots(n, type, discretizations, data.cont, data.disc)
}

#monitor.time(n=10, type="randplanar", repetitions=10)
monitor.time(n=25, type="randplanar", repetitions=10)
#monitor.time(n=50, type="randplanar", repetitions=10)
#monitor.time(n=100, type="randplanar", repetitions=10)


#setwd("~/eclipse/workspaces/Networks/SpatialMeasures");source("src/evaluation/time.R")
#TODO plot: change the 0.5 label to a zero, maybe remove the decimal parts
#TODO fix the log thingy for n=25
