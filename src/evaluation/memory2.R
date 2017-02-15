############################################################################
# Alternative script to memory.R for evaluating memory usage, and to try 
# solving the bug met with R's debugger. Without success, though.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/evaluation/memory2.R")
############################################################################
source("src/evaluation/common.R")




data.folder <- "data"
eval.folder <- file.path(data.folder,"eval")




MEM_FILE <- "data/profiling.txt"	# R profiler output file (temporary)
MEM_SLEEP <- 2 						# sleep between to tries, in seconds
TRY_LIMIT <- 100					# max number of tries
DURATION_LIMIT <- 300 				# max duration for trying, in seconds

fake.zero <- 0.35					# used when plotting with a log scale
axis.brk <- 0.45
leg.offset <- 0.09




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
load.disc.table <- function(n=5, type="randplanar", iteration=1, g)
{	tlog(2,"Loading the discretization table")
	
	# init file names
	it.folder <- file.path(eval.folder,type,paste0("n=",n),paste0("it=",iteration))
	disc.file <- file.path(it.folder,"discretizations.txt")
	
	# get the table
	if(file.exists(disc.file))
		disc.table <- as.matrix(read.table(file=disc.file,header=TRUE))
	else
		stop("Discretization table not found. You must first execute the time.R script.")
	
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
# returns: the list of tables containing the memory usage.
############################################################################
process.continuous.straightness <- function(n=5, type="randplanar", iteration=1, g)
{	tlog("Processing the continuous average straightness")
	it.folder <- file.path(eval.folder,type,paste0("n=",n),paste0("it=",iteration))
	graph.file <- file.path(eval.folder,type,paste0("n=",n),paste0("it=",iteration),paste0("disc=0",".graphml"))
	
	# for the whole graph
	table.file <- file.path(it.folder,"continuous-graph.txt")
	if(file.exists(table.file))
		graph.table <- as.matrix(read.table(file=table.file,header=TRUE))
	else
		stop("Whole graph table not found (\"",table.file,"\"). You must first execute the time.R script.")
	# check if memory usage was already measured
#	if("Memory" %in% colnames(graph.table) && !is.na(graph.table[1,"Memory"]))
#		tlog(2,"The memory usage for the graph is already available (",graph.table[1,"Memory"]," MB): nothing to do here")
#	else
	{	tlog(2,"Processing average over the whole graph")
		cmd <- paste0("Rscript --vanilla ",file.path(getwd(),"src","evaluation","memory-aux.R")," ",getwd()," ",graph.file," ",0," ",0)
		memuse <- NA
		again <- TRUE
		while(again)
		{	again <- FALSE
			invisible(system(command=cmd, intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL))
			mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[[1]],
					error=function(e){again <<- TRUE})
			if(again)
				tlog(4,"Error while trying to process memory usage. Trying again.")
			else
				memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
		}
		tlog(4,"Continuous average straightness memory usage: ",memuse," MB")
		
		# record the data as a text file
		tlog(4,"Record results as file \"",table.file,"\"")
		if("Memory" %in% colnames(graph.table))
			graph.table[1,"Memory"] <- memuse
		else
		{	graph.table <- cbind(graph.table,memuse)
			colnames(graph.table)[ncol(graph.table)] <- "Memory"
		}
		write.table(x=graph.table,file=table.file,row.names=FALSE,col.names=TRUE)
	}
	
	# process each node in the graph
	table.file <- file.path(it.folder,"continuous-nodes.txt")
	if(file.exists(table.file))
		nodes.table <- as.matrix(read.table(file=table.file,header=TRUE))
	else
		stop("Individual nodes table not found (\"",table.file,"\"). You must first execute the time.R script.")
	tlog(2,"Processing continuous average for each individual node")
	memuses <- c()
	for(i in 1:vcount(g))
	{	# check if memory usage was already measured
#		if("Memory" %in% colnames(nodes.table) && !is.na(nodes.table[i,"Memory"]))
#		{	tlog(4,"The memory usage for the node is already available (",nodes.table[i,"Memory"]," MB): nothing to do here")
#			memuses <- c(memuses,nodes.table[i,"Memory"])
#		}
#		else
		{	cmd <- paste0("Rscript --vanilla ",file.path(getwd(),"src","evaluation","memory-aux.R")," ",getwd()," ",graph.file," ",0," ",i)
			memuse <- NA
			again <- TRUE
			while(again)
			{	again <- FALSE
				invisible(system(command=cmd, intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL))
				mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[[1]],
						error=function(e){again <<- TRUE})
				if(again)
					tlog(4,"Error while trying to process memory usage. Trying again.")
				else
					memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
			}
			memuses <- c(memuses,memuse)
			tlog(4,"Continuous average straightness memory usage for node ",i,": ",memuse," MB")
		}
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
# type: randplanar (planar random graph) or erdosrenyi (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# g: the original graph.
#
# returns: the list of lists of tables containing the memory usage.
############################################################################
process.discrete.straightness <- function(n=5, type="randplanar", iteration=1, g, cont.tables)
{	tlog("Processing the discrete approximation of the average straightness")
	it.folder <- file.path(eval.folder,type,paste0("n=",n),paste0("it=",iteration))
	graph.file <- file.path(eval.folder,type,paste0("n=",n),paste0("it=",iteration),paste0("disc=0",".graphml"))
	
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
		
		# for the whole graph
		table.file <- file.path(it.folder,paste0("disc=",d,"-graph.txt"))
		if(file.exists(table.file))
			graph.table <- as.matrix(read.table(file=table.file,header=TRUE))
		else
			stop("Whole graph table not found (\"",table.file,"\"). You must first execute the time.R script.")
		# check if memory usage was already measured
#		if("Memory" %in% colnames(graph.table) && !is.na(graph.table[1,"Memory"]))
#			tlog(4,"The memory usage for the graph is already available (",graph.table[1,"Memory"]," MB): nothing to do here")
#		else
		{	tlog(4,"Processing average over the whole graph")
			memuse <- NA
			if(d>0)
			{	cmd <- paste0("Rscript --vanilla ",file.path(getwd(),"src","evaluation","memory-aux.R")," ",getwd()," ",graph.file," ",disc.table[d+1,"Granularity"]," ",0)
				again <- TRUE
				tries <- 0
				start.time <- Sys.time()
				while(again)
				{	again <- FALSE
#print(cmd);stop()					
					invisible(system(command=cmd, intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL))
					if(tries==TRY_LIMIT-1) Sys.sleep(MEM_SLEEP)
					mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[[1]],
							error=function(e){again <<- TRUE})
					end.time <- Sys.time()
					duration <- difftime(end.time,start.time,units="s")
					if(again && tries<TRY_LIMIT && duration<DURATION_LIMIT)
					{	tries <- tries + 1
						tlog(4,"Error while trying to process memory usage. Trying again (",tries,"/",TRY_LIMIT,").")
#print(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1))						
					}
					else if(again && (tries>=TRY_LIMIT || duration>=DURATION_LIMIT))
					{	tlog(4,"Error while trying to process memory usage. Too many tries or too long: giving up (",tries,"/",TRY_LIMIT," tries -- ",duration,"/",DURATION_LIMIT," seconds).")
						again <- FALSE
					}
					else
						memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
				}
			}
			tlog(6,"Discrete average straightness memory usage: ",memuse," MB")
			
			# record the data as a text file
			tlog(4,"Record results as file \"",table.file,"\"")
			if("Memory" %in% colnames(graph.table))
				graph.table[1,"Memory"] <- memuse
			else
			{	graph.table <- cbind(graph.table,memuse)
				colnames(graph.table)[ncol(graph.table)] <- "Memory"
			}
			write.table(x=graph.table,file=table.file,row.names=FALSE,col.names=TRUE)
		}
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
		{	# check if memory usage was already measured
#			if("Memory" %in% colnames(nodes.table) && !is.na(nodes.table[i,"Memory"]))
#			{	tlog(4,"The memory usage for the node is already available (",nodes.table[i,"Memory"]," MB): nothing to do here")
#				memuses <- c(memuses,nodes.table[i,"Memory"])
#			}
#			else
			{	memuse <- NA
				if(d>0)
				{	cmd <- paste0("Rscript --vanilla ",file.path(getwd(),"src","evaluation","memory-aux.R")," ",getwd()," ",graph.file," ",disc.table[d+1,"Granularity"]," ",i)
					again <- TRUE
					tries <- 0
					start.time <- Sys.time()
					while(again)
					{	again <- FALSE
						invisible(system(command=cmd, intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL))
						if(tries==TRY_LIMIT-1) Sys.sleep(MEM_SLEEP)
						mem.stats <- tryCatch(summaryRprof(MEM_FILE, memory="stats", diff=FALSE, index=1)[[1]],
								error=function(e){again <<- TRUE})
						end.time <- Sys.time()
						duration <- difftime(end.time,start.time,units="s")
						if(again && tries<TRY_LIMIT && duration<DURATION_LIMIT)
						{	tries <- tries + 1
							tlog(4,"Error while trying to process memory usage. Trying again (",tries,"/",TRY_LIMIT,").")
						}
						else if(again && (tries>=TRY_LIMIT || duration>=DURATION_LIMIT))
						{	tlog(4,"Error while trying to process memory usage. Too many tries: giving up (",tries,"/",TRY_LIMIT,").")
							again <- FALSE
						}
						else
							memuse <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
					}
				}
				memuses <- c(memuses,memuse)
				tlog(6,"Discrete average straightness memory usage for node ",i,": ",memuse," MB")
			}
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
#		   the previously processed memory usages.
############################################################################
generate.rep.plots <- function(n=5, type="randplanar", iteration=1, disc.table, cont.tables, disc.tables)
{	tlog("Generating plots and tables for the iteration ",iteration)
	it.folder <- file.path(eval.folder,type,paste0("n=",n),paste0("it=",iteration))
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
			x.lim <- c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE))
			log.axes <- ""
			mem.inset <- 0.03
			mem.leg.pos <- "topleft"
		}
		else if(xaxis=="avgseg")
		{	xlab <- "Average segmentation"
			xvals <- disc.table[,"AverageSegmentation"]
			xvals.zero <- which(xvals==0)
			xvals.positive <- which(xvals>0)
			#x.lim <- c(min(xvals[xvals.positive],na.rm=TRUE),max(xvals[xvals.positive],na.rm=TRUE))
			x.lim <- c(fake.zero,max(xvals[xvals.positive],na.rm=TRUE))
			log.axes <- "x"
			mem.inset <- c(0.03+leg.offset, 0.03)
			mem.leg.pos <- "topleft"
		} 
		else if(xaxis=="granularity")
		{	xlab <- "Granularity"
			xvals <- disc.table[,"Granularity"]
			x.lim <- c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE))
			log.axes <- ""
			mem.inset <- 0.03
			mem.leg.pos <- "topright"
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
				plot(NULL,
						xlab=xlab, ylab=ylab,
						log=log.axes,
						xaxt=if(log.axes=="x") "n" else "s",
						xlim=x.lim,
						ylim=c(min(c(graph.yvals,graph.cont.val),na.rm=TRUE),max(c(graph.yvals,graph.cont.val),na.rm=TRUE))
				)
				if(log.axes=="x")
				{	aty <- c(fake.zero, axTicks(1)[-1])
					x.labels <- c(0,aty[-1])
					axis(1,at=aty,labels=x.labels)
				}
				lines(x=x.lim,
						y=rep(graph.cont.val,2),
						col="RED"
				)
				points(x=if(log.axes=="x") xvals[xvals.positive] else xvals, 
						y=if(log.axes=="x") graph.yvals[xvals.positive] else graph.yvals,
						col="BLUE",
				)
				if(log.axes=="x") points(x=rep(fake.zero,length(xvals.zero)),
						y=graph.yvals[xvals.zero],
						col="BLUE"
					)
				if(log.axes=="x") axis.break2(1,axis.brk,style="gap")
				legend(x="bottomright",legend=c("Discrete average","Continuous average"),
						inset=mem.inset,
						fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			plot.file <- file.path(it.folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(NULL,
						xlab=xlab, ylab=ylab,
						log=log.axes,
						xaxt=if(log.axes=="x") "n" else "s",
						xlim=x.lim,
						ylim=c(min(c(nodes.yvals,nodes.cont.vals),na.rm=TRUE),max(c(nodes.yvals,nodes.cont.vals),na.rm=TRUE))
				)
				if(log.axes=="x")
				{	aty <- c(fake.zero, axTicks(1)[-1])
					x.labels <- c(0,aty[-1])
					axis(1,at=aty,labels=x.labels)
				}
				for(j in 1:length(nodes.cont.vals))
				{	lines(x=x.lim,
							y=rep(nodes.cont.vals[j],2),
							col="RED"#add.alpha("RED", 0.25)
					)
				}
				points(x=if(log.axes=="x") rep(xvals[xvals.positive],nrow(nodes.yvals)) else rep(xvals,nrow(nodes.yvals)), 
						y=if(log.axes=="x") c(t(nodes.yvals[,xvals.positive])) else c(t(nodes.yvals)),
						col="BLUE"#add.alpha("BLUE", 0.25),pch=20,
				)
				if(log.axes=="x") points(x=rep(fake.zero,length(xvals.zero)*nrow(nodes.yvals)),
						y=c(t(nodes.yvals[,xvals.zero])),
						col="BLUE"
					)
				if(log.axes=="x") axis.break2(1,axis.brk,style="gap")
				legend(x="bottomright",legend=c("Discrete average","Continuous average"),
						inset=mem.inset,
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
# type: randplanar (planar random graph) or erdosrenyi (Erdös-Rényi random
#		graph).
# iteration: number of the iteration (cf. number of repetitions).
# discretizations: list of discretization tables.
# data.cont: list of lists of tables for the continuous average straightness.
# data.disc: list of list of tables for the discrete average straightness.
############################################################################
generate.overall.plots <- function(n=10, type="randplanar", discretizations, data.cont, data.disc)
{	tlog("Generating plots and tables for all the repetitions")
	folder <- file.path(eval.folder,type,paste0("n=",n))
	
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
			x.lim <- c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE))
			log.axes <- ""
			mem.inset <- 0.03
			mem.leg.pos <- "topleft"
		}
		else if(xaxis=="avgseg")
		{	xlab <- "Average segmentation"
			xvals <- avgseg
			xvals.zero <- which(xvals==0)
			xvals.positive <- which(xvals>0)
			#x.lim <- c(min(xvals[xvals.positive],na.rm=TRUE),max(xvals[xvals.positive],na.rm=TRUE))
			x.lim <- c(fake.zero,max(xvals[xvals.positive],na.rm=TRUE))
			log.axes <- "x"
			mem.inset <- c(0.03+leg.offset, 0.03)
			mem.leg.pos <- "topleft"
		} 
		else if(xaxis=="granularity")
		{	xlab <- "Granularity"
			xvals <- granularities
			x.lim <- c(min(xvals,na.rm=TRUE),max(xvals,na.rm=TRUE))
			log.axes <- ""
			mem.inset <- 0.03
			mem.leg.pos <- "topright"
		} 
		
		# memory
		{	yaxis <- "memory"
			ylab <- "Memory (MB)"
			
			# graph plots
			plot.file <- file.path(folder,paste0("graph-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(NULL,
						xlab=xlab, ylab=ylab,
						log=log.axes,
						xaxt=if(log.axes=="x") "n" else "s",
						xlim=x.lim,
						ylim=c(min(c(graph.disc.memory,graph.cont.memory),na.rm=TRUE),max(c(graph.disc.memory,graph.cont.memory),na.rm=TRUE))
				)
				if(log.axes=="x")
				{	aty <- c(fake.zero, axTicks(1)[-1])
					x.labels <- c(0,aty[-1])
					axis(1,at=aty,labels=x.labels)
				}
				for(j in 1:length(graph.cont.memory))
				{	lines(x=x.lim,
							y=rep(graph.cont.memory[j],2),
							col="RED"#add.alpha("RED", 0.25)
					)
				}
				points(x=if(log.axes=="x") xvals[xvals.positive] else xvals, 
						y=if(log.axes=="x") graph.disc.memory[xvals.positive] else graph.disc.memory,
						col="BLUE"#add.alpha("BLUE", 0.25),pch=20,
				)
				if(log.axes=="x") points(x=rep(fake.zero,length(xvals.zero)),
						y=graph.disc.memory[xvals.zero],
						col="BLUE"
					)
				if(log.axes=="x") axis.break2(1,axis.brk,style="gap")
				legend(x="bottomright",legend=c("Discrete average","Continuous average"),
						inset=mem.inset,
						fill=c("BLUE","RED"))
			dev.off()
			
			# node plots
			plot.file <- file.path(folder,paste0("nodes-",yaxis,"-vs-",xaxis,".pdf"))
			pdf(file=plot.file)
				plot(NULL,
						xlab=xlab, ylab=ylab,
						log=log.axes,
						xaxt=if(log.axes=="x") "n" else "s",
						xlim=x.lim,
						ylim=c(min(c(nodes.disc.memory,nodes.cont.memory),na.rm=TRUE),max(c(nodes.disc.memory,nodes.cont.memory),na.rm=TRUE))
				)
				if(log.axes=="x")
				{	aty <- c(fake.zero, axTicks(1)[-1])
					x.labels <- c(0,aty[-1])
					axis(1,at=aty,labels=x.labels)
				}
				for(j in 1:length(nodes.cont.memory))
				{	lines(x=x.lim,
							y=rep(nodes.cont.memory[j],2),
							col="RED"#add.alpha("RED", 0.25)
					)
				}
				points(x=if(log.axes=="x") rep(xvals[,xvals.positive],nrow(nodes.disc.memory)) else rep(xvals,nrow(nodes.disc.memory)), 
						y=if(log.axes=="x") c(t(nodes.disc.memory[xvals.positive])) else c(t(nodes.disc.memory)),
						col="BLUE"#add.alpha("BLUE", 0.25),pch=20,
				)
				if(log.axes=="x") points(x=rep(fake.zero,length(xvals.zero)*nrow(nodes.disc.memory)),
						y=c(t(nodes.disc.memory[xvals.zero])),
						col="BLUE"
					)
				if(log.axes=="x") axis.break2(1,axis.brk,style="gap")
				legend(x="bottomright",legend=c("Discrete average","Continuous average"),
						inset=mem.inset,
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
# type: type of graph, randplanar (planar random graph) or erdosrenyi 
#		(Erdös-Rényi random graph).
# repetitions: number of instances of the graph to generate and process.
############################################################################
monitor.memory <- function(n=5, type="randplanar", repetitions=10)
{	# process the specified number of repetitions
	data.disc <- list()
	data.cont <- list()
	discretizations <- list()
	for(r in 1:repetitions)
	{	# retrieve or create the graph
		g <- init.graph(n, type, iteration=r, folder=eval.folder)
		
		# retrieve the discretization table
		disc.table <- load.disc.table(n, type, iteration=r, g)
		
		# deal with the continuous version
		cont.tables <- process.continuous.straightness(n, type, iteration=r,g)
		
		# deal with the discrete version
		disc.tables <- process.discrete.straightness(n, type, iteration=r, g, cont.tables)
		
		# generate plots
		data.disc[[r]] <- generate.rep.plots(n, type, iteration=r, disc.table, cont.tables, disc.tables)
		
		data.cont[[r]] <- cont.tables
		discretizations[[r]] <- disc.table
	}
	
	# generate plots using the repetitions
	generate.overall.plots(n, type, discretizations, data.cont, data.disc)
}

monitor.memory(n=10, type="randplanar", repetitions=10)
#monitor.memory(n=25, type="randplanar", repetitions=10)
#monitor.memory(n=50, type="randplanar", repetitions=10)
#monitor.memory(n=100, type="randplanar", repetitions=10)
