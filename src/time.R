############################################################################
# Generate the time-related plots from the article cited in the project
# readme file.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
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
		display.model(g, large=TRUE, filename=graph.filename, out.folder=graph.folder, export=TRUE, formats=c("pdf",NA))
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
		disc.table <- read.table(file=disc.file,header=TRUE)
	
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
				graph.file <- file.path(it.folder,paste0("disc=",nrow(disc.table)-1,".graphml"))
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
	tlog(2,"Processing average over the whole graph")
	start.time <- Sys.time()
	pus.value <- mean.straightness.graph(graph=g)
	end.time <- Sys.time()
	pus.duration <- difftime(end.time,start.time,units="s")
	tlog(4,"Continuous average straightness: ",pus.value," - Duration: ",pus.duration," s")
	gc()
	
	# record the data as a text file
	table.file <- file.path(it.folder,"continuous-graph.txt")
	tlog(4,"Record results as file \"",table.file,"\"")
	data <- matrix(c(pus.value,pus.duration),ncol=2)
	colnames(data) <- c("Straightness","Duration")
	write.table(x=data,file=table.file,row.names=FALSE,col.names=TRUE)
	
	# for each node in the graph
	tlog(2,"Processing average for each individual node")
	pus.values <- c()
	pus.durations <- c()
	for(i in 1:vcount(g))
	{	# process the straightness
		start.time <- Sys.time()
		pus.value <- mean.straightness.nodes.graph(graph=g, u=i)
		pus.values <- c(pus.values,pus.value)
		end.time <- Sys.time()
		pus.duration <- difftime(end.time,start.time,units="s")
		pus.durations <- c(pus.durations,pus.duration)
		tlog(4,"Continuous average straightness for node ",i,": ",pus.value," - Duration: ",pus.duration," s")
		gc()
	}

	# record the data as a text file
	table.file <- file.path(it.folder,"continuous-nodes.txt")
	tlog(4,"Record results as file \"",table.file,"\"")
	data <- cbind(pus.values,pus.durations)
	colnames(data) <- c("Straightness","Duration")
	write.table(x=data,file=table.file,row.names=FALSE,col.names=TRUE)
}



############################################################################
############################################################################
mymain <- function(n=10, type="PLANAR", repetitions=10)
{	gc()
	
	for(r in 1:repetitions)
	{	# retrieve or create the graph
		g <- init.graph(n, type, iteration=r)
		
		# retrieve the discretization table exists
		disc.table <- init.disc.table(n, type, iteration=r, g)
		
		# deal with the continuous version
		process.continuous.straightness(n,type,iteration=r,g)
		
		
		
		
		#mode <- "graph" # node graph
		for(mode in c("graph","node"))
		{	tlog("Processing the average straightness for node ",node)
			
			########## continuous version
			tlog(2,"Processing the continuous average")
			start.time <- Sys.time()
			if(mode=="node")
				pus <- mean.straightness.nodes.graph(graph=g, u=node)					# process node straightness
			else
				pus <- mean.straightness.graph(graph=g)									# process graph straightness
			end.time <- Sys.time()
			pus.duration <- difftime(end.time,start.time,units="s")
			tlog(4,"Average point straightness: ",pus," - Duration: ",pus.duration," s")
			gc()
			
			# record the data as a text file
			table.file <- paste("data/n=",vcount(g),"-",mode,"-continuous",".txt",sep="")
			data <- matrix(c(pus,pus.duration),ncol=3)
			colnames(data) <- c("Straightness","Duration")
			write.table(x=data,file=table.file,row.names=FALSE,col.names=TRUE)
			
			########## discrete version
			tlog(2,"Processing the discrete approximations")
			grans <- c(0,seq(from=max(E(g)$dist)/2,to=0.004,by=-0.001))#seq(from=0.10,to=0.004,by=-0.0001))
			prev.n <- 0
			est.str <- c(); est.nbr <- c(); est.duration <- c() ; used.grans <- c()
			for(d in 1:length(grans))
			{	tlog(4,"Iteration "	,d,"/",length(grans)," granularity: ",grans[d])
				start.time <- Sys.time()
				g2 <- add.intermediate.nodes(g, granularity=grans[d])			# create additional nodes
				if(vcount(g2)!=prev.n)											# check that the number of nodes is at least different compared to the previous graph
				{	nbr <- vcount(g2)											# total number of nodes
					if(mode=="node")
						str <- mean.straightness.nodes(graph=g2,v=node)[1,1]	# process approximate node straightness
					else
						str <- mean.straightness.nodes(graph=g2, v=NA)[1]		# process nodal approximate graph straightness
					end.time <- Sys.time()
					duration <- difftime(end.time,start.time,units="s")
					est.nbr <- c(est.nbr,nbr)
					est.str <- c(est.str,str)
					used.grans <- c(used.grans,grans[d])
					est.duration <- c(est.duration,duration)
					prev.n <- nbr
					tlog(6,"Number of nodes: ",nbr," - Duration: ",duration," s - Straightness: ",str," (Difference: ",abs(str-pus),")")
				}
				g2 <- NULL
				gc()
			}
			
			# record the data as a text file
			table.file <- paste("data/n=",vcount(g),"-",mode,"-discrete",".txt",sep="")
			data <- cbind(used.grans,est.nbr,est.str,est.duration)
			colnames(data) <- c("Granularity","Nodes","Straightness","Duration")
			write.table(x=data,file=table.file,row.names=FALSE,col.names=TRUE)
			
			# generate straightness plot
			pdf(file=paste("data/n=",vcount(g),"-",mode,"-straightness",".pdf",sep=""))		# open PDF file
			plot(x=est.nbr, y=est.str,														# plot approximations
					xlab="Number of nodes", ylab="Straightness",
					col="BLUE" 
					,ylim=c(min(c(est.str,pus)),max(c(est.str,pus)))
			)
			lines(x=c(min(est.nbr),max(est.nbr)),y=rep(pus,2),col="RED")					# plot exact value
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
			# generate computational cost plot
			pdf(file=paste("data/n=",vcount(g),"-",mode,"-time",".pdf",sep=""))				# open PDF file
			plot(x=est.nbr, y=est.duration,													# plot approximations
					xlab="Number of nodes", ylab="Time (s)",
					col="BLUE"
					,ylim=c(min(c(est.duration,pus.duration)),max(c(est.duration,pus.duration)))
			)
			lines(x=c(min(est.nbr),max(est.nbr)),y=rep(pus.duration,2),col="RED")			# plot exact value
			legend(x="bottomright",legend=c("Approximation","Exact value"),
					fill=c("BLUE","RED"))
			dev.off()
			
			gc()
		}
	}
}

#mymain(n=10,type="PLANAR",repetitions=3)
