############################################################################
# Reproduce the experimental process described in the article described in
# the project readme file.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/main.R")
############################################################################
source("src/common/misc.R")
source("src/common/plot.R")
source("src/common/transformations.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")



mem.file <- "data/profiling.txt"
planar <- TRUE

mymain <- function(){
	
#for(n in c(10,25,50,100,250,500))
for(n in c(10))
{   tlog("++++++++++++++++++++++ Processing a network of size n=",n)
	gc()
	
############################################################################
# init the graph
tlog("Initializing the graph")
graph.file <- paste("n=",n,"-graph",sep="")						# graph file name
if(file.exists(graph.file))
	g <- read.graph(graph.file,format="graphml")
else
{	# create a planar graph
	if(planar)
	{	g <- graph.empty(n=n, directed=FALSE)									# create empty graph
		V(g)$x <- runif(vcount(g),min=-1,max=1)									# setup the node spatial positions
		V(g)$y <- runif(vcount(g),min=-1,max=1)
		g <- connect.triangulation(g)											# use triangulation to init the links
	}
	# create an Erdös-Rényi graph
	else
	{	g <- erdos.renyi.game(n=n,p.or.m=0.1,directd=FALSE)
		V(g)$x <- runif(vcount(g),min=-1,max=1)									# setup the node spatial positions
		V(g)$y <- runif(vcount(g),min=-1,max=1)
	}
	g <- distances.as.weights(g)											# add inter-node distances as link attributes
	V(g)$label <- 1:vcount(g)
	display.model(g, large=TRUE, filename=graph.file, out.folder="data/", export=TRUE, formats=c("pdf",NA))
}
node <- sample(1:vcount(g),1)												# randomly draw a node for later processing	



############################################################################
#mode <- "graph" # node graph
for(mode in c("graph","node"))
{	if(mode=="node")
		tlog("Processing the average straightness for node",node)
	else
		tlog("Processing the average straightness for the whole graph")
	
	tlog(2,"Processing the continuous average")
	again <- TRUE
	while(again)
	{	again <- FALSE
		Rprof(mem.file, memory.profiling=TRUE, interval=0.0002)
		start.time <- Sys.time()
		if(mode=="node")
			pus <- mean.straightness.nodes.graph(graph=g, u=node)					# process node straightness
		else
			pus <- mean.straightness.graph(graph=g)									# process graph straightness
		end.time <- Sys.time()
		pus.duration <- difftime(end.time,start.time,units="s")
		Sys.sleep(1)
		Rprof(NULL)
		mem.stats <- tryCatch(summaryRprof(mem.file, memory="stats", diff=FALSE, index=1)[["\"mymain\""]],
				error=function(e)
				{	#print(e)
					again<<-TRUE
				})
		if(again)
			tlog(4,"Error while trying to process memory usage. Trying again.")
		else
			pus.mem <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20		 	# used memory, in MB
		gc()
	}
	tlog(4,"Average point straightness: ",pus," - Duration: ",pus.duration," s - Memory: ",pus.mem," MB")
	
	tlog(2,"Processing the discrete approximations")
	grans <- c(0,seq(from=max(E(g)$dist)/2,to=0.004,by=-0.001))#seq(from=0.10,to=0.004,by=-0.0001))
	prev.n <- 0
	est.str <- c(); est.nbr <- c(); est.duration <- c() ; used.grans <- c() ; est.mem <- c()
	i <- 1
	for(d in 1:length(grans))
	{	tlog(4,"Iteration "	,d,"/",length(grans)," granularity: ",grans[d])
		again <- TRUE
		while(again)
		{	again <- FALSE
			Rprof(mem.file, memory.profiling=TRUE, interval=0.0002)
			start.time <- Sys.time()
			g2 <- add.intermediate.nodes(g, granularity=grans[d])			# create additional nodes
			if(vcount(g2)!=prev.n)											# check that the number of nodes is at least different compared to the previous graph
			{	prev.n <- vcount(g2)
				nbr <- vcount(g2)											# total number of nodes
#if(d==1)
#	Sys.sleep(1)							# otherwise, too fast for Rprof
				if(mode=="node")
					str <- mean.straightness.nodes(graph=g2,v=node)[1,1]	# process approximate node straightness
				else
					str <- mean.straightness.nodes(graph=g2, v=NA)[1]		# process nodal approximate graph straightness
				end.time <- Sys.time()
				duration <- difftime(end.time,start.time,units="s")
				Rprof(NULL)
tlog(6,"xxxxxxxx")
				mem.stats <- tryCatch(summaryRprof(mem.file, memory="stats", diff=FALSE, index=1)[["\"mymain\""]],
						error=function(e)
						{	#print(e)
							again<<-TRUE
						})
				if(again)
					tlog(6,"Error while trying to process memory usage. Trying again.")
				else
				{	mem <- (mem.stats[1]*8 + mem.stats[3]*8 + mem.stats[5]*56)/2^20
					est.nbr <- c(est.nbr,nbr)
					est.str <- c(est.str,str)
					used.grans <- c(used.grans,grans[d])
					est.duration <- c(est.duration,duration)
					est.mem <- c(est.mem,mem)
					tlog(6,"Number of nodes: ",nbr," - Duration: ",duration," s - Memory: ",mem," MB - Straightness: ",str," (Difference: ",abs(str-pus),")")
				}
				g2 <- NULL
				gc()
			}
			else
			{	Rprof(NULL)
				g2 <- NULL
				gc()
			}
		}
	}
	
	# record the data as a text file
	table.file <- paste("data/n=",vcount(g),"-",mode,"-continuous",".txt",sep="")
	data <- matrix(c(pus,pus.duration,pus.mem),ncol=3)
	colnames(data) <- c("Straightness","Duration","Memory")
	write.table(x=data,file=table.file,row.names=FALSE,col.names=TRUE)
	table.file <- paste("data/n=",vcount(g),"-",mode,"-discrete",".txt",sep="")
	data <- cbind(used.grans,est.nbr,est.str,est.duration,est.mem)
	colnames(data) <- c("Granularity","Nodes","Straightness","Duration","Memory")
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
	# generate memory cost plot
	pdf(file=paste("data/n=",vcount(g),"-",mode,"-mem",".pdf",sep=""))				# open PDF file
	plot(x=est.nbr, y=est.mem,														# plot approximations
			xlab="Number of nodes", ylab="Max Memory (MB)",
			col="BLUE"
			,ylim=c(min(c(est.mem,pus.mem)),max(c(est.mem,pus.mem)))
	)
	lines(x=c(min(est.nbr),max(est.nbr)),y=rep(pus.mem,2),col="RED")				# plot exact value
	legend(x="bottomright",legend=c("Approximation","Exact value"),
			fill=c("BLUE","RED"))
	dev.off()
	
	gc()
}
}

#TODO check if the tolerance is really necessary
#TODO check how it works when some points are aligned (use addpoints function or just build a simple 3-node graph)
#TODO compare completely numerical and partially numerical integrations
}

# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/main.R")
# Rprof(mem.file, memory.profiling=TRUE)
# mymain();Rprof(NULL);
# mem.stats <- summaryRprof(mem.file, memory="stats", diff=FALSE, index=1) #tseries stats both
# mem <- max(sapply(mem.stats,function(v) v[1]*8+v[3]*8+v[5]*56)/2^20)
