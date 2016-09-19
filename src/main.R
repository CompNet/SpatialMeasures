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


mymain <- function(){
	
mem.file <- "data/profiling.txt"
	
#for(n in c(10,25,50,100,250,500))
for(n in c(10))
{   cat("++++++++++++++++++++++ Processing a network of size n=",n,"\n",sep="")
############################################################################
# init the graph
cat("Initializing the graph\n")
#g <- graph.empty(n=n, directed=FALSE)									# create empty graph
#V(g)$x <- runif(vcount(g),min=-1,max=1)									# setup the node spatial positions
#V(g)$y <- runif(vcount(g),min=-1,max=1)
#g <- connect.triangulation(g)											# use triangulation to init the links
	g <- erdos.renyi.game(n=n,p.or.m=0.1,directd=FALSE)
	V(g)$x <- runif(vcount(g),min=-1,max=1)									# setup the node spatial positions
	V(g)$y <- runif(vcount(g),min=-1,max=1)	
g <- distances.as.weights(g)											# add inter-node distances as link attributes
V(g)$label <- 1:vcount(g)
graph.file <- paste("n=",vcount(g),"-graph",sep="")						# display network
display.model(g, large=TRUE, filename=graph.file, out.folder="data/", export=TRUE, formats=c("pdf",NA))
node <- sample(1:vcount(g),1)	



############################################################################
# average straightness for a given node
#mode <- "graph" # node graph
for(mode in c("graph","node"))
{	if(mode=="node")
		cat("Processing the average straightness for node",node,"\n")
	else
		cat("Processing the average straightness for the whole graph\n")
	
	cat("..Processing the continuous average\n")
	Rprof(mem.file, memory.profiling=TRUE)
	Sys.sleep(time=1)															# otherwise, too fast for Rprof
	start.time <- Sys.time()
	if(mode=="node")
		pus <- mean.straightness.nodes.graph(graph=g, u=node)					# process node straightness
	else
		pus <- mean.straightness.graph(graph=g)									# process graph straightness
	end.time <- Sys.time()
	Rprof(NULL);
	mem.stats <- summaryRprof(mem.file, memory="stats", diff=FALSE)
	pus.mem <- max(sapply(mem.stats,function(v) v[1]*8+v[3]*8+v[5]*56)/2^20) 	# used memory, in MB 
	pus.duration <- difftime(end.time,start.time,units="s")
	cat("....Average point straightness:",pus," - Duration: ",pus.duration," s - Memory:",pus.mem," MB\n")
	
	cat("..Processing the discrete approximations\n")
	grans <- c(0,seq(from=max(E(g)$dist)/2,to=0.004,by=-0.001))#seq(from=0.10,to=0.004,by=-0.0001))
	prev.n <- 0
	est.str <- c(); est.nbr <- c(); est.duration <- c() ; used.grans <- c() ; est.mem <- c()
	i <- 1
	for(d in 1:length(grans))
	{	cat("....Iteration ",d,"/",length(grans)," granularity: ",grans[d],"\n",sep="")
		Rprof(mem.file, memory.profiling=TRUE)
		Sys.sleep(time=1)												# otherwise, too fast for Rprof
		start.time <- Sys.time()
		g2 <- add.intermediate.nodes(g, granularity=grans[d])			# create additional nodes
		if(vcount(g2)!=prev.n)											# check that the number of nodes is at least different compared to the previous graph
		{	prev.n <- vcount(g2)
#			display.model(g2, large=FALSE, filename="graph", out.folder="data/", export=FALSE, formats=NA)
			nbr <- vcount(g2)											# total number of nodes
			if(mode=="node")
				str <- mean.straightness.nodes(graph=g2,v=node)[1,1]	# process approximate node straightness
			else
				str <- mean.straightness.nodes(graph=g2, v=NA)[1]		# process nodal approximate graph straightness
			end.time <- Sys.time()
			duration <- difftime(end.time,start.time,units="s")
			Rprof(NULL);
			mem.stats <- summaryRprof(mem.file, memory="stats", diff=FALSE)
			mem <- max(sapply(mem.stats,function(v) v[1]*8+v[3]*8+v[5]*56)/2^20) 
			est.nbr <- c(est.nbr,nbr)
			est.str <- c(est.str,str)
			used.grans <- c(used.grans,grans[d])
			cat("......Number of nodes: ",nbr," - Duration: ",duration," s - Memory: ",mem," MB - Straightness: ",str," (Error: ",abs(str-pus),")","\n",sep="")
			est.duration <- c(est.duration,duration)
			g2 <- NULL; gc();
			est.mem <- c(est.mem,mem)
		}
	}
	
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
	# record the data as a text file
	table.file <- paste("data/n=",vcount(g),"-",mode,"-continuous",".txt",sep="")
	data <- matrix(c(pus,pus.duration,pus.mem),ncol=3)
	colnames(data) <- c("Straightness","Duration","Memory")
	write.table(x=data,file=table.file,row.names=FALSE,col.names=TRUE)
	table.file <- paste("data/n=",vcount(g),"-",mode,"-discrete",".txt",sep="")
	data <- cbind(used.grans,est.nbr,est.str,est.duration,est.mem)
	colnames(data) <- c("Granularity","Nodes","Straightness","Duration","Memory")
	write.table(x=data,file=table.file,row.names=FALSE,col.names=TRUE)
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
# mem.stats <- summaryRprof(mem.file, memory="stats", diff=FALSE) #tseries stats both
# mem <- max(sapply(mem.stats,function(v) v[1]*8+v[3]*8+v[5]*56)/2^20)
