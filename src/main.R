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




############################################################################
# Test unilateral point straightness
############################################################################
# init
g <- graph.empty(n=50, directed=FALSE)									# create empty graph
V(g)$x <- runif(vcount(g),min=-1,max=1)									# setup the node spatial positions
V(g)$y <- runif(vcount(g),min=-1,max=1)
V(g)$x[1] <- 0															# put the network center at (0,0)
V(g)$y[1] <- 0
source("src/raw-generation.R")											# use triangulation to init the links
g <- connect.triangulation(g)
#g <- add.edges(g,edges=c(1,2,2,3))
g <- distances.as.weights(g)											# add inter-node distances as link attributes
V(g)$label <- 1:vcount(g)
source("src/model/plot.R")												# display network
graph.file <- paste("n=",vcount(g),"-graph",sep="")
display.model(g, large=TRUE, filename=graph.file, out.folder="data/", export=TRUE, formats=c("pdf",NA))

# processing
center <- 1																# node id of the center node
for(mode in c("unilateral","bilateral"))
{	cat("++++ Processing mode=",mode,"\n",sep="")
	
	# exact values
	start.time <- Sys.time()
	#g2 <- g; V(g2)$type <- "original";V(g2)$label <- 1:vcount(g2)												# just used once for the figure
	if(mode=="unilateral")
		pus <- point.unilateral.straightness(g,center)					# process and display straightness
	else
		pus <- point.bilateral.straightness(g)							# process and display straightness
	#display.model(g2, large=FALSE, filename="graph", out.folder="data/", export=FALSE, formats=c("pdf",NA))	# just used once for the figure
	end.time <- Sys.time()
	pus.duration <- end.time - start.time
	cat("Average point straightness:",pus," time needed: ",pus.duration,"\n")
	
	# approximate values
	unilat <- NA
	if(mode=="unilateral")
		unilat <- center
	grans <- seq(from=max(E(g)$dist)/2,to=0.004,by=-0.001)#seq(from=0.10,to=0.004,by=-0.0001)
#	grans <- seq(from=0.10,to=0.004,by=-0.0001)#seq(from=0.10,to=0.004,by=-0.001)
	prev.n <- vcount(g)
	est.str <- c(); est.nbr <- c(); est.duration <- c()
	i <- 1
	for(d in 1:length(grans))
	{	cat("Iteration ",d,"/",length(grans)," granularity: ",grans[d],"\n",sep="")
		start.time <- Sys.time()
		g2 <- add.intermediate.nodes(g, granularity=grans[d])				# create additional nodes
		if(vcount(g2)!=prev.n)
		{	prev.n <- vcount(g2)
#			display.model(g2, large=FALSE, filename="graph", out.folder="data/", export=FALSE, formats=NA)
			nbr <- vcount(g2)													# total number of nodes
			str <- mean.node.straightness(g2,unilateral=unilat)[1]				# process approximate straightness
			end.time <- Sys.time()
			duration <- end.time - start.time
			est.nbr <- c(est.nbr,nbr)
			est.str <- c(est.str,str)
			cat("nbr=",nbr," duration=",duration," str=",str," (err=",abs(str-pus),")","\n",sep="")
			est.duration <- c(est.duration,duration)
			gc()
			#print(duration)
		}
	}
	#print(cbind(grans,est.nbr,est.str,est.duration))
	#start.time <- Sys.time()
	#str <- mean.node.straightness(g,unilateral=TRUE)[1]					# adding the original graph
	#end.time <- Sys.time()
	#duration <- end.time - start.time
	#est.nbr <- c(vcount(g),est.nbr)
	#est.str <- c(str,est.str)
	#est.duration <- c(duration, est.duration)
	
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
			col="BLUE")
	lines(x=c(min(est.nbr),max(est.nbr)),y=rep(pus.duration,2),col="RED")			# plot exact value
	legend(x="bottomright",legend=c("Approximation","Exact value"),
			fill=c("BLUE","RED"))
	dev.off()
}


### comparison on many graphs
#source("src/raw-generation.R")
#
#im <- 50
#exact.val <- rep(NA,im)
#approx.val <- rep(NA,im)
#semi.val <- rep(NA,im) #TODO
#for(i in 1:im)
#{	cat("Iteration ",i,"/",im,"\n",sep="")
#	
#	# init graph
##	g <- graph.empty(n=10, directed=FALSE)
##	V(g)$x <- runif(vcount(g),min=-1,max=1); V(g)$x[1] <- 0
##	V(g)$y <- runif(vcount(g),min=-1,max=1); V(g)$y[1] <- 0
##	g <- connect.triangulation(g)
##	g <- distances.as.weights(g)
#	
#	g <- graph.empty(n=3, directed=FALSE)
#	V(g)$x <- runif(vcount(g),min=-1,max=1); V(g)$x[1] <- 0
#	V(g)$y <- runif(vcount(g),min=-1,max=1); V(g)$y[1] <- 0
#	g <- add.edges(g,c(1,2,2,3))
#	g <- distances.as.weights(g)
#	
#	# exact value
#	exact.val[i] <- point.bilateral.straightness(g)
#	# approximate value
#	g2 <- add.intermediate.nodes(g, granularity=0.004)#0.0005)
#	approx.val[i] <- mean.node.straightness(g2, unilateral=NA)[1]
#	gc()
#	# semi-approximate value
##	semi.val[i] <- mean.node.straightness2(g2, unilateral=NA)[1]
#	gc()
#	
#	cat("....Exact: ",exact.val[i]," semi-approx: ",semi.val[i]," approximated: ",approx.val[i],"\n",sep="")
#}
#
## generate straightness plot
#idx <- order(exact.val)
#pdf(file=paste("data/","comparison-straightness",".pdf",sep=""))
#plot(x=1:im, y=exact.val[idx],
#	xlab="Test id", ylab="Straightness",
#	col="RED" 
#	,ylim=c(min(c(exact.val,approx.val)),max(c(exact.val,approx.val)))
#)
#points(x=1:im, y=semi.val[idx],col="GREEN")
#points(x=1:im, y=approx.val[idx],col="BLUE")
#legend(x="bottomright",legend=c("Exact value","Semi-Approximation","Approximation"),
#		fill=c("RED","GREEN","BLUE"))
#dev.off()


# TODO in the perf plots, instead of node numbers, use the discrete step and total network length
# number of segments ?
# the idea is to have something comparable independently from the network size and scale.
# so that the reader can think: ok, for the same precision, I need to split the graph in this small segments
