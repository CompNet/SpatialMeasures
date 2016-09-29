############################################################################
# Script used to produce the figure of the paper which illustrate the graph
# distance between two non-node points. Some manual process was involved, the
# script doesn't directly generate the figure.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/figures/proof.R")
############################################################################
library("igraph")

source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")
source("src/misc/generation.R")




# init the graph
g <- graph.empty(7,directed=FALSE)
#g <- add.edges(g,c(1,7,7,2,2,3,3,5,4,5,4,6,6,1,3,6,5,6))
g <- add.edges(g,c(1,7,7,2,2,3,3,4,4,5,5,6,6,1,3,6,6,4)) # G'
V(g)$x <- c(
	0,		# v1 = v_1
	-0.5,	# v2 = u_1
	0.75,	# v3
	2,		# v4 = u_2
	3,		# v5 = v_2
	1.25,	# v6
#	-0.4342404	# p1
#	-0.3212319	# p1
#	-0.1365337	# p1
	-0.2		# p1'
)
V(g)$y <- c(
	0,		# v1 = v_1
	3,		# v2 = u_1
	3.5,	# v3
	3.5,	# v4 = u_2
	2,		# v5 = v_2
	1.5,	# v6
#	2.605442		# p1
#	1.927391		# p1
#	0.8192024		# p1
	1.2			# p1'
)
g <- distances.as.weights(g)
plot(g,rescale=FALSE,xlim=c(min(V(g)$x),max(V(g)$x)),ylim=c(min(V(g)$y),max(V(g)$y)))

# get the graph distances
sp <- shortest.paths(graph=g, weights=E(g)$dist)

# process lambda_u and lambda_v
lu <- (sp[4,1]-sp[4,2]+sp[1,2])/2
qux <- V(g)$x[2]+lu/sp[1,2]*(V(g)$x[1]-V(g)$x[2])
quy <- V(g)$y[2]+lu/sp[1,2]*(V(g)$y[1]-V(g)$y[2])
cat("lamnda_u=",lu," qu=(",qux,",",quy,")\n",sep="")
points(x=qux,y=quy,col="RED")
lv <- (sp[5,1]-sp[5,2]+sp[1,2])/2
qvx <- V(g)$x[2]+lv/sp[1,2]*(V(g)$x[1]-V(g)$x[2])
qvy <- V(g)$y[2]+lv/sp[1,2]*(V(g)$y[1]-V(g)$y[2])
cat("lamnda_v=",lv," qv=(",qvx,",",qvy,")\n",sep="")
points(x=qvx,y=qvy,col="BLUE")

# process lambda_2
l2 <- (sp[7,5]-sp[7,4]+sp[4,5])/2
q2x <- V(g)$x[4]+l2/sp[4,5]*(V(g)$x[5]-V(g)$x[4])
q2y <- V(g)$y[4]+l2/sp[4,5]*(V(g)$y[5]-V(g)$y[4])
cat("lambda_2=",l2," q2=(",q2x,",",q2y,")\n",sep="")
points(x=q2x,y=q2y,col="GREEN")

# process point coordinates p depending on their relative positions
lp1 <- 2.210879
p1x <- V(g)$x[2]+lp1/sp[1,2]*(V(g)$x[1]-V(g)$x[2])
p1y <- V(g)$y[2]+lp1/sp[1,2]*(V(g)$y[1]-V(g)$y[2])
cat("ell_p_1=",lp1," p1=(",p1x,",",p1y,")\n",sep="")
lp2 <- 1.2
p2x <- V(g)$x[4]+lp2/sp[4,5]*(V(g)$x[5]-V(g)$x[4])
p2y <- V(g)$y[4]+lp2/sp[4,5]*(V(g)$y[5]-V(g)$y[4])
cat("ell_p_2=",lp2," p2=(",p2x,",",p2y,")\n",sep="")

print(get.shortest.paths(g,from=7, to=4, weights=E(g)$dist)$vpath)
print(get.shortest.paths(g,from=7, to=5, weights=E(g)$dist)$vpath)
