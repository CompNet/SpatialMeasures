############################################################################
# Script used to produce some of the paper figures from the experiments section.
# This one focuses on regular orbitele graphs (i.e. spider webs).
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/figures/experiments2.R")
############################################################################
source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")
source("src/misc/generation.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")



########################################
# generate a square graph
tlog("Generate a square graph")
g <- produce.square.graph(n=14,area=20)

# process several variants of the average straightness
tlog("Process the average straightness of each link")
link.str <- mean.straightness.links.graph(graph=g)
tlog("Process the average straightness of each node")
node.str <- mean.straightness.nodes.graph(graph=g)

# record them 
tlog("Record these results")
write.table(x=link.str,file="data/figures/graph2-link-str.txt",row.names=FALSE,col.names=FALSE)
write.table(x=node.str,file="data/figures/graph2-node-str.txt",row.names=FALSE,col.names=FALSE)

# plot them
tlog("Generate the corresponding plots")
myplot.graph(g, node.str=NA, link.str=link.str, large=TRUE, filename="graph2-link-str", out.folder="data/figures", export=FALSE, formats="pdf")
myplot.graph(g, node.str=node.str, link.str=NA, large=TRUE, filename="graph2-node-str", out.folder="data/figures", export=FALSE, formats="pdf")



######################################
# generate a radio-concentric graph
tlog("Generate a radio-concentric graph")
g2 <- produce.radiocentric.graph(r=8,s=10,area=20)

# process several variants of the average straightness
tlog("Process the average straightness of each link")
link.str2 <- mean.straightness.links.graph(graph=g2)
tlog("Process the average straightness of each node")
node.str2 <- mean.straightness.nodes.graph(graph=g2)

# record them 
tlog("Record these results")
write.table(x=link.str2,file="data/figures/graph3-link-str.txt",row.names=FALSE,col.names=FALSE)
write.table(x=node.str2,file="data/figures/graph3-node-str.txt",row.names=FALSE,col.names=FALSE)

# plot them
tlog("Generate the corresponding plots")
myplot.graph(g2, node.str=NA, link.str=link.str2, large=TRUE, filename="graph3-link-str", out.folder="data/figures", export=FALSE, formats="pdf")
myplot.graph(g2, node.str=node.str2, link.str=NA, large=TRUE, filename="graph3-node-str", out.folder="data/figures", export=FALSE, formats="pdf")



#######################################
# generate an orbitele graph
tlog("Generate spider graph")
g3 <- produce.orbitele.graph(r=8,s=7,area=20)

# process several variants of the average straightness
tlog("Process the average straightness of each link")
link.str3 <- mean.straightness.links.graph(graph=g3)
tlog("Process the average straightness of each node")
node.str3 <- mean.straightness.nodes.graph(graph=g3)

# record them 
tlog("Record these results")
write.table(x=link.str3,file="data/figures/graph4-link-str.txt",row.names=FALSE,col.names=FALSE)
write.table(x=node.str3,file="data/figures/graph4-node-str.txt",row.names=FALSE,col.names=FALSE)

# plot them
tlog("Generate the corresponding plots")
myplot.graph(g3, node.str=NA, link.str=link.str3, large=TRUE, filename="graph4-link-str", out.folder="data/figures", export=FALSE, formats="pdf")
myplot.graph(g3, node.str=node.str3, link.str=NA, large=TRUE, filename="graph4-node-str", out.folder="data/figures", export=FALSE, formats="pdf")
