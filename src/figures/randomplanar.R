############################################################################
# Script used to produce some of the paper figures from the experiments section.
# This one focuses on a random planar graph.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/figures/randomplanar.R")
############################################################################
source("src/misc/log.R")
source("src/misc/plot.R")
source("src/misc/transformations.R")
source("src/misc/generation.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")



#######################################
# load the original graph
tlog("Load the original graph")
g <- read.graph("data/figures/graph1.graphml",format="graphml")

# process several variants of the average straightness
tlog("Process the average straightness of each link")
link.str <- mean.straightness.links.graph(graph=g)
tlog("Process the average straightness of each node")
node.str <- mean.straightness.nodes.graph(graph=g)

# plot them
tlog("Generate the corresponding plots")
myplot.graph(g, node.str=NA, link.str=link.str, large=TRUE, filename="graph1-original-link-str", out.folder="data/figures", export=FALSE, formats="pdf")
myplot.graph(g, node.str=node.str, link.str=NA, large=TRUE, filename="graph1-original-node-str", out.folder="data/figures", export=FALSE, formats="pdf")

#######################################
# remove a few nodes from the graph
tlog("Remove nodes 16-18 from the graph")
g2 <- delete.vertices(g,16:18)

# process the average straightness again
tlog("Process the average straightness of each link of the modified graph")
link.str2 <- mean.straightness.links.graph(graph=g2)
tlog("Process the average straightness of each node of the modified graph")
node.str2 <- mean.straightness.nodes.graph(graph=g2)

# plot them again
tlog("Generate the corresponding plots")
myplot.graph(g2, node.str=NA, link.str=link.str2, large=TRUE, filename="graph1-modified-link-str", out.folder="data/figures", export=FALSE, formats="pdf")
myplot.graph(g2, node.str=node.str2, link.str=NA, large=TRUE, filename="graph1-modified-node-str", out.folder="data/figures", export=FALSE, formats="pdf")

