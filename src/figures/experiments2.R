############################################################################
# Script used to produce some of the paper figures from the experiments section.
#
# Vincent Labatut 09/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/figures/experiments.R")
############################################################################
source("src/common/misc.R")
source("src/common/plot.R")
source("src/common/transformations.R")
source("src/common/generation.R")

source("src/straightness/continuous.R")
source("src/straightness/discrete.R")



#######################################
# generate a square graph
g <- produce.square.graph(n=15)

# process several variants of the average straightness
tlog("Process the average straightness of each link")
link.str <- mean.straightness.links.graph(graph=g)
tlog("Process the average straightness of each node")
node.str <- mean.straightness.nodes.graph(graph=g)

# plot them
tlog("Generate the corresponding plots")
myplot.graph(g, node.str=NA, link.str=link.str, large=TRUE, filename="graph2-link-str", out.folder="data/figures", export=FALSE, formats="pdf")
myplot.graph(g, node.str=node.str, link.str=NA, large=TRUE, filename="graph2-node-str", out.folder="data/figures", export=FALSE, formats="pdf")

#######################################
# generate a radio-concentric graph
g2 <- produce.radiocentric.graph(r=8,s=10)

# process several variants of the average straightness
tlog("Process the average straightness of each link")
link.str <- mean.straightness.links.graph(graph=g2)
tlog("Process the average straightness of each node")
node.str <- mean.straightness.nodes.graph(graph=g2)

# plot them
tlog("Generate the corresponding plots")
myplot.graph(g2, node.str=NA, link.str=link.str, large=TRUE, filename="graph3-link-str", out.folder="data/figures", export=FALSE, formats="pdf")
myplot.graph(g2, node.str=node.str, link.str=NA, large=TRUE, filename="graph3-node-str", out.folder="data/figures", export=FALSE, formats="pdf")


