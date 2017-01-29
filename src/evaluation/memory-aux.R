############################################################################
# Script used by memory2.R to monitor the memory. It allows to start a fresh
# new R session each time we process a network. Indeed, the R garbage collector
# does not seem to really clean the memory when invoked. This allows to force 
# clearing the memory.
# 
# Vincent Labatut 10/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/evaluation/memory-aux.R")
# Rscript --vanilla memory-aux.R workspace graph_file discretization_coeff node_id
############################################################################

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
print(args)

# load the necessary libs
library("igraph")

# parameter constants
PAR_WP    <- 1
PAR_GRAPH <- 2
PAR_DISC  <- 3
PAR_NODE  <- 4

# memory monitor constants
MEM_FILE <- "data/profiling.txt"	# R profiler output file (temporary)
MEM_INTER <- 0.00002				# R profiler update rate
MEM_SLEEP <- 5 						# sleep between to tries, in seconds

# use/get parameter values
setwd(args[PAR_WP])
node <- as.integer(args[PAR_NODE])
disc <- as.integer(args[PAR_DISC])
g <- read.graph(args[PAR_GRAPH],format="graphml")

# load the necessary scripts
source("src/evaluation/common.R")

# process it
temp <- function(g, node, disc)
{	Rprof(MEM_FILE, memory.profiling=TRUE, interval=MEM_INTER)
		if(disc==0)
		{	# whole graph
			if(node==0)
			{	mean.straightness.graph(graph=g)
				Sys.sleep(MEM_SLEEP)
			}
			# specific node
			else
				mean.straightness.nodes.graph(graph=g, u=node)
			
		}
		else #if(args[PAR_MODE]=="DISCRETE")
		{	#Sys.sleep(MEM_SLEEP)
			g2 <- add.intermediate.nodes(g, granularity=disc)
			# whole graph
			if(node==0)
			{	mean.straightness.nodes(graph=g2, v=NA)[1]
				#Sys.sleep(MEM_SLEEP)
			}
			# specific node
			else
				mean.straightness.nodes(graph=g2,v=node)[1,1]
		}
	#Sys.sleep(MEM_SLEEP)
	Rprof(NULL)
}
temp(g,node,disc)
#invisible(T)
