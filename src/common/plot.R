############################################################################
# Functions related to the plotting of the generated networks, including
# additional information such as link color to represent link addition/deletion.
#
# Note: this script contains some functions taken from the Web.
#
# Vincent Labatut 12/2015
#
# setwd("~/eclipse/workspaces/Networks/SpiderNet")
# source("src/common/plot.R")
############################################################################



############################################################################################
# Adds an alpha chanel to an existing color.
# Source code by Markus Gesmann, retrieved from the following URL:
#	https://gist.github.com/mages/5339689#file-add-alpha-r
############################################################################################
add.alpha <- function(col, alpha=1)
{	if(missing(col))
		stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb)/255, 2, 
			function(x) 
				rgb(x[1], x[2], x[3], alpha=alpha))  
}



############################################################################################
# Displays/exports a generated graph. The graph is exported as both an image (PDF and/or PNG) and 
# various data files. The plotting takes some additional data into account, such as added/removed
# links, polygons, etc.
#
# g: the graph to process.
# node.str: the straightness values associated to the links.
# link.str: the straightness values associated to the links.
# large: whether to use (graphically) large nodes containing their id, or small ones.
# filename: the base name of the various files to generate.
# out.folder: folder containing the generated files.
# export: whether to record only image files (FALSE) or also data files (TRUE).
# formats: formats of the produced plot files (NA for onscreen plotting).
# ...: other parameters, passed to the regular plot function.
############################################################################################
myplot.graph <- function(g, node.str=NA, link.str=NA, large=TRUE, filename=NA, out.folder=NA, export=TRUE, formats=c(NA,"png","pdf"), ...)
{	# plot parameters
	vertex.sizes <- 1
	vertex.label <- ""
	vertex.color <- "BEIGE"
	if("type" %in% list.vertex.attributes(g))
	{	vertex.sizes[V(g)$type=="original"] <- 15
		vertex.label[V(g)$type=="original"] <- V(g)$label[V(g)$type=="original"]
		vertex.sizes[V(g)$type=="extra"] <- 5
		vertex.label[V(g)$type=="extra"] <- ""
	}
#vertex.label <- 1:vcount(g)
	# set node colors
	if(!all(is.na(node.str)))
	{	cscale <- colorRamp(c("BLUE","CYAN","YELLOW","RED"))
		vertex.color <- apply(cscale(node.str), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
	}
	
	# possibly create output folder
	if(!is.na(filename) & !is.na(out.folder))
		dir.create(path=out.folder,showWarnings=FALSE,recursive=TRUE)
	
	# plot parameters
	mn <- 0
	mx <- max(c(V(g)$x,V(g)$y))
	mn <- min(c(V(g)$x,V(g)$y))	#TODO added for model debugging
	lm <- c(mn,mx)
#lm <- c(-1,+1)	#TODO added for model debugging
#	print(lm)
	
	# set link colors
	if(all(is.na(link.str)))
	{	link.cols <- rep("BLACK",ecount(g))
		link.widths <- rep(1,ecount(g))
	}
	else
	{	cscale <- colorRamp(c("BLUE","CYAN","YELLOW","RED"))
		link.cols <- apply(cscale(link.str), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
		link.widths <- 4*link.str
	}
		
	# generate plots
	for(format in formats)
	{	if(!is.na(format) & !is.na(filename) & !is.na(out.folder))
		{	if(format=="pdf")
				pdf(file=file.path(out.folder,paste0(filename,".pdf")))
			else if(format=="png")
				png(filename=file.path(out.folder,paste0(filename,".png")),
						width=800,height=800,units="px",pointsize=20,bg="white"
				)
		}
		
		if(is.na(format) | (!is.na(filename) & !is.na(out.folder)))
		{	if(large) 
				plot(g,main=g$title,
						vertex.color=vertex.color,
						vertex.label.color="BLACK",vertex.label.cex=0.5,
						edge.color=link.cols,edge.width=link.widths,
						rescale=FALSE,axes=TRUE,asp=1,xlim=lm,ylim=lm,
						...
				) 
			else 
				plot(g,main=g$title,
						vertex.size=vertex.sizes,vertex.color=vertex.color,
						vertex.label=vertex.label,vertex.label.cex=0.5,vertex.label.color="BLACK",
						edge.color=link.cols,edge.width=link.widths,
						rescale=FALSE,axes=TRUE,asp=1,xlim=lm,ylim=lm,
						...
				)
		}
		
		# add the legend
#		if(!all(is.na(link.str)) || !all(is.na(node.str)))
#			legend.gradient(pnts=, cols=c("BLUE","CYAN","YELLOW","RED"), 
#					limits=c(0,1), title="Straightness"
#			)
		
		# add the plot title
		if(!is.null(g$name))
			title(g$name, 
#					cex.main=0.5
			)
		
		if(!is.na(format) & !is.na(filename) & !is.na(out.folder))
			dev.off()
	}
	
	# produce network files
	if(!is.na(filename) & !is.na(out.folder) & export)
	{	write.graph(graph=g, file=file.path(out.folder,paste0(filename,".graphml")), format="graphml")
		write.graph(graph=g, file=file.path(out.folder,paste0(filename,".net")), format="pajek")
		write.graph(graph=g, file=file.path(out.folder,paste0(filename,".edgelist")), format="edgelist")
	}
}
