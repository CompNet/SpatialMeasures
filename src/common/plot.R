############################################################################
# Functions related to the plotting of the generated networks, including
# additional information such as link color to represent link addition/deletion.
#
# Note: this script contains some functions taken from the Web.
#
# Vincent Labatut 12/2015
#
# setwd("~/eclipse/workspaces/Networks/SpiderNet")
# source("src/model/plot.R")
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
# large: whether to use (graphically) large nodes containing their id, or small ones.
# filename: the base name of the various files to generate.
# out.folder: folder containing the generated files.
# export: whether to record only image files (FALSE) or also data files (TRUE).
# formats: formats of the produced plot files (NA for onscreen plotting).
# ...: other parameters, passed to the regular plot function.
############################################################################################
display.model <- function(g, large=TRUE, filename=NA, out.folder=NA, export=TRUE, formats=c(NA,"png","pdf"), ...)
{	# plot parameters
	vertex.sizes <- 1
	vertex.label <- ""
	if("type" %in% list.vertex.attributes(g))
	{	vertex.sizes[V(g)$type=="original"] <- 15
		vertex.label[V(g)$type=="original"] <- V(g)$label[V(g)$type=="original"]
		vertex.sizes[V(g)$type=="extra"] <- 5
		vertex.label[V(g)$type=="extra"] <- ""
	}
vertex.label <- 1:vcount(g) #TODO remove this after debug
	
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
	link.cols <- rep("BLACK",ecount(g))
	link.cols[E(g)$added] <- "GREEN"
	link.widths <- rep(1,ecount(g))
	link.widths[E(g)$added] <- 6
	
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
						rescale=FALSE,axes=TRUE,asp=1,xlim=lm,ylim=lm,
						edge.color=link.cols,edge.width=link.widths,
						...
				) 
			else 
				plot(g,main=g$title,
						vertex.size=vertex.sizes,vertex.label=vertex.label,vertex.label.cex=0.5,
						edge.color=link.cols,edge.width=link.widths,
						rescale=FALSE,axes=TRUE,asp=1,xlim=lm,ylim=lm,
						...
				)
		}
		
		# add removed links
		if(!is.null(g$removed))
		{	print(g$removed)
			segments(V(g)$x[g$removed[,1]],V(g)$y[g$removed[,1]],
					V(g)$x[g$removed[,2]],V(g)$y[g$removed[,2]],
					col="RED",lty="dashed",lwd=4
			)
		}
		
		# add arrows to represent move vectors
		idx <- (1:vcount(g))[V(g)$sx!=0 |  V(g)$sy!=0]
		if(length(idx>0))
			arrows(V(g)$x[idx],V(g)$y[idx],
					V(g)$x[idx]+V(g)$sx[idx],V(g)$y[idx]+V(g)$sy[idx],
					length=0.1,col="BLUE",lwd=2
			)
		idx <- (1:vcount(g))[V(g)$ax!=0 |  V(g)$ay!=0]
		if(length(idx>0))
			arrows(V(g)$x[idx],V(g)$y[idx],
					V(g)$x[idx]+V(g)$ax[idx],V(g)$y[idx]+V(g)$ay[idx],
					length=0.1,col="RED",lwd=2
			)
		idx <- (1:vcount(g))[V(g)$cx!=0 |  V(g)$cy!=0]
		if(length(idx>0))
			arrows(V(g)$x[idx],V(g)$y[idx],
					V(g)$x[idx]+V(g)$cx[idx],V(g)$y[idx]+V(g)$cy[idx],
					length=0.1,col="GREEN",lwd=2
			)
		
		# add virtual links (possibly produced during cover processing)
		if(!is.null(g$virt.links) && nrow(g$virt.links)>0)
		{	for(i in 1:nrow(g$virt.links))
			{	n1 = g$virt.links[i, 1]
				n2 = g$virt.links[i, 2]
				
				segments(V(g)$x[n1], V(g)$y[n1], 
							V(g)$x[n2], V(g)$y[n2],
							col="GREY",lty="dashed",lwd=4
				)
			}
		}
		
		# add polygons (possibly produced during cover processing)
		if(!is.null(g$polygons) && length(g$polygons)>0)
		{	for(i in 1:length(g$polygons))
			{	for(j in 1:length(g$polygons[[i]]))
				{	# get the triangle
					triangle <- g$polygons[[i]][[j]]
					# get its three vertices
					a <- triangle[1]
					b <- triangle[2]
					c <- triangle[3]
					# get their positions
					ax <- V(g)$x[a]; ay <- V(g)$y[a]
					bx <- V(g)$x[b]; by <- V(g)$y[b]
					cx <- V(g)$x[c]; cy <- V(g)$y[c]
					# draw the polygon
					polygon(c(ax, bx, cx), c(ay,by,cy),
							density=NA,
							col=add.alpha(i, 0.2),
							lty="blank", lwd=1
					)
				}
			}
		}
		
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
