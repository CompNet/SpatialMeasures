############################################################################
# Functions related to the plotting of the generated networks, including
# additional information such as link color to represent link addition/deletion.
#
# Note: this script contains some functions taken from the Web.
#
# Vincent Labatut 12/2015
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("d:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/misc/plot.R")
############################################################################
library(plotrix)




############################################################################
# Generate a gradient legend, to be added to an existing plot.
# 
# Source code by John Colby, retrieved from the following URL:
#	http://stackoverflow.com/a/9314880/1254730
############################################################################
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title="") 
{	scale = (length(lut)-1)/(max-min)
	
#	dev.new(width=1.75, height=5)
	plot(c(0,10), c(min,max), type="n", bty="n", xaxt="n", xlab=title, yaxt="n", ylab="", main="")
	axis(2, ticks, las=1)
	for (i in 1:(length(lut)-1))
	{	y <- (i-1)/scale + min
		rect(0,y,10,y+1/scale, col=lut[i], border=NA)
	}
}
#pdf(file="data/legend.pdf",width=1.75,height=5)
#color.bar(colorRampPalette(c("BLUE","CYAN","YELLOW","RED"))(100), min=0, max=1, title="Straightness")
#dev.off()



############################################################################################
# Adds an alpha chanel to an existing color.
#
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
# autoscale: TRUE to let igraph automatically scale the graph representation.
# ...: other parameters, passed to the regular plot function.
############################################################################################
myplot.graph <- function(g, node.str=NA, link.str=NA, large=TRUE, filename=NA, out.folder=NA, export=TRUE, formats=c(NA,"png","pdf"), autoscale=FALSE, ...)
{	# possibly create output folder
	if(!is.na(filename) & !is.na(out.folder))
		dir.create(path=out.folder,showWarnings=FALSE,recursive=TRUE)
	
	# plot parameters
	vertex.sizes <- rep(1,vcount(g))
	vertex.label <- rep("",vcount(g))
	vertex.color <- rep("BEIGE",vcount(g))
	vertex.shape <- rep("circle",vcount(g))
	vertex.frame.color <- rep("BLACK",vcount(g))
	
	# adjust in function of the vertex types
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
	
	# mark certain nodes
	if("marked" %in% list.vertex.attributes(g))
	{	#vertex.color[V(g)$marked] <- "WHITE"
		vertex.shape[V(g)$marked] <- "square"
		vertex.frame.color[V(g)$marked] <- "PURPLE"
	}
	
	# set link colors
	if(all(is.na(link.str)))
	{	link.cols <- rep("BLACK",ecount(g))
		link.widths <- rep(1,ecount(g))
	}
	else
	{	cscale <- colorRamp(c("BLUE","CYAN","YELLOW","RED"))
		link.cols <- apply(cscale(link.str), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
		#link.widths <- 4*link.str
		link.widths <- rep(3,ecount(g))
	}

	# mark certain links
	if("marked" %in% list.edge.attributes(g))
	{	link.widths[E(g)$marked] <- link.widths[E(g)$marked]*3
		link.cols[E(g)$marked] <- "PURPLE"
	}
	
	# setup plot dimensions
	mn <- 0
	mx <- max(c(V(g)$x,V(g)$y))
	mn <- min(c(V(g)$x,V(g)$y))	#TODO added for model debugging
	lm <- c(mn,mx)
#lm <- c(-1,+1)	#TODO added for model debugging
#	print(lm)
	
	# generate plots
	for(format in formats)
	{	if(!is.na(format) & !is.na(filename) & !is.na(out.folder))
		{	format <- tolower(format)
			if(format=="pdf")
				pdf(file=file.path(out.folder,paste0(filename,".pdf")))
			else if(format=="png")
				png(filename=file.path(out.folder,paste0(filename,".png")),
						width=800,height=800,units="px",pointsize=20,bg="white"
				)
		}
		
		if(is.na(format) | (!is.na(filename) & !is.na(out.folder)))
		{	if(large) 
			{	if(autoscale)
				{	plot(g,main=g$title,
						vertex.color=vertex.color,
						vertex.label.color="BLACK",vertex.label.cex=0.5,
						vertex.frame.color=vertex.frame.color,vertex.shape=vertex.shape,
						edge.color=link.cols,edge.width=link.widths,
						rescale=TRUE,
						axes=FALSE,
						...
					)
				}
				else
				{	plot(g,main=g$title,
						vertex.color=vertex.color,
						vertex.label.color="BLACK",vertex.label.cex=0.5,
						vertex.frame.color=vertex.frame.color,vertex.shape=vertex.shape,
						edge.color=link.cols,edge.width=link.widths,
						rescale=FALSE,asp=1,xlim=lm,ylim=lm,
						axes=FALSE,
						...
					)
				}
			}
			else 
			{	if(autoscale)
				{	plot(g,main=g$title,
						vertex.size=vertex.sizes,vertex.color=vertex.color,
						vertex.label=vertex.label,vertex.label.cex=0.5,vertex.label.color="BLACK",
						vertex.frame.color=vertex.frame.color,vertex.shape=vertex.shape,
						edge.color=link.cols,edge.width=link.widths,
						rescale=TRUE,
						axes=FALSE,
						...
					)
				}
				else
				{	plot(g,main=g$title,
						vertex.size=vertex.sizes,vertex.color=vertex.color,
						vertex.label=vertex.label,vertex.label.cex=0.5,vertex.label.color="BLACK",
						vertex.frame.color=vertex.frame.color,vertex.shape=vertex.shape,
						edge.color=link.cols,edge.width=link.widths,
						rescale=FALSE,asp=1,xlim=lm,ylim=lm,
						axes=FALSE,
						...
					)
				}
			}
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




############################################################################
# There is a problem with plotrix function axis.break when inserting a gap
# in a log axis: it draws on top of the existing plot without removing
# the appropriate parts. This function fixes that (see the commented lines).
# 
# See the original documentation to know what parameters the function takes.
############################################################################
axis.break2 <- function (axis = 1, breakpos = NULL, pos = NA, bgcol = "white", breakcol = "black", style = "slash", brw = 0.02) 
{	figxy <- par("usr")
	xaxl <- par("xlog")
	yaxl <- par("ylog")
	xw <- (figxy[2] - figxy[1]) * brw
	yw <- (figxy[4] - figxy[3]) * brw
	if (!is.na(pos)) 
		figxy <- rep(pos, 4)
	if (is.null(breakpos)) 
		breakpos <- ifelse(axis%%2, figxy[1] + xw * 2, figxy[3] + 
						yw * 2)
# >> VL no need to normalize the coordinates 
#    if (xaxl && (axis == 1 || axis == 3)) 
#        breakpos <- log10(breakpos)
#    if (yaxl && (axis == 2 || axis == 4)) 
#        breakpos <- log10(breakpos)
# << VL no need to normalize the coordinates 
	switch(axis, br <- c(breakpos - xw/2, figxy[3] - yw/2, breakpos + 
							xw/2, figxy[3] + yw/2), br <- c(figxy[1] - xw/2, breakpos - 
							yw/2, figxy[1] + xw/2, breakpos + yw/2), br <- c(breakpos - 
							xw/2, figxy[4] - yw/2, breakpos + xw/2, figxy[4] + yw/2), 
			br <- c(figxy[2] - xw/2, breakpos - yw/2, figxy[2] + 
							xw/2, breakpos + yw/2), stop("Improper axis specification."))
	old.xpd <- par("xpd")
	par(xpd = TRUE)
	if (xaxl) 
		br[c(1, 3)] <- 10^br[c(1, 3)]
	if (yaxl) 
		br[c(2, 4)] <- 10^br[c(2, 4)]
	if (style == "gap") {
		if (xaxl) {
			figxy[1] <- 10^figxy[1]
			figxy[2] <- 10^figxy[2]
		}
		if (yaxl) {
			figxy[3] <- 10^figxy[3]
			figxy[4] <- 10^figxy[4]
		}
		if (axis == 1 || axis == 3) {
			rect(breakpos, figxy[3], breakpos + xw, figxy[4], 
					col = bgcol, border = bgcol)
			xbegin <- c(breakpos, breakpos + xw)
			ybegin <- c(figxy[3], figxy[3])
			xend <- c(breakpos, breakpos + xw)
			yend <- c(figxy[4], figxy[4])
# >> VL no need to normalize the coordinates 
#            if (xaxl) {
#                xbegin <- 10^xbegin
#                xend <- 10^xend
#            }
# << VL no need to normalize the coordinates 
		}
		else {
			rect(figxy[1], breakpos, figxy[2], breakpos + yw, 
					col = bgcol, border = bgcol)
			xbegin <- c(figxy[1], figxy[1])
			ybegin <- c(breakpos, breakpos + yw)
			xend <- c(figxy[2], figxy[2])
			yend <- c(breakpos, breakpos + yw)
			if (xaxl) {
				xbegin <- 10^xbegin
				xend <- 10^xend
			}
		}
		par(xpd = TRUE)
	}
	else {
		rect(br[1], br[2], br[3], br[4], col = bgcol, border = bgcol)
		if (style == "slash") {
			if (axis == 1 || axis == 3) {
				xbegin <- c(breakpos - xw, breakpos)
				xend <- c(breakpos, breakpos + xw)
				ybegin <- c(br[2], br[2])
				yend <- c(br[4], br[4])
				if (xaxl) {
					xbegin <- 10^xbegin
					xend <- 10^xend
				}
			}
			else {
				xbegin <- c(br[1], br[1])
				xend <- c(br[3], br[3])
				ybegin <- c(breakpos - yw, breakpos)
				yend <- c(breakpos, breakpos + yw)
				if (yaxl) {
					ybegin <- 10^ybegin
					yend <- 10^yend
				}
			}
		}
		else {
			if (axis == 1 || axis == 3) {
				xbegin <- c(breakpos - xw/2, breakpos - xw/4, 
						breakpos + xw/4)
				xend <- c(breakpos - xw/4, breakpos + xw/4, breakpos + 
								xw/2)
				ybegin <- c(ifelse(yaxl, 10^figxy[3 + (axis == 
													3)], figxy[3 + (axis == 3)]), br[4], br[2])
				yend <- c(br[4], br[2], ifelse(yaxl, 10^figxy[3 + 
												(axis == 3)], figxy[3 + (axis == 3)]))
				if (xaxl) {
					xbegin <- 10^xbegin
					xend <- 10^xend
				}
			}
			else {
				xbegin <- c(ifelse(xaxl, 10^figxy[1 + (axis == 
													4)], figxy[1 + (axis == 4)]), br[1], br[3])
				xend <- c(br[1], br[3], ifelse(xaxl, 10^figxy[1 + 
												(axis == 4)], figxy[1 + (axis == 4)]))
				ybegin <- c(breakpos - yw/2, breakpos - yw/4, 
						breakpos + yw/4)
				yend <- c(breakpos - yw/4, breakpos + yw/4, breakpos + 
								yw/2)
				if (yaxl) {
					ybegin <- 10^ybegin
					yend <- 10^yend
				}
			}
		}
	}
	segments(xbegin, ybegin, xend, yend, col = breakcol, lty = 1)
	par(xpd = FALSE)
}




############################################################################################
# Draws two histograms on the same plot.
#
# x1, x2: the raw series whose histograms are desired.
# breaks: number of breaks (processed over both series).
# x.label: text label of the x axis.
# series.names: names of the series (two values text vector).
# leg.pos: position of the legend ("topleft", "bottomright", etc.)
############################################################################################
multi.hist <- function(x1, x2, breaks=10, x.label, series.names, leg.pos)
{	# compute common breaks
	data <- c(x1,x2)
	mn <- min(data)
	mx <- max(data)
	bk <- seq(from=mn,to=mx,by=(mx-mn)/breaks)
	
	# process bar sizes
	histo1 <- hist(x1, breaks=bk, plot=FALSE)
	histo2 <- hist(x2, breaks=bk, plot=FALSE)
	
	# process y limits
	y.lim <- c(0, max(c(histo1$counts,histo2$counts)))
	
	# set colors
	c1 <- add.alpha("BLUE", 0.25)
	c2 <- add.alpha("RED", 0.25)
	
	# create an empty plot
	plot(NULL, 
			xlab=x.label, ylab="Frequency", main="",
			ylim=y.lim, xlim=c(mn,mx),
	)
	
	# add the histograms
	plot(histo1, 
			col=c1, border=NA,
			add=TRUE
	)
	plot(histo2, 
			col=c2, border=NA, 
			add=TRUE
	)
	
	# add the legend
	legend(x=leg.pos,legend=series.names,
			inset=0.03,
			fill=c("BLUE","RED"))
}



############################################################################################
# Generates a bar plot comparing two centrality measures. The reference measure is used as
# a baseline, and to order the nodes on the x axis. The comparison measure is used to process
# the ranking difference with the reference measure, and the result appears as the bar heights.
#
# ref.vals: values for the reference measure (a numerical vector).
# comp.vals: values for the comparison measure (same).
# ref.measure: name of the reference measure.
# comp.measure: name of the comparison measure.
# alpha: parameter used to build the matrix (or possibly NA if none was defined).
# folder: folder in which to generate the plots.
# formats: format of the generated file ("PDF", "PNG", or both).
############################################################################################
rank.diff.barplot <- function(disc.vals, cont.vals, out.folder=NA, formats=c("pdf", "png"))
{	disc.measure <- "Discrete Approximation"
	cont.measure <- "Continuous average Straightness" 
	
	# set up node order
	disc.rk <- rank(disc.vals,ties.method="min")
	cont.rk <- rank(cont.vals,ties.method="min")
	diff <- cont.rk - disc.rk
	idx <- order(disc.vals, decreasing=TRUE)
	
	# setup file name
	plot.path <- file.path(out.folder,"rank-comparison")
	
	# record data
	data.file <- paste0(plot.path,".txt")
	data <- cbind(disc.rk,cont.rk,diff)
	colnames(data) <- c("DiscRank","ContRank","Difference")
	write.table(x=data,file=data.file,row.names=FALSE,col.names=TRUE)
	
	# record plot
	for(format in formats)
	{	if(!is.na(format))
		{	format <- tolower(format)
			plot.file <- paste0(plot.path,".",format)
			if(format=="pdf")
				pdf(file=plot.file,bg="white")
			else if(format=="png")
				png(filename=plot.file,width=800,height=800,units="px",pointsize=20,bg="white")
		}
		barplot(diff[idx], 
			col="RED",
			border="BLACK",
#			main="Rank changes between the continuous and discrete average Straightness values",
#			ylim=c(-length(disc.vals),length(disc.vals)),
			xlab="Nodes ordered by decreasing discrete average Straightness",
			ylab="Rank changes obtained with the continuous average Straightness"
		)
		
		if(!is.na(format))
			dev.off()
	}
}
