############################################################################
# Extracts road networks from OpenStreetMap.
# https://www.openstreetmap.org
#
# Vincent Labatut 01/2017
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/misc/extraction.R")
############################################################################
library("tools")
library("igraph")
library("osmar")
# https://cran.r-project.org/web/packages/osmar/index.html
# version 1.1-7




data.folder <- "data"
urban.folder <- file.path(data.folder,"urban")




# list the cities and their boxes
cities <- list(
#	      avignon=c(  4.7669, 43.9240,   4.8455, 43.9619),
#	avignon-small=c(  4.7969, 43.9420,   4.8203, 43.9543),
#	         sfax=c( 10.6272, 34.6507,  10.8614, 34.8640),
#	       beijin=c(116.2635, 39.8233, 116.4887, 39.9897),
#	     istanbul=c( 28.9208, 40.9846,  29.0671, 41.0721),
	troisrivieres=c(-72.6351, 46.3104, -72.4902, 46.4083),
		    tokyo=c(139.5895, 35.5378, 139.9136, 35.8590),
	      newyork=c(-74.0465, 40.5389, -73.7787, 40.9083)
)




# process each city
for(c in 1:length(cities))
{	name <- names(cities)[c]
	cat("Processing city ",name,"\n",sep="")
	
	osm.file <- file.path(urban.folder,paste0(name,".osm"))
#	src <- osmsource_api()
	src <- osmsource_file(osm.file)

	box <- corner_bbox(cities[[name]][1], cities[[name]][2], cities[[name]][3], cities[[name]][4])	#large
	city <- get_osm(box, source=src)
	
	# extract the graph of streets
	street.ids <- find(city, way(tags(k=="highway")))		# get the ids of the highways
	#city.sub <- subset(city, way_ids=street.ids)			# get the subset of the data containing these streets
	#street.ids <- find(city.sub, way(tags(k=="name")))		# in this subset, find the named streets
	nal.ids <- find_down(city, way(street.ids))				# get the ids of the ways and nodes constituting the named streets
	city.sub <- subset(city, ids=nal.ids)					# get the subset of the data containing these streets and nodes
	
	# convert to an igraph object (with spatial coordinates)
	g <- as_igraph(city.sub)
	g <- as.undirected(graph=g,mode="collapse")
	idx <- match(V(g)$name,city.sub$nodes$attrs[,"id"])
	V(g)$x <- city.sub$nodes$attrs[idx,"lon"]*1000
	V(g)$y <- city.sub$nodes$attrs[idx,"lat"]*1000
	
	# filter the nodes located out of the box
	idx <- which(V(g)$x<cities[[name]][1] || V(g)$y<cities[[name]][3] || V(g)$x>cities[[name]][3] || V(g)$y>cities[[name]][4])
	print(length(idx))
	g <- delete_vertices(graph=g, v=idx)

	plot(g,vertex.label=NA,vertex.size=1)
	
	# export as graphml
	net.file <- paste0(file_path_sans_ext(osm.file),".graphml")
	write.graph(g, net.file, format="graphml")
}
