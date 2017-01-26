############################################################################
# Extracts road networks from OpenStreetMap.
# https://www.openstreetmap.org
# In order to obtain some clean graphs, this script was completed with some 
# manual processing: convert to Gephi, remove tendrils, then programmatically
# keep only the giant component (to remove all remaining noise).
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




data.folder <- "data"
urban.folder <- file.path(data.folder,"urban")




# list the cities and their boxes
cities <- list(
 		  abidjan=c( -4.0314,   5.3125,  -4.0100,   5.3359),	#TODO
#	 alicesprings=c(133.8220, -23.7248, 133.9050, -23.6639),
#	      avignon=c(  4.7669,  43.9240,   4.8455,  43.9619),
#	       beijin=c(116.2635,  39.8233, 116.4887,  39.9897),
 		 bordeaux=c( -0.6805,  44.7757,  -0.4968,  44.8924),	#TODO
 			dakar=c(-17.5338,  14.6434, -17.4185,  14.8104),	#TODO
 		 hongkong=c(114.1150,  22.1912, 114.2640,  22.2961),	#TODO
#	     istanbul=c( 28.9208,  40.9846,  29.0671,  41.0721),
#	   karlskrona=c( 15.5604,  56.1488,  15.6016,  56.1704),
 		   lisbon=c( -9.1810,  38.6973,  -9.0938,  38.7511),	#TODO
 		ljubljana=c( 14.4455,  46.0116,  14.5790,  46.0857),	#TODO
 		manhattan=c(-74.0224,  40.6994, -73.9219,  40.8510)		#TODO
##	      newyork=c(-74.0465,  40.5389, -73.7787,  40.9083),
#	         sfax=c( 10.6272,  34.6507,  10.8614,  34.8640),
#	     soustons=c( -1.4313,  43.7396,  -1.3158,  43.7793),
#		    tokyo=c(139.5895,  35.5378, 139.9136,  35.8590),
#	troisrivieres=c(-72.6351,  46.3104, -72.4902,  46.4083),
)




# process each city
for(c in 1:length(cities))
{	name <- names(cities)[c]
	cat("Processing city ",name,"\n",sep="")
	city.folder <- file.path(urban.folder,name)
	
	osm.file <- file.path(city.folder,"data.osm")
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
	V(g)$x <- city.sub$nodes$attrs[idx,"lon"]*1000			# multiply by 1000 so that the network can be browsed in Gephi
	V(g)$y <- city.sub$nodes$attrs[idx,"lat"]*1000
	
	# filter the nodes located out of the box
	idx <- which(V(g)$x<cities[[name]][1] || V(g)$y<cities[[name]][3] || V(g)$x>cities[[name]][3] || V(g)$y>cities[[name]][4])
	print(length(idx))
	g <- delete_vertices(graph=g, v=idx)

	plot(g,vertex.label=NA,vertex.size=1)
	
	# export as graphml
	#net.file <- paste0(file_path_sans_ext(osm.file),".graphml")
	net.file <- file.path(city.folder,"graph.graphml")
	write.graph(g, net.file, format="graphml")
}


# normalization (finally not necessary)
#c <- 1
#name <- names(cities)[c]
#city.folder <- file.path(urban.folder,name)
#net.file <- file.path(city.folder,"graph.graphml")
#g <- read.graph(net.file,format="graphml")
#temp <- c(V(g)$x,V(g)$y)
#V(g)$x <- (V(g)$x - min(temp)) / (max(temp) - min(temp)) * 100
#V(g)$y <- (V(g)$y - min(temp)) / (max(temp) - min(temp)) * 100
#net.file <- file.path(city.folder,"graph0.graphml")
#write.graph(g,net.file,format="graphml")
