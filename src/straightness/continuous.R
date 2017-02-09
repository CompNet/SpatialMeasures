############################################################################
# Processing the continuous average straightness, based on integrations.
# There are several variants, described in the article (see the readme file).
#
# Vincent Labatut 05/2016
#
# setwd("~/eclipse/workspaces/Networks/SpatialMeasures")
# setwd("c:/eclipse/workspaces/Networks/SpatialMeasures")
# source("src/straightness/continuous.R")
############################################################################
library("igraph")
#library("pracma") # test for another integration library (was not faster than the default method)




############################################################################
# Checks if the specified floats are equal, for a given tolerance. Useful for
# certain rounding bugs.
#
# x,y: the values to compare.
# tolerance: the tolerance (optional, default=10^-10).
############################################################################
TOLERANCE <- 1e-10
tol.eq <- function(x,y,tolerance=TOLERANCE)
{	abs(x-y)<tolerance
}



############################################################################
# Checks if the three specified points are aligned. 
#
# Note: uses function tol.eq to compare the points positions.
#
# x1,y1: coordinates of the first point.
# x2,y2: coordinates of the second point.
# x3,y3: coordinates of the third point.
#
# returns: TRUE iff the points are aligned (colinear).
############################################################################
check.alignment <- function(x1, y1, x2, y2, x3, y3)
{	disp <- FALSE
	
	result <- FALSE
	if(disp) cat("........(",x1,",",y1,") (",x2,",",y2,") (",x3,",",y3,")\n",sep="")
	if(disp) cat("........x1==x2: ",tol.eq(x1,x2)," x1==x3: ",tol.eq(x1,x3),"\n",sep="")	
	
	# possible the same points
	if(tol.eq(x1,x2) && tol.eq(y1,y2) 
			|| tol.eq(x1,x3) && tol.eq(y1,y3)
			|| tol.eq(x2,x3) && tol.eq(y2,y3))
		result <- TRUE
	
	# possibly vertical
	else if(tol.eq(x1,x2))
		result <- tol.eq(x1,x3)
	
	# not vertical
	else if(!tol.eq(x1,x3))
	{	slope1 <- (y1-y2)/(x1-x2)
		slope2 <- (y1-y3)/(x1-x3)
		result <- tol.eq(slope1,slope2)
		if(disp) cat("........slope1: ",(y1-y2)/(x1-x2)," slope2: ",(y1-y3)/(x1-x3)," slope1==slope2: ",abs(slope1-slope2)<1e-10,"\n",sep="")
	}
	
	return(result)
}



############################################################################
# Auxiliary functions processing the lambda_u and lambda_v break-even distances 
# for the two specified links.
# 
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# u2,v2: end nodes of the second link.
#
# returns: a vector containing lambda_u and lambda_v.
############################################################################
aux.process.lambdauv <- function(g.dist, u1, v1, u2, v2)
{	disp <- FALSE
	
	# process Euclidean distances
	#ed.u1v1 <- get.dist(u1,v1,e.dist)
	ed.u1v1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
	
	if(disp) cat("..........(u1,v1)=(",u1,",",v1,") vs (u2,v2)=(",u2,",",v2,") e.dist[u1,v1]=",ed.u1v1,"\n",sep="")
	if(disp) cat("..........g.dist[u2,v1]=",g.dist[u2,v1]," -- g.dist[u2,u1]=",g.dist[u2,u1]," g.dist[v2,v1]=",g.dist[v2,v1]," -- g.dist[v2,u1]=",g.dist[v2,u1],"\n",sep="")
	
	# apply the article formulas
	lambdau <- (g.dist[u2,v1] - g.dist[u2,u1] + ed.u1v1)/2
	lambdav <- (g.dist[v2,v1] - g.dist[v2,u1] + ed.u1v1)/2
	if(disp) cat("..........lambdau:",lambdau," lambdav:",lambdav,"\n",sep="")
	
	# must sometimes correct some rounding problems
	if(lambdau<0)
		lambdau <- 0
	else if(lambdau>ed.u1v1)
		result <- ed.u1v1
	if(lambdav<0)
		lambdav <- 0
	else if(lambdav>ed.u1v1)
		lambdav <- ed.u1v1
	
	result <- c(lambdau,lambdav)
	return(result)
}



############################################################################
# Auxiliary functions processing the lambda_2 threshold for specified links
# and point.
# 
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# u2,v2: end nodes of the second link.
# lambdau, lambdav: pre-processed break-even distances.
#
# returns: the lambda_2 break-even distance.
############################################################################
aux.process.lambda2 <- function(g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
{	disp <- FALSE
	
	# process Euclidean distances
	#ed.u1v1 <- get.dist(u1,v1,e.dist)
	ed.u1v1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
	#ed.u2v2 <- get.dist(u2,v2,e.dist)
	ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
	
	if(disp) cat("..........(u1,v1)=(",u1,",",v1,") vs (u2,v2)=(",u2,",",v2,") e.dist[u1,v1]=",ed.u1v1," e.dist[u2,v2]=",ed.u2v2,"\n",sep="")
	if(disp) cat("..........g.dist[u1,u2]=",g.dist[u1,u2]," -- g.dist[u1,v2]=",g.dist[u1,v2]," g.dist[v1,u2]=",g.dist[v1,u2]," -- g.dist[v1,v2]=",g.dist[v1,v2],"\n",sep="")
	
	# apply the article formulas
	result <- NA
	if(ellp1<=lambdau && ellp1<=lambdav)
		result <- (g.dist[u1,v2] - g.dist[u1,u2] + ed.u2v2)/2
	else if(ellp1<=lambdau && ellp1>lambdav)
		result <- (g.dist[v1,v2] - g.dist[u1,u2] + ed.u1v1 + ed.u2v2 - 2*ellp1)/2
	else if(ellp1>lambdau && ellp1<=lambdav)
		result <- (g.dist[u1,v2] - g.dist[v1,u2] - ed.u1v1 + ed.u2v2 + 2*ellp1)/2
	else if(ellp1>lambdau && ellp1>lambdav)
		result <- (g.dist[v1,v2] - g.dist[v1,u2] + ed.u2v2)/2
	if(disp) cat("..........lambda2:",result,"\n",sep="")
	
	# must sometimes correct some rounding problems
	if(result<0)
		result <- 0
	else if(result>ed.u2v2)
		result <- ed.u2v2
	
	return(result)
}



############################################################################
# Auxilary function used during the processing of the straightness, corresponds
# to the antiderivative of the f function from the article.
#
# a,b,c,d,e,f: constants used as parameters of the integrated function.
# ell: variable of the function. 
# 
# returns: value of the antiderivative needed to process the straightness.
############################################################################
aux.process.straightness.antiderivative <- function(a, b, c, d, e, f, ell)
{	disp <- FALSE
	
	# debug mode
	if(disp)
	{	cat("........a=",a,"; b=",b,"; c=",c,"; d=",d,"; e=",e,"; f=",f,"; ell=",sprintf("%.25f",ell),"\n",sep="")
	
		temp11 <- sqrt(b^2 + d^2)
		temp12 <- sqrt(b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2)
		temp1311 <- b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f
		temp1312 <- (a^2 + c^2)*f^2
		temp1313 <- (e + f*ell)
		cat("........    ","temp1311: ",temp1311," temp1312: ",temp1312," temp1313: ",temp1313,"\n",sep="")
		temp131 <- (temp1311+temp1312)^(3/2)*temp1313
		temp13 <- log(abs(temp131))
		temp1 <- temp11 * temp12 * temp13
		temp2 <- ((-(b^2*e) + a*b*f + d*(-(d*e) + c*f))
					* log(abs(2*(a*b + c*d + b^2*ell + d^2*ell + sqrt(b^2 + d^2)*sqrt(a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2)))))
		temp31 <- b^2 + d^2
		temp32 <- a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2
		temp33 <- b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2
		cat("........    ","temp31: ",temp31," temp32: ",temp32," temp33: ",temp33,"\n",sep="")
		temp341 <- -(a*b*e) - c*d*e + a^2*f + c^2*f - b^2*e*ell - d^2*e*ell + a*b*f*ell + c*d*f*ell
		temp342 <- b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2
		temp343 <- a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2
		temp34 <- (2*f^3*(temp341 + sqrt(temp342)*sqrt(temp343)))
		cat("........    ","temp341: ",temp341," temp342: ",temp342," temp343: ",temp343," sqrt(temp342)*sqrt(temp343): ",sqrt(temp342)*sqrt(temp343)," temp34: ",temp34,"\n",sep="")
		temp3 <- (sqrt(temp31)
					* (f*sqrt(temp32)
						- sqrt(temp33)
						* log(abs(temp34))))
		result <- (temp1 + temp2 + temp3) / (sqrt(b^2 + d^2)*f^2)
		cat("........    ","temp11: ",temp11," temp12: ",temp12," temp131: ",temp131," temp13: ",temp13," temp1: ",temp1," temp2: ",temp2," temp3: ",temp3," result: ",result,"\n",sep="")
		cat("........    ","result: ",result,"\n",sep="")
#		stop()
	}
	
	# normal mode
	else
	{	result <- 
			((sqrt(b^2 + d^2)
					* sqrt(b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2)
					* log(abs((b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2)^(3/2)*(e + f*ell)))
					+ (-(b^2*e) + a*b*f + d*(-(d*e) + c*f))
					* log(abs(2*(a*b + c*d + b^2*ell + d^2*ell + sqrt(b^2 + d^2)*sqrt(a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2))))
					+ sqrt(b^2 + d^2)
					* (f*sqrt(a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2)
						- sqrt(b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2)
						* log(abs(2*f^3*(-(a*b*e) - c*d*e + a^2*f + c^2*f - b^2*e*ell - d^2*e*ell + a*b*f*ell + c*d*f*ell 
													+ sqrt(b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2)
													* sqrt(a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2))))))
				/(sqrt(b^2 + d^2)*f^2))
	}
	
	return(result)
}



############################################################################
# Auxilary function used during the processing of the straightness, corresponds
# to the double integral of the f function performed in the article.
# 
# a0,b0,c0,d0,e0,f0,g1,h1,i1,g2,h2,i2: constants used as parameters of the integrated 
#								 function.
# fellp2: functions used to process ellp2.
# ellp2up: upper bound of ellp2.
# lb,ub: lower and upper bound of the integral to process.
#
# returns: value of the double integral needed to process the straightness.
############################################################################
aux.process.straightness.integral <- function(a0,b0, c0, d0,e0, f0, g1,h1, i1, g2,h2, i2, fellp2, ell2pup, lb, ub)
{	disp <- FALSE
	if(disp) cat("......lb: ",sprintf("%.23f",lb)," ub: ",sprintf("%.23f",ub),"\n",sep="")
	
	result <- 0
	if(!tol.eq(lb,ub))
	{	fun <- Vectorize(function(x)
				{	ell2mid <- fellp2(x)
					if(disp) cat("......for x=",x,", from 0 through ",ell2mid," to ",ell2pup,"\n",sep="")
					if(tol.eq(ell2mid,0))
						part1 <- 0
					else
					{	part11 <- aux.process.straightness.antiderivative(a0+b0*x,c0,d0+e0*x,f0,g1+h1*x,i1,ell2mid)
						part12 <- aux.process.straightness.antiderivative(a0+b0*x,c0,d0+e0*x,f0,g1+h1*x,i1,0)
						part1 <- part11 - part12
					}
					if(tol.eq(ell2pup,ell2mid))
						part2 <- 0
					else
					{	part21 <- aux.process.straightness.antiderivative(a0+b0*x,c0,d0+e0*x,f0,g2+h2*x,i2,ell2pup)
						part22 <- aux.process.straightness.antiderivative(a0+b0*x,c0,d0+e0*x,f0,g2+h2*x,i2,ell2mid)
						part2 <- part21 - part22
					}
					res <- part1 + part2
					if(disp) cat("......res = ",part1," + ",part2," = ",res,"\n",sep="")
					return(res)
				})
		
		# default integration method
		error.flag <- FALSE
		intres <- tryCatch(integrate(f=fun,lower=lb,upper=ub,abs.tol=1e-15), error=function(e) error.flag <<- TRUE)
		if(error.flag)
		{	cat("......",a0,",",b0,",",c0,",",d0,",",e0,",",f0,",",g1,",",h1,",",i1,",",g2,",",h2,",",i2,",fellp2,",ell2pup,",",lb,",",ub,"\n",sep="")
			print(fellp2)
#stop("Error during the numerical integration.)
		}
		else
		{	result <- intres$value
			if(disp) cat("......",result," vs. ",ell2pup*(ub-lb)," (",intres$abs.error,") \n",sep="")
		}
		# test using the pracma package (turned out not to be faster than base, so eventually not used)
		#result <- integral(fun=fun, xmin=lb, xmax=ub, method="Kronrod",vectorized=TRUE)
	}
	
	return(result)
}



############################################################################
# Processes the straightness between two points (not necessarily nodes).
# 
# graph: considered graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# u2,v2: end nodes of the second link.
# ellp2: position of p_2 on the second link.
# lambdau, lambdav: pre-processed break-even distances.
# lambda2: pre-processed break-even distance.
# 
# returns: straightness value between p_1 and p_2. 
############################################################################
aux.straightness.point.point <- function(graph, g.dist, u1, v1, ellp1, u2, v2, ellp2, lambdau, lambdav, lambda2)
{	# get the node coordinates
	#xu1 <- V(graph)$x[u1]
	#yu1 <- V(graph)$y[u1]
	#xu2 <- V(graph)$x[u2]
	#yu2 <- V(graph)$y[u2]
	#xv1 <- V(graph)$x[v1]
	#yv1 <- V(graph)$y[v1]
	#xv2 <- V(graph)$x[v2]
	#yv2 <- V(graph)$y[v2]
	# note: although supposedly equivalent, the below version is much faster than the above one
	xu1 <- vertex_attr(graph, name="x", index=u1)
	yu1 <- vertex_attr(graph, name="y", index=u1)
	xu2 <- vertex_attr(graph, name="x", index=u2)
	yu2 <- vertex_attr(graph, name="y", index=u2)
	xv1 <- vertex_attr(graph, name="x", index=v1)
	yv1 <- vertex_attr(graph, name="y", index=v1)
	xv2 <- vertex_attr(graph, name="x", index=v2)
	yv2 <- vertex_attr(graph, name="y", index=v2)
	
	# process Euclidean distances
	#ed.u1v1 <- get.dist(u1,v1,e.dist)
	ed.u1v1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
	#ed.u2v2 <- get.dist(u2,v2,e.dist)
	ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
	ed.p1p2 <- sqrt((xu2 + ellp2/ed.u2v2*(xv2-xu2) - xu1 - ellp1/ed.u1v1*(xv1-xu1))^2 
				+ (yu2 + ellp2/ed.u2v2*(yv2-yu2) - yu1 - ellp1/ed.u1v1*(yv1-yu1))^2)
	
	# process the straightness
	result <- NA
	# if the points are disconnected
	if(is.infinite(g.dist[u1,u2]))
		result <- 0
	# if p1 and p2 are lying on the same link
	else if(setequal(c(u1,v1),c(u2,v2)) 
			|| ellp1==0 && u1 %in% c(u2,v2)
			|| ellp1==ed.u1v1 && v1 %in% c(u2,v2)
			|| ellp2==0 && u2 %in% c(u1,v1)
			|| ellp2==ed.u2v2 && v2 %in% c(u1,v1))
		result <- 1
	# general case
	else
	{	# u1u2
		if(ellp1<=lambdau && ellp1<=lambdav && ellp2<=lambda2
			|| ellp1<=lambdau && ellp1>lambdav && ellp2<=lambda2)
			result <- ed.p1p2 / (ellp1 + g.dist[u1,u2] + ellp2)
		# u1v2
		else if(ellp1<=lambdau && ellp1<=lambdav && ellp2>lambda2
			|| ellp1>lambdau && ellp1<=lambdav && ellp2>lambda2)
			result <- ed.p1p2 / (ellp1 + g.dist[u1,v2] + ed.u2v2 - ellp2)
		# v1u2
		else if(ellp1>lambdau && ellp1<=lambdav && ellp2<=lambda2
			|| ellp1>lambdau && ellp1>lambdav && ellp2<=lambda2)
			result <- ed.p1p2 / (ed.u1v1 - ellp1 + g.dist[v1,u2] + ellp2)
		# v1v2
		else if(ellp1<=lambdau && ellp1>lambdav && ellp2>lambda2
			|| ellp1>lambdau && ellp1>lambdav && ellp2>lambda2)
			result <- ed.p1p2 / (ed.u1v1 - ellp1 + g.dist[v1,v2] + ed.u2v2 - ellp2)
	}
		
	return(result)
}



############################################################################
# Processes the total straightness between a point and a link (i.e. all the 
# points constituting this link).
# 
# graph: considered graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# u2,v2: end nodes of the second link.
# lambdau, lambdav: pre-processed break-even distances.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# 
# returns: total straightness between p_1 and (u_2,v_2). 
############################################################################
aux.total.straightness.point.link <- function(graph, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav, use.primitive=TRUE)
{	disp <- FALSE
	
	# process the appropriate lambda2
	lambda2 <- aux.process.lambda2(g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
	
	# get the node coordinates
	#xu1 <- V(graph)$x[u1]
	#yu1 <- V(graph)$y[u1]
	#xu2 <- V(graph)$x[u2]
	#yu2 <- V(graph)$y[u2]
	#xv1 <- V(graph)$x[v1]
	#yv1 <- V(graph)$y[v1]
	#xv2 <- V(graph)$x[v2]
	#yv2 <- V(graph)$y[v2]
	# note: although supposedly equivalent, the below version is much faster than the above one
	xu1 <- vertex_attr(graph, name="x", index=u1)
	yu1 <- vertex_attr(graph, name="y", index=u1)
	xu2 <- vertex_attr(graph, name="x", index=u2)
	yu2 <- vertex_attr(graph, name="y", index=u2)
	xv1 <- vertex_attr(graph, name="x", index=v1)
	yv1 <- vertex_attr(graph, name="y", index=v1)
	xv2 <- vertex_attr(graph, name="x", index=v2)
	yv2 <- vertex_attr(graph, name="y", index=v2)
	
	# process Euclidean distances
	#ed.u1v1 <- get.dist(u1,v1,e.dist)
	ed.u1v1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
	#ed.u2v2 <- get.dist(u2,v2,e.dist)
	ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
	
	# process the straightness
	result <- NA
	# if the points are disconnected
	if(is.infinite(g.dist[u1,u2]))
		result <- 0
	# if u_2 and v_2 are at the same position
	else if(tol.eq(xu2,xv2) && tol.eq(yu2,yv2)) 
		result <- 0
	# if p_1 lies on link (u_2,v_2)
	else if(setequal(c(u1,v1),c(u2,v2))
			|| ellp1==0 && u1 %in% c(u2,v2)
			|| ellp1==ed.u1v1 && v1 %in% c(u2,v2))
		result <- ed.u2v2
	# if p_1 is aligned with (u_2,v_2) and optimally connected to this link
	else if(tol.eq(ellp1,0) && check.alignment(xu1,yu1,xu2,yu2,xv2,yv2) && tol.eq(ed.u2v2,g.dist[u1,u2])
			|| tol.eq(ellp1,ed.u1v1) && check.alignment(xv1,yv1,xu2,yu2,xv2,yv2) && tol.eq(ed.u2v2,g.dist[v1,u2]))
			# this can be the case only if p_1 is u_1 or v_1
		result <- ed.u2v2
	# general case
	else
	{	# if the point is aligned with the link, we sometimes get some rounding problems (could not find why, exactly): 
		# in this case we need to use the numerical integration instead of the primitive
		if(use.primitive && check.alignment(xu1+ellp1/ed.u1v1*(xv1-xu1),yu1+ellp1/ed.u1v1*(yv1-yu1),xu2,yu2,xv2,yv2))
		{	if(disp)
			{	cat("......Point aligned with the link (and not optimally connected) >> no primitive\n")
				cat("......",xu1+ellp1/ed.u1v1*(xv1-xu1),",",yu1+ellp1/ed.u1v1*(yv1-yu1),",",xu2,",",yu2,",",xv2,",",yv2,"\n",sep="")
			}
			use.primitive <- FALSE
		}
		
		# integrate using the antiderivative
		if(use.primitive)
		{	if(disp) cat("......(u1,v1)=(",u1,",",v1,") - ellp1=",ellp1," - (u2,v2)=(",u2,",",v2,") - lambda2=",lambda2," - e.dist[u2,v2]:",ed.u2v2,"\n",sep="")
			if(disp) cat("......lambda_u=",lambdau," - lambda_v=",lambdav,"\n",sep="")
			
			# set the common parameters
			a <- xu2 - xu1 - ellp1*(xv1-xu1)/ed.u1v1
			b <- (xv2-xu2)/ed.u2v2
			c <- yu2 - yu1 - ellp1*(yv1-yu1)/ed.u1v1
			d <- (yv2-yu2)/ed.u2v2
			
			# set the variable parameters
			e1 <-  ellp1 + g.dist[u1,u2]	; e2 <- ellp1 + g.dist[u1,v2] + ed.u2v2	; e3 <- ed.u1v1 - ellp1 + g.dist[v1,u2]	; e4 <- ed.u1v1 - ellp1 + g.dist[v1,v2] + ed.u2v2
			f1 <- 1							; f2 <- -1								; f3 <- 1								; f4 <- -1
			
			# process the integral
			part1 <- 0
			part2 <- 0
start.time <- Sys.time()			
			if(lambdau <= lambdav)
			{	if(ellp1 <= lambdau)
				{	if(disp) cat("......case ellp1 <= lambdau, lambdav\n",sep="")
					if(!tol.eq(lambda2,0))
						part1 <- aux.process.straightness.antiderivative(a,b,c,d,e1,f1,lambda2) - aux.process.straightness.antiderivative(a,b,c,d,e1,f1,0)
					if(!tol.eq(lambda2,ed.u2v2))
						part2 <- aux.process.straightness.antiderivative(a,b,c,d,e2,f2,ed.u2v2) - aux.process.straightness.antiderivative(a,b,c,d,e2,f2,lambda2)
				}
				else if(ellp1 <= lambdav)
				{	if(disp) cat("......case lambdau < ellp1 <= lambdav\n",sep="")
					if(!tol.eq(lambda2,0))
						part1 <- aux.process.straightness.antiderivative(a,b,c,d,e3,f3,lambda2) - aux.process.straightness.antiderivative(a,b,c,d,e3,f3,0)
					if(!tol.eq(lambda2,ed.u2v2))
						part2 <- aux.process.straightness.antiderivative(a,b,c,d,e2,f2,ed.u2v2) - aux.process.straightness.antiderivative(a,b,c,d,e2,f2,lambda2)
				}
				else
				{	if(disp) cat("......case lambdau, lambdav < ellp1\n",sep="")
					if(!tol.eq(lambda2,0))
						part1 <- aux.process.straightness.antiderivative(a,b,c,d,e3,f3,lambda2) - aux.process.straightness.antiderivative(a,b,c,d,e3,f3,0)
					if(!tol.eq(lambda2,ed.u2v2))
						part2 <- aux.process.straightness.antiderivative(a,b,c,d,e4,f4,ed.u2v2) - aux.process.straightness.antiderivative(a,b,c,d,e4,f4,lambda2)
				}
			}
			else
			{	if(ellp1 <= lambdav)
				{	if(disp) cat("......case ellp1 <= lambdau, lambdav\n",sep="")
					if(!tol.eq(lambda2,0))
						part1 <- aux.process.straightness.antiderivative(a,b,c,d,e1,f1,lambda2) - aux.process.straightness.antiderivative(a,b,c,d,e1,f1,0)
					if(!tol.eq(lambda2,ed.u2v2))
						part2 <- aux.process.straightness.antiderivative(a,b,c,d,e2,f2,ed.u2v2) - aux.process.straightness.antiderivative(a,b,c,d,e2,f2,lambda2)
				}
				else if(ellp1 <= lambdau)
				{	if(disp) cat("......case lambdau < ellp1 <= lambdau\n",sep="")
					if(!tol.eq(lambda2,0))
						part1 <- aux.process.straightness.antiderivative(a,b,c,d,e1,f1,lambda2) - aux.process.straightness.antiderivative(a,b,c,d,e1,f1,0)
					if(!tol.eq(lambda2,ed.u2v2))
						part2 <- aux.process.straightness.antiderivative(a,b,c,d,e4,f4,ed.u2v2) - aux.process.straightness.antiderivative(a,b,c,d,e4,f4,lambda2)
				}
				else
				{	if(disp) cat("......case lambdau, lambdav < ellp1\n",sep="")
					if(!tol.eq(lambda2,0))
						part1 <- aux.process.straightness.antiderivative(a,b,c,d,e3,f3,lambda2) - aux.process.straightness.antiderivative(a,b,c,d,e3,f3,0)
					if(!tol.eq(lambda2,ed.u2v2))
						part2 <- aux.process.straightness.antiderivative(a,b,c,d,e4,f4,ed.u2v2) - aux.process.straightness.antiderivative(a,b,c,d,e4,f4,lambda2)
				}
			}
end.time <- Sys.time()
duration <- difftime(end.time,start.time,units="s")
#cat(duration,"")	
			
			# combine to get the result
			if(disp) cat("......part1:",part1," (max:",lambda2,") part2:",part2," (max:",ed.u2v2-lambda2,")\n",sep="")
			result <- part1 + part2
#z<-0;while(z<10e1){z<-z+1;result<-0}
		}
		
		# numerical integration (much slower)
		else
		{	# define the function to integrate
start.time <- Sys.time()
			fun <- Vectorize(function(x)
			{	res <- aux.straightness.point.point(graph, g.dist, u1, v1, ellp1, u2, v2, ellp2=x, lambdau, lambdav, lambda2)
				return(res)
			})
			# perform the approximate integration
			result <- integrate(f=fun,lower=0,upper=ed.u2v2,abs.tol=1e-15)$value
end.time <- Sys.time()
duration <- difftime(end.time,start.time,units="s")
#cat(duration,"")	
		}
	}
	
	if(disp) cat("......total.straightness.point.link:",result,"\n")
	return(result)
}



############################################################################
# Processes the average straightness between a point and a link (i.e. all the 
# points constituting this link).
# 
# graph: considered graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# u2,v2: end nodes of the second link.
# lambdau, lambdav: pre-processed break-even distances.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# 
# returns: average straightness between p_1 and (u_2,v_2). 
############################################################################
aux.mean.straightness.point.link <- function(graph, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav, use.primitive=TRUE)
{	disp <- FALSE
	
	# process Euclidean distances
	#ed.u2v2 <- get.dist(u2,v2,e.dist)
	ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
	
	# get the total straightness between the point and the link
	total <- aux.total.straightness.point.link(graph, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav, use.primitive)
	if(disp) cat("....u1:",u1," v1:",v1," ellp1:",ellp1," -- u2:",u2," v2:",v2," -- Total:",total," e.dist[u2,v2]:",ed.u2v2,"\n",sep="")
	
	# normalize to get the mean
	result <- total / ed.u2v2
	
	return(result)
}



############################################################################
# Processes the average straightness between each specified node and a link 
# (i.e. all the points constituting this link).
# 
# graph: considered graph.
# u: vector of nodes of interest.
# e: link of interest.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# g.dist: graph distance between all pairs of nodes. If missing, it is processed 
#		  by the function (processing once before calling this function can 
#		  therefore speed up the processing, if there are several calls).
# 
# returns: vector containing the average straightness between each node of u
#          and the link e=(u_2,v_2). 
############################################################################
mean.straightness.nodes.link <- function(graph, u=1:vcount(graph), e, use.primitive=TRUE, g.dist)
{	disp <- FALSE
	
	# init the result vector
	result <- c()
	
	# init the distances
	if(missing(g.dist))
	{	eatt <- list.edge.attributes(graph)
		if("dist" %in% eatt)
			g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
		else
			stop("You need to either provide the graph distances (using parameter g.dist) or a graph possessing an edge attribute called \"dist\"")
	}
	
	# get the second link
	el <- get.edgelist(graph)
	u2 <- el[e,1]
	v2 <- el[e,2]
	
	# process each specified node
	for(i in 1:length(u))
	{	if(disp) cat("..Processing node ",u[i]," vs. link (",u2,",",v2,")\n",sep="")
		
		# if the node is an isolate, no need to go further
		if(degree(graph,u[i])==0)
			str <- 0
		
		# otherwise, we process the mean straightness for the node
		else
		{	# get the relative position
			tmp <- as_ids(neighbors(graph,v=u[i]))[1]
			if(tmp<u[i])
			{	u1 <- tmp
				v1 <- u[i]
				#ellp1 <- get.dist(u1,v1,e.dist)
				ellp1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
			}
			else
			{	u1 <- u[i]
				v1 <- tmp
				ellp1 <- 0
			}
			
			# process the lambdas
			temp <- aux.process.lambdauv(g.dist, u1, v1, u2, v2)
			lambdau <- temp[1]
			lambdav <- temp[2]
			
			# get the mean straightness between the point and the link
			str <- aux.mean.straightness.point.link(graph, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav, use.primitive)
			if(disp) cat("..Avg straightness:",str,"\n",sep="")
		}
		
		# add to the result vector
		result <- c(result,str)
	}
	
	return(result)
}



############################################################################
# Processes the average straightness between a point and the graph (i.e. all 
# the points constituting the graph).
# 
# graph: considered graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# ellp1: position of p_1 on the first link.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# 
# returns: average straightness between p_1 and all the points constituting the graph. 
############################################################################
aux.mean.straightness.point.graph <- function(graph, g.dist, u1, v1, ellp1, use.primitive=TRUE)
{	disp <- FALSE
	
	total.str <- 0
	total.lgt <- 0
	
	# process each link
	el <- get.edgelist(graph)
	for(i in 1:nrow(el))
	{	# get the second link
		u2 <- el[i,1]
		v2 <- el[i,2]
		if(disp) cat("..(",u1,",",v1," ;",ellp1,") vs (",u2,",",v2,") (link #",i,"/",nrow(el),")\n",sep="")
		
		# get its length
		#ed.u2v2 <- get.dist(u2,v2,e.dist)
		ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
		
		# process the lambdas
		temp <- aux.process.lambdauv(g.dist, u1, v1, u2, v2)
		lambdau <- temp[1]
		lambdav <- temp[2]
		
		# get the mean straightness between the point and the link
		start.time <- Sys.time()
			part.str <- aux.total.straightness.point.link(graph, g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav, use.primitive)
		end.time <- Sys.time()
		duration <- difftime(end.time,start.time,units="s")
#if(is.infinite(part.str))
#	stop("..Infinite total straightness in aux.total.straightness.point.link")
		if(disp) cat("part.str: ",part.str," e.dist[u2,v2]:",ed.u2v2," (duration: ",duration,"s)\n",sep="")
		total.str <- total.str + part.str
		total.lgt <- total.lgt + ed.u2v2
	}
	
	# normalize to get the mean
	result <- total.str / total.lgt
	
	return(result)
}



############################################################################
# Processes the average straightness between each one of the specified nodes 
# and the graph (i.e. all the points constituting the graph). 
# 
# This can be interpreted as a centrality measure corresponding to the 
# average accessibility of the node from (or to) the rest of the graph.
# 
# graph: considered graph.
# u: vector of nodes of interest.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# g.dist: graph distance between all pairs of nodes. If missing, it is processed 
#		  by the function (processing once before calling this function can 
#		  therefore speed up the processing, if there are several calls).
# 
# returns: vector containing the average straightness between each node of u
#          and the points constituting the graph.
############################################################################
mean.straightness.nodes.graph <- function(graph, u=1:vcount(graph), use.primitive=TRUE, g.dist)
{	disp <- FALSE
	
	# init the result vector
	result <- c()
	
	# init the distances
	if(missing(g.dist))
	{	if(disp) cat("The graph distance is missing: we must process it\n",sep="")
		eatt <- list.edge.attributes(graph)
		if("dist" %in% eatt)
			g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
		else
			stop("You need to either provide the graph distances (using parameter g.dist) or a graph possessing an edge attribute called \"dist\"")
	}
	
	# process each specified node
	for(i in 1:length(u))
	{	if(disp) cat("Processing node ",u[i],"\n",sep="")
		
		# if the node is an isolate, no need to go further
		if(degree(graph,u[i])==0)
			str <- 0
		
		# otherwise, we process the mean straightness for the node
		else
		{	# get the relative position
			tmp <- as_ids(neighbors(graph,u[i]))[1]
			if(tmp<u[i])
			{	u1 <- tmp
				v1 <- u[i]
				#ellp1 <- get.dist(u1,v1,e.dist)
				ellp1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
			}
			else
			{	u1 <- u[i]
				v1 <- tmp
				ellp1 <- 0
			}
			
			# get the mean straightness between the point and the graph
			start.time <- Sys.time()
				str <- aux.mean.straightness.point.graph(graph, g.dist, u1, v1, ellp1, use.primitive)
			end.time <- Sys.time()
			duration <- difftime(end.time,start.time,units="s")
			if(disp) cat("Avg Straightness: ",str," (duration: ",duration,"s)\n",sep="")
#if(is.infinite(str))
#	stop("Error: infinite straightness")
		}
		
		# add to the result vector
		result <- c(result,str)
	}
	
	return(result)
}



############################################################################
# Processes the total straightness between two links (i.e. all the pairs of points 
# such that each point of the pair is located on a different link).
# 
# graph: considered graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# u2,v2: end nodes of the second link.
# lambdau, lambdav: pre-processed break-even distances.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# 
# returns: total straightness between (u_1,v_1) and (u_2,v_2). 
############################################################################
aux.total.straightness.link.link <- function(graph, g.dist, u1, v1, u2, v2, lambdau, lambdav, use.primitive=TRUE)
{	disp <- FALSE
	if(disp) cat("..Processing links (",u1,",",v1,") and (",u2,",",v2,")\n",sep="")
	
	# get the node coordinates
	#xu1 <- V(graph)$x[u1]
	#yu1 <- V(graph)$y[u1]
	#xu2 <- V(graph)$x[u2]
	#yu2 <- V(graph)$y[u2]
	#xv1 <- V(graph)$x[v1]
	#yv1 <- V(graph)$y[v1]
	#xv2 <- V(graph)$x[v2]
	#yv2 <- V(graph)$y[v2]
	# note: although supposedly equivalent, the below version is much faster than the above one
	xu1 <- vertex_attr(graph, name="x", index=u1)
	yu1 <- vertex_attr(graph, name="y", index=u1)
	xu2 <- vertex_attr(graph, name="x", index=u2)
	yu2 <- vertex_attr(graph, name="y", index=u2)
	xv1 <- vertex_attr(graph, name="x", index=v1)
	yv1 <- vertex_attr(graph, name="y", index=v1)
	xv2 <- vertex_attr(graph, name="x", index=v2)
	yv2 <- vertex_attr(graph, name="y", index=v2)
	
	# process Euclidean distances
	#ed.u1v1 <- get.dist(u1,v1,e.dist)
	ed.u1v1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
	#ed.u2v2 <- get.dist(u2,v2,e.dist)
	ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
	
	# process the straightness
	result <- NA
	# if the links are disconnected
	if(is.infinite(g.dist[u1,u2]))
		result <- 0
	# if the links are the same
	else if(setequal(c(u1,v1),c(u2,v2)))
		result <- (ed.u1v1^2)/2
	# if u_1 and v_1 are at the same position
	else if(tol.eq(xu1,xv1) && tol.eq(yu1,yv1))
		result <- ed.u2v2
	# if u_2 and v_2 are at the same position
	else if(tol.eq(xu2,xv2) && tol.eq(yu2,yv2)) 
		result <- ed.u1v1
	# if both links are aligned and connected by optimal paths
	else if(check.alignment(xu1,yu1,xv1,yv1,xu2,yu2) && check.alignment(xu1,yu1,xv1,yv1,xv2,yv2) && tol.eq(ed.u2v2,g.dist[u1,u2]))
		result <- ed.u1v1*ed.u2v2
	# general case
	else
	{	# if the point is aligned with the link, we sometimes get some rounding problems (could not find why, exactly): 
		# in this case we need to use the numerical integration
		if(use.primitive && check.alignment(xu1,yu1,xv1,yv1,xu2,yu2) && check.alignment(xu1,yu1,xv1,yv1,xv2,yv2))
		{	if(disp) cat("......Both links are aligned (and not optimally connected) >> no primitive\n")
			use.primitive <- FALSE
		}
		
		# use the pre-processed antiderivative of function f (faster)
		if(use.primitive)
		{	# check if the links have a common end-node
			#if(length(intersect(c(u1,v1),c(u2,v2)))>0)
			#	cat("WARNING: Links (",u1,",",v1,") and (",u2,",",v2,") have a common end-node\n")
			
			# possibly round the lambdas
			if(tol.eq(lambdau,0))
				lambdau <- 0
			else if(tol.eq(lambdau,ed.u1v1))
				lambdau <- ed.u1v1
			if(tol.eq(lambdav,0))
				lambdav <- 0
			else if(tol.eq(lambdav,ed.u1v1))
				lambdav <- ed.u1v1
			
			# function used to process lambda2
			fellp2 <- function(ellp1) aux.process.lambda2(g.dist, u1, v1, ellp1, u2, v2, lambdau, lambdav)
			
			# common constants
			a0 <- xu2 - xu1
			b0 <- (xu1 - xv1) / ed.u1v1
			c0 <- (xv2 - xu2) / ed.u2v2
			d0 <- yu2 - yu1
			e0 <- (yu1 - yv1) / ed.u1v1
			f0 <- (yv2 - yu2) / ed.u2v2
			
			# specific cases
			gu1u2 <- g.dist[u1,u2]						; hu1u2 <- 1	; iu1u2 <- 1
			gu1v2 <- g.dist[u1,v2] + ed.u2v2			; hu1v2 <- 1	; iu1v2 <- -1
			gv1u2 <- g.dist[v1,u2] + ed.u1v1			; hv1u2 <- -1	; iv1u2 <- 1
			gv1v2 <- g.dist[v1,v2] + ed.u1v1 + ed.u2v2	; hv1v2 <- -1	; iv1v2 <- -1
			
			part1 <- 0 ; part2 <- 0 ; part3 <- 0
			if(lambdau <= lambdav)
			{	# from 0 to lambda_u
				part1 <- aux.process.straightness.integral(a0,b0, c0, d0,e0, f0, g1=gu1u2,h1=hu1u2, i1=iu1u2, g2=gu1v2,h2=hu1v2, i2=iu1v2, fellp2=fellp2, ell2pup=ed.u2v2, lb=0,       ub=lambdau)
				# from lambda_u to lambda_v
				part2 <- aux.process.straightness.integral(a0,b0, c0, d0,e0, f0, g1=gv1u2,h1=hv1u2, i1=iv1u2, g2=gu1v2,h2=hu1v2, i2=iu1v2, fellp2=fellp2, ell2pup=ed.u2v2, lb=lambdau, ub=lambdav)
				# from lambda_v to ||u1v1||
				part3 <- aux.process.straightness.integral(a0,b0, c0, d0,e0, f0, g1=gv1u2,h1=hv1u2, i1=iv1u2, g2=gv1v2,h2=hv1v2, i2=iv1v2, fellp2=fellp2, ell2pup=ed.u2v2, lb=lambdav, ub=ed.u1v1)
			}
			else
			{	# from 0 to lambda_v
				part1 <- aux.process.straightness.integral(a0,b0, c0, d0,e0, f0, g1=gu1u2,h1=hu1u2, i1=iu1u2, g2=gu1v2,h2=hu1v2, i2=iu1v2, fellp2=fellp2, ell2pup=ed.u2v2, lb=0,       ub=lambdav)
				# from lambda_v to lambda_u
				part2 <- aux.process.straightness.integral(a0,b0, c0, d0,e0, f0, g1=gu1u2,h1=hu1u2, i1=iu1u2, g2=gv1v2,h2=hv1v2, i2=iv1v2, fellp2=fellp2, ell2pup=ed.u2v2, lb=lambdav, ub=lambdau)
				# from lambda_u to ||u1v1||
				part3 <- aux.process.straightness.integral(a0,b0, c0, d0,e0, f0, g1=gv1u2,h1=hv1u2, i1=iv1u2, g2=gv1v2,h2=hv1v2, i2=iv1v2, fellp2=fellp2, ell2pup=ed.u2v2, lb=lambdau, ub=ed.u1v1)
			}
			result <- part1 + part2 + part3
		}
		# use an approximate numerical integration instead of the antiderivative
		else
		{	# define the function to integrate
			fun <- Vectorize(function(x)
					{	res <- aux.total.straightness.point.link(graph, g.dist, u1, v1, ellp1=x, u2, v2, lambdau, lambdav, use.primitive)
						return(res)
					})
			
			# define the integral bounds
			b1 <- min(lambdau,lambdav)
			b2 <- max(lambdau,lambdav)
			
			# perform the approximate integration
			res1 <- integrate(f=fun,lower=0,upper=b1,abs.tol=1e-15)$value
			res2 <- integrate(f=fun,lower=b1,upper=b2,abs.tol=1e-15)$value
			res3 <- integrate(f=fun,lower=b2,upper=ed.u1v1,abs.tol=1e-15)$value
			result <- res1 + res2 + res3
		}
	}
	
	if(disp) cat("..Result:",result,"\n",sep="")
	return(result)
}
	


############################################################################
# Processes the average straightness between two links (i.e. all the pairs of 
# points such that each point of the pair is located on a different link).
# 
# graph: considered graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the first link.
# u2,v2: end nodes of the second link.
# lambdau, lambdav: pre-processed break-even distances.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# 
# returns: average straightness between (u_1,v_1) and (u_2,v_2). 
############################################################################
aux.mean.straightness.link.link <- function(graph, g.dist, u1, v1, u2, v2, lambdau, lambdav, use.primitive=TRUE)
{	result <- NA
	
	# if the links are the same
	if(setequal(c(u1,v1),c(u2,v2)))
		result <- 1
	# general case
	else
	{	# get the total straightness between the links
		total <- aux.total.straightness.link.link(graph, g.dist, u1, v1, u2, v2, lambdau, lambdav, use.primitive)
		# normalize to get the mean
		#ed.u1v1 <- get.dist(u1,v1,e.dist)
		ed.u1v1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
		#ed.u2v2 <- get.dist(u2,v2,e.dist)
		ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
		result <- total / (ed.u2v2*ed.u1v1)
	}
	
	return(result)
}



############################################################################
# Processes the average straightness between two links (i.e. all the pairs of 
# points such that each point of the pair is located on a different link).
# 
# graph: considered graph.
# e1: id of the first link. 
# e2: id of the second link. 
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# g.dist: graph distance between all pairs of nodes. If missing, it is processed 
#		  by the function (processing once before calling this function can 
#		  therefore speed up the processing, if there are several calls).
# 
# returns: average straightness between e_1=(u_1,v_1) and e_2=(u_2,v_2). 
############################################################################
mean.straightness.link.link <- function(graph, e1, e2, use.primitive=TRUE, g.dist)
{	disp <- FALSE
	
	result <- NA
	el <- get.edgelist(graph)
	
	# init the distances
	if(missing(g.dist))
	{	eatt <- list.edge.attributes(graph)
		if("dist" %in% eatt)
			g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
		else
			stop("You need to either provide the graph distances (using parameter g.dist) or a graph possessing an edge attribute called \"dist\"")
	}
	
	# get the first link end-nodes
	u1 <- el[e1,1]
	v1 <- el[e1,2]
	# get the second link end-nodes
	u2 <- el[e2,1]
	v2 <- el[e2,2]
	cat("Processing link ",e1,"=(",u1,",",v1,") vs. ",e2,"=(",u2,",",v2,")\n",sep="")
	
	# process the lambdas
	temp <- aux.process.lambdauv(g.dist, u1, v1, u2, v2)
	lambdau <- temp[1]
	lambdav <- temp[2]
	
	# process the straightness
	result <- aux.mean.straightness.link.link(graph, g.dist, u1, v1, u2, v2, lambdau, lambdav, use.primitive)
	return(result)
}



############################################################################
# Processes the average straightness between a link and the whole graph (i.e. 
# all the pairs of points such that one point is on the link and the other is
# on the graph).
# 
# graph: considered graph.
# g.dist: pre-processed graph distances for all pairs of nodes in the graph.
# u1,v1: end nodes of the link of interest.
# u2,v2: end nodes of the second link.
# exclude.self: whether or not the straightness between the link and itself 
#				should be considered in the average.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# 
# returns: average straightness between (u_1,v_1) and the rest of the graph.
############################################################################
aux.mean.straightness.link.graph <- function(graph, g.dist, u1, v1, exclude.self=FALSE, use.primitive=TRUE)
{	el <- get.edgelist(graph)
	
	# get the length of the first link
	#ed.u1v1 <- get.dist(u1,v1,e.dist)
	ed.u1v1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
	
	# process all links (possibly excluding (u_1,v_1)
	total.str <- 0
	total.length <- 0
	for(i in 1:nrow(el))
	{	# get the second link end-nodes
		u2 <- el[i,1]
		v2 <- el[i,2]
		
		# get its length
		#ed.u2v2 <- get.dist(u2,v2,e.dist)
		ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
		
		# process the lambdas
		temp <- aux.process.lambdauv(g.dist, u1, v1, u2, v2)
		lambdau <- temp[1]
		lambdav <- temp[2]
		
		# check if the link is the same than (u_1,v_1)
		if(!setequal(c(u1,v1),c(u2,v2)) || !exclude.self)
		{	# process the total straightness
			str <- aux.total.straightness.link.link(graph, g.dist, u1, v1, u2, v2, lambdau, lambdav, use.primitive)
			# add to result
			total.str <- total.str + str
			total.length <- total.length + ed.u1v1*ed.u2v2
		}
	}
	
	# process the normalizing value
	if(exclude.self)
		denom <- total.length
	else
		denom <- total.length - (ed.u1v1^2)/2
	
	# process the average
	result <- total.str / denom
	return(result)
}



############################################################################
# Processes the average straightness between each one of the specified links 
# and the graph (i.e. all the points constituting the graph). 
# 
# This can be interpreted as a centrality measure corresponding to the 
# average accessibility of the link from (or to) the rest of the graph.
# 
# graph: considered graph.
# e: vector of nodes of interest.
# exclude.self: whether or not the straightness between one link and itself 
#				should be considered in the average.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# g.dist: graph distance between all pairs of nodes. If missing, it is processed 
#		  by the function (processing once before calling this function can 
#		  therefore speed up the processing, if there are several calls).
# 
# returns: vector containing the average straightness between each node of u
#          and the points constituting the graph.
############################################################################
mean.straightness.links.graph <- function(graph, e=1:ecount(graph), exclude.self=FALSE, use.primitive=TRUE, g.dist)
{	# init the result vector
	result <- c()
	el <- get.edgelist(graph)
	
	# init the distances
	if(missing(g.dist))
	{	eatt <- list.edge.attributes(graph)
		if("dist" %in% eatt)
			g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
		else
			stop("You need to either provide the graph distances (using parameter g.dist) or a graph possessing an edge attribute called \"dist\"")
	}
	
	# process averages based on the previous totals
	for(i in e)
	{	# get the link end-nodes
		u1 <- el[i,1]
		v1 <- el[i,2]
		
		# get the total straightness between the point and the link
		str <- aux.mean.straightness.link.graph(graph, g.dist, u1, v1, exclude.self, use.primitive)
		
		# add to the result vector
		result <- c(result,str)
	}
	
	return(result)
}



############################################################################
# Processes the average straightness between each one pair of points constituting
# the graph of interest. 
# 
# graph: considered graph.
# exclude.self: whether or not the straightness between one link and itself 
#				should be considered in the average.
# use.primitive: whether to use the pre-processed primitive (faster), or to
# 				 to integrate numerically.
# g.dist: graph distance between all pairs of nodes. If missing, it is processed 
#		  by the function (processing once before calling this function can 
#		  therefore speed up the processing, if there are several calls).
# 
# returns: a single value representing the average straightness for the whole
#		   graph.
############################################################################
mean.straightness.graph <- function(graph, exclude.self=FALSE, use.primitive=TRUE, g.dist)
{	result <- NA
	el <- get.edgelist(graph)
	
	# init the distances
	if(missing(g.dist))
	{	eatt <- list.edge.attributes(graph)
		if("dist" %in% eatt)
			g.dist <- shortest.paths(graph=graph, weights=E(graph)$dist)
		else
			stop("You need to either provide the graph distances (using parameter g.dist) or a graph possessing an edge attribute called \"dist\"")
	}
	
	# process all pairs of links
	total.str <- 0
	total.length <- 0
	for(i in 1:nrow(el))
	{	# get the first link end-nodes
		u1 <- el[i,1]
		v1 <- el[i,2]
		
		# get its length
		#ed.u1v1 <- get.dist(u1,v1,e.dist)
		ed.u1v1 <- g.dist[u1,v1]	# equivalent (but faster) to the above line, since the nodes are neighbors
		
		for(j in i:nrow(el))
		{	# check if the links are the same one
			if(i!=j || !exclude.self)
			{	# get the second link end-nodes
				u2 <- el[j,1]
				v2 <- el[j,2]
				
				# get its length
				#ed.u2v2 <- get.dist(u2,v2,e.dist)
				ed.u2v2 <- g.dist[u2,v2]	# equivalent (but faster) to the above line, since the nodes are neighbors
				
				# process the lambdas
				temp <- aux.process.lambdauv(g.dist, u1, v1, u2, v2)
				lambdau <- temp[1]
				lambdav <- temp[2]
				
				# process straightness
				str <- aux.total.straightness.link.link(graph, g.dist, u1, v1, u2, v2, lambdau, lambdav, use.primitive)
				# add to result
				total.str <- total.str + str
				total.length <- total.length + ed.u1v1*ed.u2v2
			}
		}
	}
	
	# possibly substract superfluous values
	if(exclude.self)
		denom <- total.length
	else
	{	denom <- total.length
		#temp <- sapply(1:nrow(el), function(i) (get.dist(el[i,1],el[i,2],e.dist)^2)/2)
		temp <- sapply(1:nrow(el), function(i) (g.dist[el[i,1],el[i,2]]^2)/2)	# equivalent (but faster) to the above line, since the nodes are neighbors
		denom <- denom - sum(temp)
	}
	
	# process the average
	result <- total.str / denom
	return(result)
}
