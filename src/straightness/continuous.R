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


source("src/common/misc.R")
source("src/common/transformations.R")


############################################################################
# Auxilary function used during the processing of the exact straightness.
#
# a,b,c,d,e,f: constants used as parameters of the integrated function.
# ell: variable of the function. 
# 
# returns: value of the antiderivative needed to process the exact straightness.
############################################################################
point.unilateral.straightness.antiderivative <- function(a,b,c,d,e,f,ell)
{	
#	cat("a=",a,"; b=",b,"; c=",c,"; d=",d,"; e=",e,"; f=",f,"; ell=",sprintf("%.25f",ell),"\n",sep="")
#	
#	# convert to higher precision
#	a <- mpfr(a, 64)
#	b <- mpfr(b, 64)
#	c <- mpfr(c, 64)
#	d <- mpfr(d, 64)
#	e <- mpfr(e, 64)
#	f <- mpfr(f, 64)
#	ell <- mpfr(ell, 64)
#	
#	# process the value
#	temp11 <- sqrt(b^2 + d^2)
#	temp12 <- sqrt(b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2)
#	temp1311 <- b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f
#	temp1312 <- (a^2 + c^2)*f^2
#	temp1313 <- (e + f*ell)
#	cat("\t","temp1311: ",as.numeric(temp1311)," temp1312: ",as.numeric(temp1312)," temp1313: ",as.numeric(temp1313),"\n",sep="")
##	cat("\t","temp1311: ",temp1311," temp1312: ",temp1312," temp1313: ",temp1313,"\n",sep="")
#	temp131 <- (temp1311+temp1312)^(3/2)*temp1313
#	temp13 <- log(abs(temp131))
#	temp1 <- temp11 * temp12 * temp13
#	temp2 <- ((-(b^2*e) + a*b*f + d*(-(d*e) + c*f))
#				* log(abs(2*(a*b + c*d + b^2*ell + d^2*ell + sqrt(b^2 + d^2)*sqrt(a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2)))))
#	temp3log <- (2*f^3*(-(a*b*e) - c*d*e + a^2*f + c^2*f - b^2*e*ell - d^2*e*ell + a*b*f*ell + c*d*f*ell 
#								+ sqrt(b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2)
#								* sqrt(a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2)))
#	cat("temp3log: ",as.numeric(temp3log),"\n",sep="")
#	temp3 <- (sqrt(b^2 + d^2)
#				* (f*sqrt(a^2 + c^2 + 2*a*b*ell + 2*c*d*ell + (b^2 + d^2)*ell^2)
#					- sqrt(b^2*e^2 + d^2*e^2 - 2*a*b*e*f - 2*c*d*e*f + (a^2 + c^2)*f^2)
#					* log(abs(temp3log))))
#	result <- (temp1 + temp2 + temp3) / (sqrt(b^2 + d^2)*f^2)
#	cat("\t"," temp11: ",as.numeric(temp11)," temp12: ",as.numeric(temp12)," temp131: ",as.numeric(temp131)," temp13: ",as.numeric(temp13)," temp1: ",as.numeric(temp1)," temp2: ",as.numeric(temp2)," temp3: ",as.numeric(temp3)," result: ",as.numeric(result),"\n",sep="")
##	cat("\t"," temp11: ",temp11," temp12: ",temp12," temp131: ",temp131," temp13: ",temp13," temp1: ",temp1," temp2: ",temp2," temp3: ",temp3," result: ",result,"\n",sep="")
#	result <- as.numeric(result)
	
	result <- 
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
	
#	cat("\tresult: ",result,"\n",sep="")
	return(result)
}


############################################################################
# Auxilary function used during the processing of the exact straightness.
#
# xa,ya: coordinates of A, the first end point of the considered link.
# xb,yb: coordinates of B, the second end point of the considered link.
# xc,yc: coordinates of C, the node considered as the network center.
# dgca: geodesic distance between C and A.
# dgcb: geodesic distance between C and B.
# ddab: Euclidean distance between A and B.
# 
# returns: total unilateral straightness over the specified link.
############################################################################
point.unilateral.straightness.link <- function(xa, ya, xb, yb, xc, yc, dgca, dgcb, ddab)
{	# process break-even point
	ellp <- (dgcb - dgca + ddab) / 2
#	cat("..ellp=",ellp,"\n",sep="")
	
	# common
	a <- xa - xc
	b <- (xb - xa) / ddab
	c <- ya - yc
	d <- (yb - ya) / ddab
##if(pp!=0 && pp!=1)	# to generate the article figure representing all str break-even points
#{	xpp <- xa + pp*b
#	ypp <- ya + pp*d
#	na <- which(V(g2)$x==xa & V(g2)$y==ya)[1]
#	nb <- which(V(g2)$x==xb & V(g2)$y==yb)[1]
#	g2 <<- delete.edges(g2,E(g2)[na %--% nb])
#	g2 <<- add.vertices(g2,1)
#	npp <- vcount(g2)
#	print(na);print(nb);print(npp)
#	V(g2)$type[npp] <<- "extra"
#	V(g2)$x[npp] <<- xpp
#	V(g2)$y[npp] <<- ypp
#	g2 <<- add.edges(g2,c(na,npp,npp,nb))
#}

	# from 0 to p'
	if(abs(ellp)>1e-10)
	{	e <- dgca
		f <- 1
		part1 <- point.unilateral.straightness.antiderivative(a,b,c,d,e,f,ellp) - point.unilateral.straightness.antiderivative(a,b,c,d,e,f,0)
	}
	else
		part1 <- 0
	
	# from p' to 1
	if(abs(ellp-ddab)>1e-10)
	{	e <- dgcb + ddab
		f <- -1
		part2 <- point.unilateral.straightness.antiderivative(a,b,c,d,e,f,ddab) - point.unilateral.straightness.antiderivative(a,b,c,d,e,f,ellp)
	}
	else
		part2 <- 0
	
#	cat("..part1=",part1," part2=",part2,"\n",sep="")
	result <- part1 + part2
	return(result)
}


############################################################################
# Processes the average point unilateral straightness, i.e. the average straightness 
# between the center of the network and any point constituting the network, 
# including those which are not nodes. In other words, we consider paths stopping on 
# a link.
#
# g: the graph to process.
# center: id of the node considered as the center.
# 
# returns: a numerical value (the mean straightness).
############################################################################
point.unilateral.straightness <- function(g, center=1)
{	# init
	el <- get.edgelist(g,names=FALSE)
	g <- distances.as.weights(g)
	sp <- shortest.paths(graph=g, weights=E(g)$dist)
	xc <- V(g)$x[center]
	yc <- V(g)$y[center]
	
	# process point straightness for all links in the network
	total.straightness <- 0
	total.length <- 0
	for(e in 1:nrow(el))
	{	# get the required data
		na <- el[e,1]			# get node A
		nb <- el[e,2]			# get node B
		xa <- V(g)$x[na]		# get their spatial position
		ya <- V(g)$y[na]
		xb <- V(g)$x[nb]
		yb <- V(g)$y[nb]
		dgca <- sp[center,na] 	# geodesic distance between the center and A
		dgcb <- sp[center,nb]	# same between the center and B
		#TODO should check if these geodesic distances are defined (disconnected graph?)
		ddab <- E(g)$dist[e]	# Euclidean distance between A and B
		
#		cat("Link (",na,",",nb,"):\n",sep="")
		
		# center must be processed separately (otherwise: undefined primitive)
		if(na==center || nb==center)
		{	# in this case (center and one of its neighbors), the straightness
			# is by definition one over the whole link, so integrating it is trivial
			pus <- ddab
		}
		# general case
		else
		{	pus <- point.unilateral.straightness.link(xa,ya,xb,yb,xc,yc,dgca,dgcb,ddab)
		}
		
		# update both variables
#		cat("..",pus," vs. ",ddab,"\n",sep="")
		total.straightness <- total.straightness + pus
		total.length <- total.length + ddab
	}
	
	# divide to get average 
	result <- total.straightness / total.length
#	cat("Total: ",total.straightness,"/",total.length," = ",result,"\n",sep="")
	
	return(result)
}


############################################################################
# Auxilary function used during the processing of the exact straightness.
# 
# a,b,c,d,e,f,g,h,i: constants used as parameters of the integrated function.
# fell2a,fell2b: functions used to process ell2.
# lb,ub: lower and upper bound of the integral to process.
#
# returns: value of the definite integral needed to process the exact straightness.
############################################################################
point.bilateral.straightness.integral <- function(a,b, c, d,e, f, g,h, i, fell2a, fell2b, lb, ub)
{	#cat("......lb: ",sprintf("%.23f",lb)," ub: ",sprintf("%.23f",ub),"\n",sep="")
	
	if(abs(lb-ub)<1e-10)
		result <- 0
	else
	{	fun <- Vectorize(function(x)
				{	f2lb <- fell2a(x)
					f2ub <- fell2b(x)
#					cat("........for x=",x,", from ",f2lb," to ",f2ub,"\n",sep="")
					if(abs(f2ub-f2lb)<1e-10)
						res <- 0
					else
					{	part1 <- point.unilateral.straightness.antiderivative(a+b*x,c,d+e*x,f,g+h*x,i,f2ub)
						part2 <- point.unilateral.straightness.antiderivative(a+b*x,c,d+e*x,f,g+h*x,i,f2lb)
						res <- part1 - part2
#						cat("..........part1: ",part1," part2: ",part2,"\n",sep="")
					}
					return(res)
				})
		
		result <- integrate(f=fun,lower=lb,upper=ub)$value
		
		fun2 <- Vectorize(function(x) fell2b(x)-fell2a(x))
		maxval <- integrate(f=fun2,lower=lb,upper=ub)$value
#		cat("........",result," vs. ",maxval,"\n",sep="")
	}
	
	return(result)
}
point.bilateral.straightness.integral.bis <- function(a,b, c, d,e, f, g1,h1, i1, g2,h2, i2, fell2, ell2up, lb, ub)
{	#cat("......lb: ",sprintf("%.23f",lb)," ub: ",sprintf("%.23f",ub),"\n",sep="")
	
	if(abs(lb-ub)<1e-10)
		result <- 0
	else
	{	fun <- Vectorize(function(x)
				{	ell2mid <- fell2(x)
#					cat("........for x=",x,", from 0 through ",ell2mid," to ",ell2up,"\n",sep="")
					if(ell2mid<1e-10)
						part1 <- 0
					else
					{	part11 <- point.unilateral.straightness.antiderivative(a+b*x,c,d+e*x,f,g1+h1*x,i1,ell2mid)
						part12 <- point.unilateral.straightness.antiderivative(a+b*x,c,d+e*x,f,g1+h1*x,i1,0)
						part1 <- part11 - part12
#						cat("..........part1 = ",part11," - ",part12," = ",part1,"\n",sep="")
#						cat(".......... max1 = ",ell2mid," >> ",part1<=ell2mid,"\n",sep="")
						
#						fun2 <- function(y) sqrt((a+b*x+c*y)^2 + (d+e*x+f*y)^2)/(g1+h1*x+i1*y)
#						verif <- integrate(f=fun2,lower=0,upper=ell2mid)$value
#						cat("..........verif = ",verif," >> ",abs(part1-verif)<1e-10,"\n",sep="")
#						curve(fun2,0,ell2mid,ylim=c(0,4));stop()
					}
					if(abs(ell2up-ell2mid)<1e-10)
						part2 <- 0
					else
					{	part21 <- point.unilateral.straightness.antiderivative(a+b*x,c,d+e*x,f,g2+h2*x,i2,ell2up)
						part22 <- point.unilateral.straightness.antiderivative(a+b*x,c,d+e*x,f,g2+h2*x,i2,ell2mid)
						part2 <- part21 - part22
#						cat("..........part2 = ",part21," - ",part22," = ",part2,"\n",sep="")
#						cat(".......... max2 = ",ell2up-ell2mid," >> ",part2<=(ell2up-ell2mid),"\n",sep="")
						
#						fun21 <- function(y) sqrt((a+b*x+c*y)^2 + (d+e*x+f*y)^2)
#						fun22 <- function(y) g2+h2*x+i2*y
#						fun2 <- function(y) fun21(y)/fun22(y)
#						verif <- integrate(f=fun2,lower=ell2mid,upper=ell2up)$value
#						cat("..........verif = ",verif," >> ",abs(part2-verif)<1e-10,"\n",sep="")
						
#						curve(fun2,ell2mid,ell2up,ylim=c(0,4),col="RED")
#						curve(fun21,ell2mid,ell2up,add=TRUE,col="GREEN")
#						curve(fun22,ell2mid,ell2up,add=TRUE,col="BLUE")
#						stop()
					}
					res <- part1 + part2
#					cat("..........res = ",part1," + ",part2," = ",res,"\n",sep="")
					return(res)
				})
		
		intres <- integrate(f=fun,lower=lb,upper=ub,abs.tol=1e-15)
		result <- intres$value
		
#		cat("........",result," vs. ",ell2up*(ub-lb)," (",intres$abs.error,") \n",sep="")
	}
	
	return(result)
}


############################################################################
# Auxilary function used during the processing of the exact straightness.
#
# xa1,ya1: coordinates of A1, the first end point of the first considered link.
# xb1,yb1: coordinates of B1, the second end point of the first considered link.
# xa2,ya2: coordinates of A2, the first end point of the second considered link.
# xb2,yb2: coordinates of B2, the second end point of the second considered link.
# dga1a2: geodesic distance between A1 and A2.
# dga1b2: geodesic distance between A1 and B2.
# dgb1a2: geodesic distance between B1 and A2.
# dgb1b2: geodesic distance between B1 and B2.
# dda1b1: Euclidean distance between A1 and B1.
# dda2b2: Euclidean distance between A2 and B2.
# 
# returns: total unilateral straightness over the specified link.
############################################################################
point.bilateral.straightness.link <- function(xa1, ya1, xb1, yb1, xa2, ya2, xb2, yb2, dga1a2, dga1b2, dgb1a2, dgb1b2, dda1b1, dda2b2)
{	TOLERANCE <- 1e-10
	
	# process break-even points
	ell1p <- (dda1b1 + dgb1a2 - dga1a2) / 2
	if(abs(ell1p)<TOLERANCE)
		ell1p <- 0
	else if(abs(dda1b1-ell1p)<TOLERANCE)
		ell1p <- dda1b1
	ell1pp <- (dda1b1 + dgb1b2 - dga1b2) / 2
	if(abs(ell1pp)<TOLERANCE)
		ell1pp <- 0
	else if(abs(dda1b1-ell1pp)<TOLERANCE)
		ell1pp <- dda1b1
	fell2p1 <- function(ell1)
		{	res <- (dga1b2 - dga1a2 + dda2b2) / 2
			if(abs(res)<TOLERANCE)
				res <- 0
			else if(abs(dda2b2-res)<TOLERANCE)
				res <- dda2b2
			return(res)
		}
	fell2p2 <- function(ell1)
		{	res <- (dgb1b2 - dga1a2 + dda1b1 + dda2b2 - 2*ell1) / 2
			if(abs(res)<TOLERANCE)
				res <- 0
			else if(abs(dda2b2-res)<TOLERANCE)
				res <- dda2b2
			return(res)
		}
	fell2p3 <- function(ell1)
		{	res <- (dga1b2 - dgb1a2 - dda1b1 + dda2b2 + 2*ell1) / 2
			if(abs(res)<TOLERANCE)
				res <- 0
			else if(abs(dda2b2-res)<TOLERANCE)
				res <- dda2b2
			return(res)
		}
	fell2p4 <- function(ell1)
		{	res <- (dgb1b2 - dgb1a2 + dda2b2) / 2
			if(abs(res)<TOLERANCE)
				res <- 0
			else if(abs(dda2b2-res)<TOLERANCE)
				res <- dda2b2
			return(res)
		}
	
	# debug
#	cat("(",xa1,",",ya1,")-[",dda1b1,"]-(",xb1,",",yb1,")          (",xa2,",",ya2,")-[",dda2b2,"]-(",xb2,",",yb2,")\n",sep="")
#	cat("dga1a2: ",dga1a2," dga1b2: ",dga1b2," dgb1a2: ",dgb1a2," dgb1b2: ",dgb1b2,"\n",sep="")
#	cat("ell1p: ",ell1p," ell1pp: ",ell1pp,"\n",sep="")
	
	# common constants
	a <- xa2 - xa1
	b <- (xa1 - xb1) / dda1b1
	c <- (xb2 - xa2) / dda2b2
	d <- ya2 - ya1
	e <- (ya1 - yb1) / dda1b1
	f <- (yb2 - ya2) / dda2b2
	
	# specific cases
	ga1a2 <- dga1a2						; ha1a2 <- 1	; ia1a2 <- 1
	ga1b2 <- dga1b2 + dda2b2			; ha1b2 <- 1	; ia1b2 <- -1
	gb1a2 <- dgb1a2 + dda1b1			; hb1a2 <- -1	; ib1a2 <- 1
	gb1b2 <- dgb1b2 + dda1b1 + dda2b2	; hb1b2 <- -1	; ib1b2 <- -1
	
	if(ell1p <= ell1pp)
	{	#cat(".from 0 through ell1p=",ell1p," through ell1pp=",ell1pp," to ",dda1b1,"\n",sep="")
		# from 0 to l1'
#		part1 <- (point.bilateral.straightness.integral(a,b, c, d,e, f, g=ga1a2,h=ha1a2, i=ia1a2, fell2a=function(x)0, fell2b=fell2p1, lb=0, ub=ell1p)
#				+ point.bilateral.straightness.integral(a,b, c, d,e, f, g=ga1b2,h=ha1b2, i=ia1b2, fell2a=fell2p1, fell2b=function(x)dda2b2, lb=0, ub=ell1p))
		part1 <- point.bilateral.straightness.integral.bis(a,b, c, d,e, f, g1=ga1a2,h1=ha1a2, i1=ia1a2, g2=ga1b2,h2=ha1b2, i2=ia1b2, fell2=fell2p1, ell2up=dda2b2, lb=0, ub=ell1p)
		max1 <- ell1p*dda2b2
		
		# from l1' to l1''
#		part2 <- (point.bilateral.straightness.integral(a,b, c, d,e, f, g=gb1a2,h=hb1a2, i=ib1a2, fell2a=function(x)0, fell2b=fell2p3, lb=ell1p, ub=ell1pp)
#				+ point.bilateral.straightness.integral(a,b, c, d,e, f, g=ga1b2,h=ha1b2, i=ia1b2, fell2a=fell2p3, fell2b=function(x)dda2b2, lb=ell1p, ub=ell1pp))
		part2 <- point.bilateral.straightness.integral.bis(a,b, c, d,e, f, g1=gb1a2,h1=hb1a2, i1=ib1a2, g2=ga1b2,h2=ha1b2, i2=ia1b2, fell2=fell2p3, ell2up=dda2b2, lb=ell1p, ub=ell1pp)
		max2 <- (ell1pp - ell1p)*dda2b2
		
		# from l1'' to A1B1
#		part3 <- (point.bilateral.straightness.integral(a,b, c, d,e, f, g=gb1a2,h=hb1a2, i=ib1a2, fell2a=function(x)0, fell2b=fell2p4, lb=ell1pp, ub=dda1b1)
#				+ point.bilateral.straightness.integral(a,b, c, d,e, f, g=gb1b2,h=hb1b2, i=ib1b2, fell2a=fell2p4, fell2b=function(x)dda2b2, lb=ell1pp, ub=dda1b1))
		part3 <- point.bilateral.straightness.integral.bis(a,b, c, d,e, f, g1=gb1a2,h1=hb1a2, i1=ib1a2, g2=gb1b2,h2=hb1b2, i2=ib1b2, fell2=fell2p4, ell2up=dda2b2, lb=ell1pp, ub=dda1b1)
		max3 <- (dda1b1 - ell1pp)*dda2b2
	}
	else
	{	#cat(".from 0 through ell1pp=",ell1pp," through ell1p=",ell1p," to ",dda1b1,"\n",sep="")
		# from 0 to l1''
#		part1 <- (point.bilateral.straightness.integral(a,b, c, d,e, f, g=ga1a2,h=ha1a2, i=ia1a2, fell2a=function(x)0, fell2b=fell2p1, lb=0, ub=ell1pp)
#				+ point.bilateral.straightness.integral(a,b, c, d,e, f, g=ga1b2,h=ha1b2, i=ia1b2, fell2a=fell2p1, fell2b=function(x)dda2b2, lb=0, ub=ell1pp))
		part1 <- point.bilateral.straightness.integral.bis(a,b, c, d,e, f, g1=ga1a2,h1=ha1a2, i1=ia1a2, g2=ga1b2,h2=ha1b2, i2=ia1b2, fell2=fell2p1, ell2up=dda2b2, lb=0, ub=ell1pp)
		max1 <- ell1pp*dda2b2
		
		# from l1'' to l1'
#		part2 <- (point.bilateral.straightness.integral(a,b, c, d,e, f, g=ga1a2,h=ha1a2, i=ia1a2, fell2a=function(x)0, fell2b=fell2p2, lb=ell1pp, ub=ell1p)
#				+ point.bilateral.straightness.integral(a,b, c, d,e, f, g=gb1b2,h=hb1b2, i=ib1b2, fell2a=fell2p2, fell2b=function(x)dda2b2, lb=ell1pp, ub=ell1p))
		part2 <- point.bilateral.straightness.integral.bis(a,b, c, d,e, f, g1=ga1a2,h1=ha1a2, i1=ia1a2, g2=gb1b2,h2=hb1b2, i2=ib1b2, fell2=fell2p2, ell2up=dda2b2, lb=ell1pp, ub=ell1p)
		max2 <- (ell1p - ell1pp)*dda2b2
		
		# from l1' to A1B1
#		part3 <- (point.bilateral.straightness.integral(a,b, c, d,e, f, g=gb1a2,h=hb1a2, i=ib1a2, fell2a=function(x)0, fell2b=fell2p4, lb=ell1p, ub=dda1b1)
#				+ point.bilateral.straightness.integral(a,b, c, d,e, f, g=gb1b2,h=hb1b2, i=ib1b2, fell2a=fell2p4, fell2b=function(x)dda2b2, lb=ell1p, ub=dda1b1))
		part3 <- point.bilateral.straightness.integral.bis(a,b, c, d,e, f, g1=gb1a2,h1=hb1a2, i1=ib1a2, g2=gb1b2,h2=hb1b2, i2=ib1b2, fell2=fell2p4, ell2up=dda2b2, lb=ell1p, ub=dda1b1)
		max3 <- (dda1b1 - ell1p)*dda2b2
	}
	
#	cat("..part1=",part1," part2=",part2," part3=",part3,"\n",sep="")
#	cat(".. max1=",max1, "  max2=", max2,"  max3=", max3,"\n",sep="")
	result <- part1 + part2 + part3
	return(result)
}


############################################################################
# Processes the average point bilateral straightness, i.e. the average straightness 
# between any pair of points constituting the network (including, but not limited to
# nodes). In other words, we also consider paths stopping anywhere on a link.
#
# g: the graph to process.
# 
# returns: a numerical value (the mean straightness).
############################################################################
point.bilateral.straightness <- function(g)
{	# init
	el <- get.edgelist(g,names=FALSE)
	g <- distances.as.weights(g)
	sp <- shortest.paths(graph=g, weights=E(g)$dist)
	
	# process point straightness for all pairs of links in the network
	total.straightness <- 0
	total.length <- 0
	for(e1 in 1:nrow(el))
	{	# get the required data
		na1 <- el[e1,1]			# get node A1
		nb1 <- el[e1,2]			# get node B1
		xa1 <- V(g)$x[na1]		# get their spatial position
		ya1 <- V(g)$y[na1]
		xb1 <- V(g)$x[nb1]
		yb1 <- V(g)$y[nb1]
		dda1b1 <- E(g)$dist[e1]	# Euclidean distance between A1 and B1
		
		for(e2 in e1:nrow(el))
		{	# get the required data
			na2 <- el[e2,1]			# get node A2
			nb2 <- el[e2,2]			# get node B2
			xa2 <- V(g)$x[na2]		# get their spatial position
			ya2 <- V(g)$y[na2]
			xb2 <- V(g)$x[nb2]
			yb2 <- V(g)$y[nb2]
			dda2b2 <- E(g)$dist[e2]	# Euclidean distance between A2 and B2
			
			# debug
#			cat("Link (",na1,",",nb1,") vs. (",na2,",",nb2,"):\n",sep="")
			
			# geodesic distances
			dga1a2 <- sp[na1,na2] 	# geodesic distance between A1 and A2
			dga1b2 <- sp[na1,nb2]	# geodesic distance between A1 and B2
			dgb1a2 <- sp[nb1,na2] 	# geodesic distance between B1 and A2
			dgb1b2 <- sp[nb1,nb2]	# geodesic distance between B1 and B2
			
			# the case "one link over the same link" must be processed separately (otherwise: undefined primitive)
			if(e1==e2)
			{	# in this case, the straightness is by definition one 
				# for all pairs of points over both links, so integrating it is trivial
				# we divide by 2 because the graph is undirected, so we don't want to count paths twice
				pbs <- dda1b1 * dda2b2 / 2
				length <- dda1b1 * dda2b2 / 2
			}
			# general case
			else
			{	pbs <- point.bilateral.straightness.link(xa1, ya1, xb1, yb1, xa2, ya2, xb2, yb2, dga1a2, dga1b2, dgb1a2, dgb1b2, dda1b1, dda2b2)
				# test using python instead, to integrate
				#pbs <- python.call(
				#		py.foo="point_bilateral_straightness_link", 
				#		xa1, ya1, xb1, yb1, xa2, ya2, xb2, yb2, dga1a2, dga1b2, dgb1a2, dgb1b2, dda1b1, dda2b2
				#)
				length <- dda1b1*dda2b2
			}
			
			# update both variables
#			cat("..",pbs," vs. ",dda1b1*dda2b2,"\n",sep="")
#			if(pbs>dda1b1*dda2b2)
#				cat("..ERRRRRRRROR\n")
			total.straightness <- total.straightness + pbs
			total.length <- total.length + length
		}
	}
	
	# divide to get average 
	result <- total.straightness / total.length
#	cat("Total: ",total.straightness,"/",total.length," = ",result,"\n",sep="")
	
	return(result)
}
