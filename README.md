# SpatialMeasures v.0.1
=======
*Continuous average Straightness for spatial graphs*

* Copyright 2016-17 Vincent Labatut 

SpatialMeasures is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. For source availability and license information see `licence.txt`

* Lab site: http://lia.univ-avignon.fr/
* GitHub repo: https://github.com/CompNet/SpatialMeasures
* Contact: Vincent Labatut <vincent.labatut@univ-avignon.fr>

-----------------------------------------------------------------------

# Description
This set of R scripts was designed to process several variants of the *Straightness* (aka. *Directness* and probably others): the ratio of the Euclidean to the graph distance. It is a measure designed to study spatial graphs, i.e. graphs embedded in an Euclidean space (nodes have spatial positions, links have spatial length, etc.).

First, this toolbox can process the Straightness using the traditional approach, i.e. considering only paths connecting two nodes. It can process the Straightness between two specific nodes, or the Straightness averaged over all the pairs of nodes in the graph.

Second, this toolbox can also handle the Straightness through a *continuous* approach (by opposition to the *discrete* traditional approach), and incidentally this is the point of the article [Lbt'16]. The Straightness is generalized to deal with point-to-point paths, i.e. ones connecting a location on a link to another such position (by opposition to nodes, which are necessarily located on the ends of links). Based on this, we can process the Straightness by integration over the link. This allows deriving the following continuous average variants:
* Average Straightness between a point and an link (or rather: all the points on this link);
* Average Straightness between a point and the rest of the graph (or rather: all the points constituting this graph);
* Average Straightness between two links (or rather: all pairs of points located on these links, each point being on a different one);
* Average Straightness between a link and the rest of the graph (or rather: all the points located on the link on one side, and all the points constituting the graph on the other side);
* Average Straightness over the graph (or rather: between all the pairs of points constituting the graph).

Besides the functions used to process the average measures themselves, the scripts also allow to replicate the experiment conducted in my article [Lbt'16].


# Data
Most of the data is generated randomly and recorded in the `data` folder (see the *Organization* section for more details).

We also experiment on a few publicly available real-world networks, also contained in the data folder `data`.

**<To be completed>**


# Organization
Here are the folders composing the project:
* Folder `src`: contains the source code (R scripts).
  * Folder `evaluation`: scripts used to evaluate the time and memory usage (see [Lbt'16]).
  * Folder `figures`: scripts used to produce some of the figures presented in the article.
  * Folder `misc`: scripts used for generating, modifying and plotting the graphs.
  * Folder  `straightness`: scripts used to process the discrete and continuous average variants of the Straightness.
  * Script `main.R`: shows an example of how to use process the Straightness on your own graphs.  
* Folder `data`: contains the files used by our scripts, and generated by our scripts.
  * Folder `figures`: files produced by the scripts in `src/figures`, corresponding to graph representations appearing in the paper.  
  * Folder `time`: files related to the temporal evaluation of Straightness processing (these are plots and graphs).
  * Folder `memory`: same thing for the memory aspect of performance.

**<To be edited>**


# Installation
1. Install the [`R` language](https://www.r-project.org/)
2. Install the following R packages:
   * [`geometry`](https://cran.r-project.org/web/packages/geometry/index.html) (tested with version 0.3.6)
   * [`igraph`](http://igraph.org/r/) (tested with version 1.0.1)
3. Download this project from GitHub and unzip the archive.


# Use
In order to replicate the experiments from the article, perform the following operations:
1. Open the `R` console.
2. Set the current directory as the working directory, using `setwd("<my directory>")`.
3. Run one of the scripts `src/figures/exp2.R` to generate the article graph figures, `src/evaluation/time.R` or `src/evaluation/memory.R` to produce the plots related to performance.

If you just want to apply the measures to your own graph, then need to use the following functions, as shown in `src/main.R`:
* xxxxx
* xxxxx
**<To be completed>**

Note the specified graph must be an `igraph` object, and its nodes must be described by their position in a 2D space taking.
For this purpose, the graph must contain two nodal attributes called `x` and `y`. See the `igraph` documentation to know how to define such attributes.

Also, note you only need the `R` script `src/straightness/continuous.R` to process the variants of the continuous average straightness on your graphs. The rest of the scripts is just there for plotting, testing, evaluating comparing, etc.
  

# Extension
You may want to apply the approach described in our paper to average other spatial measures. You would need either to process analytically the integral on the considered measure, as described for the first integration of the Straightness in the article. Alternatively, you may use numerical integration, as described for the second integration of the Straightness in the article. In any case, you can use the functional decomposition from `continuous.R`, but you will have to adapt most of it to the processing of your specific measure (i.e. our source code is not generic). 


# Dependencies
* [`igraph`](http://igraph.org/r/) package: used to build and handle graphs.
* [`geometry`](https://cran.r-project.org/web/packages/geometry/index.html) package: used to triangulate and generate planar random graphs.


# To-do List
* <To be completed>


# References
* **[Lbt'16]** Labatut, V. Continuous Average Straightness in Spatial Graphs, Submitted to Journal of Complex Networls, 2016.
**<URL goes here>**
