# SpatialMeasures v.0.2
=======
*Continuous average Straightness for spatial graphs*

* Copyright 2016-17 Vincent Labatut 

SpatialMeasures is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. For source availability and license information see `licence.txt`

* Lab site: http://lia.univ-avignon.fr/
* GitHub repo: https://github.com/CompNet/SpatialMeasures
* Contact: Vincent Labatut <vincent.labatut@univ-avignon.fr>

-----------------------------------------------------------------------

# Description
This set of R scripts was designed to process several variants of the *Straightness* (aka. *Directness* and probably other names): the ratio of the Euclidean to the graph distance. It is a measure designed to study spatial graphs, i.e. graphs embedded in an Euclidean space (nodes have spatial positions, links have spatial length, etc.).

First, this toolbox can process the Straightness using the traditional approach, i.e. considering only paths connecting two nodes. It can process the Straightness between two specific nodes, or the Straightness averaged over certain pairs of nodes in the graph (including all pairs).

Second, this toolbox can also compute the average Straightness through a *continuous* approach (by opposition to the *discrete* traditional approach), and incidentally this is the point of the article [Lbt'17]. The Straightness is generalized to deal with point-to-point paths, i.e. paths connecting a location on a link to another such position (by opposition to nodes, which are necessarily located on the ends of links). Based on this, we can average the Straightness by integration over the link. This allows deriving the following continuous average variants:
* Average Straightness between a point and an link (or rather: all the points on this link);
* Average Straightness between a point and the rest of the graph (or rather: all the points constituting this graph);
* Average Straightness between two links (or rather: all pairs of points located on these links, each point being on a different one);
* Average Straightness between a link and the rest of the graph (or rather: all the points located on the link on one side, and all the points constituting the graph on the other side);
* Average Straightness over the graph (or rather: between all the pairs of points constituting the graph).

Besides the functions used to process the average measures themselves, the scripts also allow to replicate the different experiments conducted in the article [Lbt'17].


# Data
Most of the data is generated randomly and recorded in the `data` folder (see the *Organization* section for more details).

We also experiment on a few real-world networks. These are all road networks retrieved from [OpenStreetMap](https://www.openstreetmap.org), using the script `src/misc/extraction.R`. These networks are available on [Figshare](https://doi.org/10.6084/m9.figshare.4721407), as explained in the *Installation* section.


# Organization
Here are the folders composing the project:
* Folder `src`: contains the source code (R scripts).
  * Folder `evaluation`: scripts used to evaluate the time and memory usage (see [Lbt'17]).
  * Folder `figures`: scripts used to produce some of the figures presented in the article, and process certain of their properties.
  * Folder `misc`: scripts used for extracting, generating, modifying and plotting the graphs.
  * Folder `straightness`: scripts used to process the discrete and continuous average variants of the average Straightness.
  * Script `main.R`: shows an example of how to process the Straightness on your own graphs through this toolkit.  
* Folder `data`: contains the files used and generated by our scripts.
  * Folder `figures`: files produced by the scripts in `src/figures`, corresponding to graph representations appearing in the paper.  
  * Folder `eval`: files related to the evaluation of the time processing and memory usage needed when processing the average Straightness.
  * Folder `urban`: plots and data generated when processing the Straightness of real-world road networks.


# Installation
1. Install the [`R` language](https://www.r-project.org/)
2. Install the following R packages:
   * [`geometry`](https://cran.r-project.org/web/packages/geometry/index.html): required (tested with version 0.3.6).
   * [`igraph`](http://igraph.org/r/): required (tested with version 1.0.1).
   * [`plotrix`](https://cran.r-project.org/web/packages/plotrix/): required (tested with version 3.5-11).
   * [`osmar`](https://cran.r-project.org/web/packages/osmar/index.html): only if you want to extract new road networks via the `src/misc/extraction.R script` (tested with version 1.1-7)
   * [`splancs`](https://cran.r-project.org/web/packages/splancs/index.html): only if you want to generate spiderweb-like networks via the `src/figures/regular.R script` (tested with version 2.01-39)
3. Download this project from GitHub and unzip the archive.
4. If you want to apply the scripts on the real-world road networks from the article [Lbt'17], download the `graph` archive from [Figshare](https://doi.org/10.6084/m9.figshare.4721407) and unzip it so that it overwrites the `data/urban` folder. 


# Use
In order to replicate the experiments from the article, perform the following operations:

1. Open the `R` console.
2. Set the current projetct directory as the working directory, using `setwd("my/path/to/the/project/SpatialMeasures")`.
3. Run one of the following scripts:
   * `src/figures/regular.R`, `src/figures/random.R` or `src/figures/urban.R` to generate the graph figures from the article [Lbt'17].
   * `src/evaluation/time.R`, `src/evaluation/memory.R` or `src/evaluation/urban.R` to produce the plots/tables related to performance.
   * `src/evaluation/urban.R`

If you just want to apply the measures to your own graph, then you need to use the following functions, as shown in `src/main.R`:
* `mean.straightness.nodes.link(g, u, e)`: average Straightness between node `u` and link `e`.
* `mean.straightness.nodes.graph(g, u)`: average Straightness between node `u` and the rest of graph `g`.
* `mean.straightness.link.link(g, e1, e2)`: average Straightness between links `e1` and `e2`.
* `mean.straightness.links.graph(g, e)`: average Straightness between link `e` and the rest of graph `g`.
* `mean.straightness.graph(g)`: average Straightness between each pair of points constituting graph `g`.

Note the specified graph must be an `igraph` object, and its nodes must be described by their position in a 2D space taking.
For this purpose, the graph must contain two nodal attributes called `x` and `y`. See the `igraph` documentation to know how to define such attributes.

Also, note you only need the `R` script `src/straightness/continuous.R` to process the variants of the continuous average straightness on your graphs. The rest of the scripts is just there for plotting, testing, evaluating, comparing, etc. You need to uncomment the function `get.dist` in `continuous.R`, though.
  

# Extension
You may want to apply the approach described in our paper to average other spatial measures. You would need either to process analytically the integral on the considered measure, as described for the first integration of the Straightness in the article. Alternatively, you may use numerical integration, as described for the second integration of the Straightness in the article. In any case, you can use the functional decomposition from `continuous.R`, but you will have to adapt most of it to the processing of your specific measure (i.e. our source code is not generic). 


# Dependencies
* [`igraph`](http://igraph.org/r/) package: used to build and handle graphs.
* [`geometry`](https://cran.r-project.org/web/packages/geometry/index.html) package: used to triangulate and generate planar random graphs.
* [`plotrix`](https://cran.r-project.org/web/packages/plotrix/) package: used to generate certain plots.
* [`osmar`](https://cran.r-project.org/web/packages/osmar/index.html) package: used to extract real-world road networks.
* [`splancs`](https://cran.r-project.org/web/packages/splancs/index.html): used when generating spiderweb-like networks.


# To-do List
* Check if tolerance is now actually needed (after dealing with the specific cases).


# References
* **[Lbt'17]** Labatut, V. Continuous Average Straightness in Spatial Graphs, Submitted to [Journal of Complex Networks](https://academic.oup.com/comnet), 2017.
**<URL goes here>**
