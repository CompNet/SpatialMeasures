# SpatialMeasures v.0.1
=======
*Continuous average measures for spatial graphs*

* Copyright 2016 Vincent Labatut 

SpatialMeasures is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. For source availability and license information see `licence.txt`

* Lab site: http://lia.univ-avignon.fr/
* GitHub repo: https://github.com/CompNet/SpatialMeasures
* Contact: Vincent Labatut <vincent.labatut@univ-avignon.fr>

**Note:** this is a development version. Please, use the last available release instead.

-----------------------------------------------------------------------

# Description
This set of R scripts was designed to process continuous averages of several classic spatial measures (which are usually processed in a discrete way in the literature). 
For now, it can handle the following measures:
* Straightness (aka. Directness and probably others): the ratio of the Euclidean to the graph distance.
* Tortuosity (aka. Circuity, Directness, Route Factor, Detour Index, and probably others): the ratio of the graph to the Euclidean distance (i.e. the reciprocal of the Straightness).  

Besides the functions used to process the average measures themselves, the scripts also allow to replicate the experiment conducted in my article [L'16]. 


# Data
<To be completed>


# Organization
Here are the folders composing the project:
<To be edited>
* Folder `src`: contains the source code (R scripts).
* Folder `in`: contains the files used by our scripts, i.e. the inputs.
  * Folder `pils`: results of the correlation clustering method (or any other graph partitioning method), for each considered parameter set (year, policy, etc). 
  * Folder `raw`: the raw data extracted from the VoteWatch website.
    * Folder `aggregated`: this folder contains several CSV files build from the original data:
      * `all-votes.csv`: concatenation of all vote outcomes for all documents and all MEPS. Can be considered as a compact representation of the data contained in the folder `votes_by_document`.
      * `mep-details.csv`: list of the MEPs having voted at least once in the considered term, with their details.
      * `mep-loyalty.csv`: same thing than `allvotes.csv`, but for the loyalty (i.e. whether or not the MEP voted like the majority of the MEPs in his political group).
      * `policy-freq.csv`: list of the topics considered during the term, with the corresponding number of documents.
      * `vote-details.csv`: list of the voted texts with their details.
    * `original`: this folder contains a collection of CSV files, each one describing the outcome of the vote session relatively to one specific document.
* Folder `out`: contains the file produced by our scripts. See the *Use* section for more details.


# Installation
1. Install the [`R` language](https://www.r-project.org/)
2. Install the following R packages:
   * [`igraph`](http://igraph.org/r/) (tested with version 1.0.1)
   * [`geometry`](https://cran.r-project.org/web/packages/geometry/index.html) (tested with versions 0.3.5 and 0.3.6)
3. Download (and possibly unzip) this project from GitHub.


# Use
In order to replicate the experiments from the article, perform the following operations:

1. Open the `R` console.
2. Set the current directory as the working directory, using `setwd("<my directory>")`.
3. Run the main script `code/main.R`.

The script will produce the following files in the folder `output_files`:
* `agreement`: histograms representing the distributions of agreement and rebellion indices. Each subfolder corresponds to a specific topic.
* `community_algorithms_csv`: Performances obtained by the partitioning algorithms (for both community detection and correlation clustering). Each subfolder corresponds to a specific topic.
  * `xxxx_cluster_information.csv`: table containing several variants of the imbalance measure, for the considered algorithms.
* `community_algorithms_results`: Comparison of the partitions detected by the various algorithms considered, and distribution of the cluster/community sizes. Each subfolder corresponds to a specific topic.
  * `xxxx_cluster_comparison.csv`: table comparing the partitions detected by the community detection algorithms, in terms of Rand index and other measures.
  * `xxxx_ils_cluster_comparison.csv`: like `xxxx_cluster_comparison.csv`, except we compare the partition of community detection algorithms with that of the ILS.
  * `xxxx_yyyy_distribution.pdf`: histogram of the community (or cluster) sizes detected by algorithm `yyyy`.
* `graphs`: the networks extracted from the vote data. Each subfolder corresponds to a specific topic.
  * `xxxx_complete_graph.graphml`: network at the `Graphml` format, with all the information: nodes, edges, nodal attributes (including communities), weights, etc. 
  * `xxxx_edges_Gephi.csv`: only the links, with their weights (i.e. vote similarity). 
  * `xxxx_graph.g`: network at the `g` format (for ILS). 
  * `xxxx_net_measures.csv`: table containing some stats on the network (number of links, etc.).
  * `xxxx_nodes_Gephi.csv`: list of nodes (i.e. MEPs), with details.
  * `xxxx_infomap_community.txt`: membership vector generated by the application of the InfoMap algorithm.
  * `xxxx_multilevel_community.txt`: membership vector generated by the application of the Multi-Level algorithm.
  * `xxxx_fastgreedy_community.txt`: membership vector generated by the application of the FastGreedy algorithm.
  * `xxxx_walktrap_community.txt`: membership vector generated by the application of the Walktrap algorithm.

If you just want to use the measures, then call the following functions:

* xxxxx
* xxxxx
* <To be completed>

Note the specified graph must be an `igraph` object, and its nodes must be described by their position in a 2D space taking.
For this purpose, the graph must contain two nodal attributes called `x` and `y`. See the `igraph` documentation to know how to
define such attributes.  

# Extension
You may want to apply the approach described in our paper to average other spatial measures. This is possible
with this source code, provided you define... <instructions here>


# Dependencies
* [`igraph`](http://igraph.org/r/) package: used to build and handle graphs.
* [`geometry`](https://cran.r-project.org/web/packages/geometry/index.html) package: used to triangulate and generate planar random graphs.


# To-do List
* <To be completed>


# References
* **[L'16]** Labatut, V. <Title goes here>, Submitted to <Journal goes here>, 2016.
<URL goes here>
