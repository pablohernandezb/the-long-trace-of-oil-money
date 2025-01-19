#------------------------------------------------------------------------------
#  Title: main.R													  
#  Created by: Pablo Hernandez Borges 										  
#  Created on: 18 Apr 2021												  
#  Last modified on: 22 Apr 2021											  
#  Last modified by: Pablo Hernandez Borges								  
#  Purpose: This .R file is the main file to study corruption networks in
#           Venezuela using the dataset from the webpage Chavismo Inc.
#------------------------------------------------------------------------------



#==================================================================#
######===================    CONFIGURATION    ================######
#==================================================================#
print(R.version.string, side=1,line=4,adj=1)

# Check and set a working directory for your session
getwd()				# what is your current working directory? 

# Setting the working directory.
setwd("C:/Users/pabherna/OneDrive - Texas Tech University/Spring 2021/POLS 5367 - International Political Economy/Research Paper/data")

#==================================================================#
####========   Installing and Loading Packages in R   =========##### 
#==================================================================#

install.packages("statnet")
# statnet is a suite of packages and installs a number of packages belonging to the suite.
# the most important packages that come with it include "network" and "sna"

install.packages("igraph")
#igraph is another widely used network analysis tool with different functions than statnet

install.packages("rgl") # we'll use this in 3D visualization.

#==================================================================#
####== Setting up the working space and loading the libraries ==####
#==================================================================#

# Loading libraries
set.seed(47401)         # Set a seed for reproducible results
library(statnet)				
library(rgl)
search() 			          

#======================   Using help in R   =======================#

help(dir)								# screen on the function "dir"
?dir        						# another way, same thing

help.search('dir')     	# return topics with the term
??dir   								# another way, same thing. 

help.start()    				# browser window, full site search

#====================================================================#
####====== Loading an edgelist and vertex attributes into R ======####
#====================================================================#

### Relationships Legend ###

# International trials	1
# Business connection	2
# Company connection	3
# Complaints for corruption	4
# Designates in charge	5
# Enemies	6
# Facilitators	7
# Family 	8
# Friends	9
# Human Rights violation	10
# Integrates company	11
# Occupied functions	12
# Sanctioned by	13
# Student colleagues	14
# Suscribed contract	15


# Read in the data:
nodes <- read.csv("chavismoinc_nodes.csv", header=T, as.is=T)

# Changing nodes to 0
nodes$country_size[nodes$country_size == -1] <- 0

edges <- read.csv("chavismoinc_edges.csv", header=T, as.is=T)

# Examine the data:
head(nodes) 
head(edges) 

# Convert the data into the network format used by the statnet family.
# We can generate a 'network' object from an edgelist, 
# an adjacency matrix, or an incidence matrix. (same is true in igraph) 
?edgeset.constructors

#==================================================================#

# Both the "network" and the "igraph" packages can be used to convert
# these nodes and edges into network format.

#==================================================================#
####=============   Using the "network" package   ==============####
#==================================================================#

# Load the network package if not loaded already (should be loaded because it's part of statnet).
library(network)

# Remember to set the ignore.eval to F for weighted networks.
corrupt_net_1 <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                       loops=F, multiple=F, ignore.eval = F, directed = F)
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
corrupt_net_1

# You can access the edges, vertices, and the network matrix using:
corrupt_net_1[,]    # access the adjacency matrix

corrupt_net_1 %n% "net.name" <- "Corruption Network in Venezuela" #  set the network attribute "net.name" to "Media Network"
corrupt_net_1 %n% "net.name" # Check whether we've been able to set it properly.

corrupt_net_1 %v% "country"    # The node (v for vertex) attribute "country"
corrupt_net_1 %e% "type"     # The edge attribute "type"
corrupt_net_1 %v% "col_intl_trial"     # The edge attribute "type"

# these edge and node attributes came from the objects "edges" and "nodes" that
# we used when creating the network object. With the syntax above we can inspect
# any one we want. We can also set a new attribute or assign new value to an
# existing attribute using the `<-` (asssignment) operation.

# For example, we can set a node attribute called "col" that we will
# later use for visualizing node color. Right hand side chooses color
# based on the "media.type" node attribute.

# Plotting by Actors who has an International Trial.
corrupt_net_1 %v% "col_intl_trial" <- c("dodgerblue","firebrick1")[nodes$intl_trial+1]
plot(corrupt_net_1, vertex.cex=(corrupt_net_1 %v% "country_size")/5, vertex.col = (corrupt_net_1 %v% "col_intl_trial"))

# Now, plot the network:
plot(corrupt_net_1, vertex.cex=(corrupt_net_1 %v% "country_size")/5, vertex.col = (corrupt_net_1 %v% "col"))

# Plotting by Actors with a Cutpoint (whas International Trials Strong and Weak)

cutpoints <- rep.int(0, 1584)
for(i in 1:length(cutpoints)) {
  if(i<858) {
    if(cutpoints_strong[i]==FALSE && cutpoints_weak[i]==FALSE) {
      next
    } else if (cutpoints_strong[i]==FALSE && cutpoints_weak[i]==TRUE) {
      cutpoints[i] <- 1
    } else if (cutpoints_strong[i]==TRUE && cutpoints_weak[i]==TRUE){
      cutpoints[i] <- 3
    } else if (cutpoints_strong[i]==TRUE && cutpoints_weak[i]==FALSE){
      cutpoints[i] <- 2
    }
  } else if (i==858) {
    next
  } else {
    if(cutpoints_strong[i+1]==FALSE && cutpoints_weak[i+1]==FALSE) {
      next
    } else if (cutpoints_strong[i+1]==FALSE && cutpoints_weak[i+1]==TRUE) {
      cutpoints[i] <- 1
    } else if (cutpoints_strong[i+1]==TRUE && cutpoints_weak[i+1]==TRUE){
      cutpoints[i] <- 3
    } else if (cutpoints_strong[i+1]==TRUE && cutpoints_weak[i+1]==FALSE){
      cutpoints[i] <- 2
    }
  }
}

set.vertex.attribute(corrupt_net_1, "Cutpoints", cutpoints)
set.vertex.attribute(corrupt_net_1, "Intl_trial", Intl_Trial)

corrupt_net_1 %v% "col_cutpoints" <- c("dodgerblue","goldenrod2","orangered","firebrick1")[cutpoints+1]
plot(corrupt_net_1, vertex.cex=(((corrupt_net_1 %v% "Cutpoints")+1)/3), vertex.col = (corrupt_net_1 %v% "col_cutpoints"))

plot(corrupt_net_1, vertex.cex=(corrupt_net_1 %v% "country_size")/4, vertex.col = (corrupt_net_1 %v% "col_intl_trial"), mode = "kamadakawai")
plot(corrupt_net_1, vertex.cex=(corrupt_net_1 %v% "country_size")/4, vertex.col = (corrupt_net_1 %v% "col_cutpoints"), mode = "kamadakawai")

#==================================================================#
####================   Basic Network Measures   ================####
#==================================================================#

summary(corrupt_net_d)                                   # Get an overall summary
print(corrupt_net_d)                                     # Simple print method
network.dyadcount(corrupt_net_d)                         # How many dyads in nflo?
network.edgecount(corrupt_net_d)                         # How many edges are present?
network.size(corrupt_net_d)                              # How large is the network?
gden(corrupt_net_d, mode = "dgraph")                      # Network Density

#Degree Distribution
#Calculating In-Degree and Out-Degree to Visualize the Total Degree Distribution: What is the distribution of Connectiveness?
InDegree <- degree(corrupt_net_d, cmode="indegree")     #Computing the in-degree of each node
OutDegree <- degree(corrupt_net_d, cmode="outdegree")   #Computing the out-degree of each node

par(mar = rep(2, 4))
par(mfrow=c(2,2)) # Set up a 2x2 display
hist(InDegree, xlab="Indegree", main="In-Degree Distribution", prob=FALSE)
hist(OutDegree, xlab="Outdegree", main="Out-Degree Distribution", prob=FALSE)
hist(InDegree+OutDegree, xlab="Total Degree", main="Total Degree Distribution", prob=FALSE)
par(mfrow=c(1,1)) # Restore display

#Average Path Length 

#Walks: A walk is a sequence of nodes and ties, starting and ending with nodes, in which each node is incident with the edges
#...following and preceding it in the sequence (Wasserman and Faust 1994, p. 105).
# The beginning and ending node of a walk may be differeent, some nodes may be included more than once, and some ties may be included more than once.
#Paths: A path is a walk where all the nodes and all the ties are distinct.
#A shortest path between two nodes is refrred to as a geodesic (Wasserman and Faust 1994, p. 110)
#Average path length or the geodesic distance is the average number of steps along the shortest paths for all possible pairs of nodes.

# By default, nodes that cannot reach each other have a geodesic distance of infinity. 
# Because, Inf is the constant for infinity, we need to replace INF values to calculate the shortest path length.
# Here we replace infinity values with 0 for visualization purposes.

corrupt_Geo <- geodist(corrupt_net_1, inf.replace=0)
#AHS_Geo <- geodist(AHS_Network)                #Matrix with Infinity
corrupt_Geo

#The length of the shortest path for all pairs of nodes.
corrupt_Geo$gdist 

#The number of shortest path for all pairs of nodes.
corrupt_Geo$counts  

#Shortest Path Matrix
Geo_Dist = corrupt_Geo$gdist
hist(Geo_Dist)

###########################
#   MESO-LEVEL MEASURES   #
###########################

#Dyads
#Null-Dyads: Pairs of nodes with no arcs between them
#Asymmetric dyads: Pairs of nodes that have an arc between the two nodes going in one direction or the other, but not both
#Mutual/Symmetric Dyad: Pairs of nodes that have arcs going to and from both nodes  <--> 

#Number of Symmetric Dyads
mutuality(corrupt_net_1)

#Dyadic Ratio: Ratio of Dyads where (i,j)==(j,i) to all Dyads
grecip(corrupt_net_d, measure="dyadic")

#Edgwise Ratio: Ratio of Reciprocated Edges to All Edges
grecip(corrupt_net_d, measure="edgewise")

as.sociomatrix(corrupt_net_1)                            # Show it as a sociomatrix (mess!)
corrupt_net_1[,]                                         # Another way to make a mess
corrupt_net_1[1:10,1:10]                                 # Display the first 10 rows and 10 cols

#==================================================================#
#================   Basic transformation options   ================#
#==================================================================#

corruptmax <- symmetrize(corrupt_net_1, rule="weak")	  # symmetrize this network based on "maximum"
                                                        # i.e. i<->j iff i->j or i<-j (OR rule)
corrupt_net_2 <- network(corruptmax, directed=FALSE)		# show it as a network object
corruptmax											                        # compare to original

corruptmin <- symmetrize(corrupt_net_1, rule="strong")	# symmetrize this network based on "minimum"
                                            # i.e. i<->j iff i->j and i<-j (AND rule)
network(corruptmin, directed=FALSE)					  # show it as a network object
corruptmin					# notice how many fewer edges there are in drugmin than in original network
?symmetrize										# consult the help file!


#==================================================================#
#=================   Basic visualization options   ================#
#==================================================================#
plot(corrupt_net_2,displaylabels=T)                  # basic plot with labels, no options
                                    # It might look like a mess, zoom in to have
                                    # a better look.
plot(corrupt_net_d,displaylabels=F,mode="circle")    # The old circle graph
gplot(corrupt_net_d)                                 # Requires sna package (in statnet)
# Sometimes when using gplot, you will get an error message about setting label length
# You can safely ignore this. 
# but you have many other choices - make good use of them.
# Run the following four plots together.
# Different layout options are specified with the "mode" argument.

# start here
par(mfrow=c(2,2)) 	#start a 2 by 2 paneled image that will display all 4 plots below
plot.network(corrupt_net_2, displaylabels=TRUE, displayisolates=TRUE, arrowhead.cex=.5, 
             label.cex=.5, vertex.cex=1, edge.col=8, label.col=(1), vertex.col=2, 
             label.border=1, vertex.border=1, sub="Default Layout, isolates and labels included")
plot.network(corrupt_net_2, displaylabels=FALSE, displayisolates=FALSE, arrowhead.cex=.5, 
             label.cex=.5, vertex.cex=1, edge.col=8, label.col=1, vertex.col=4, 
             label.border=1, vertex.border=1, mode="fruchtermanreingold", 
             sub="F-R layout, no isolate or labels")
plot.network(corrupt_net_2, displaylabels=FALSE, displayisolates=TRUE, arrowhead.cex=.5, 
             label.cex=.5, vertex.cex=1, edge.col=8, label.col=1, vertex.col=(nodes$country), 
             label.border=1, vertex.border=1, mode="kamadakawai", 
             sub="K-K layout, isolates included, no labels, color gender")
plot.network(corrupt_net_2, displaylabels=FALSE, displayisolates=FALSE, arrowhead.cex=.5, 
             label.cex=.5, vertex.cex=1, edge.col=8, label.col=1, vertex.col =(nodes$country_size), 
             label.border=1, vertex.border=1, mode="kamadakawai", 
             sub="K-K layout, no isolates or labels, color ethnicity")
par(mfrow=c(1,1))                            # Restore display
title(main="Comparing Plot Options", cex=3)
# end here

# learn more about your options!    
?plot.network
?gplot.layout
# your color options, in order
palette()

# Now, for some fun
# Note - we've had a little problem with this in some versions of Mac OS. 
gplot3d(corrupt_net_2, gmode = "digraph", mode = "fruchtermanreingold", displayisolates = FALSE, 
        vertex.radius = 5, vertex.col =(nodes$country_size))
# Select the image with your mouse and drag it around to get the 3-D effect.

# For more information about gplot and gplotd3:
?gplot
?gplot3d

### another option for drawing 

# Using the Florence data that we loaded previously.
plot.sociomatrix(flo, diaglab=F, cex.lab=.4)
?plot.sociomatrix

#==================================================================#
#=====================   BASIC MEASURES   =========================#
#==================================================================#

##### Density #####
gden(corrupt_net_2)								# density
?gden                     # for more information

##### Reciprocity #####
grecip(corrupt_net_d, measure = c("dyadic"))		        #All pairs, including null
grecip(corrupt_net_d, measure = c("dyadic.nonnull"))   #Like UCINET Dyad-based recip
grecip(corrupt_net_d, measure = c("edgewise"))         #Like UCINET Arc-based"
grecip(corrupt_net_d,"edgewise.lrr")		# log of the ratio of the edgewise reciprocity to density

# What does this mean? Want to know more? Type ?grecip for more info.
dyad.census(corrupt_net_d)                 # counts of mutuals, asymmetrics, and nulls
dyad.census(corrupt_net_d)                  # only two options in a symmetric graph

##### Triad Census #####
triad.census(drug)                 # Directed triad census
directedtriads <- triad.census(drug)				# we can assign it to an object if we like.
directedtriads			# you've created a new R object called directedtriads
ls()								# Now it's in out workspace.
# Now, we will write this output into a simple comma separated values file
# that we could easily use in another software package
write.csv(directedtriads,file="directedtriads.csv", row.names=FALSE)
?write.csv
?write.table  #other write options

triad.census(flo, mode="graph")		# Undirected triad census
undirtriads <- triad.census(flo, mode="graph")	
undirtriads						
ls()

##### Transitivity #####
gtrans(corrupt_net_d, mode="digraph", measure = c("weakcensus"), use.adjacency = FALSE)			
# returns the # of transitive ordered triples. 
gtrans(corrupt_net_d, mode="digraph", measure = c("weak"), use.adjacency = FALSE)						
#"weak" is the most commonly reported transitivity measurement

# statnet doesn't have an implemented measurement like UCINET's 2 leg/3 leg, 
# though you could program it if interested and you could calculate it using the triadic 
# census and the constituent transitive ordered triples!

# For information on why we should use "use.adjacency = FALSE" and other useful info:
?gtrans

##### Centrality #####
# Degree
degree(corrupt_net_2)                # A vector that counts all edges incident to a node 
                                     # a.k.a degree centrality
deg <- degree(corrupt_net_2)				 # An alternate way, now creating an R object
ls()								                 # see your new object? 

ideg <- degree(corrupt_net_d, cmode="indegree")                    # Indegree 
odeg <- degree(corrupt_net_d, cmode="outdegree")                   # Outdegree 
all(degree(corrupt_net_d) == ideg+odeg)                            
# This should be true, as the sum of all indegrees and outdegrees
# over the network equals the total degree for each node.

# Network level centralization measures
centralization(corrupt_net_d, degree, cmode="indegree") 			# in-degree centralization
centralization(corrupt_net_d, degree, cmode="outdegree") 		# out-degree centralization

# Now you can do things with these 
# If Indegree and Outdegree were perfectly correlated, we would be looking at 
# a perfectly linear relationship with a line with slope of 1 and an intercept of 0. 
plot(ideg,odeg, type="o", col="red") # Plot ideg by odeg
# youch - that's ugly

# Plot simple histograms of the degree distribution and tile them together:
par(mfrow=c(2,2))                                             # Set up a 2x2 display
hist(ideg, xlab="Indegree", main="Indegree Distribution", prob=TRUE)
hist(odeg, xlab="Outdegree", main="Outdegree Distribution", prob=TRUE)
hist(ideg+odeg, xlab="Total Degree", main="Total Degree Distribution", prob=TRUE)
par(mfrow=c(1,1))  											# Restore display

prest <- prestige(corrupt_net_1)		# prestige without a specified mode is just indegree centrality
all(prest == ideg)        # see that this is the case

plot(ideg,prest, type="o", col="red") # Plot ideg by prest
# Now you see a linear relationship - because prestige's default measure and 
# ideg measure the exact same thing! 

# Make a plot where the node size is proportional to in-degree centrality. 
gplot(corrupt_net_d, vertex.cex=(ideg)^0.4/2, vertex.sides=50, edge.col = 8,
      label.cex=0.4,arrowhead.cex=.5, label=network.vertex.names(corrupt_net_d),
      displayisolates=F, mode = "kamadakawai")
?gplot 				#so many options!

# Try it again, with some additional parameters. For example, vertex color is
# reflective of ethnicity. This time save your image as a pdf.
pdf("myplot1.pdf")
gplot(drug, vertex.cex=(ideg)^0.4/2, vertex.sides=50, edge.col = 8, vertex.col =(ethnicity),
      label.cex=0.4,arrowhead.cex=.5, label=network.vertex.names(drug), displayisolates=F, mode = "fruchtermanreingold")
dev.off()
# Go find your new pdf file in your working directory
# Comes in handy for placing an image into Latex. R allows for saving images
# in different formats. To see a number of them:
?png

# Betweenness centrality
bet <- betweenness(corrupt_net_d, gmode="digraph")  # Betweenness for a directed graph (digraph)
bet
centralization(corrupt_net_d, betweenness, mode="digraph")
#this is going to be a mess, because betweenness is unscaled and simply divided by 100
gplot(corrupt_net_1, vertex.cex=(bet)/100, gmode="digraph", vertex.col = (corrupt_net_1 %v% "col_intl_trial"))   
#so, try this, which takes the square root of the raw betweenness score and divides it by 20
gplot(corrupt_net_d, vertex.cex=sqrt(bet)/20, gmode="digraph", vertex.col = (corrupt_net_1 %v% "col_intl_trial"))   

#You can try to avoid this and get the scaled version 
bet_scaled <- betweenness(corrupt_net_1, gmode="digraph", rescale=T)     
# Scaled betweenness for a directed graph (digraph)
bet_scaled
# but you still have issues, so you have to multiply it this time to make it easier to see
gplot(corrupt_net_1, vertex.cex=(bet_scaled*50), gmode="digraph", vertex.col =(corrupt_net_1 %v% "col_intl_trial")) 
# remember, you have many options for rescaling a value. The important thing is that you 
# keep the rank. Adding a small constant value to all betweenness scores can help make the 
# smallest nodes visible while maintaining the rank.

# In the end, these decisions about how to scale node sizes, color, etc.
# come down to making sure a plot conveys the idea it is supposed to
# and is visually appealing.

#==================================================================#
#=================   SUBGROUP IDENTIFICATIONS   ===================#
#==================================================================#

clique.census(corrupt_net_1)  # returns a list
clique.census(corrupt_net_1,  # Find maximal cliques of varying sizes
              tabulate.by.vertex=FALSE,
              enumerate=FALSE) 
clique.census(corruptmax, # Remember symmetrizing on the maximum (under the "weak" rule)?
              tabulate.by.vertex=FALSE,
              enumerate=FALSE) 

### Now, we'll return to the florence marriage relations.
?flo
gplot(flo,displaylabels=T)
is.connected(flo)                                     # Strongly connected?
is.connected(flo, connected="weak")                   # Weakly connected?
# it can't be, because there is an isolate!
geodist(flo)                                          # geodesics
reachability(flo)                                     # reachbility matrix
# note everyone but 12 is reachable.
# That's why geodist shows it as of Inf(inite) distance to all other nodes.
kcores(corrupt_net_1)                                           # k-cores (by degree)
bicomponent.dist(corrupt_net_1)                                 # bicomponents
cutpoints(corrupt_net_1)                                        # find cutpoints
# see your cutpoints
gplot(corrupt_net_1,vertex.col=2+cutpoints(corrupt_net_1,mode="graph", return.indicator=T))

#==================================================================#
#===================   COMMUNITY DETECTION   ======================#
#==================================================================#

# We will now switch to igraph for some community detection.
detach("package:statnet")
library(igraph)

edgelist <- read.csv("chavismoinc_edges.csv", sep= "\t") # Something went wrong here
# Can you guess what?
# We can check the documentation to see the standard settings
?read.csv
# or check out the edgelist we've created:
head(edgelist)
# The first line was taken as column names. Instead:
edgelist <- read.csv("chavismoinc_edges.csv", header=TRUE)

# Now we are all set to create a graph in igraph
?graph_from_data_frame
corruptnet <- graph_from_data_frame(edgelist)

# We can now inspect the object
typeof(corruptnet)
class(corruptnet)
corruptnet
is.bipartite(corruptnet)
is.connected(corruptnet, mode = "strong")
is.connected(corruptnet, mode = "weak")

# The data also contains attributes that we can attach to the graph
attributes <- read.csv("hospitalattributes.csv") # Loads the data and saves it in the dataframe format
head(attributes) # Check if data was loaded properly

# Now we can attach the attributes to our graph
V(hospitalnet)$Sex=as.character(attributes$SEX[match(V(hospitalnet)$name,attributes$id)])
V(hospitalnet)$Age=as.character(attributes$AGE[match(V(hospitalnet)$name,attributes$id)])
V(hospitalnet)$Race=as.character(attributes$RACE[match(V(hospitalnet)$name,attributes$id)])
V(hospitalnet)$Title=as.character(attributes$TITLE[match(V(hospitalnet)$name,attributes$id)])
V(hospitalnet)$Unit=as.character(attributes$UNIT[match(V(hospitalnet)$name,attributes$id)])

# Show the value of the "Unit" attribute for every node
V(hospitalnet)$Unit

# Easy way to summarize the unique values assumed by this attribute
unique(V(hospitalnet)$Unit)

# Now we can assign colors to nodes based on their attributes, in this case their unit
V(hospitalnet)$color=V(hospitalnet)$Unit 
V(hospitalnet)$color=gsub("Critical Unit","red",V(hospitalnet)$color)
V(hospitalnet)$color=gsub("Boundary Spanners","blue",V(hospitalnet)$color)
V(hospitalnet)$color=gsub("Unit 1","green",V(hospitalnet)$color)
V(hospitalnet)$color=gsub("Unit 2","yellow",V(hospitalnet)$color)

# set arbitrary size for the nodes while preserving degree order
V(hospitalnet)$size=log(degree(hospitalnet))*3

#Let's compare some different plotting algorithms
pdf("my_plot2.pdf")
par(mfrow=c(2,2)) 	#start a 2 by 2 paneled image that will display all 4 plots below

plot.igraph(hospitalnet,vertex.label=NA,layout=layout.fruchterman.reingold,edge.arrow.size=0.05,edge.curved=TRUE)
legend(x=-1.5, y=-1.1, c("Critical Unit","Boundary Spanners", "Unit 1", "Unit 2"), pch=21, 
       pt.cex=2, cex=.8, bty="n", ncol=1, pt.bg=c("red", "blue", "green", "yellow"))

plot.igraph(hospitalnet,vertex.label=NA,layout=layout_in_circle,edge.arrow.size=0.05,edge.curved=TRUE)
legend(x=-1.5, y=-1.1, c("Critical Unit","Boundary Spanners", "Unit 1", "Unit 2"), pch=21, 
       pt.cex=2, cex=.8, bty="n", ncol=1, pt.bg=c("red", "blue", "green", "yellow"))

plot.igraph(hospitalnet,vertex.label=NA,layout=layout.mds, edge.arrow.size=0.05,edge.curved=TRUE)
legend(x=-1.5, y=-1.1, c("Critical Unit","Boundary Spanners", "Unit 1", "Unit 2"), pch=21, 
       pt.cex=2, cex=.8, bty="n", ncol=1, pt.bg=c("red", "blue", "green", "yellow"))

plot.igraph(hospitalnet,vertex.label=NA,layout=layout.davidson.harel, edge.arrow.size=0.05,edge.curved=TRUE)
legend(x=-1.5, y=-1.1, c("Critical Unit","Boundary Spanners", "Unit 1", "Unit 2"), pch=21, 
       pt.cex=2, cex=.8, bty="n", ncol=1, pt.bg=c("red", "blue", "green", "yellow"))
#turn the graphic device off
dev.off()

#==================================================================#

# We can also use community detection to see if we can find groups
# in the graph that are more densely connected than other groups 

# Unfortunately, community detection support is not great in igraph.
# If you are interested in better implementations you should use Python.

# transform the network into an undirected network
corruptnetcomm <- graph_from_data_frame(edgelist, directed = F)
# check if everything was properly loaded
corruptnetcomm
# See documentation for algorithm
?cluster_edge_betweenness
# Apply community detection
cluster_edge_betweenness(corruptnetcomm)
cluster_louvain(corruptnetcomm)

# save the results in a new object
communities1 <- cluster_edge_betweenness(corruptnetcomm)
communities2 <- cluster_louvain(corruptnetcomm)
communities1
communities2
typeof(communities2)
class(communities2)
# Plot the communities for visual comparison
par(mfrow=c(1,2))
plot(communities1,corruptnetcomm, vertex.label=NA)
plot(communities2,corruptnetcomm, vertex.label=NA)

par(mfrow=c(1,1))

###########################
#   MESO-LEVEL MEASURES   #
###########################

#Dyads
#Null-Dyads: Pairs of nodes with no arcs between them
#Asymmetric dyads: Pairs of nodes that have an arc between the two nodes going in one direction or the other, but not both
#Mutual/Symmetric Dyad: Pairs of nodes that have arcs going to and from both nodes  <--> 

#Number of Symmetric Dyads
mutuality(corrupt_network_d)

#Dyadic Ratio: Ratio of Dyads where (i,j)==(j,i) to all Dyads
grecip(corrupt_network_d, measure="dyadic")

#Edgwise Ratio: Ratio of Reciprocated Edges to All Edges
grecip(corrupt_network_d, measure="edgewise")

#Directed Triad Census
#Triads can be in Four States
#Empty: A, B, C
#An Edge: A -> B, C
#A Star (2 Edges): A->B->C
#Closed: A->B->C->A

#Triad types (per Davis & Leinhardt):
#003  A, B, C, empty triad.
#012  A->B, C 
#102  A<->B, C  
#021D A<-B->C 
#021U A->B<-C 
#021C A->B->C
#111D A<->B<-C
#111U A<->B->C
#030T A->B<-C, A->C
#030C A<-B<-C, A->C.
#201  A<->B<->C.
#120D A<-B->C, A<->C.
#120U A->B<-C, A<->C.
#120C A->B->C, A<->C.
#210  A->B<->C, A<->C.
#300  A<->B<->C, A<->C, completely connected.

triad.census(corrupt_net_d)

#Hierarchy Measures: Components,Cut Points, K-Cores, and Cliques
#Components: Components are maximally connected subgraphs (Wasserman and Faust 1994, p. 109). 
#Recall that community 7 has two large components and several small dyads and triads.
#There are two types of components: strong and weak.
#Strong components are components connected through directed paths (i --> j, j --> i)
#Weak components are components connected through semi-paths (--> i <-- j --> k)
components(corrupt_net_d, connected="strong")
components(corrupt_net_d, connected="weak")

#Which node belongs to which component?
Corrupt_Comp <- component.dist(corrupt_net_d, connected="strong")
Corrupt_Comp

Corrupt_Comp$membership # The component each node belongs to
Corrupt_Comp$csize      # The size of each component
Corrput_Comp$cdist      # The distribution of component sizes

#Cut-Sets and Cut-Points: Cut-sets describe the connectivity of the graph based on the removal of nodes, while cut-points describe
#...the connectivity of the graph based on the removal of lines (Harary 1969)
#k refers to the number of nodes or lines that would need to be removed to reduce the graph to a disconnected state.

cutpoints_strong <- cutpoints(corrupt_net_d,mode="dgraph", connected="strong",return.indicator=T)
set.vertex.attribute(corrupt_net_d, "cutpoint_strong", cutpoints_strong)

cutpoints_weak <- cutpoints(corrupt_net_d,mode="dgraph", connected="weak",return.indicator=T)
set.vertex.attribute(corrupt_net_d, "cutpoint_weak", cutpoints_weak)

gplot(corrupt_net_d,vertex.col=2+cutpoints(corrupt_net_d,mode="dgraph", connected="strong",return.indicator=T))
gplot(corrupt_net_d,vertex.col=2+cutpoints(corrupt_net_d,mode="dgraph", connected="weak",return.indicator=T))
#The plot only shows subgraphs consisting of nodes with a degree of 2 or more.
#The green nodes indicate cut-ponts where the removal of the node would separate one subgraph from another.

#Let's remove one of the cutpoints and count components again.
Corrupt_Cut <- corrupt_net_d[-11,-11]
#"-11" selects all the elments in the first row/column.
#So, AHS_Cut will be AHS_Network with node 1 removed.

components(Corrupt_Cut, connected="strong")  #There are 74 strong components in AHS_Cut compared to 73 in AHS_Network

#Bi-Components: Bi-Components refer to subgraphs that require at least the removal of two nodes or two lines to transform it into a 
#...disconnected set of nodes. 
#In large highly connected networks, we frequently analyze the properties of the largest bi-component to get a better understanding
#...of the social system represented by the network.
Corrupt_BiComp <- bicomponent.dist(corrupt_net_d) 

#Identify Cohesive Subgroups
#K-Cores: A k-core is a subgraph in which each node is adjacent to at least a minimum number of, k, to the other nodes in the subgraph.
#..., while a k-plex specifies the acceptable number of lines that can be absent from each node (Wasserman and Faust 1994, p. 266). 
kcores(corrupt_net_d) 
#Show the nesting of cores
Corrupt_kc<-kcores(corrupt_net_d,cmode="indegree")
gplot(corrupt_net_d,vertex.col=rainbow(max(Corrupt_kc)+1)[Corrupt_kc+1])

#Now, showing members of the 4-core only (All Nodes Have to Have a Degree of 4)
gplot(corrupt_net_d[Corrupt_kc>3,Corrupt_kc>3],vertex.col=rainbow(max(Corrupt_kc)+1)[Corrupt_kc[Corrupt_kc>3]+1])

#Cliques:  A clique is a maximally complete subgraph of three or more nodes.
#In other words, a clique consists of a subset of nodes, all of which are adjacent to each other, and where there are no other 
#...nodes that are also adjacent to all of the members of the clique (Luce and Perry 1949)

#We need to symmetrize recover all ties between i and j.
set.network.attribute(corrupt_net_d, "directed", FALSE) 

#The clique census returns a list with several important elements 
#Let's assign that list to an object we'll call AHS_Cliques.
#The clique.comembership parameter takes values "none" (no co-membership is computed),
#"sum" (the total number of shared cliques for each pair of nodes is computed),
#bysize" (separate clique co-membership is computed for each clique size)

Corrupt_Cliques <- clique.census(corrupt_net_d, mode = "graph", clique.comembership="sum")
Corrupt_Cliques # an object that now contains the results of the clique census

#The first element of the result list is clique.count: a matrix containing the number of cliques of different 
#...sizes (size = number of nodes in the clique).
#The first column (named Agg) gives you the total  number of cliques of each size,
#The rest of the columns show the number of cliques each node participates in.

#Note that this includes cliques of sizes 1 & 2. We have those when the largest fully connected structure includes just 1 or 2 nodes.
Corrupt_Cliques$clique.count

#The second element is the clique co-membership matrix:
Corrupt_Cliques$clique.comemb

# The third element of the clique census result is a list of all found cliques:
# (Remember that a list can have another list as its element)
Corrupt_Cliques$cliques # a full list of cliques, all sizes

Corrupt_Cliques$cliques[[1]] # cliques size 1
Corrupt_Cliques$cliques[[2]] # cliques of size 2
Corrupt_Cliques$cliques[[3]] # cliques of size 3
Corrupt_Cliques$cliques[[4]] # cliques of size 4

length(Corrupt_Cliques$cliques[[1]])
length(Corrupt_Cliques$cliques[[2]])
length(Corrupt_Cliques$cliques[[3]])
length(Corrupt_Cliques$cliques[[4]])

###########################
#   NODE LEVEL MEASURES   #
###########################

#Restoring Our Directed Network
set.network.attribute(corrupt_net_d, "directed", TRUE) 

#Reachability
#An actor is "reachable" by another if there exists any set of connections by which we can trace from the source to the target actor, 
#regardless of how many other nodes fall between them (Wasserman and Faust 1994, p. 132).
#If the network is a directed network, then it possible for actor i to be able to reach actor j, but for j not to be able to reach i.
#We can classify how connected one node is to another by considering the types of paths connecting them.
#Weakly Connected: The nodes are connected by a semi-path (--> i <--- j ---> k)
#Unilaterally Connected: The nodes are connected by a path (i --> j --> k)
#Strongly Connected: The nodes are connected by a path from i to k and a path from k to i.
#Recursively Connected: The nodes are strongly connected, and the nodes along the path from i to k and from k to i are the same in reverse order.
#e.g., i <--> j <--> k 

#Interpreting the reachability matrix, the first column indicates a specific node, the second an alter (alters can occur multiple times),
#and the third column indicates the number of paths connecting the two (total is a cumulative count of the number of paths in the network).
#For example, interpreting row 2, node 2 can reach node 235 through 235 paths (470-235), whereas in the middle of the list node 343 can reach node 1 through only 1 path.
Corrupt_Reach <- reachability(corrupt_net_d) 
??reachablity #For more information on this measure

#Degree Centrality: Total, In-Degree, Out-Degree

#In-Degree Centrality: The number of nodes adjacent to node i (Wasserman and Faust 1994, p. 126). i <--
InDegree <- degree(corrupt_net_d, cmode="indegree")
InDegree <- InDegree * .15                #Scaling in-degree to avoid high in-degree nodes from crowding out the rest of the nodes

set.vertex.attribute(corrupt_net_d, "InDegree", InDegree)

#Out-Degree Centrality: The number of nodes adjacent from node i (Wasserman and Faust, p. 126). i -->
OutDegree <- degree(corrupt_net_d, cmode="outdegree")
OutDegree <- OutDegree * .5                 #Scaling in-degree to avoid high in-degree nodes from crowding out the rest of the nodes

set.vertex.attribute(corrupt_net_d, "OutDegree", OutDegree)

#Total Degree Centrality: The Total Number of Adjacent Nodes (In-Degree + Out-Degree)
TotalDegree <- OutDegree + InDegree
TotalDegree <- TotalDegree * .4

set.vertex.attribute(corrupt_net_d, "TotalDegree", TotalDegree)

#Try Sizing by the Different Degrees
set.seed(12345)
as.data.frame(corrupt_net_d) %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_edges(color = "lightgray") +
  geom_nodes(color = Color_Race, size = InDegree) +       
  #geom_nodelabel_repel (color = Race, label = Race) +#   For networks with fewer nodes, we might want to label
  theme_blank() + 
  geom_density_2d()

#Path Centralities: Closeness Centrality, Information Centrality, Betweenness Centrality

#Closeness Centrality: Closeness centrality measures the geodesic distances of node i to all other nodes.
#Functionally, this measures range from 0 to 1, and is the inverse average distance between actor i and all other actors (Wasserman and Faust 1994, p. 185)
#This measure does not work well when there are disconnected components because the distances between components cannot be summed as
#...they are technically infinite. There are several work arounds, see Acton and Jasny's alternative below.

CN_Closeness <- closeness(corrupt_net_d, gmode="digraph", cmode="directed")
CN_Closeness
hist(CN_Closeness , xlab="Closness", prob=TRUE) 

#Alternative Approach to Measuring Closesness from the Geodesic Distances Matrix from Acton's and Jasny's Statnet Tutorial
Closeness <- function(x){ # Create an alternate closeness function!
  geo <- 1/geodist(x)$gdist # Get the matrix of 1/geodesic distance
  diag(geo) <- 0 # Define self-ties as 0
  apply(geo, 1, sum) # Return sum(1/geodist) for each vertex
}

Closeness <-  Closeness(corrupt_net_d)                        #Applying the function
Closeness
hist( Closeness , xlab="Alt. Closeness", prob=TRUE)         #Better behaved!

set.vertex.attribute(corrupt_net_d, "Closeness", Closeness)

#Information Centrality: Information Centrality measures the information flowing from node i.
#In general, actors with higher information centrality are predicted to have greater control over the flow of information within a network.
#Highly information-central individuals tend to have a large number of short paths to many others within the social structure.
?infocent  #For more information

CN_Info <- infocent(corrupt_net_d, rescale=TRUE)
CN_Info
hist(CN_Info , xlab="Information Centrality", prob=TRUE) 

set.vertex.attribute(corrupt_net_d, "Information", CN_Info)

gplot(corrupt_net_d, vertex.cex=(CN_Info)*250, gmode="graph") # Use w/gplot
#As suggested by the histogram there is relatively little variation in information centrality in this graph.

#Betweenness Centrality: The basic intuition behind Betweenness Centrality is that the actor between all the other actors in the 
#...has some control over the paths in the network. 
#Functionally, Betweenness Centrality is the ratio of the sum of all shortest paths linking j and k that includes node i over 
#...all the shortest paths linking j and k (Wasserman and Faust 1994, p. 191)

CN_Betweenness <- betweenness(corrupt_net_d, gmode="digraph")  
CN_Betweenness
hist(CN_Betweenness , xlab="Betweenness Centrality", prob=TRUE) 

set.vertex.attribute(corrupt_net_d, "Betweenness_Centrality", CN_Betweenness)

gplot(corrupt_net_d, vertex.cex=sqrt(CN_Betweenness)/25, gmode="digraph") 

#Comparing Closeness and Betweenness Centralities
cor(Closeness, CN_Betweenness)                             #Correlate our adjusted measure of closeness with betweenness
plot(Closeness, CN_Betweenness)                            #Plot the bivariate relationship

#Measures of Power in Influence Networks: Bonachich and Eigenvector Centrality

#Bonachich Centrality: The intuition behind Bonachich Power Centrality is that the power of node i is recursively defined 
#...by the sum of the power of its alters. 
#The nature of the recursion involved is then controlled by the power exponent: positive values imply that vertices become 
#...more powerful as their alters become more powerful (as occurs in cooperative relations), while negative values imply 
#...that vertices become more powerful only as their alters become weaker (as occurs in competitive or antagonistic relations).
?bonpow   #For more information about the measure

#Eigenvector Centrality: Conceptually, the logic behind eigenvectory centrality is that node i's influence is proportional to the 
#...to the centraltities' of the nodes adjacent to node i. In other words, we are important because we know highly connected people.
#Mathematically, we capture this concept by calculating the values of the first eigenvector of the graph's adjacency matrix.
?evcent   #For more information.

CN_Eigen <- evcent(corrupt_net_d)
CN_Eigen
hist(CN_Eigen , xlab="Eigenvector Centrality", prob=TRUE) 

set.vertex.attribute(corrupt_net_d, "Eigen_Centrality", CN_Eigen)

gplot(corrupt_net_d, vertex.cex=CN_Eigen*10, gmode="digraph") 


###########################
#   POSITIONAL ANALYSIS   #
###########################

#Burt's (1992) measures of structural holes are supported by iGraph and ego network variants of these measures are supported by egonet
#...the egonet package is compatable with the sna package.

#You can find descriptions and code to run Burt's measures in igraph at: http://igraph.org/r/doc/constraint.html

#Brokerage: The brokerage measure included in the SNA package builds on past work on borkerage (Marsden 1982), but is a more 
#...explicitly group oriented measure. Unlike Burt's (1992) measure, the Gould-Fernandez measure requires specifying a group variable
#...based on an attribute. I use race in the example below.

#Brokerage Roles: Group-Based Concept
#w_I: Coordinator Role (Mediates Within Group Contact)
#w_O: Itinerant Broker Role (Mediates Contact between Individuals in a group to which the actor does not belong)
#b_{IO}: Representative: (Mediates incoming contact from out-group members)
#b_{OI}: Gatekeeper: (Mediates outgoing contact from in-group members)
#b_O: Liason Role: (Mediates contact between individuals of two differnt groups, neither of which the actor belongs)
#t: Total or Cumulative Brokerage (Any of the above paths)
?brokerage   #for more information

Intl_Trial <- as.vector(nodes$intl_trial)

CN_Brokerage <- brokerage(corrupt_net_d, Intl_Trial)
CN_Brokerage
hist(CN_Brokerage$cl, xlab="Cumulative Brokerage", prob=TRUE) 

CN_Brokerage_Normalized <- round(CN_Brokerage$z.nli, 2) # Normalized values
CN_Brokerage_Raw <- CN_Brokerage$raw.nli # Raw values

CN_CBrokerage <- (CN_Brokerage$cl)
gplot(corrupt_net_d, vertex.cex=CN_CBrokerage*2, gmode="digraph") 

#Structural Equivalence
#Structural equivalence: Similarity/Distance Measures Include:
#Correlation
#Euclidean Distance
#Hamming Distance
#Gamma Correlation
sedist(corrupt_net_d, mode="digraph", method="correlation")

#Cluster based on structural equivalence:
CN_Clustering <- equiv.clust(corrupt_net_d, mode="digraph",plabels=network.vertex.names(corrupt_net_d))
CN_Clustering                        #Specification of the equivalence method used
plot(CN_Clustering)                  #Plot the dendrogram
rect.hclust(CN_Clustering$cluster, h=30)

#Generating a Block Model based on the Structural Equivalence Clustering
CN_BM <- blockmodel(corrupt_net_d, CN_Clustering, h=30)
CN_BM

#Extract the block image for Visualization
bimage <- CN_BM$block.model
bimage
bimage[is.nan(bimage)] <- 1

#Visualizing the block image (with self-reflexive ties)
gplot(bimage, diag=TRUE, edge.lwd=bimage*5, vertex.cex=sqrt(table(CN_BM$block.membership))/2,
      gmode="graph", vertex.sides=50, vertex.col=gray(1-diag(bimage)))

##### QEP Models #####

corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_intl_trials <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='International trials'))

# Now, plot the network:
plot(cn_net_intl_trials, vertex.cex=(cn_net_intl_trials %v% "country_size")/5, vertex.col = (cn_net_intl_trials %v% "col"))
summary(cn_net_intl_trials)                                   # Get an overall summary
print(cn_net_intl_trials)                                     # Simple print method
network.dyadcount(cn_net_intl_trials)                         # How many dyads in nflo?
network.edgecount(cn_net_intl_trials)                         # How many edges are present?
network.size(cn_net_intl_trials)                              # How large is the network?
gden(cn_net_intl_trials, mode = "digraph")                      # Network Density

#Obtaining subgraphs by relationships
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_business <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Business connection'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_company <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Company connection'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_complains <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Complaints for corruption'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_designates <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Designates in charge'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_enemies <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Enemies'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_facilitators <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Facilitators'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_family <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Family '))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_friends <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Friends'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_hhrr <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Human Rights violation'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_integrates_comp <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Integrates company'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_occupied_func <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Occupied functions'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_sanctioned <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Sanctioned by'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_student <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Student colleagues'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)
cn_net_contract <- delete.edges(corrupt_net_d,eid=which(corrupt_net_d%e%'type'!='Suscribed contract'))
corrupt_net_d <- network(edges, vertex.attr=nodes, matrix.type="edgelist", 
                         loops=F, multiple=F, ignore.eval = F, directed = T)

# Checking for integrity

# Economic: 733
network.size(cn_net_business)
network.edgecount(cn_net_business) 
network.size(cn_net_company)
network.edgecount(cn_net_company) 
network.size(cn_net_facilitators)
network.edgecount(cn_net_facilitators) 
network.size(cn_net_integrates_comp)
network.edgecount(cn_net_integrates_comp) 
network.size(cn_net_contract)
network.edgecount(cn_net_contract) 

# Political: 1273
network.size(cn_net_complains)
network.edgecount(cn_net_complains) 
network.size(cn_net_designates)
network.edgecount(cn_net_designates) 
network.size(cn_net_hhrr)
network.edgecount(cn_net_hhrr)
network.size(cn_net_occupied_func)
network.edgecount(cn_net_occupied_func) 
network.size(cn_net_sanctioned)
network.edgecount(cn_net_sanctioned) 

# Social: 129 edges
network.size(cn_net_family)
network.edgecount(cn_net_family) 
network.size(cn_net_friends)
network.edgecount(cn_net_friends) 
network.size(cn_net_enemies)
network.edgecount(cn_net_enemies) 
network.size(cn_net_student)
network.edgecount(cn_net_student) 

gden(cn_net_intl_trials, mode = "digraph") 
gden(cn_net_business, mode = "digraph")
gden(cn_net_company, mode = "digraph") 
gden(cn_net_facilitators, mode = "digraph") 
gden(cn_net_integrates_comp, mode = "digraph") 
gden(cn_net_contract, mode = "digraph") 
gden(cn_net_complains, mode = "digraph") 
gden(cn_net_designates, mode = "digraph") 
gden(cn_net_hhrr, mode = "digraph")
gden(cn_net_occupied_func, mode = "digraph") 
gden(cn_net_sanctioned, mode = "digraph") 
gden(cn_net_family, mode = "digraph") 
gden(cn_net_friends, mode = "digraph") 
gden(cn_net_enemies, mode = "digraph")
gden(cn_net_student, mode = "digraph")

list_full_g <- list(cn_net_complains,cn_net_designates,cn_net_hhrr,cn_net_occupied_func,cn_net_sanctioned,cn_net_business,cn_net_company,cn_net_facilitators,cn_net_integrates_comp,cn_net_contract,cn_net_family,cn_net_friends,cn_net_hhrr,cn_net_enemies,cn_net_student)

for (g in list_full_g) {
  
print(paste(centralization(g, degree, cmode="indegree"),centralization(g, degree, cmode="outdegree"),gtrans(g, mode="digraph", measure = c("weak"),use.adjacency = FALSE), gtrans(g, mode="digraph", measure = c("weakcensus"), use.adjacency = FALSE), sep=" "))

}

# Creating variables as matrices

# Dependent Variable
cn_m_intl_trials <- as.matrix(cn_net_intl_trials) 

# Independent Variables
cn_m_business <- as.matrix(cn_net_business)
cn_m_company <- as.matrix(cn_net_company)
cn_m_complains <- as.matrix(cn_net_complains)
cn_m_designates <- as.matrix(cn_net_designates)
cn_m_enemies <- as.matrix(cn_net_enemies)
cn_m_facilitators <- as.matrix(cn_net_facilitators)
cn_m_family <- as.matrix(cn_net_family)
cn_m_friends <- as.matrix(cn_net_friends)
cn_m_hhrr <- as.matrix(cn_net_hhrr)
cn_m_integrages_comp <- as.matrix(cn_net_integrates_comp)
cn_m_occupied_func <- as.matrix(cn_net_occupied_func)
cn_m_sanctioned <- as.matrix(cn_net_sanctioned)
cn_m_student <- as.matrix(cn_net_student)
cn_m_contract <- as.matrix(cn_net_contract)

# Creating lists with independent variables by dimension
list_economy <- list(cn_m_business,cn_m_company,cn_m_facilitators,cn_m_integrages_comp,cn_m_contract)
list_social <- list(cn_m_family,cn_m_friends,cn_m_enemies,cn_m_student)
list_political <- list(cn_m_complains,cn_m_designates,cn_m_hhrr,cn_m_occupied_func,cn_m_sanctioned)
list_full <- list(cn_m_complains,cn_m_designates,cn_m_hhrr,cn_m_occupied_func,cn_m_sanctioned,cn_m_business,cn_m_company,cn_m_facilitators,cn_m_integrages_comp,cn_m_contract,cn_m_complains,cn_m_designates,cn_m_hhrr,cn_m_occupied_func,cn_m_sanctioned)

# Test
m.test <- netlogit(cn_m_intl_trials,list(cn_m_business),reps=1)

mean(m.test$dist[,2] >= m.test$tstat[2])
mean(m.test$dist[,2] <= m.test$tstat[2])
mean(abs(m.test$dist[,2]) <= abs(m.test$tstat[2]))

plot(density(m.test$dist[,2]))
abline(v=m.test$tstat[2])

# Models

m.full <- netlogit(cn_m_intl_trials,list_full,reps=100)
m.econ <- netlogit(cn_m_intl_trials,list_economy,reps=100)
m.pol <- netlogit(cn_m_intl_trials,list_political,reps=100)
m.soc <- netlogit(cn_m_intl_trials,list_personal,reps=100)

std_coef <- function (model) {
  sigma2 <- sum(model$residuals ^ 2) / model$df.residual
  Rinv <- backsolve(model$qr$qr, diag(model$rank), 0)
  sqrt(rowSums(Rinv ^ 2) * sigma2)
}

std_coef(m.pol)

m.pol$qr$qr

mean(m.test$dist[,2] >= m.test$tstat[2])
mean(m.test$dist[,2] <= m.test$tstat[2])
mean(abs(m.test$dist[,2]) <= abs(m.test$tstat[2]))

plot(density(m.test$dist[,2]))
abline(v=m.test$tstat[2])


#### Saving Nodes data ####

Closeness[1585] <- 0
CN_Betweenness[1585] <- 0
CN_Eigen[1585] <- 0
CN_Info[1585] <- 0
CN_CBrokerage[1585] <- 0
cutpoints_strong[1585] <- FALSE
cutpoints_weak[1585] <- FALSE
InDegree[1585] <- 0
OutDegree[1585] <- 0
TotalDegree[1585] <- 0

df <- data.frame(unlist(Closeness),unlist(CN_Betweenness),unlist(CN_CBrokerage),unlist(CN_Eigen),unlist(CN_Info),unlist(cutpoints_strong),unlist(cutpoints_weak),unlist(InDegree),unlist(OutDegree),unlist(TotalDegree))
write.csv(df,"output_nodes_analysis.csv", row.names = FALSE)

df_cn_brokerage_norm <- as.data.frame(CN_Brokerage_Normalized)
write.csv(df_cn_brokerage_norm,"output_brokerage.csv", row.names = TRUE)

# Saving the Workspace

save.image(file = "workspace.RData")

###### Documentation ######

# Network Help File
# https://rdrr.io/github/statnet/network/man/network.html

# Colors in R
# http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

# Preparing Network Data in R
# https://www.mjdenny.com/Preparing_Network_Data_In_R.html

# Descriptive Statistics for Networks (R Labs, Duke Network Analysis Center)
# https://dnac.ssri.duke.edu/r-labs/2017/02_descriptive_statistics.php

# SIENA (2) Lab
# https://dnac.ssri.duke.edu/r-labs/2018/siena_intuition.php

# Analyzing Brokerage in Networks Using R
# https://rpubs.com/pjmurphy/320424

# Using QAP Logistic Regression
# https://rpubs.com/pjmurphy/338798