#
###  http://www.kateto.net/polnet2015
#

install.packages("igraph")
install.packages("network") 
install.packages("sna")
install.packages("ndtv")

library("igraph")
library("network") 
library("sna")
library("ndtv")

setwd("/home/vimal/Desktop/")
nodes <- read.csv("Vertex.csv", header=T, as.is=T)
links <- read.csv("Edge.csv", header=T, as.is=T)
net <- graph.data.frame(links, nodes, directed=T)


### Fixed node size
#V(net)$size <- 10
### Node size according to  degree of node
deg <- igraph::degree(net, mode="all")
V(net)$size <- deg *0.3
### Node color according to  degree of node
colrs <- rainbow(length(V(net))+1)
V(net)$color <- colrs[V(net)$size]
#V(net)$color <- 'lightblue'



l <- layout.lgl


plot(net, layout=l,vertex.label.cex=0.3,edge.arrow.size=.2,edge.color="grey", vertex.frame.color="#ffffff"
     , vertex.label.color="black", main="EC network")
### Plot Color Legend
plot(rep(1,length(V(net))),col=rainbow(length(V(net))),pch=19,cex=sort(V(net)$size),xlab="Degree of node")
###########################

V(net)$vertex <- 2

l <- layout.circle(net)


plot(net, layout=l,vertex.label.cex=0.8,edge.arrow.size=.2,edge.color="orange",vertex.color="orange", vertex.frame.color="#ffffff"
     , vertex.label.color="black", main="Full network")


l <- layout.random(net)
l <- layout.lgl
l <- layout.fruchterman.reingold(net, repulserad=vcount(net)^3, area=vcount(net)^2.4)
l <- layout.circle(net)

l <- layout.fruchterman.reingold(net)

a <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)







plot(net, layout=a,vertex.label.cex=0.3,edge.arrow.size=.2,edge.color="orange",vertex.color="orange", vertex.frame.color="#ffffff",
     vertex.label=V(net)$media, vertex.label.color="black", main="Full network TF family")




cut.off <- mean(links$weight) 
net.sp <- delete.edges(net, E(net)[weight<cut.off])
#l <- layout.fruchterman.reingold(net.sp, repulserad=vcount(net)^2.1)
l <- layout.sphere(net.sp)
plot(net.sp, layout=l,vertex.label=V(net)$media,vertex.label.cex=0.8,edge.arrow.size=.2,edge.color="orange",vertex.color="orange", vertex.frame.color="#ffffff",vertex.label.color="black", main="Key part of network") 
