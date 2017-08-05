library(igraph)
set.seed(449)

biogrid <- read.delim("./data/BIOGRID.txt",stringsAsFactors = F)
names(biogrid)
attach(biogrid)

HSnet <- graph.data.frame(
        data.frame(Entrez.Gene.Interactor.A,Entrez.Gene.Interactor.B),
        directed=F)
#plot(HSnet)

A <- get.adjacency(HSnet)
# multiple edges
A <- A[1:15,1:15]

g1<-graph.adjacency(A)

png("Adjacency.png", width=800,height = 600)

col=c("yellow", "slategray", "tan1","lightgray", "tan3",
      "lightyellow" ,"lightskyblue", "hotpink1","plum2", "wheat1",
      "lightsalmon", "orange", "lightblue", "lightgoldenrod", "lightpink")

coords <- layout_(g1, as_star())
# plot(g1,rescale=F, # axes=TRUE, layout=coords,ylim=c(-2,2),xlim=c(-2,2),
#      asp = 0,
#      vertex.size = 5, vertex.label.cex = 0.8)
plot(g1, vertex.size = 10, vertex.label.cex = 0.8, asp = 0, axes = FALSE,
     layout = coords, vertex.color= col)
dev.off()
# the following is FALSE if the graph is not simple
is.simple(HSnet)

# remove multiple edges and self-loops
HSnet <- simplify(HSnet, remove.multiple = TRUE, remove.loops = TRUE,
                  edge.attr.comb = getIgraphOpt("edge.attr.comb"))
is.simple(HSnet)
A <- get.adjacency(HSnet)
# only single edges now
A <- A[1:15,1:15]
print(A)
g2<-graph.adjacency(A);

png("Adjacency3.png", width=800,height = 600)

col=c("yellow", "slategray", "tan1","lightgray", "tan3",
      "lightyellow" ,"lightskyblue", "hotpink1","plum2", "wheat1",
      "lightsalmon", "orange", "lightblue", "lightgoldenrod", "lightpink")

coords <- layout_(g2, as_star())
# plot(g2,rescale=F, # axes=TRUE, layout=coords,ylim=c(-2,2),xlim=c(-2,2),
#      asp = 0,
#      vertex.size = 5, vertex.label.cex = 0.8)
plot(g2, vertex.size = 10, vertex.label.cex = 0.8, asp = 0, axes = FALSE,
     layout = coords, vertex.color= col)
dev.off()

# for this application we remove nodes of very high degree;
# these are usually house-keeping
# genes that are necessary to keep a cell alive,
# but are usually not specific to a particular disease.
overly.attached.proteins <- which(igraph::degree(HSnet)>1000)
HSnet <- igraph::delete.vertices(HSnet, overly.attached.proteins )
# the following is TRUE if the graph is connected.
igraph::is.connected(HSnet)

# read the gene-id table
gene.table <- read.delim("./data/gene-id-table.txt")
names(gene.table)

# read the scores for Autism
gene.score<-read.csv("./data/gene-score.csv",stringsAsFactors=F)
attach(gene.score)
names(gene.score)
# display the scores
unique(Score)
# identify the genes that have significant scores
signif.scores<-c("3","1S","1","2S","2","3S")
signif.genes<-Gene.Symbol[which(Score %in% signif.scores)]
signif.EIDs <- gene.table[which(gene.table[,1] %in% signif.genes),2]
# Now use the protein interaction network HSnet, created previously,
# to determine the genes that are present in the network and
# known to cause Autism
geneEIDs <- as.numeric(V(HSnet)$name)
HSnetN<-HSnet
V(HSnetN)$name<-1:length(V(HSnet))
signif.ids<-which(geneEIDs %in% signif.EIDs)
length(signif.ids)

source('steiner_tree.R')
# Identify the Steiner Tree and note the time this function call takes
system.time(HS.stree <- steiner_tree(terminals=signif.ids, graph=HSnetN))
# Output the overlap between significant Autism and vertices in the Steiner Tree
length(intersect(signif.ids,V(HS.stree)$name))
labels<-gene.table[as.numeric(V(HS.stree)$name),1]
labels<-as.character(labels)
# identify the genes that have significant scores, and assign the color “red” to them
colors<-rep("skyblue",length(V(HS.stree)))
colors[which(as.numeric(V(HS.stree)$name) %in% signif.ids)] = "red"
# assign colors to the vertices of the tree
V(HS.stree)$color = colors
# plot and save to file
pdf("ASD_interactome.pdf",width=12, height=12)

system.time(plot(HS.stree,vertex.label=labels,vertex.size=5,vertex.label.cex=0.8))
dev.off()

png("ASD_interactome.png", width=800,height = 600)
system.time(plot(HS.stree,vertex.label=labels,vertex.size=5,vertex.label.cex=0.8))
plot(HS.stree, vertex.label=labels,vertex.size=5,vertex.label.cex=0.8,
     asp = -1, axes = FALSE)

dev.off()

library(sna, quietly=TRUE)
# a function that computes the connectivity scores for a network
# here the scores are diameters of the network and average geodesic
# distance between any two nodes
c.scores<-function(graph) {
        n<-length(V(graph))
        sp<-shortest.paths(graph)
        neighbors<-sum(sp==1)/2
        neighbors2<-sum(sp==2)/2
        return(c(2*neighbors/(n*(n-1)),2*neighbors2/(n*(n-1))))
}
clus<-clusters(HSnetN, mode=c("weak"))
connected.ids<-which(clus$membership==1)
length(connected.ids)
# Generate N randomly chosen subnetworks. Note: this will take a while if N is set large.
N<-50
strees<-list(N)
effs<-numeric(N)
nei<-numeric(N)
nei2<-numeric(N)
for (i in 1:N){
        new.ids<-sample(x=connected.ids,size=length(signif.ids))
        strees[[i]] <- steiner_tree(terminals=new.ids, graph=HSnetN)
        effs[i]<-efficiency(get.adjacency(strees[[i]],sparse=F))
        cs<-c.scores(strees[[i]])
        nei[i]<-cs[1]
        nei2[i]<-cs[2]
}

# print the efficiencies and connectivity scores for each of the N random graphs
effs
nei
nei2
# Finally, print the efficiency score and connectivity scores for the
# Autism Interactome
efficiency(get.adjacency(HS.stree,sparse=F))
c.scores(HS.stree)

# compute the betweeness centrality scores for each node
betweeness_centrality_scores = igraph::betweenness(HS.stree)
# now identify only those NOT already known to be significant
significant_centrality = c()
count = 0
for (i in 1:length(betweeness_centrality_scores)){
        if (!(as.numeric(names(betweeness_centrality_scores[i])) %in% signif.ids)) {
                significant_centrality = c(significant_centrality,
                                           betweeness_centrality_scores[i])
        }
}
# sort
significant_centrality = sort(significant_centrality, decreasing=TRUE)
length(significant_centrality)

