
library(gplots)

#Generates a random matrix

data <- matrix(rnorm(500), 100, 5, dimnames=list(paste("g", 1:100, sep=""), paste("t", 1:5, sep="")))

#export the matrix as csv

write.csv(data, "df.csv", row.names=TRUE)

#read the csv matrix

a <- read.csv("df.csv")

# remove the first column containing gene names and convert to matrix

a_numeric <- as.matrix(a[, -1]) 

head(a_numeric)

#hierarchical clustering

heatmap.2(a_numeric)



#Stepwise Approach with Tree Cutting



## Row- and column-wise clustering 

hr <- hclust(as.dist(1-cor(t(a_numeric), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(a_numeric, method="spearman")), method="complete") 

## Tree cutting

mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 

## Plot heatmap 

mycol <- colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
heatmap.2(a_numeric, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 

#Principal Component Analysis (PCA)

pca <- prcomp(a_numeric, scale=T)
summary(pca) # Prints variance summary for all principal components

plot(pca$x, pch=20, col="blue", type="n") # To plot dots, drop type="n"
text(pca$x, rownames(pca$x), cex=0.8)



#Clustering
#Scaling
dim(a_numeric)

## Scaling
a_numericscaled <- t(scale(t(a_numeric))) # Centers and scales a_numeric row-wise
apply(a_numericscaled, 1, sd)



#Distance Matrices
#Euclidean distance matrix

dist(a_numeric[1:4,], method = "euclidean")

#Correlation-based distance matrix
#Correlation matrix

c <- cor(t(a_numeric), method="pearson") 
as.matrix(c)[1:4,1:4]

#Correlation-based distance matrix

d <- as.dist(1-c)
as.matrix(d)[1:4,1:4]

#Hierarchical Clustering with hclust

hr <- hclust(d, method = "complete", members=NULL)
names(hr)

par(mfrow = c(1, 2)); plot(hr, hang = 0.1); plot(hr, hang = -1) 


#Tree plotting I

plot(as.dendrogram(hr), edgePar=list(col=3, lwd=4), horiz=T) 

#Tree plotting II

library(ape) 
plot.phylo(as.phylo(hr), type="p", edge.col=4, edge.width=2, 
           show.node.label=TRUE, no.margin=TRUE)

#Tree Cutting

hr

## Print row labels in the order they appear in the tree
hr$labels[hr$order] 

#Tree cutting with cutree

mycl <- cutree(hr, h=max(hr$height)/2)
mycl[hr$labels[hr$order]] 

#Heatmaps

library(gplots)
heatmap.2(a_numeric, col=redgreen(75))


#All in one step: clustering and heatmap plotting

library(pheatmap); library("RColorBrewer")
pheatmap(a_numeric, color=brewer.pal(9,"Blues"))


#Customizing heatmaps

hc <- hclust(as.dist(1-cor(a_numeric, method="spearman")), method="complete")
mycol <- colorpanel(40, "darkblue", "yellow", "white")
heatmap.2(a_numeric, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(mycl))

#K-Means Clustering with PAM

library(cluster)
pamy <- pam(d, 4)
(kmcol <- pamy$clustering)


heatmap.2(a_numeric, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(kmcol))


#K-Means Fuzzy Clustering

library(cluster)
fannyy <- fanny(d, k=4, memb.exp = 1.5)
round(fannyy$membership, 2)[1:4,]

fannyy$clustering 

## Returns multiple cluster memberships for coefficient above a certain 
## value (here >0.1)
fannyyMA <- round(fannyy$membership, 2) > 0.10 
apply(fannyyMA, 1, function(x) paste(which(x), collapse="_"))



#Principal Component Analysis (PCA)

library(scatterplot3d)
pca <- prcomp(a_numeric, scale=TRUE)
names(pca)

summary(pca) # Prints variance summary for all principal components.

scatterplot3d(pca$x[,1:3], pch=20, color="blue") 






