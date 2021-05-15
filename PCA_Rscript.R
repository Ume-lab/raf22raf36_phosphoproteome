
 data = read.csv('pca.csv',head=T,row.names=1)
 
 pca = prcomp(data, scale=T)
 pca
 
plot(pca$x[,1], pca$x[,2], pch=19, cex=0.6, xlab="PC1", ylab="PC2")
text(pca$x[,1], pca$x[,2], rownames(data), cex=0.6, pos=4)
