# CSVデータの読み込み
 data = read.csv('pca.csv',head=T,row.names=1)
 
# データの標準化(割合であったり、回数であったりするためデータの単純比較ができないので)
 pca = prcomp(data, scale=T)
# 主成分分析結果の表示
 pca
 
# 主成分分析の結果をプロットします。
plot(pca$x[,1], pca$x[,2], pch=19, cex=0.6, xlab="PC1", ylab="PC2")
text(pca$x[,1], pca$x[,2], rownames(data), cex=0.6, pos=4)
