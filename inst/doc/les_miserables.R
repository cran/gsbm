## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(igraph)
library(gsbm)
library(missSBM)
library(RColorBrewer)

## ----load graph----------------------------------------------------------
data(les_miserables)
A<- les_miserables$A
names <- les_miserables$names
net <- graph_from_adjacency_matrix(A, mode = "undirected")
V(net)$name <- names
V(net)$color <- "gray80"
deg <- degree(net, mode="all")
V(net)$size <- deg
plot(net, vertex.label.cex = 0.4)

## ----classical_SBM-------------------------------------------------------
vBlocks <- 1:10
collection_sbm <- missSBM::estimate(prepare_data(A), vBlocks = vBlocks, sampling = "node")
L_missSBM <- collection_sbm$bestModel$fittedSBM$connectProb
colo <- round(collection_sbm$bestModel$fittedSBM$blocks)
colo <- sapply(1:nrow(A), function(i) which.max(colo[i,]))
pal3 <- brewer.pal(10, "Set3")
V(net)$color <- pal3[colo]
V(net)$label <- NA
V(net)$size <- deg
plot(net)

## ----robust_SBM----------------------------------------------------------
lambda1 <- 4
lambda2 <- 5
res <- gsbm_mcgd(A, lambda1 = lambda1, lambda2 = lambda2)
outliers <- names[which(colSums(res$S)>0)]
sv <- svd(res$L)
pc <- sv$u[,1:4]
rownames(pc) <- names
pc <- pc[setdiff(names, outliers),]
com <- kmeans(pc, centers=4, nstart=50)
com$cluster
colo2 <- 1:nrow(A)
names(colo2) <- names
comu <- com$cluster
comu[which(comu==4)] <- 6
colo2[setdiff(names, outliers)] <- pal3[comu]
colo2[outliers] <- "red"
labels <- names(A)
names(labels) <- names(A)
labels[setdiff(names, outliers)] <- NA
V(net)$label <- NA
V(net)$color <- colo2
V(net)$size <- deg
E(net)$arrow.size <- 5
plot(net, vertex.label.dist=20)

