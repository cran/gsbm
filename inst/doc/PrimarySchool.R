## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load data----------------------------------------------------------------
library(Matrix)
library(gsbm)
library(igraph)
library(RColorBrewer)
library(combinat)

data(PrimarySchool)
A <- PrimarySchool$A
class <- PrimarySchool$class
n <- dim(A)[1]
class.names <- levels(class)
class.effectifs <- table(class)

## ----compute frequency of interaction-----------------------------------------
# Compute frequency of interactions
K = length(class.names)
Q = matrix(rep(0, K*K), ncol = K, nrow = K)
for (k1 in 1:K){
  com1 = class.names[k1]
  for (k2 in 1:K){
    com2 = class.names[k2]
    Q[k1, k2] <- mean(A[class == com1, class == com2], na.rm = TRUE)
    Q[k2, k1] <- Q[k1, k2]
  }
}
Q <- round(Q, 2)
colnames(Q) <- class.names
rownames(Q) <- class.names
print(Q)

## ----Plot network-------------------------------------------------------------
# Vertex labels
v.label <- rep(NA, n)

# Graph without NA's
A.noNA <- matrix(0, ncol = n, nrow = n)
A.noNA[!is.na(A)] <- A[!is.na(A)]
g <- graph_from_adjacency_matrix(A.noNA, mode = "undirected")

# Layout
W<- matrix(0, n, n)
for (i in 1:n){
  for (j in 1:n){
      if ((!is.na(A[i,j])) && A[i,j]==1){
        if (class[i]==class[j] && class[i]!= "teachers"){
          W[i,j] <- 2
        }
        else{
          W[i,j] <- 1
        }
      }
  }
}
gw <- graph_from_adjacency_matrix(W, mode = "undirected")
lay <- layout_with_fr(gw)
pal <- brewer.pal(11, "Paired")

# Plot network
plot(g, vertex.label = v.label, vertex.size = (2 + colSums(A, na.rm =T)/5), vertex.color = class, layout = lay, palette = pal, edge.width = 0.3, edge.curved=0.1)

legend("bottomleft", legend = class.names, col = pal, pch = 1, bty = 'n', cex = 0.5, ncol = 1)

## ----run GSBM-----------------------------------------------------------------
# Investigate these outliers
T1<-Sys.time()
res <- gsbm_mcgd(A, 8.5, 8, maxit = 100)
T2<-Sys.time()
Tdiff= difftime(T1, T2)
print(paste0("running time : ", Tdiff))
outliers <- which(colSums(res$S)>0)
no <- length(outliers)
outliers.class <- class[outliers]
outliers.Q <- matrix(0,nrow = no, ncol = K)
for (o in 1:no){
  for (k in 1:K){
    comk = class.names[k]
    outliers.Q[o,k] <- mean(A[outliers[o],class == comk], na.rm = T)
  }
  print(paste0("node ", outliers[o], " is in class ", outliers.class[o]))
  print(outliers.Q[o,], 2)
  print("")
}

## ----recover the classes------------------------------------------------------
# Compute svd and run k-means
L.U <- svd(res$L[-outliers, -outliers], nu = 10)$u 
est_class <- kmeans(L.U, 10, nstart = 100)$cluster
length(est_class)

# Recover the class among the permutations of labels
equ_label = combinat::permn(1:10)
error <- rep(0, length(equ_label))
for (i in 1:length(equ_label)){
  error[i] <- sum(class[-outliers] != (class.names[unlist(equ_label[i])])[est_class])
}
error[which.min(error)]
class_est <- class.names[unlist(equ_label[which.min(error)])][est_class]
class_est_out <- rep(NA, n)
class_est_out[outliers] <- "outlier"
class_est_out[-outliers] <- class_est
class_est_out <- as.factor(class_est_out)

# Plot network
plot(g, vertex.label = v.label, vertex.size = (2 + colSums(A, na.rm =T)/5), vertex.color = class_est_out, layout = lay, palette = pal, edge.width = 0.3, edge.curved=0.1)

class_est_out <- as.factor(class_est_out)



