library(wordspace)
library(cluster)
library(wordspace)

dem <- function(X){
  D <- philentropy::distance(t(X), method = "euclidean")
  tt <- matrix(NA,nrow = ncol(X), ncol = 1)
  for (i in 1:ncol(X)) {
    tt[i] <- D[i,][order(D[i,])[7]]    
  }
  return(tt)
}

euclidean <- function(a, b) sqrt(sum((a - b)^2))

make.similarity <- function(my.data,tt) {
  N <- ncol(my.data)
  S <- matrix(rep(NA,N^2), ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      S[i,j] <- exp((-euclidean(my.data[,i], my.data[,j])^2/(tt[i]*tt[j])))
    }
    S[i,i] <- 0
  }
  return(S)
}

graph_laplacian <- function(W, normalized = TRUE)
{
  stopifnot(nrow(W) == ncol(W)) 
  
  g = rowSums(W) # degrees of vertices
  n = nrow(W)
  
  if(normalized)
  {
    D_half = diag(1 / sqrt(g) )
    return( diag(n) - D_half %*% W %*% D_half ) # I-D^-1/2WD^-1/2
  }
  else
  {
    return( diag(g) - W )
  }
}

set.seed(2023)
t1 <- dem(std.dna)
t2 <- dem(std.rna)
t3 <- dem(std.tmt)


S_rna = make.similarity(std.rna, t2)
LS_rna = graph_laplacian(S_rna,normalized = T)
ei_rna = eigen(LS_rna, symmetric = TRUE)

S_tmt = make.similarity(std.tmt, t3)
LS_tmt = graph_laplacian(S_tmt,normalized = T)
ei_tmt = eigen(LS_tmt, symmetric = TRUE)

S_dna = make.similarity(std.dna, t1)
LS_dna= graph_laplacian(S_dna,normalized = T)
ei_dna = eigen(LS_dna, symmetric = TRUE)


silhouette_score <- function(k){
  n <- 79
  ev_rna <- (ei_rna$vectors[,(n-k+1):n]) 
  ev_tmt <- (ei_tmt$vectors[,(n-k+1):n])
  ev_dna <- (ei_dna$vectors[,(n-k+1):n])
  
  L_mod1 <- LS_rna+LS_dna+LS_tmt-
    (ev_rna%*%t(ev_rna)+ev_dna%*%t(ev_dna)+ ev_tmt%*%t(ev_tmt))
  
  ei_mod1 = eigen(L_mod1, symmetric = TRUE)
  
  
  ev_mod <- ei_mod1$vectors[,(n-k+1):(n)]
  n_ev_mod <- normalize.rows(ev_mod)
  km <- kmeans(n_ev_mod, centers=k, nstart=10)
  ss <- silhouette(km$cluster, dist(n_ev_mod))
  mean(ss[, 3])
  
}

set.seed(2023)
library(wordspace)
k <- 2:8
avg_sil <- sapply(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
title("ROSMAP: Optimal number of clusters")
