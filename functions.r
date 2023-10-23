Euclidean <- function(vec1, vec2, n=5) {
sqr.sum = 0
for (i in seq(n)) {sqr.sum <- sqr.sum + (vec1[i]-vec2[i])^2}
dis <- sqrt(sqr.sum)
return(dis)
}

shrinkage.vector <- function(med) {
N <- nrow(med)
dis.list <- vector()
for (i in seq(N*N)){
row = (i-1) %/% N +1
col = (i-1) %% N +1 
print(paste(row, "-", col))
dis.list[i] <- Euclidean(med[row,], med[col,])
}
dis.list <- unlist(dis.list)
A <- matrix(dis.list, nrow=N, byrow=TRUE)
A.inv <- solve(A)
a2 <- A.inv %*% rep(1,N)
return(a2)
}
