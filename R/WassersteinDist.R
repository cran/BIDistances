WassersteinDist=function(Data,p=1,InverseWeighting=FALSE){
#INPUT:
  #Data [1:n,1:d]   data matrix of n cases and d feautures
# Output
# DistanceMatrix [1:n,1:n]    matrix of distances, symmetric
#author MCT, 12,2023
  #transponieren, damit wir die zeilen als variablen ansehen koennen fuer die eigentlich wasserstein1d gilt
  Data=t(Data)
  if (!requireNamespace('transport', quietly = TRUE)) {
    message(
      'Subordinate package (transport) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return("Subordinate package (transport) is missing.
                Please install the package which is defined in 'Suggests'."
    )
  }
  num_cols <- ncol(Data)

distance_matrix <- matrix(NA, nrow=num_cols, ncol=num_cols)

dimnames(distance_matrix) <- list(colnames(Data), colnames(Data))

compute_wasserstein <- function(col1, col2,p=1,InverseWeighting) {
  if(isTRUE(InverseWeighting)){
    weights1 <- rep(1/length(col1), length(col1))
    weights2 <- rep(1/length(col2), length(col2))
    return(transport::wasserstein1d(a = col1, b = col2, wa = weights1, wb = weights2,p = p))
  }else{
    return(transport::wasserstein1d(a = col1, b = col2,p = p))
  }
}
# Compute the distance for each pair of columns
for (i in 1:num_cols) {
  for (j in 1:num_cols) {
    if (i > j) {
      distance_matrix[i, j] <- compute_wasserstein(Data[, i], Data[, j],p=p,InverseWeighting=InverseWeighting)
    } else {
      distance_matrix[i, j] <- 0  # Distance to self is zero
    }
  }
}
#Zu tun damit aus unteren diagonal matrix eine komplette matrix wird
return(as.matrix(as.dist(distance_matrix)))
}
