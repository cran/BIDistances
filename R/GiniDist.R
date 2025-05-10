GiniDist=function(Data){
  #INPUT:
  #Data [1:n,1:d]   data matrix of n cases and d feautures
  # Output
  # DistanceMatrix [1:n,1:n]    matrix of distances, symmetric
  #author MCT, 12,2023
  #transponieren, damit wir die zeilen als variablen ansehen koennen fuer die eigentlich wasserstein1d gilt
  Data=t(Data)
  if (!requireNamespace('ineq', quietly = TRUE)) {
    message(
      'Subordinate package (ineq) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return("Subordinate package (transineqport) is missing.
                Please install the package which is defined in 'Suggests'."
    )
  }
  num_cols <- ncol(Data)

  distance_matrix <- matrix(NA, nrow=num_cols, ncol=num_cols)

  dimnames(distance_matrix) <- list(colnames(Data), colnames(Data))

  gini_coefficient <- function(distribution) {
    return(ineq::ineq(distribution, type = "Gini"))
  }

  compute_gini <- function(col1, col2) {
    gini1 <- gini_coefficient(col1)
    gini2 <- gini_coefficient(col2)
    return(gini_distance <- abs(gini1 - gini2))
  }
  # Compute the distance for each pair of columns
  for (i in 1:num_cols) {
    for (j in 1:num_cols) {
      if (i > j) {
        distance_matrix[i, j] <- compute_gini(Data[, i], Data[, j])
      } else {
        distance_matrix[i, j] <- 0  # Distance to self is zero
      }
    }
  }
  #Zu tun damit aus unteren diagonal matrix eine komplette matrix wird
  return(as.matrix(as.dist(distance_matrix)))
}
