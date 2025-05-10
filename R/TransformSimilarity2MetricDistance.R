TransformSimilarity2MetricDistance=function(Similarity){
  #1986Legendre_Metric_and_Euclidean_properties_of_dissimilarity_c.pdf
  #page10
  #siehe auch [Bock, 1974, pp.77-79]
  if(!is.matrix(Similarity)){
    warning('Matrix Input expected, calling as.matrix')
    Similarity=as.matrix(Similarity)
  }
  n=nrow(Similarity)
  if(any(diag(Similarity)!=1)){
    stop('Transformation ist not allowed because diagonal does not equal one.')
  }
  ab = range(Similarity)
  if(length(ab) != 2){
    warning("Range could not be successfuly computed.")
    return(NULL)
  }
  if(!is.vector(ab) | !is.numeric(ab)){
    warning("Range could not be successfuly computed.")
    return(NULL)
  }
  if(any(is.na(ab))){
    warning("Range could not be successfuly computed.")
    return(NULL)
  }
  if(ab[1] < 0 | ab[2] > 1){
    stop('Transformation ist not allowed because eiter negative values or values higher one given.')
  }
  Distance=sqrt((1-Similarity))
  return(Distance)
}
