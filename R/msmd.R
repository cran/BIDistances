msmd = function(Values1, Values2, ParameterC){
  #
  # https://github.com/JanaHolznigenkemper/msm_speedup/tree/main/src/MSMDistances
  #
  # V = msmd(required variable setting)
  # V = msmd(require + exemplary optional variable setting)
  #
  # DESCRIPTION
  # Time Warp Edit Distance for two (different) time series
  #
  # INPUT
  # Values1[1:n]        Numeric vector with values of the first time series
  # Values2[1:m]        Numeric vector with values of the second time series
  # Time1[1:n]          Numeric vector with time stamps of the first time series
  # Time2[1:m]          Numeric vector with time stamps of the second time series
  # Nu                  Numeric: Elasticity parameter - nu >=0 needed for
  #                     distance measure
  # Lambda              Numeric: Penalty for deletion operation
  # Degree              Integer: Degree of the p norm for local cost
  #
  # OPTIONAL
  #
  #
  # OUTPUT
  # List with elements
  # Value         Distance measure
  #
  # DETAILS
  # abc
  #
  # Authors: Marteau, Pierre-Francois 2009, translated from C to R by Quirin Stier, 09/2023
  #

  if(!is.numeric(Values1) | !is.vector(Values1)){
    stop("twed: Values1 must be a numeric vector.")
  }

  if(!is.numeric(Values2) | !is.vector(Values2)){
    stop("twed: Values2 must be a numeric vector.")
  }


  N1          = length(Values1)
  N2          = length(Values2)
  #DP          = matrix(0, nrow = N1, ncol = N2)
  #Initialize DP Matrix and set first row and column to infinity
  #DP[1, 2:N2] = Inf
  #DP[2:N1, 1] = Inf

  tmpArray = rep(Inf, N2+1)

  tmpVar = 0
  for(i in 2:N1){
    tmpVari = Values1[i]
    for(j in 2:N2){
      tmpVarj     = Values2[j]
      d1          = tmpVar + abs(tmpVari - tmpVarj)
      d2          = tmpArray[j]     + hlp_fun(tmpVari, Values1[i - 1], tmpVarj)
      d3          = tmpArray[j - 1] + hlp_fun(tmpVarj, tmpVari,        Values2[j - 1])
      tmpVar      = tmpArray[j]
      tmpArray[j] = min(d1, d2, d3)
    }
    tmpVar = Inf
  }
  return(tmpArray[N2])
}

hlp_fun = function(new_point, x, y){
  if(new_point < min(x, y) || new_point > max(x, y)) {
    return(min(abs(new_point - x), abs(new_point - y)));
  }else{
    return(0)
  }
}


