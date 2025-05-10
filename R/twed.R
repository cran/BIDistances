twed = function(Values1, Values2, Time1, Time2, Nu = 1, Lambda = 1, Degree = 2){
  # 
  # https://github.com/jzumer/pytwed/tree/master
  # 
  # V = twed(required variable setting)
  # V = twed(require + exemplary optional variable setting)
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
  # TWED                  TWED distance between time series Values1 (Time1) and
  #                       Values2 (Time2)
  # DPMatrix[1:n, 1:m]    Numeric matrix
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
  
  if(!is.numeric(Time1) | !is.vector(Time1)){
    stop("twed: Time1 must be a numeric vector.")
  }
  
  if(!is.numeric(Time2) | !is.vector(Time2)){
    stop("twed: Time2 must be a numeric vector.")
  }
  
  if(Nu < 0){
    stop("twed: Nu must be >= 0.")
  }
  
  if(Degree < 0){
    stop("twed: Degree must be in [0,Inf].")
  }
  
  N1          = length(Values1)
  N2          = length(Values2)
  DP          = matrix(0, nrow = N1, ncol = N2)
  # Initialize DP Matrix and set first row and column to infinity
  DP[1, 2:N2] = Inf
  DP[2:N1, 1] = Inf
  
  for(i in 2:N1){
    for(j in 2:N2){
      # Deletion in Values1
      del_a = DP[i-1, j  ] + Dlp(Values1[i - 1], Values1[i], p = Degree) + Nu * (Time1[i] - Time1[i - 1]) + Lambda
      # Deletion in Values2
      del_b = DP[i,   j-1] + Dlp(Values2[j - 1], Values2[j], p = Degree) + Nu * (Time2[j] - Time2[j - 1]) + Lambda
      # Keep data points in both time series
      MyMatch = DP[i - 1, j - 1] +
        Dlp(Values1[i], Values2[j], p = Degree) +
        Dlp(Values1[i - 1], Values2[j - 1], p = Degree) +
        Nu * (abs(Time1[i] - Time2[j]) + abs(Time1[i - 1] - Time2[j - 1]))
      # Choose the operation with the minimal cost and update DP Matrix
      DP[i, j] = min(del_a, del_b, MyMatch)
    }
  }
  
  distance = DP[N1, N2]
  
  return(list("TWED" = distance, "DPMatrix" = DP))
}

Dlp = function(A, B, p=2){
  cost = sum((abs(A-B))**p)**(1/p)
  return(cost)
}



