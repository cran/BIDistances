InnerVariableDifferences = function (Variable){
  # all inner differences di = xi-yi for all xi,yi in Variable
  #
  # INPUT
  # Variable(1:n)    matrix of size nxd without NaNs in columns where defined = 1
  #
  # OUTPUT
  # Di     matrix of all inner differences
  # mt
  # zurueckrechnen auf unquadrierte Differenzen

  N  = length(Variable)
  DM = InnerVariableDifferencesRcpp(Variable, N)

  return(DM)
}
