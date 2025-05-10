VariablePrecision= function(Variable){
  # [MinAbsDiff,MinAbsNZDiff,MinExpo] = VariablePrecision(Variable)
  # find he smallest (NonZero) Difference of any Value in Variable
  #
  # INPUT
  # Variable(1:n,1:d)    numeric data or
  #
  # OUTPUT
  # MinAbsDiff(1:d)        Absolute value of the smallest Difference of any Value in Variable
  # MinAbsNZDiff(1:d)      Absolute value of the smallest NonZero Difference of any Value in Variable
  # MinExpo(1:d)           The smalest Nachkommastelle in Variable NOT JET IMPLEMENTED
  
  zeros <-function (n,m=n,o=0) {
    # zeros(n)     returns an n-by-n matrix of 0s. 
    # zeros(n,1)   returns a vector of 0s 
    # zeros(n,m)   returns an n-by-m matrix of zeros
    # zeros(n,m,o) returns an 3D matrix  of zeros
    
    # ALU
    
    if (m==1) { # vector wird zurueckgegeben
      return(c(1:n)*0) ;   
    }else{      # return n-by-m matrix of ones.         
      if  (o==0){
        return(matrix(0,n,m));
      }else{   #3D matrix
        nullen = rep(0, m*n*o);  # soviele nullen wie in die 3D matrix pasen
        return(array(nullen,c(n,m,o)));
        
      } # end  if  (o==0)
    } # end if (m==1) 
  } # end  function  zeros
  
  n = dim(Variable)[1]
  d = dim(Variable)[2]
  
  if(is.null(d)){
    Diffs        = InnerVariableDifferences(Variable)
    MinAbsDiff   = min(abs(as.vector(Diffs)))
    MinAbsNZDiff = min(Diffs[Diffs > 0])
    MinExpo      = NaN
  }else{
    # d>1
    MinAbsDiff   = zeros(d, 1)
    MinAbsNZDiff = zeros(d, 1)
    MinExpo      = zeros(d, 1)
    
    for(i in 1:d) {
      Diffs           = InnerVariableDifferences(Variable[, i])
      MinAbsDiff[i]   = min(abs(as.vector(Diffs)))
      MinAbsNZDiff[i] = min(Diffs[Diffs > 0])
      MinExpo[i]      = NaN
    }# for i
  }# if d = 1
  
  return(list("MinAbsDiff"   = MinAbsDiff,
              "MinAbsNZDiff" = MinAbsNZDiff,
              "MinExpo"      = MinExpo))
}