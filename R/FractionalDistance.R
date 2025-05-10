FractionalDistance=function(Data,p){
  #file:///C:/Subversion/PRO/Research/WissenAusDaten2014/70DistanceDistributions/06OtherPublications/2001Aggarwal_On_the_Surprising_Behavior_of_Distance_Metric_in_High_Dimensional_Space.pdf
  n=nrow(Data)
  Distance=matrix(0,n,n)
  dimnames(Distance) <- list(rownames(Data), rownames(Data))

  for(i in 1:n){
    for(j in 1:n){
      if(i>j){
        Distance[i,j]=sum((abs(Data[i,]-Data[j,]))^p)^(1/p)
      }
    }
  }

  #Zu tun damit aus unteren diagonal matrix eine komplette matrix wird
  return(as.matrix(as.dist(Distance)))
}
