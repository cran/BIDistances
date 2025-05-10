DistanceDistributions=function(Data, DistanceMethods=c('bhjattacharyya', 'bray',
                                                       'canberra', 'chord',
                                                       'divergence', 'euclidean',
                                                       'minkowski', 'geodesic',
                                                       'hellinger', 'kullback',
                                                       'manhattan', 'maximum',
                                                       'soergel', 'wave',
                                                       'whittaker'),
                               CosineNonParallel = TRUE, CorrelationDist = TRUE,
                               Mahalanobis = FALSE, Podani = FALSE,
                               PlotIt = FALSE, PlotSampleSize = 5e3){

  #chebychev=maximum =L_inf norm
  #manhatten = L_1 norm
  #[Thrun, 2021]  Thrun, M. C.: The Exploitation of Distance Distributions for Clustering, International Journal of Computational Intelligence and Applications, Vol. 20(3), pp. 2150016, DOI: \doi{10.1142/S1469026821500164}, 2021.

  print('Calculating Distance Measures')
  requireNamespace('parallelDist')
  if(isTRUE(Mahalanobis)){
    DistanceMethods=c('mahalanobis',DistanceMethods)
  }
  if(isTRUE(Podani)){
    DistanceMethods=c('podani',DistanceMethods)
  }
  Matrix=sapply(X = DistanceMethods,FUN = function(method,x){
    print(paste('Computing',method))
    DD=matrix(NaN,nrow(x),nrow(x))
    try({
      if(method=='minkowski'){
        DD=as.matrix(parallelDist::parDist(x = x,method = 'minkowski',p=0.1))
      }else{
        DD=as.matrix(parallelDist::parDist(x = x,method = method))
      }
    })
    return(DD[upper.tri(DD,diag = F)])
  },Data)
  Matrix=as.matrix(Matrix)

  if(any(DistanceMethods=='minkowski')){
    indp=c(0.5,4,10,100)#0.1 haben wir oben
    print(paste('Computing minkowski',indp))
    ind=which(colnames(Matrix)=='minkowski')
    Minkowski=Matrix[,ind]
    Matrix=Matrix[,-ind]
    for(i in indp){
      DD=as.matrix(parallelDist::parDist(x = Data,method = 'minkowski',p=i))
      Minkowski=cbind(Minkowski,DD[upper.tri(DD,diag = F)])
    }

    names=paste0('minkowski_p_',c(0.1,indp))
    colnames(Minkowski)=names
    Matrix=cbind(Matrix,Minkowski)
  }


  if(CosineNonParallel){
    print(paste('Computing','cosine'))
    cosineM=BIDistances::DistanceMatrix(Data,'cosine')
    cosine=cosineM[upper.tri(cosineM,diag = F)]
    Matrix=cbind(Matrix,cosine)
  }
  if(CorrelationDist){
    print(paste('Computing','pearson'))
    pearsonM=1-cor(t(Data),method = 'pearson')
    pearson =pearsonM [upper.tri(pearsonM ,diag = F)]
    Matrix=cbind(Matrix,pearson)
  }
  #Liste=c(Liste,Liste[[6]]^2)

  #MDplot(Matrix,c(DistanceMethods),Scaling = 'CompleteRobust',SampleSize=1e+05)
  n=nrow(Matrix)
  ind=apply(Matrix,2,function(x) sum(is.finite(x)))
  Matrix=Matrix[,ind==n,drop=FALSE]

  NUniquepervar=apply(Matrix,MARGIN = 2,function(x) {
    x=x[is.finite(x)]
    return(length(unique(x)))
  })
  Matrix=Matrix[,NUniquepervar>1e2,drop=FALSE]

  requireNamespace('e1071')
  kurt=apply(Matrix,MARGIN = 2,function(x) {
    x=x[is.finite(x)]
    return(e1071::kurtosis(x))
  })
  NonRelevant=Matrix[,kurt>=10,drop=FALSE]
  Matrix=Matrix[,kurt<10,drop=FALSE]


  requireNamespace('diptest')
  diptestVal=apply(Matrix,2,function(x) {
    x=x[is.finite(x)]
    if(length(x)>0){
      y=suppressWarnings(diptest::dip.test(x))
      if(identical(y, numeric(0))) y=0
      return(suppressWarnings(y$p.value))
    }else{
      return(1)
    }
  })

  NonRelevant2=Matrix[,diptestVal>0.05,drop=FALSE]
  Matrix=Matrix[,diptestVal<0.05,drop=FALSE]
  if(dim(Matrix)[2]==0){
		Matrix=cbind(NonRelevant2,NonRelevant)
		print('No appropriate structurs in distance distributions found. Data is either density-based or other distances have to investigated.')
		obj=suppressWarnings(MDplot(Matrix,SampleSize = PlotSampleSize,Ordering="Columnwise",Scaling="Percentalize"))+ggtitle('Inappropriate Distances')
		return(list(DistanceMatrix=NULL,DistanceChoice='Not Found',OrderedDistances=Matrix,ggobject=obj))
  }
  #kurt=-DataVisualizations::RobustNormalization(kurt,Centered = F,Capped = T,pmin = 0.001,pmax = 0.999)
  #requireNamespace('modes')
  bimodalitycoef=apply(Matrix,2,function(x) {
    x=x[is.finite(x)]
    if(length(x)>0){
    y=suppressWarnings(DataVisualizations::BimodalityAmplitude(x))
      if(identical(y, numeric(0))) y=0
      return(y)
    }else{
      return(0)
    }
  })
  # iqr_vals=apply(Matrix,2,function(x) {
  #   x=x[is.finite(x)]
  #   if(length(x)>0){
  #     y=quantile(x,c(0.25,0.75),na.rm = T)
  #     iqr=y[2]-y[1]
  #     return(iqr)
  #   }else{
  #     return(0)
  #   }
  # })
  # iqr_vals=iqr_vals/max(iqr_vals,na.rm = T)#DatabionicSwarm::RobustNormalization(iqr_vals,Centered = F,Capped = T)
  # arr=(iqr_vals+3*bimodalitycoef)/4

  ind=order(bimodalitycoef,decreasing = T)
  Matrix=cbind(Matrix[,ind,drop=FALSE],NonRelevant2,NonRelevant)

  if(all(bimodalitycoef<0.3)){
    Matrix=cbind(Matrix,NonRelevant2,NonRelevant)
    print('No appropriate structurs in distance distributions found. Data is either density-based or other distances have to investigated.')
    obj=suppressWarnings(MDplot(Matrix,SampleSize = PlotSampleSize,Ordering="Columnwise",Scaling="Percentalize"))+ggtitle('Inappropriate Distances')
    return(list(DistanceMatrix=NULL,DistanceChoice='Not Found',OrderedDistances=Matrix,ggobject=obj))
  }

  obj=suppressWarnings(MDplot(Matrix,SampleSize = PlotSampleSize,Ordering="Columnwise",Scaling="Percentalize"))+ggtitle('Most-Relevant Distances (left) to Inappropriate Distances (right)')
  if(PlotIt){
    print('Plotting Distance Distributions')
    print(suppressWarnings(obj))
  }

  DistanceChoice=colnames(Matrix)[1]

  switch (DistanceChoice,
          cosine = {
            Distance=cosineM
          },
          pearson={
            Distance=pearsonM
          },
          minkowski_p_0.1={
            Distance=as.matrix(parallelDist::parDist(x = Data,method = 'minkowski',p=0.1))
          },
          minkowski_p_0.5={
            Distance=as.matrix(parallelDist::parDist(x = Data,method = 'minkowski',p=0.5))
          },
          minkowski_p_4={
            Distance=as.matrix(parallelDist::parDist(x = Data,method = 'minkowski',p=0.4))
          },
          minkowski_p_10={
            Distance=as.matrix(parallelDist::parDist(x = Data,method = 'minkowski',p=10))
          },
          minkowski_p_100={
            Distance=as.matrix(parallelDist::parDist(x = Data,method = 'minkowski',p=100))
          },
          {
            Distance=as.matrix(parallelDist::parDist(x = Data,method =DistanceChoice))
          }

  )
  return(list(DistanceMatrix=Distance,DistanceChoice=DistanceChoice,OrderedDistances=Matrix,ggobject=obj))
}
