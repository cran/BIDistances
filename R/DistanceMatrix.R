DistanceMatrix = function(X,method='euclidean',dim=2,outputisvector=FALSE){
  # Dmatrix = DistanceMatrix(X,distance)
  # computes the distance between objects in the data matrix, X,
  # using the method specified by distance, which can be any of the following character strings
  #
  # INPUT
  # X[d,n]         Daten bestehend aus d Datensaetzen/Werten/Zeilen von n Vektoren/Variablen/Spalten ohne NaN
  #                Distanz wird jeweils zwischen zwei Zeilen berechnet
  #
  # Optional
  # method          method specified by distance string:
  #                 'euclidean','sqEuclidean','mahalanobis','pearson','spearman','kendall','cityblock=manhatten','cosine','chebychev'=max(abs(x-y)),'jaccard','minkowski','manhattan','binary', 'canberra'=sum abs(x-y)/sum(abs(x)-abs(y)), 'maximum', 'braycur'=sum abs(x -y)/abs(x+y), "cosine"
  # dim             if method="minkowski", or wasserstein, choose scalar, The ppth root of the sum of the ppth powers of the differences of the components. For wasserstein default is 1
  # outputisvector  Falls Vector der Distanzen benoetigt wird, wird Paket pracma fuer squareform geladen
  #
  # OUTPUT
  # Dmatrix       Distance-Matrix: Pairwise distance between pairs of objects oder Vektor(outputisvector=TRUE)

  # Autor: MT
  # 1. Editor:
  #
  # EXAMPLE
  #
  # Nota
  # needs package stats, vegan, sphet, biotools
  if (!is.matrix(X)) {
    X = as.matrix(X)
    warning('Data was not in matrix datatype!')
  }
  if (any(is.nan(X), na.rm = TRUE))
    stop('NaNs in Data found!')

  m = nrow(X) #Zeilenanzahl
  n = ncol(X) #Spaltenanzahl
  if (m < n) {
    #warning('Pruefe Anzahl von Spalten und Zeilen der Input-Matrix!')
    #warning('Vektorraum unterbestimmt: Definition eines gueltigen Distanzmasses zweifelhaft, weil jede Variable, also jeder Vektor, zu Wenig Werte also Zeilen besitzt!')
  }
  method=tolower(method)
  #the same
  if(method=="cityblock") method="manhatten"

  switch(method,
         euclidean = {
             Dmatrix = as.matrix(parallelDist::parallelDist(X))
           },
         sqeuclidean = {
             Dmatrix = as.matrix(parallelDist::parallelDist(X))
            Dmatrix = (Dmatrix) ^ 2
           },
         binary = {
           Dmatrix = as.matrix(dist(X, method = "binary"))
           },
         braycur = {
           #requireRpackage('sphet')
           if(!requireNamespace("sphet")){
             Dmatrix = sphet::distance(coord = X, measure = "braycur",
                                       output = FALSE, type = "distance")
             }
           },
         # cityblock = {
         # Dmatrix = as.matrix(dist(X, method = "manhattan"))
   # },
   wasserstein = {
     #Dmatrix = as.matrix(arules::dissimilarity(X, method = "pearson"))
     Dmatrix=WassersteinDist(X,p=dim,InverseWeighting = FALSE)
   },
   pearsond = {
     #Dmatrix = as.matrix(arules::dissimilarity(X, method = "pearson"))
     Dmatrix=1-cor(t(X),method = 'pearson')
   },
   spearmand = {
     Dmatrix=1-cor(t(X),method = 'spearman')
   },
   kendalld = {
     Dmatrix=1-cor(t(X),method = 'kendall')
   },
    pearsonm = {
      Dmatrix=TransformSimilarity2MetricDistance((cor(t(X),method = 'pearson')+1)/2)#1-cor(t(X),method = 'pearson')
    },
   spearmanm = {
     Dmatrix=TransformSimilarity2MetricDistance((cor(t(X),method = 'spearman')+1)/2)#1-cor(t(X),method = 'spearman')
   },
   kendallm = {
     Dmatrix=TransformSimilarity2MetricDistance((cor(t(X),method = 'kendall')+1)/2)##1-cor(t(X),method = 'kendall')
   },
   # maximum = {
     # Dmatrix = as.matrix(dist(X, method = "maximum"))
   # },
   # canberra = {
     # Dmatrix = as.matrix(dist(X, method = "canberra"))
   # },
   fractional = {
    Dmatrix = FractionalDistance(X, dim)
   },
   fagerdissimilarity = {
     if(!requireNamespace("parallelDist")){
       DD = as.matrix(parallelDist::parDist(X, method = 'fager'))
       Dmatrix=-DD
     }
  },
    cosine = {
      #requireRpackage("lsa")
	    #requireNamespace('lsa')
      #sim = lsa::cosine(X)
      #Dmatrix = max(sim, na.rm = T) - sim
      Dmatrix=CosinusDistance(X)
      #Dmatrix=cosine(X)
    },
    chebychev = {
      #requireRpackage('sphet')
      #Dmatrix=chebyshev(X)
      if(!requireNamespace("sphet")){
        Dmatrix = sphet::distance(coord = X, measure = "chebyshev",
                                  output = FALSE, type = "distance")
      }
    },
   # jaccard = {
    #  Dmatrix = jaccard(X)
    #},
    #mahalanobis = {
    #  Sx = cov(X)
     # Dmatrix = Mahalanobis(X, Sx)
    #},
    minkowski = {
      if(!requireNamespace("parallelDist")){
        Dmatrix = as.matrix(parallelDist::parDist(X, method = method,p=dim))
      }
      #Dmatrix = as.matrix(dist(X, method = "minkowski", p = dim))
    },
   # manhattan = {
   #   Dmatrix = as.matrix(dist(X, method = "manhattan"))
   # },
    {#all other not mentioned above
      #if(!requireNamespace("parallelDist")){
        Dmatrix = as.matrix(parallelDist::parDist(X, method = method))
      #}
      #return("Error! Bitte den String der Methodenauswahl ueberpruefen")
    }
  )
  if(!outputisvector) {
    return(Dmatrix)
  } else{
    #if(!require(pracma)){install.packages("pracma")}else{library(pracma)}  intern implementiert
    #return(squareform(Dmatrix))}
    return(Dmatrix[upper.tri(Dmatrix)])
  }
}
