Tfidf_dist=function(FeatureMatrix_Gene2GoTerm,tf_fun=mean){
  # V <- Tfidf_dist(FeatureMatrix_Gene2GoTerm)
  # TfidfDistances <- V$dist
  # TfidfWeights   <- V$TfidfWeights
  #
  # Computes the term frequency inverse document frequency (tfidf) distance
  # for a FeatureMatrix_Gene2GoTerm.
  # In case of genes with annotated GOterms from gene ontology
  # genes can be interpreted as documents and GOterms as terms. 
  #
  # INPUT:
  # FeatureMatrix_Gene2GoTerm[1:n, 1:m] 	matrix of 1:n genes as rows and 1:m GoTerms as columns
  #
  #                                 FeatureMatrix_Gene2GoTerm[i,j] > 0 iff GOterm j 
  #																	is relevant for gene i. The FeatureMatrix_Gene2GoTerm[i,j] > 1 if the specific gene
  #                                 is annotated by in a specific GO-Term with more than one evidence code 
  #																	FeatureMatrix_Gene2GoTerm[i,j] is the frequency of
  #																	term js occurance in document i.
  # tf_fun                          Default is the mean of non zero entries as defined in [Stier & Thrun, 2023, Informatica] as the augmented frequency. The alternative definition would use sum()                     
  # OUTPUT:
  # TfidfDistances[1:n, 1:n]	      dist object of tdfidf distances between the documents = absolute difference of TfidfWeights
  # TfidfWeights[1:n]               the term frequence inverse document frequency weights used for the distance   
  #
  #Details
  # GO-Terms represent documents in the IR sense and the Description terms are the genes of a set. In [Stier & Thrun, 2023, Informatica] only manually curated genes were used.
  
  #Author: Michael Thrun
  
  #wie haeufig kommt term im document vor?
  ##d.h. wie haeufig kommt ein gene in GO-termen vor?
  fin = rowSums(FeatureMatrix_Gene2GoTerm)
  maxind = which.max(fin)
  #maximum occurance of any gene overall go terms it is annotated in
  maxfin = max(fin)
  
  #augmented frequency, to prevent a bias towards longer documents (GoTerms with genes that have a lot of evidence codes), e.g. 
  # raw frequency would be sume
  #augmented frequency of occurrence of a gene G in a given GO-Term T 
  augmented_freq=apply(FeatureMatrix_Gene2GoTerm,1,function(x) tf_fun(x[x!=0]))
  #Normalization by the maximal occurrence of that gene in any Go-Term:
  maxOcc=apply(FeatureMatrix_Gene2GoTerm,1,max)
  #Term frequency, tf is in general the relative frequency of term t (Gene) within document d (GoTerm), 
  tf = rep(NaN, length(fin))
  for (i in 1:dim(FeatureMatrix_Gene2GoTerm)[1]) {
    tf[i] =  augmented_freq[i]/maxOcc[i]
  }
  
  #Compute Inverse document frequency (IDF)
  idf = rep(NaN, length(fin))
  #idf logarithmically counts the number of GO terms in which any gene of the set is annotated
  #divided by be the number of GO terms to which the specific gene is annotated.
  for (i in 1:dim(FeatureMatrix_Gene2GoTerm)[1]) {
    idf[i] = log(1 + maxfin / fin[i])
  }
  tfidf=tf*idf

  Distance = matrix(NaN, length(idf), length(idf))
  for (i in 1:length(idf)) {
    for (j in 1:length(idf)) {
      Distance[i, j] = abs(tfidf[i] - tfidf[j])
    }
  }
  return(list(dist=as.dist(Distance),TfidfWeights=tfidf))
  #Term=term
  #gene=document
}