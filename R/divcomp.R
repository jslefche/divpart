#' @title divcomp
#' 
#' @description 
#' \code{divcomp} Calculate diversity of order q
#' 
#' @author Jon Lefcheck
#' 
#' @param mat, a sample (rows) - by - species (columns) abundance matrix
#' @param dissim, a species-by-species dissimilarity matrix. If absent, uses taxonomic dissimilarities (0, 1)
#' @param q, order of diversity. 0 = species richness, 1 = Shannon, 2 = Simpson
#' 
#' @return Returns a vector of diversity values corresponding to each sample

divcomp = function(mat, dissim = NULL, q = 0) {
  
  # Relativize abundance matrix
  p = as.matrix(apply(mat, 1, function(x) x/sum(x, na.rm = TRUE)))
  
  p[, apply(p, 2, function(x) all(is.na(x)))] = 0
  
  # Create taxonomic dissimilarity matrix if dissim is empty
  if(is.null(dissim)) {
    
    dissim = matrix(1, nrow = nrow(p), ncol = nrow(p))
    
    diag(dissim) = 0
    
  } else {
    
    dissim = as.matrix(dissim)
    
  }
  
  if(ncol(p) != nrow(dissim)) stop("Number of rows must match between mat and groups")
  
  # Get similarity matrix from dissimilarity matrix
  Z = 1 - dissim

  # Calculate expected (average) similarity matrix and store in a matrix Zp
  Zp = matrix(0, nrow(p), ncol(p))
  
  for(k in 1:ncol(p)) { 
    
    for(i in 1:nrow(p)) { 
      
      for(j in 1:nrow(p)) { 
        
        Zp[i, k] = Zp[i, k] + Z[i, j] * p[j, k]
        
      } } } 

  # Compute diversity for each sample
  ret = rep(0, ncol(p))
  
  for(k in 1:length(ret)) {
    
    if(q == 1) {
      
      s = rep(0, nrow(Zp)) 
      
      for(i in 1:nrow(Zp)) {
        
        s[i] = Zp[i, k] ^ p[i, k]
        
      }
      
      ret[k] = 1/prod(s)
      
    } else {
      
      ret[k] = sum(p[, k] * (Zp[, k] ^ (q - 1)), na.rm = TRUE) ^ (1 / (1 - q))
      
    }
    
  }
  
  names(ret) = rownames(mat)
  
  # Return results
  return(ret)
  
}