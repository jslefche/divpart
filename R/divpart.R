#' @title divpart
#' 
#' @description 
#' \code{divcomp} Hierarchical diversity partitioning
#' 
#' @author Jon Lefcheck
#' 
#' @param mat, a sample (rows) - by - species (columns) abundance matrix
#' @param groups, a sample (rows) - by - groups (columns) matrix corresponding to the hierarhical sampling levels
#' @param dissim, a species-by-species dissimilarity matrix. If absent, uses taxonomic dissimilarities (0, 1)
#' @param q, order of diversity. 0 = species richness, 1 = Shannon, 2 = Simpson
#' 
#' @return Returns a data.frame of alpha, beta, and gamma diversity values for each level of the hierarchy

divpart = function(mat, groups, dissim = NULL, q = 0) {
  
  if(nrow(mat) != nrow(groups)) stop("Number of rows must match between mat and groups")
  
  # Organize groups from most levels to fewest
  groups. = sapply(colnames(groups), function(x) length(unique(groups[, x])))
  
  groups = groups[, names(rev(sort(groups.)))]
  
  # Split groups into a list
  groups.l = lapply(apply(groups, 2, list), unlist)
  
  # Split each grouping level and apply analysis
  do.call(rbind, lapply(length(groups.l):2, function(i) {
      
    # Calculate local diversity for each group in the level below the current group
    alpha = if(i - 1 == 1) {
        
        divcomp(mat, dissim = dissim, q = q)
        
      } else {
        
        sapply(unique(groups.l[[i - 1]]), function(j) {
      
          # Subset data for group
          mat. = mat[j == groups.l[[i - 1]], ]
          
          # Summarize at the group level
          mat. = t(as.matrix(colSums(mat.)))
            
          # Calculate local diversity
          divcomp(mat., dissim = dissim, q = q) 
          
        } )
        
        }
      
    # Get mean alpha diversity
    mean.alpha = mean(alpha, na.rm = TRUE)
      
    # Calculate gamma diversity across groups for the current level
    gamma = if(i - 1 == 1) {
      
      mat. =  t(as.matrix(colSums(mat)))
      
      divcomp(mat., dissim = dissim, q = q)
      
    } else {
      
      sapply(unique(groups.l[[i]]), function(j) {
      
      # Subset data for group
      mat. = mat[j == groups.l[[i]], ]
      
      # Summarize at the group level
      mat. = t(as.matrix(colSums(mat.)))
      
      # Calculate local diversity
      divcomp(mat., dissim = dissim, q = q) 
      
      } )
      
    }
    
    # Get additive beta diversity
    beta = gamma - mean.alpha
    
    # Return in data.frame
    dat = data.frame(
      level = i,
      name = colnames(groups)[i - 1],
      alpha = mean.alpha,
      gamma = gamma,
      beta = beta
    )
    
    rownames(dat) = NULL
    
    return(dat)
    
  } ) )
  
}