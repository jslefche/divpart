#' @title divpart
#'
#' @description
#' \code{divpart} Hierarchical diversity partitioning
#'
#' @author Jon Lefcheck
#'
#' @param mat, a sample (rows) - by - species (columns) abundance matrix
#' @param groups, a sample (rows) - by - groups (columns) matrix corresponding to the hierarhical sampling levels
#' @param dissim, a species-by-species dissimilarity matrix. If absent, uses taxonomic dissimilarities (0, 1)
#' @param q, order of diversity. 0 = species richness, 1 = Shannon, 2 = Simpson
#'
#' @return Returns a data.frame of alpha, beta, and gamma diversity values for each level of the hierarchy

divpart = function(mat, groups = NULL, dissim = NULL, q = 0) {

  if(nrow(mat) != nrow(groups)) stop("Number of rows must match between mat and groups")

  # Create groups at lowest level if groups are absent
  if(is.null(groups)) groups = matrix(1:nrow(mat), nrow = nrow(mat))

  # Create groups at lowest level if not already present in groups
  if(!any(sapply(groups, function(x) length(unique(x))) == nrow(mat))) {

    groups = cbind(groups, matrix(1:nrow(mat), nrow = nrow(mat)))

    colnames(groups)[ncol(groups)] = "ID"

  }

  # Organize groups from most levels to fewest
  groups. = sapply(colnames(groups), function(x) length(unique(groups[, x])))

  groups = groups[, names(rev(sort(groups.)))]

  # Split groups into a list
  groups.l = lapply(apply(groups, 2, list), unlist)

  # Split each grouping level and apply analysis
  do.call(rbind, lapply(length(groups.l):2, function(i) {

    # Calculate local diversity for each group in the level below the current group
    alpha = sapply(unique(groups.l[[i - 1]]), function(j) {

          # Subset data for group
          mat. = mat[j == groups.l[[i - 1]], , drop = FALSE]

          # Summarize at the group level
          if(nrow(mat.) > 1) mat. = t(as.matrix(colSums(mat.)))

          # Calculate local diversity
          divcomp(mat., dissim = dissim, q = q)

          } )

    # Get mean alpha diversity
    mean.alpha = sapply(unique(groups.l[[i]]), function(j) {

      mean(alpha[groups[groups[, i] == j, i - 1]], na.rm = TRUE)

    } )

    # Calculate gamma diversity across groups for the current level
    gamma = sapply(unique(groups.l[[i]]), function(j) {

      # Subset data for group
      mat. = mat[j == groups.l[[i]], ]

      # Summarize at the group level
      if(nrow(mat.) > 1) mat. = t(as.matrix(colSums(mat.)))

      # Calculate local diversity
      divcomp(mat., dissim = dissim, q = q)

      } )

    # Get additive beta diversity
    beta.add = gamma - mean.alpha

    # Get multiplicative beta diversity
    beta.mult = gamma / mean.alpha

    # Return in data.frame
    dat = data.frame(
      level = colnames(groups)[i],
      name = names(gamma),
      alpha = mean.alpha,
      gamma = gamma,
      beta.add = beta.add,
      beta.mult = beta.mult
    )

    rownames(dat) = NULL

    # Get plot-level mean alpha
    if(i == 2) {

      # Get plot-level diversity
      alpha. = mean(divcomp(mat, dissim = dissim, q = q), na.rm = T)

      gamma. = sum(colSums(mat) > 0)

      beta.add. = gamma. - alpha.

      beta.mult. = gamma. / alpha.

      # Add to dat
      dat = rbind(dat,
            data.frame(
              level = colnames(groups)[i - 1],
              name = paste0("1:", nrow(mat)),
              alpha = alpha.,
              gamma = NA, #gamma.,
              beta.add = NA, #beta.add.,
              beta.mult = NA #beta.mult.)
            )

    }

    return(dat)

  } ) )

}
