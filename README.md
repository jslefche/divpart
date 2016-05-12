# divpart: Hierarchical diversity partitioning using effective numbers in R

### 
This package calculates alpha, beta, and gamma diversity for various indices of diversity (species richness, Shannon entropy, Simpson diversity) at successive levels of a hierarchical sampling scheme. 

Diversity indices can incorporate species dissimilarities, such as those from functional traits or phylogenetic trees, and thus can also provide an estimate of functional and phylogenetic turnover across spatial scales.

### Examples
```
# Install package
library(devtools)
install_github("jslefche/divpart")

library(divpart)

# Create sample-by-species abundance matrix
set.seed(9) 

mat = matrix(rpois(100, 1), nrow = 20, ncol = 10)

rownames(mat) = paste0("plot", 1:20)
colnames(mat) = paste0("species", 1:10)

# Create matrix of grouping variables
groups = data.frame(plot = 1:nrow(mat), subsite = rep(letters[1:4], 5), site = rep(letters[1:2], each = 10), region = "A")

# Run diversity partitioning
divpart(mat, groups, q = 1)

```
