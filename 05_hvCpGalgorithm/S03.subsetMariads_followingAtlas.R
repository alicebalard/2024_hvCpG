
## Step 1: Extract column counts from Atla data
## (done in S02 as it needs loads of ram and time!)
n_matrices <- length(my_list_mat)
ncols_list1 <- sapply(my_list_mat, ncol)

## Results:
n_matrices <- 
ncols_list1 <-

set.seed(123)  # for reproducibility

# Step 2: Randomly select matrices from Maria's datasets
selected_indices <- sample(seq_along(list2), n_matrices)

# Step 3: For each selected matrix, randomly subset columns to match list1
list3 <- mapply(function(mat2, target_ncol) {
  actual_ncol <- ncol(mat2)
  if (actual_ncol < target_ncol) {
    stop("Matrix in list2 has fewer columns than required.")
  }
  cols <- sample(actual_ncol, target_ncol)
  mat2[, cols, drop = FALSE]
}, list2[selected_indices], ncols_list1, SIMPLIFY = FALSE)
