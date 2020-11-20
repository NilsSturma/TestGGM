# create matric containing all permutations of 1:n
permutations <- function(n){
  return(matrix(unlist(combinat::permn(1:n)), ncol=n, byrow=TRUE))
}

compute_S <- function(n,i,L){
  K = floor((n-1)/L)
  integer_set = setdiff(1:n, i)
  len = length(integer_set)
  len = len - (len %% L)
  integer_set = integer_set[1:len]
  S = matrix(integer_set, nrow=K, ncol=3, byrow=TRUE)
  return(S)
}