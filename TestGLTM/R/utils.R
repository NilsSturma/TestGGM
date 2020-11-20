# create matric containing all permutations of 1:n
permutations <- function(n){
  return(matrix(unlist(combinat::permn(1:n)), ncol=n, byrow=TRUE))
}
