

# Create noise

shift_unif = function(x,a,b){
  a + (b-a) * x
}

n = 4
a = -0.5
b = 0.5
A = matrix(shift_unif(runif(n**2), a, b), ncol=n)   # elements in (-1,1)
noise = t(A) %*% A # alwas pos. semidefinit
