##Solution

## 1
cov.mc = function(n){
  if(n < 1) stop('n should be at least 1')
  u = runif(n)
  e = exp(u)
  UE = mean(u*e)
  U = mean(u)
  E = mean(e)
  UE - U*E
}
## 2
theta = 3/2 - exp(1)/2
set.seed(3)
theta.hat = cov.mc(5000)
abs(theta - theta.hat)

## 3
sampleN = function(k){
  N = rep(NA, k)
  for(j in 1:k){
    n = 0 # starting value for N. Will be updated
    v = 1 # auxiliary variable used to multiply
    while(TRUE){
      u = runif(1)
      f = prod(v,u)
      if(f >= exp(-3)){
        v = f
        n = n + 1
      }else break
    }
    N[j] = n
  }
  N
}
## 4
set.seed(3)
Nsamp = sampleN(20000)
mean(Nsamp == 3)
