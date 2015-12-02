library(copula)

beta.kernel.copula.surface = function (u,v,bx,by,p) {
  s = seq(1/p, len=(p-1), by=1/p)
  mat = matrix(0,nrow = p-1, ncol = p-1)
  for (i in 1:(p-1)) {
    a = s[i]
    for (j in 1:(p-1)) {
      b = s[j]
      cat("a=", a, " | b=", b, "\n")
      mat[i,j] = sum(dbeta(a,u/bx,(1-u)/bx) *
                       dbeta(b,v/by,(1-v)/by)) / length(u)
    } 
  }
return(data.matrix(mat)) }
