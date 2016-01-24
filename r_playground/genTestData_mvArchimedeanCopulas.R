## A script which produces test inputs/outputs for the Matlab implementation
## of Frank, Gumbel, and Clayton multivariate copula density functions

rm(list = ls())
cat("\014")

## hardcoded, but for now not sure how to handle this in R :/
setwd("/home/kiran/phd_code/copula/r_playground/")

## Generate Clayton Copula data
library(copula)
library(utils)

alpha <- 4;
clayton.cop <- claytonCopula(alpha, dim=3)
frank.cop <- frankCopula(alpha, dim=3)
gumbel.cop <- gumbelCopula(alpha, dim=3)

## generate input to unit hypercube
uu <- seq(0.1,0.99,0.1)
u <- expand.grid(uu,uu,uu)

sz <- nrow(u)

claytonpdf <- numeric(sz)
frankpdf   <- numeric(sz)
gumbelpdf  <- numeric(sz)

for(ii in 1:sz)
{
  claytonpdf[ii] <- dCopula(unlist(u[ii,]),clayton.cop)
  frankpdf[ii] <- dCopula(unlist(u[ii,]),frank.cop)
  gumbelpdf[ii] <- dCopula(unlist(u[ii,]),gumbel.cop)
}

write.table(u, file = "testfiles/mvArchimedeanCopula_input.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(claytonpdf, file = "testfiles/claytonPdf3D_output.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(frankpdf, file = "testfiles/frankPdf3D_output.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(gumbelpdf, file = "testfiles/gumbelPdf3D_output.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)