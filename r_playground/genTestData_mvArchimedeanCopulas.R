#**************************************************************************
#* 
#* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
#*
#* This program is free software: you can redistribute it and/or modify
#* it under the terms of the GNU General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#* GNU General Public License for more details.
#*
#* You should have received a copy of the GNU General Public License
#* along with this program.  If not, see <http://www.gnu.org/licenses/>.

## A script which produces test inputs/outputs for the Matlab implementation
## of Frank, Gumbel, and Clayton multivariate copula density functions

rm(list = ls())
cat("\014")

## hardcoded, but for now not sure how to handle this in R :/
setwd("/home/kiran/phd_code/copula/r_playground/")

## unload any existing loaded copula packages
#detach("package:copula", unload=TRUE)
library(copula, lib.loc="/home/kiran/R_sources/install")  ## load our personal copula library, where we
                                                          ## made small print-outs to understand better what was
                                                          ## going on with some of the code :)
library(utils)


## Generate data
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

##dCopula(c(0.1,0.1,0.1),gumbel.cop)

write.table(u, file = "testfiles/mvArchimedeanCopula_input.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(claytonpdf, file = "testfiles/claytonPdf3D_output.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(frankpdf, file = "testfiles/frankPdf3D_output.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(gumbelpdf, file = "testfiles/gumbelPdf3D_output.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)