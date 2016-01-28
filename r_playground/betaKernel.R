#******************************************************************************
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
