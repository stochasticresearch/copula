#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
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
#**************************************************************************

import math
import numpy as np

import random

"""
continuingRV.py shows the effect of "continuing" a random variable, as defined by
Denuit and Lambert in the paper "Constraints on Concordance Measures in bivariate
Discrete Data," published in the Journal of Multivariate Analysis (93) 2005.  
"""

def random_discrete_rv(x, cdf_domain):
    """
    Generates a random discrete distribution and returns the PDF as a list of tuple's
    
    Inputs:
    x -- the domain of the discrete random variable to be generated
    cdf_domain   -- the desired points over which to calculate the CDF
    
    Outputs:
    y -- A dictionary that looks as follows:
         {'pdf': [(x_1,y_1),(x_2,y_2), ... (x_n,y_n)],
          'cdf': [(x_1,y_1),(x_2,y_2), ... (x_n,y_n)]}
    """
    
    # calculate the PDF
    pdf = []
    total = 0
    for ii in range(len(x)):
        val = round(random.uniform(0,1.0/len(x)),2)
        total = total + val
        # ensure the total probability sums to 1
        if(ii+1==len(x)):
            val = 1-(total-val)
        
        tup = (x[ii],val)
        pdf.append(tup)
    unzipped = zip(*pdf)
    
    # calculate the CDF over the domain of interest
    # the formula for a discrete CDF is given by: F(x) = P(X<=x) = sum(P(X=x_i); x_i<x)
    cdf = []
    for ii in range(len(cdf_domain)):
        xx = cdf_domain[ii]
        domain_vals = filter(lambda f: f<=xx, x)
        cdf_val = 0
        for jj in domain_vals:
            cdf_val = cdf_val + unzipped[1][unzipped[0].index(jj)]
        cdf.append((xx,cdf_val))
        
    y = {}
    y['pdf'] = pdf
    y['cdf'] = cdf
    
    return y

def continue_rv(discrete_rv, domain):
    """
    Returns the PDF and CDF of the continued random variable for visualization purposes.
    
    Assumptions:
      1.) Currently, we continue the RV w/ a uniform RV for simplicity, later we can expand
          on this.
      2.) The domain requested (input) must be the same as teh domain over which the provided
          CDF must be defined
    
    Inputs:
    discrete_pdf -- the discrete random variable inputted in the same output format as 
                    the function random_discrete_rv
    domain       -- the desired points over which to calculate the PDF and CDF
    
    Outputs:
    y -- A dictionary that looks as follows:
         {'pdf': [(x_1,y_1),(x_2,y_2), ... (x_n,y_n)],
          'cdf': [(x_1,y_1),(x_2,y_2), ... (x_n,y_n)]}
    """
    
    discrete_pdf_unzipped = zip(*discrete_rv['pdf'])
    discrete_cdf_unzipped = zip(*discrete_rv['cdf'])
    
    # calculate the continued pdf & cdf
      # f*(s) = f([s]+1)   ** See Denuit and Lambert for more details
    pdf = []
    for ii in range(len(domain)):
        s = domain[ii]
        s_bracket = int(s)
        try:
            idx = discrete_pdf_unzipped[0].index(s_bracket+1)
        except ValueError:
            if(s_bracket > len(discrete_pdf_unzipped[0])-1):
                idx = len(discrete_pdf_unzipped[0])-1
            else:
                idx = 0
        yy_pdf = discrete_pdf_unzipped[1][idx]
        
        pdf.append((s,yy_pdf))
    
    # calculate the continued cdf
      # F*(s) = F([s]) + (s-[s])*f([s]+1)   ** See Denuit and Lambert for more details
    cdf = []
    # TODO: we can calculate the CDF in the same loop as the pdf calculation
    for ii in range(len(domain)):
        s = domain[ii]
        s_bracket = int(s)
        
        try:
            idx = discrete_pdf_unzipped[0].index(s_bracket+1)
        except ValueError:
            if(s_bracket > len(discrete_pdf_unzipped[0])-1):
                idx = len(discrete_pdf_unzipped[0])-1
            else:
                idx = 0
    
        # add the (s-[s])*f([s]+1) part
        yy = (s-s_bracket)*(discrete_pdf_unzipped[1][idx])
        
        # add the F([s]) part
        try:
            idx = discrete_cdf_unzipped[0].index(s_bracket)
        except ValueError:
            if(s_bracket > len(discrete_cdf_unzipped[0])-1):
                idx = len(discrete_cdf_unzipped[0])-1
            else:
                idx = 0
        yy = yy + (discrete_cdf_unzipped[1][idx])
        
        cdf.append((s,yy))
    
    y = {}
    y['pdf'] = pdf
    y['cdf'] = cdf
    
    return y

if __name__=='__main__':
    x = [1,2,3,4,5]
    step = 0.1
    cdf_domain = np.arange(0,5+step, step)
    y = random_discrete_rv(x, cdf_domain)
    y_unzipped_cdf = zip(*y['cdf'])
    y_unzipped_pdf = zip(*y['pdf'])
    cy = continue_rv(y, cdf_domain)
    cy_unzipped_cdf = zip(*cy['cdf'])
    cy_unzipped_pdf = zip(*cy['pdf'])
    
    # plot for visualization
    import matplotlib.pyplot as plt
    plt.subplot(2,1,1)
    plt.plot(y_unzipped_cdf[0], y_unzipped_cdf[1], label='F(X)')
    plt.plot(cy_unzipped_cdf[0], cy_unzipped_cdf[1], label='F(X*)')
    plt.title('CDF\'s')
    plt.legend(loc=2)
    plt.grid(True)
    plt.subplot(2,1,2)
    plt.stem(y_unzipped_pdf[0], y_unzipped_pdf[1], linefmt='b-', markerfmt='ro', label='f(x)')
    plt.stem(cy_unzipped_pdf[0], cy_unzipped_pdf[1], linefmt='g-', markerfmt='bo', label='f(x*)')
    plt.legend(loc=2)
    plt.grid(True)
    plt.title('PDF\'s')
    
    
    plt.show()