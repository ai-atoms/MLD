# eosfit.py fits E(V) data to a Birch-Murnaghan equation of state. 
# Current version: 2.0
#
# Copyright (C) 2012 Kurt Lejaeghere <Kurt.Lejaeghere@UGent.be>, Center for
# Molecular Modeling (CMM), Ghent University, Ghent, Belgium
#
# eosfit.py is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# eosfit.py is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for 
# more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with eosfit.py; if not, see <http://www.gnu.org/licenses/>.

# The following code is based on the source code of eos.py from the Atomic 
# Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>.

# Python and numpy are required to use this script.
import os
import numpy as np
from sys import argv

def fit_bulk_and_a0(cname):

 infile = open("%s_volume_energ"%cname,"r")
 #infile = argv[1]

 data = np.loadtxt(infile)
 volume, bulk_modulus, energy0, bulk_deriv, residuals = BM(data)

 echarge = 1.60217733e-19

 outstr = '''\
 Equation Of State parameters - least squares fit of a Birch Murnaghan curve

 %.5f \t %.5f \t %.3f

   V0 \t \t  B0 \t \t  BP
 [A^3/at] \t [GPa] \t \t [--] 

 1-R^2: %f

 ''' % (volume, (bulk_modulus * echarge * 1.0e21), bulk_deriv, residuals[0])

 outfile = open('data'+'.eosout', 'w')
 outfile.write(outstr)

 outfile.close()
 
 return (volume, bulk_modulus * echarge * 1.0e21, energy0)


def BM(energies):

    fitdata = np.polyfit(energies[:,0]**(-2./3.), energies[:,1], 3, full=True)
    ssr = fitdata[1]
    sst = np.sum((energies[:,1] - np.average(energies[:,1]))**2.)
    residuals0 = ssr/sst
    deriv0 = np.poly1d(fitdata[0])
    deriv1 = np.polyder(deriv0, 1)
    deriv2 = np.polyder(deriv1, 1)
    deriv3 = np.polyder(deriv2, 1)

    volume0 = 0
    x = 0
    for x in np.roots(deriv1):
        if x > 0 and deriv2(x) > 0:
            volume0 = x**(-3./2.)
            break

    if volume0 == 0:
        print('Error: No minimum could be found')
        print(('The min value is %22.7f'%min(energies[:,1])))
        print(('The max value is %22.7f'%max(energies[:,1])))
        exit()
    
    derivV2 = 4./9. * x**5. * deriv2(x)
    derivV3 = (-20./9. * x**(13./2.) * deriv2(x) -
        8./27. * x**(15./2.) * deriv3(x))
    bulk_modulus0 = derivV2 / x**(3./2.)
    bulk_deriv0 = -1 - x**(-3./2.) * derivV3 / derivV2


#Generate input script for GNUplot for fitting
    prefix="bulk"
    fi ="%s_input_eos.dat"%prefix
    finput=open(fi,"w")
    for i in range(len(energies[:,0])):
      finput.write(" %12.7f  %12.7f   %12.7f\n"%(energies[i,0],deriv0(energies[i,0]**(-2./3.)),energies[i,1]))
    finput.close()
    ndiscrete=1000
    minval=min(energies[:,0])
    maxval=max(energies[:,0])
    step=(maxval-minval)/float(ndiscrete-1)
    funct="%s_fit.dat"%prefix
    fout=open(funct,"w")

    for i in range(ndiscrete):
      #debug print i, minval + float(i)*step,deriv0(minval + float(i)*step),minval, step
      x=minval + step*float(i)
      fout.write('%12.7f %15.8f\n'%(x,deriv0(x**(-2./3.))))
    fout.close()
    os.system("rm -f fit.log")
    fn = "%s_fit.gp"%prefix
    fout = open(fn,"w")
    fout.write('set term gif\n')
    fout.write('set output "%s_fit.gif"\n'%prefix)
    fout.write("plot '%s_input_eos.dat' u 1:3 pt 5 title  'DATA', '%s_input_eos.dat'  u 1:2   title 'INTER' , '%s_fit.dat' u 1:2 with line title 'FUNCTION'"%(prefix, prefix,prefix))
    fout.close()
    os.system("gnuplot '%s_fit.gp' > gnuplot.report_bulk 2>&1"%prefix)

    energy0=deriv0(volume0**(-2.0/3.0))
    
    return volume0, bulk_modulus0, energy0, bulk_deriv0, residuals0

usage = '''\
Use: python eosfit.py filename
    calculates the Birch-Murnaghan equation of state from a given file, 
    containing in its columns the volumes in A^3/atom and energies in eV/atom,
    respectively
    output is printed in filename.eosout
--help gives an overview of all options
'''

#numpy.polyfit manual:
#numpy.polyfit(x, y, deg, rcond=None, full=False, w=None, cov=False)[source]
#Least squares polynomial fit.
#
#Fit a polynomial p(x) = p[0] * x**deg + ... + p[deg] of degree deg to points (x, y). 
#Returns a vector of coefficients p that minimises the squared error.
#
#Parameters:	
#x ....:  array_like, shape (M,), x-coordinates of the M sample points (x[i], y[i]).
#y ....:  array_like, shape (M,) or (M, K) y-coordinates of the sample points. Several 
#         data sets of sample points sharing the same x-coordinates can be fitted at once by passing 
#         in a 2D-array that contains one dataset per column.
#deg ..:  int, Degree of the fitting polynomial
#rcond :  float, optional Relative condition number of the fit. 
#         Singular values smaller than this relative to the largest singular value will be ignored. 
#         The default value is len(x)*eps, where eps is the relative precision of 
#         the float type, about 2e-16 in most cases.
#full .:  bool, optional
#         Switch determining nature of return value. 
#         When it is False (the default) just the coefficients are returned, 
#         when True diagnostic information from the singular value decomposition is also returned.
#w ... :  array_like, shape (M,), optional
#         weights to apply to the y-coordinates of the sample points.
#cov ..:  bool, optional
#         Return the estimate and the covariance matrix of the 
#         estimate If full is True, then cov is not returned.
#Returns.:	
#p ......: ndarray, shape (M,) or (M, K)
#          Polynomial coefficients, highest power first. 
#          If y was 2-D, the coefficients for k-th data set are in p[:,k].
#          residuals, rank, singular_values, rcond : :
#          Present only if full = True. Residuals of the least-squares fit, 
#          the effective rank of the scaled Vandermonde coefficient matrix, its singular values, and the 
#          specified value of rcond. For more details, see linalg.lstsq.
#V ......: ndarray, shape (M,M) or (M,M,K)
#          Present only if full = False and cov`=True. The covariance matrix of the polynomial coefficient estimates. 
#          The diagonal of this matrix are the variance estimates for each coefficient. 
#          If y is a 2-D array, then the covariance matrix for the `k-th data set are in V[:,:,k]
#Warns...:	
#RankWarning :   The rank of the coefficient matrix in the least-squares fit is deficient. 
#                The warning is only raised if full = False.
#                The warnings can be turned off by
#
#
