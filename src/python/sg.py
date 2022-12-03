from math import *
from numpy import *
import random
import collections

from point import Point

class SG:

  def __init__(self, num_points, pol_degree, diff_order=0):

    """ calculates filter coefficients for symmetric savitzky-golay filter.
        see: http://www.nrbook.com/a/bookcpdf/c14-8.pdf
  
        num_points   means that 2*num_points+1 values contribute to the
                     smoother.
  
        pol_degree   is degree of fitting polynomial
  
        diff_order   is degree of implicit differentiation.
                     0 means that filter results in smoothing of function
                     1 means that filter results in smoothing the first 
                                                 derivative of function.
                     and so on ...
  
    """
  
    # setup interpolation matrix
    # ... you might use other interpolation points
    # and maybe other functions than monomials ....
  
    x = arange(-num_points, num_points+1, dtype=int)
    monom = lambda x, deg : pow(x, deg)

    A = zeros((2*num_points+1, pol_degree+1), dtype=float)
    for i in range(2*num_points+1):
      for j in range(pol_degree+1):
        A[i,j] = monom(x[i], j)
          
    # calculate diff_order-th row of inv(A^T A)
    ATA = dot(A.transpose(), A)
    rhs = zeros((pol_degree+1,), float)
    rhs[diff_order] = (-1)**diff_order
    wvec = linalg.solve(ATA, rhs)
  
    # calculate filter-coefficients
    self.coeff = dot(A, wvec)

#   print "Coefficients: ", self.coeff
    self.n = len(self.coeff)
    # empty deque on init
#   self.rb = collections.deque([],self.n)
    # full deque on init (with zeros)
    self.rb = collections.deque([0.0]*self.n,self.n)

  def sgf(self,x):

    self.rb.appendleft(x)

    N = (size(self.coeff)-1)//2
    res = convolve(self.rb, self.coeff)
    res = res[N:-N]

#   returns the whole buffer
#   return res
#   returns just the central value
    return res[(len(self.coeff)-1)//2]

  def smooth(self,signal):
    
    """ applies coefficients calculated to signal """
    
    N = (size(self.coeff)-1)//2
#   print "N = ",N
#   print "signal = ",signal
#   print "len(signal) = ",len(signal)
    res = convolve(signal, self.coeff) if len(signal) > 2 else [0.0]*len(signal)

    return res[N:-N]

def applySGFilter(gazepoints,window,degree,difforder):

  sgx = SG(window,degree,difforder)
  sgy = SG(window,degree,difforder)

  y = []
  x = []
  t = []
  err = []
  points = []

  for point in gazepoints:
    x.append(point.at(0))
    y.append(point.at(1))
    t.append(point.gettimestamp())
    err.append(point.getStatus())

# print "len(x) = ",len(x)
# print "len(y) = ",len(y)
# print "len(t) = ",len(t)

  dxdt = sgx.smooth(x)
  dydt = sgy.smooth(y)

# print "len(dxdt) = ",len(dxdt)
# print "len(dydt) = ",len(dydt)

  for i in range(len(dydt)):
    points.append(Point(dxdt[i],dydt[i],t[i],err[i]))

  return points
