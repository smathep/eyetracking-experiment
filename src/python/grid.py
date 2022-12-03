#!/usr/bin/env python3
# Grid.py
# This class represents a grid (matrix) in which we'll deposit fixation counts

# system includes
import sys
from math import *
from numpy import *

from point import Point

class Grid:
  def __init__(self, m, n):
    self.A = zeros((m,n),int)
    # how many tiles in our grid
    self.xtiles = m
    self.ytiles = n
		
  # index into the matrix of this Grid
  def at(self, i, j):
    return self.A[i,j]

  def fill(self,points,w,h):
    for pt in points:
      ix = (pt.at(0) * float(w))//(int(w)//self.xtiles)
      iy = (pt.at(1) * float(h))//(int(h)//self.ytiles)
      self.A[ix,iy] += 1
