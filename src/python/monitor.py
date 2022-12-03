#!/usr/bin/env python3
# Lagrange.py

# system includes
import sys
import os
import math
import numpy as np
from scipy import spatial

# local includes
#from point import Point
#from fixation import Fixation
#from grid import Grid
#from aoi import AOI

class Monitor:

  #        n: number of target points
  #        S: the actual (known) locations of the n target points
  def __init__(self,w=1920,h=1200,screen=22,dist=21.65,n=5):
    self.n = n
    self.S = np.matrix(np.zeros((self.n,2),dtype=float))

    self.w = float(w)
    self.h = float(h)

    r = math.sqrt(self.w*self.w + self.h*self.h)
    self.dpi = r/float(screen)

    self.D = float(dist)
#   fov = 2*math.degrees(math.atan(math.radians(screen/(2*self.D))))
    fov = 2*math.degrees(math.atan2(screen,2*self.D))
    fovx = 2*math.degrees(math.atan2(float(w)/self.dpi,2*self.D))
    fovy = 2*math.degrees(math.atan2(float(h)/self.dpi,2*self.D))


#   print "screen subtends %f (%f x %f) degrees" % \
#                   (float(fov),float(fovx),float(fovy))
#   print "screen aspect = %f" % float(float(w)/float(h))
#   print "sampling period = %f (ms)" % float(period * 1000.0)
#   print "filter window = %f (ms)" % float(dt * 1000.0)
#   print "self.dpi = %f" % self.dpi

    # get pixels for degree visual angle (10 deg offser for targets)
    degs = 10.0
#   dpix = 2*self.D*math.tan(math.radians(degs/2.0)) * self.dpi
    dpix = self.deg2pix(degs)
    dx = dpix/float(w)
    dy = dpix/float(h)

    # in (0,0) top-left normalized coordinates

    self.S[0,0] = 0.5		# center
    self.S[0,1] = 0.5

    self.S[1,0] = 0.5 + dx	# bottom_right
    self.S[1,1] = 0.5 + dy

    self.S[2,0] = 0.5 - dx	# bottom_left
    self.S[2,1] = 0.5 + dy

    self.S[3,0] = 0.5 - dx	# top_left
    self.S[3,1] = 0.5 - dy
 
    self.S[4,0] = 0.5 + dx	# top_right
    self.S[4,1] = 0.5 - dy

    self.ndict = {'center':0,\
                  'bottom_right':1,\
                  'bottom_left':2,\
                  'top_left':3,\
                  'top_right':4}

  def deg2pix(self,degs):
    dpix = 2*self.D*math.tan(math.radians(degs/2.0)) * self.dpi
    return dpix

  def pix2deg(self,dpix):
    r = math.sqrt(self.w*self.w + self.h*self.h)
    # pixel distances coming in are normalized
    dpix = dpix * r
    degs = 2*math.degrees(math.atan2(dpix/self.dpi,2*self.D))
    return degs

  def map2ctr(self,x,y):
    # map [0,1] -> [-1,1]
    x = x * 2.0 - 1.0
    y = y * 2.0 - 1.0
    return x,y
