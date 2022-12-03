#!/usr/bin/env python3
#
#  butterworth.py
#  
# From: http://www.exstrom.com/journal/sigproc/:
# 
# Hollos, Stefan and Hollos, J. Richard, "Recursive Digital Filters:
#   A Concise Guise", Exstrom Laboratories, LLC, Longmont, CO, USA
#   April, 2014
#   ISBN 9781887187244 (ebook)
#   URL: http://www.abrazol.com/books/filter1/

import os, sys, math
import numpy as np

from point import Point

class Butterworth:
  def __init__(self,n,fs,fc):

    self.n  = n/2
    a  = math.tan(math.pi * fc / fs)
    a2 = a * a

    self.A  = [0.0]*self.n
    self.d1 = [0.0]*self.n
    self.d2 = [0.0]*self.n
    self.w0 = [0.0]*self.n
    self.w1 = [0.0]*self.n
    self.w2 = [0.0]*self.n

    r  = 0.0

    for i in range(self.n):
      r = math.sin(math.pi * (2.0 * float(i) + 1.0)/(4.0 * float(self.n)))
      fs = a2 + (2.0 * a * r) + 1.0
      self.A[i] = a2 / fs
      self.d1[i] = 2.0 * (1.0 - a2)/fs
      self.d2[i] = -(a2 - (2.0 * a * r) + 1.0)/fs

  def bwlpf(self,x):

    """ input x is just a value
        bwlpf is meant to process a stream of x values
        it updates w0[], w1[], w2[], the filter output memory, as it sees
        incoming values of x
        for this reason, each stream of data must use its own instance
        of the Butterworth filter
    """
    for i in range(self.n):
      self.w0[i] = self.d1[i] * self.w1[i] + self.d2[i] * self.w2[i] + x
      x = self.A[i] * (self.w0[i] + 2.0 * self.w1[i] + self.w2[i])
      self.w2[i] = self.w1[i]
      self.w1[i] = self.w0[i]

    return x

  def dump(self):
    print("A: ", self.A)
    print("d1: ", self.d1)
    print("d2: ", self.d2)
    print("w0: ", self.w0)
    print("w1: ", self.w1)
    print("w2: ", self.w2)

def applyBWFilter(gazepoints,degree,herz,cutoff):

  bwx = Butterworth(degree,herz,cutoff)
  bwy = Butterworth(degree,herz,cutoff)

  points = []

  for point in gazepoints:
    x = bwx.bwlpf(point.at(0))
    y = bwy.bwlpf(point.at(1))
    t = point.gettimestamp()
    err = point.getStatus()
    points.append(Point(x, y, t, err))

  return points
