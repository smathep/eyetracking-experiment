#!/usr/bin/env python3

import os, sys, math
import numpy as np

from scipy import signal

from bw import Butterworth
from point import Point

import bw

def main(argv):

  herz = 30
  sfdegree = 6
  sfcutoff = 1.15

  bwx = Butterworth(sfdegree,herz,sfcutoff)
  bwx.dump()

  n = 30

  x = []
  data = []
  sdata = []
  for i in range(n):
    data.append(Point(10,0,0,'None'))
    x.append(10)
    print("Point: ", data[i].at(0))
    print("x: ", x[i])

  sdata = bw.applyBWFilter(data,sfdegree,herz,sfcutoff)
  for i in range(n):
    print("Point: ", sdata[i].at(0))

  b, a = signal.butter(sfdegree,1.0/herz,'low',False,'ba')
  print("b, a = ", b, a)

  y = signal.filtfilt(b, a, x, padlen=n-1)
  for i in range(n):
    print("y: ", y[i])

if __name__ == "__main__":
  main(sys.argv[1:])
