#!/usr/bin/env python3
# Point.py
# This class represents a gazepoint

# system includes
import sys

class Point:
  def __init__(self, x, y, t, err):
    self.error = err
    self.timestamp = t

    self.coord = []
    self.coord.append(x)
    self.coord.append(y)
		
  # index into the coords of this Point
  def at(self, k):
    return self.coord[k]

  # set coords
  def set(self, x, y):
    self.coord[0] = x
    self.coord[1] = y

  def getStatus(self):
    return self.error
	
  # a gazepoint is valid if it's normalized coordinates are in the
  # range [0,1] and both eyes are present
  def valid(self):
    if self.error == "None" and \
       self.coord[0] > 0 and self.coord[1]  > 0 and \
       self.coord[0] < 1 and self.coord[1] < 1:
      return True
    else:
      return False
	
  def gettimestamp(self):
    return self.timestamp
