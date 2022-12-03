#!/usr/bin/env python3
# Fixation.py
# this file contains the class definition
# for a Fixation object

# system includes
import sys

# local includes
from point import Point

# fixation subclasses point
class Fixation(Point):

  # constructor
  def __init__(self, x, y, t, d):
    Point.__init__(self,x,y,t,"None")
    self.number = 0
    self.duration = d
    self.percentDuration = 0.

  # getter/setter methods below this line
  def setNumber(self, n):
    self.number = n

  def setDuration(self, d):
    self.duration = d

  def getNumber(self):
    return self.number

  def getDuration(self):
    return self.duration

  def getPercentDuration(self):
    return self.percentDuration

  def normalizeDuration(self,mindur,maxdur):
    range = maxdur - mindur
    if(abs(range) < 0.00001):
      if(maxdur < 0.00001):
        # assume 0 duration
        self.percentDuration = 0.0
      else:
        # assume 1 fixation only
        self.percentDuration = 1.0
    else:
      self.percentDuration = (self.duration - mindur)/(range)
