#!/usr/bin/env python3
# Scanpath.py
# This class encapsulates the analysis program

# system includes
import sys
import os
import math
#import xml.etree.ElementTree as ET
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import tostring
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy import signal
from scipy import integrate
import xmlpp
import pywt

# local includes
from point import Point
from fixation import Fixation
from grid import Grid
from aoi import AOI
from monitor import Monitor

import bw
import sg

def modmax(d):

  # compute signal modulus
  m = [0.0]*len(d)
  for i in range(len(d)):
    m[i] = math.fabs(d[i])

  # compute thresholded (modmax) signal: if value is larger than both
  # neighbours, and strictly larger than either of them, then it is a
  # local maximum
  t = [0.0]*len(d)
  for i in range(len(d)):
    ll = m[i-1] if i >= 1 else m[i]
    oo = m[i]
    rr = m[i+1] if i < len(d)-2 else m[i]
    if (ll <= oo and oo >= rr) and (ll < oo or oo > rr):
      # compute magnitude
      t[i] = math.sqrt(d[i]**2)
    else:
      t[i] = 0.0

  return t

class Run:
  def __init__(self,subj,runn,date):
    self.subjectName = subj
    self.runNumber = runn
    self.dateTime = date

class Scanpath:
  def __init__(self):
    self.normal = [1280.0, 1024.0]
    self.refresh_rate = 50
    self.fileName = None
    self.expName = None
    self.run = Run("1","1","January")
    self.fileExtension = None
    self.filter = None
    self.pupild = []		# pupil diameter (raw, corrected for PFE)
    self.pupils = []		# pupil diameter (smoothed, corrected for PFE)
    self.pupildwt = []
    self.cA1 = []
    self.cA2 = []
    self.cA3 = []
    self.cA4 = []
    self.cA = []
    self.cD4 = []
    self.cD3 = []
    self.cD2 = []
    self.cD1 = []
    self.cD = []
    self.cD_LH = []
    self.cD4t = []
    self.cD3t = []
    self.cD2t = []
    self.cD1t = []
    self.cDt = []
    self.cD_LHt = []
    self.bpcpd = []
    self.pcpd = []
    self.gazepoints = []
    self.smthpoints = []
    self.velocity = []
    self.magnitude = []
    self.acceleration = []
    self.fixations = []
    # a jaggest list (list of lists; each sublist is a list of indeces)
    # this should be same length as list of fixations
    # each element is a list of indeces into smthpoints (and velocity)
    #   that comprise the fixation
    self.fixpoints = []
    self.fixpoints_vx = []
    self.fixpoints_vx2 = []
    self.fixpoints_vy = []
    self.fixpoints_vy2 = []
    self.msac_rates = []
    self.K = []
    self.ICA = 0.0
    self.ICA_LH = 0.0
    # similar to microsaccades, need sacpoints for list of list of indeces
    # pertaining to saccade sequences
    self.sacpoints = []

  def parseFile(self, fileName, w, h, screen, dist, herz):
    self.fileName = fileName

    # get text after last '.' symbol and force it to lowercase
    self.fileExtension = fileName.split('.')[-1].lower()

    if self.fileExtension == 'xml':
      # currently, Mirametrix stimulus dimensions are NOT
      # recorded in the input files, so we have to hardcode
      # it here.
      self.normal = [float(w), float(h)]
      self.refresh_rate = float(herz)
      self.parseXML(w,h,screen,dist)
    elif self.fileExtension == 'csv':
      self.normal = [float(w), float(h)]
      self.refresh_rate = float(herz)
      self.parseCSV(w,h,screen,dist)
    elif self.fileExtension == 'raw':
      self.normal = [float(w), float(h)]
      self.refresh_rate = float(herz)
      self.parseRAW(w,h,screen,dist)
    elif self.fileExtension == 'fxd':
      self.normal = [float(w), float(h)]
      self.refresh_rate = float(herz)
      self.parseFXD(w,h,screen,dist)
    elif self.fileExtension == 'txt':
      self.normal = [float(w), float(h)]
      self.refresh_rate = float(herz)
      self.parseClearViewEFD(w,h,screen,dist)
    else:
      print('Unsupported file extension! (%s)' % (self.fileExtension))

  # parse fxd (x y d) fixation file
  def parseFXD(self,w,h,screen,dist):
    try:
      f = open(self.fileName,'r')
    except IOError:
      print("Can't open file: " + self.fileName)
      return

    # read lines, throwing away first one (header)
  # linelist = f.readlines()
    # linelist = f.read().split('\r')
    linelist = f.read().splitlines()
    # no header in raw files
  # linelist = linelist[1:]

    for line in linelist:
      entry = line.split(' ')

      # coords should already be normalized
      x = float(entry[0])/self.normal[0]
      y = float(entry[1])/self.normal[1]
      d = float(entry[2])

      x, y = self.clamp(x, y)

      self.fixations.append(Fixation(x,y,0,d))

    # compute mean fixation duration
    durations = [0.0]*len(self.fixations)
    for i in range(len(self.fixations)):
      durations[i] = self.fixations[i].getDuration()
    avedur = np.mean(durations)
    mindur = min(durations)
    maxdur = max(durations)

        #print "Read: %d fixations, %f sec mean duration, [%f,%f]" % \
         #         (len(self.fixations), avedur, mindur, maxdur)

    # normalize
    for i in range(len(self.fixations)):
      self.fixations[i].normalizeDuration(mindur,maxdur)

  # parse a raw (x y t) source file
  def parseRAW(self,w,h,screen,dist):
    try:
      f = open(self.fileName,'r')
    except IOError:
      print("Can't open file: " + self.fileName)
      return

    # read lines, throwing away first one (header)
  # linelist = f.readlines()
    # linelist = f.read().split('\r')
    linelist = f.read().splitlines()
    # no header in raw files
  # linelist = linelist[1:]

    firstline = True
    st = 0.0
    for line in linelist:
      entry = line.split(' ')

      # coords should already be normalized
      x = float(entry[0])
      y = float(entry[1])
      d = float(entry[2])
      t = float(entry[3])
      err = 'None'

      # record first time stamp
      if firstline:
        firstline = False
        st = t

      # make time relative to start, already in seconds
      t = (t - st)
      # make time relative to start, convert to seconds from milliseconds
#     t = (t - st)/1000.0
      # make time relative to start, convert to seconds from microseconds
#     t = (t - st)/1000000.0

      x, y = self.clamp(x, y)
      self.gazepoints.append(Point(x,y,t,err))

      # this file has pupil diameter!

      # correct the pupil diameter via PFE (pupil foreshortening effect)
      monitor = Monitor(w,h,screen,dist)
      # assume normalized gaze coordinates (0,1) with screen center (.5,.5)
      cx, cy = 0.5, 0.5
      dx, dy = x - cx, y - cy
      dpix = math.sqrt(dx*dx + dy*dy)
      degs = monitor.pix2deg(dpix)
      pfe = math.sqrt(math.cos(math.radians(degs)))
#     this is the more aggressive variant but maxes out at 0.94 at center
#     when it should be 1.0
#     pfe = math.sqrt(0.992 * math.cos((math.radians(degs + 5.3)/1.121)))

      # nan values for pupild really messed up Butterworth filtering
      if not math.isnan(d):
        self.pupild.append(Point(d,d/pfe,t,err))

  # parse a raw (x y t) source file
  def parseCSV(self,w,h,screen,dist):
    try:
      f = open(self.fileName,'r')
    except IOError:
      print("Can't open file: " + self.fileName)
      return

    # read lines, throwing away first one (header)
  # linelist = f.readlines()
    # linelist = f.read().split('\r')
    linelist = f.read().splitlines()
    # eat header in csv files
    linelist = linelist[1:]

    firstline = True
    st = 0.0
    for line in linelist:
      entry = line.split(',')

#     print("entry: ",entry)

      # coords should already be normalized
      x = float(entry[0])
      y = float(entry[1])
      d = float(entry[2])
      t = float(entry[3])
      err = 'None'

      # record first time stamp
      if firstline:
        firstline = False
        st = t

      # make time relative to start, already in seconds
      t = (t - st)
      # make time relative to start, convert to seconds from milliseconds
#     t = (t - st)/1000.0
      # make time relative to start, convert to seconds from microseconds
#     t = (t - st)/1000000.0

      x, y = self.clamp(x, y)
      self.gazepoints.append(Point(x,y,t,err))

      # this file has pupil diameter!

      # correct the pupil diameter via PFE (pupil foreshortening effect)
      monitor = Monitor(w,h,screen,dist)
      # assume normalized gaze coordinates (0,1) with screen center (.5,.5)
      cx, cy = 0.5, 0.5
      dx, dy = x - cx, y - cy
      dpix = math.sqrt(dx*dx + dy*dy)
      degs = monitor.pix2deg(dpix)
      pfe = math.sqrt(math.cos(math.radians(degs)))
#     this is the more aggressive variant but maxes out at 0.94 at center
#     when it should be 1.0
#     pfe = math.sqrt(0.992 * math.cos((math.radians(degs + 5.3)/1.121)))

      # nan values for pupild really messed up Butterworth filtering
      if not math.isnan(d):
        self.pupild.append(Point(d,d/pfe,t,err))

  # print "Read: ", len(self.gazepoints), "points"

  # parse a mirametrix source file
  # (first part is old code used to parse Alicann's redone XML file)
  def parseXML(self,w,h,screen,dist):
      #print "Parsing: ",self.fileName
    tree = ET.parse(self.fileName)
      #print "tree = %s" % tree

    # should be something like Experiment, {}
    root = tree.getroot()
      #print "root.tag = %s, root.attrib = %s" % (root.tag, root.attrib)

    for rec in root:
      gazeattr = rec.attrib

      bogus = False

      t = float(gazeattr['TIME'])
      xl = float(gazeattr['LPOGX'].split()[0])
      yl = float(gazeattr['LPOGY'].split()[0])
      xr = float(gazeattr['RPOGX'].split()[0])
      yr = float(gazeattr['RPOGY'].split()[0])
      errl = gazeattr['LPOGV']
      errr = gazeattr['RPOGV']

      if(errl == "0" and errr == "0"):
        bogus = True

      if(not bogus):
        x = (xl+xr)/2.0
        y = (yl+yr)/2.0

        x, y = self.clamp(x, y)

        self.gazepoints.append(Point(x,y,t,"None"))

    # print "Read: ", len(self.gazepoints), "points"

  # dump out XML
  def dumpXML(self,fileName,w,h):
#   print "Dumping ", fileName

    # buld tree structure
    exp = ET.Element('Experiment',filename=self.expName)
#   print exp.tag, exp.attrib

    run = ET.SubElement(exp,'Run')
    run.attrib['subjectName'] = self.run.subjectName
    run.attrib['runNumber'] = self.run.runNumber
    run.attrib['dateTime'] = self.run.dateTime

#   print run.tag, run.attrib

#   for pt in self.gazepoints:
#     data = ET.SubElement(run,'GazeData')
#     data.attrib['Id'] = "1"
#     data.attrib['Time'] = "%f" % (pt.gettimestamp())
#     data.attrib['Error'] = "%s" % (pt.getStatus())
#     data.attrib['X'] = "%f" % (pt.at(0))
#     data.attrib['Y'] = "%f" % (pt.at(1))

    for fx in self.fixations:
      data = ET.SubElement(run,'FixationData')
      data.attrib['Id'] = "1"
      data.attrib['Time'] = "%f" % (fx.gettimestamp())
      data.attrib['FX'] = "%f" % (fx.at(0) * float(w))
      data.attrib['FY'] = "%f" % (fx.at(1) * float(h))
      data.attrib['Duration'] = "%f" % (fx.getDuration())
      data.attrib['PercentDuration'] = "%f" % (fx.getPercentDuration())

    # wrap experiment in ElementTree instance, output
    # this is machine-readable output, not very human-readable
#   tree = ET.ElementTree(exp)
##  tree.write(fileName)

    # this is more human-readable but quite slow
    outfile = open(fileName,'w')
#   print xmlpp.pprint(tostring(exp))
    xmlpp.pprint(tostring(exp),output=outfile)
    outfile.close()

  def dumpDAT(self,fileName,w,h):

    outfile = open(fileName,'w')
    for pt in self.gazepoints:
      x = pt.at(0) * float(w)
      y = pt.at(1) * float(h)
      str = "%f %f %f\n" % (pt.gettimestamp(), x, y)
      outfile.write(str)
    outfile.close()

  def getpcpd(self,fileName,w,h,bT,eT):

    # reset running mean and counter
    ud,ucd = 0., 0.
    k = 0

    print("bT, eT = ", bT, eT)

    # first compute mean pupil diameter, over initial T time period
    for i in range(len(self.pupils)):

      # (smoothed) data we are classifying
      d, cd = self.pupils[i].at(0), self.pupils[i].at(1)

      if self.pupils[i].gettimestamp() < bT:
        # compute running mean
        # from Brown, Robert Grover, "Introduction to Random Signal
        #   Analysis and Kalman Filtering", John Wiley & Sons, New York, NY
        #   1983 [p.182] TK5102.5 .B696
        ud = float(k)/float(k+1) * ud + 1.0/float(k+1) * d
        ucd = float(k)/float(k+1) * ucd + 1.0/float(k+1) * cd
        k += 1

    # reset running mean and counter
    upcpd, ucpcpd = 0., 0.
    k = 0

    # second compute the PCPD w.r.t. the mean pupil diameter computed above
    for i in range(len(self.pupils)):

      if self.pupils[i].gettimestamp() < eT:
        # (smoothed) data we are classifying
        d, cd = self.pupils[i].at(0), self.pupils[i].at(1)

        # compute running mean of PCPD
        pcpd, cpcpd = d - ud, cd - ucd
        upcpd = float(k)/float(k+1) * upcpd + 1.0/float(k+1) * pcpd
        ucpcpd = float(k)/float(k+1) * ucpcpd + 1.0/float(k+1) * cpcpd
        k += 1

        t = self.pupils[i].gettimestamp()
        self.pcpd.append(Point(pcpd,upcpd,t,'None'))

    outfile = open(fileName,'w')
    str = "%f %f %f\n" % (ud,upcpd,ucpcpd)
    outfile.write(str)
    outfile.close()

  def getbpcpd(self,fileName,w,h,bscanpath):

    # reset running mean and counter
    ud, ucd = 0., 0.
    k = 0

    # first compute mean pupil diameter, over initial T time period of baseline
    for i in range(len(bscanpath.pupils)):

      # (smoothed) data we are classifying
      d, cd = bscanpath.pupils[i].at(0), bscanpath.pupils[i].at(1)

      # compute running mean
      # from Brown, Robert Grover, "Introduction to Random Signal
      #   Analysis and Kalman Filtering", John Wiley & Sons, New York, NY
      #   1983 [p.182] TK5102.5 .B696
      ud = float(k)/float(k+1) * ud + 1.0/float(k+1) * d
      ucd = float(k)/float(k+1) * ucd + 1.0/float(k+1) * cd
      k += 1

    print("^^^^^^^^^^^^^^^^ BASELINE PUPIL DIAMETER MEAN: ", ud)

    # reset running mean and counter
    upcpd, ucpcpd = 0., 0.
    k = 0

    # second compute the PCPD w.r.t. the mean pupil diameter computed above
    for i in range(len(self.pupils)):

      # (smoothed) data we are classifying
      d, cd = self.pupils[i].at(0), self.pupils[i].at(1)

      # compute running mean of PCPD
      pcpd, cpcpd = d - ud, cd - ucd
      upcpd = float(k)/float(k+1) * upcpd + 1.0/float(k+1) * pcpd
      ucpcpd = float(k)/float(k+1) * ucpcpd + 1.0/float(k+1) * cpcpd

      k += 1

      t = self.pupils[i].gettimestamp()
      self.bpcpd.append(Point(pcpd,cpcpd,t,'None'))

    print("^^^^^^^^^^^^^^^^ CHANGE IN PUPIL DIAMETER MEAN: ", upcpd)
    print("^^^^^^^^^^^^^^^^ PERCENT CHANGE IN PUPIL DIAMETER MEAN: ", upcpd/ud * 100.0)
    print("^^^^^^^^^^^^^^^^ CHANGE IN PFE PUPIL DIAMETER MEAN: ", ucpcpd)
    print("^^^^^^^^^^^^^^^^ PERCENT CHANGE IN PFE PUPIL DIAMETER MEAN: ", ucpcpd/ucd * 100.0)

    outfile = open(fileName,'w')
    str = "%f %f %f\n" % (ud,upcpd,ucpcpd)
    outfile.write(str)
    outfile.close()

  def pdwtLH(self,fileName,w,h,herz,sfdegree,sfcutoff,lof,hif):
    d = []
    for pt in self.pupild:
      d.append(pt.at(0))

    # lof, hif are low frequency band, high frequency band, resp.
    # maxlevel = \floor{\log_2{\left(\frac{n_d}{n_f-1}\right)}}
    w = pywt.Wavelet('sym16')
    maxlevel = pywt.dwt_max_level(len(d), filter_len=w.dec_len)
    print("max DWT level: ", maxlevel)

    hif = 1
    lof = 1 if maxlevel//2 < 1 else maxlevel//2

    print("lo DWT level: ", lof)
    print("hi DWT level: ", hif)

    # get the wavelet detail coefficients for high and low frequenices
    cD_H = pywt.downcoef('d',d,'sym16','per',level=hif)
    cD_L = pywt.downcoef('d',d,'sym16','per',level=lof)

    # normalize
    cD_H[:] = [x / math.sqrt(2**hif) for x in cD_H]
    cD_L[:] = [x / math.sqrt(2**lof) for x in cD_L]

    cD_LH = cD_L

    # obtain the LH ratio (HF:LH in this case)
    #
    # note that I am artificially extending the signal by using the 2*i
    # index into the original data's timeline
#   print "len(cD_L): ", len(cD_L)
#   print "len(cD_H): ", len(cD_H)
    # length of lower frequency detail coefficients will be shorter than
    # the high frequency octaves, hence we should iterate over the low
    # frequency octave and then use those inndices to index to the high
    # frequency octave
    for i in range(len(cD_L)):
#     print "cD_L: ", cD_L[i]
#     print "cD_H: ", cD_H[((2**lof)/(2**hif))*i]
      # does low:high ratio really make sense?
      cD_LH[i] = cD_L[i] / cD_H[((2**lof)//(2**hif))*i]
      # why not high:low signal ratio
#     cD_LH[i] = cD_H[((2**lof)/(2**hif))*i] / cD_L[i]

    self.cD_LH = cD_LH

    # the rest of the algorithm proceeds like the original IPA:
    # find the modmax spikes in the LH signal and threshold
    cD_LHm = modmax(cD_LH)
    cD_LHt = cD_LHm
    lambda_univ = np.std(cD_LHm) * math.sqrt(2.0*np.log2(len(cD_LHm)))
#   cD_LHt = pywt.threshold(cD_LHm,lambda_univ,mode="hard")
    cD_LHt = pywt.threshold(cD_LHm,lambda_univ,mode="less")

    self.cD_LHt = cD_LHt

    outfile = open(fileName,'w')
    for i in range(len(cD_LH)):
      str = "%f %f\n" % (self.pupild[(2**lof)*i].gettimestamp(), cD_LH[i])
      outfile.write(str)
    outfile.close()

    return lof, hif

  def pIPALH(self,fileName,w,h,herz,sfdegree,sfcutoff):

    # this function computes Sandra P. Marshall's Index of Cognitive Activity
    # using the DWT of the pupil diamter: pwdt (above) must be run first!

    ts = self.pupild[0].gettimestamp()
    te = self.pupild[-1].gettimestamp()
    tt = te - ts
#   print "ts: ", ts
#   print "te: ", te
#   print "total time: ", tt
    ctr = 0
    for i in range(len(self.cD_LHt)):
      if math.fabs(self.cD_LHt[i]) > 0:
        ctr += 1
    self.ICA_LH = float(ctr)/tt if tt > 0 else 0

    outfile = open(fileName,'w')
    str = "%f\n" % (self.ICA_LH)
    outfile.write(str)
    outfile.close()

    return

    # below is stuff I tried with power spectra analysis...it seemed promising
    # at first but couldn't quite get it to work

    # see: https://pycwt.readthedocs.io/en/latest/tutorial.html#time-series-spectral-analysis-using-wavelets
    # trying to get spectral power from LF:HF ratio
    # compute the normalized wavelet power spectra
    # these guys say to use thresholding according to:
    #  cD_LH[i] = sng(cD_LH[i])max(0,|cD_LH[i]| - lambda)
    # where lambda is universal threshold
    # https://pdfs.semanticscholar.org/cc41/7e0ed55a2b1129f68ded1e3c4990b6ff24d9.pdf
    cD_LH = self.cD_LH
    lambda_univ = np.std(cD_LH) * math.sqrt(2.0*np.log2(len(cD_LH)))
    # in soft thresholding, data values with absolute value less than param
    # are replaced with substitute. Data values with absolute value greater
    # or equal to the thresholding value are shrunk toward zero by value;
    # in other words, the new value is
    # data/np.abs(data) * np.maximum(np.abs(data) - value, 0).
#   cD_LHt = pywt.threshold(cD_LH,lambda_univ,mode="soft")
    # in less thresholding, the data is replaced with substitute where data
    # is above the thresholding value. Lesser data values pass untouched.
    cD_LHt = pywt.threshold(cD_LH,lambda_univ,mode="less")

#   power = (np.abs(self.cD_LH)) ** 2
    power = (np.abs(cD_LHt)) ** 2
    print("Power: ", power)
    # integrate
    integral = integrate.simps(power)
    print("integral: ", integral)
    fWc = 1000.0	# sampling rate (1000 Hz)
    print('1.0/fWc * integral = ', 1.0/fWc * integral)
    # see: https://www2.mps.mpg.de/solar-system-school/lectures/numerical_methods/span.pdf
    rms = math.sqrt(1.0/fWc * integral)
    print("rms: ", rms)

    outfile = open(fileName,'w')
 #  str = "%f\n" % (1.0/fWc * integral)
    str = "%f\n" % (rms)
    outfile.write(str)
    outfile.close()

  def pdwtlevel(self,fileName,w,h,herz,sfdegree,sfcutoff,dwtlvl):
    d = []
    for pt in self.pupild:
      d.append(pt.at(0))

    cA = pywt.downcoef('a',d,'sym16','per',level=dwtlvl)
    cD = pywt.downcoef('d',d,'sym16','per',level=dwtlvl)

    cA[:] = [x / math.sqrt(2**dwtlvl) for x in cA]
    cD[:] = [x / math.sqrt(2**dwtlvl) for x in cD]

    self.cA = cA
    self.cD = cD

    cDm = modmax(cD)
    cDt = cDm
    lambda_univ = np.std(cDm) * math.sqrt(2.0*np.log2(len(cDm)))
    cDt = pywt.threshold(cDm,lambda_univ,mode="hard")

    self.cDt = cDt

    # note that I am artificially extending the signal by using the 2*i
    # index into the original data's timeline
    for i in range(len(cA)):
      self.pupildwt.append(Point(cA[i],cD[i],\
                           self.pupild[(2**dwtlvl)*i].gettimestamp(),'None'))

    outfile = open(fileName,'w')
    for i in range(len(cD)):
      str = "%f %f\n" % (self.pupild[(2**dwtlvl)*i].gettimestamp(), cD[i])
      outfile.write(str)
    outfile.close()

  def pIPAlevel(self,fileName,w,h,herz,sfdegree,sfcutoff):

    # this function computes Sandra P. Marshall's Index of Cognitive Activity
    # using the DWT of the pupil diamter: pwdt (above) must be run first!

    ts = self.pupild[0].gettimestamp()
    te = self.pupild[-1].gettimestamp()
    tt = te - ts
#   print "ts: ", ts
#   print "te: ", te
#   print "total time: ", tt
    ctr = 0
    for i in range(len(self.cDt)):
      if math.fabs(self.cDt[i]) > 0:
        ctr += 1
    self.ICA = float(ctr)/tt if tt > 0 else 0

    outfile = open(fileName,'w')
    str = "%f\n" % (self.ICA)
    outfile.write(str)
    outfile.close()

  def pdwt(self,fileName,w,h,herz,sfdegree,sfcutoff):
    d = []
    for pt in self.pupild:
      d.append(pt.at(0))
    # note: len(cA) == len(cD) == floor((len(d) + wavelet.dec_len - 1)/2)
#   print "len(d): ", len(d)
    # pywt.wavelist('db') returns:
    # db1..db20
    # filter length is twice the index
    # default mode is 'symmetric'
    # other modes: 'zero', 'constant', 'periodic', 'smooth', 'periodization'
    # 'per' mode gives n/2 lengthed signal regardless of wavelet size
#   (cA, cD) = pywt.dwt(d,'db1')
    # for 60 Hz, according to Marshall
#   (cA, cD) = pywt.dwt(d,'db4','per')
    # for 250 Hz, according to Marshall
#   (cA, cD) = pywt.dwt(d,'db16','per')
#   (cA1, cD1) = pywt.dwt(d,'sym16','per')
    try:
      (cA1, cD1) = pywt.wavedec(d,'sym16','per',level=1)
    except ValueError:
      return
#   print "len(cA): ", len(cA1)
#   print "len(cD): ", len(cD1)
#   (cA2, cD2, cD1) = pywt.wavedec(d,'db16','per',level=2)
    try:
      (cA2, cD2, cD1) = pywt.wavedec(d,'sym16','per',level=2)
    except ValueError:
      return
#   print "len(cA2): ", len(cA2)
#   print "len(cD2): ", len(cD2)
#   print "len(cD1): ", len(cD1)
#   (cA3, cD3, cD2, cD1) = pywt.wavedec(d,'db16','per',level=3)
    try:
      (cA3, cD3, cD2, cD1) = pywt.wavedec(d,'sym16','per',level=3)
    except ValueError:
      return
#   print "len(cA3): ", len(cA3)
#   print "len(cD3): ", len(cD3)
#   print "len(cD2): ", len(cD2)
#   print "len(cD1): ", len(cD1)
    try:
      (cA4, cD4, cD3, cD2, cD1) = pywt.wavedec(d,'sym16','per',level=4)
    except ValueError:
      return
#   print "len(cA4): ", len(cA4)
#   print "len(cD4): ", len(cD4)
#   print "len(cD3): ", len(cD3)
#   print "len(cD2): ", len(cD2)
#   print "len(cD1): ", len(cD1)

    # normalize
    cA1[:] = [x / math.sqrt(2.0) for x in cA1]
    cA2[:] = [x / math.sqrt(4.0) for x in cA2]
    cA3[:] = [x / math.sqrt(8.0) for x in cA3]
    cA4[:] = [x / math.sqrt(16.0) for x in cA4]

    cD1[:] = [x / math.sqrt(2.0) for x in cD1]
    cD2[:] = [x / math.sqrt(4.0) for x in cD2]
    cD3[:] = [x / math.sqrt(8.0) for x in cD3]
    cD4[:] = [x / math.sqrt(16.0) for x in cD4]

    self.cA1 = cA1
    self.cA2 = cA2
    self.cA3 = cA3
    self.cA4 = cA4

    self.cD4 = cD4
    self.cD3 = cD3
    self.cD2 = cD2
    self.cD1 = cD1

    # threshold using universal threshold
    # \lambda_{univ} = \hat{\sigma}\sqrt(2\log{n})
    # where \hat{\sigma} is the standard deviation of the noise
#   cD = pywt.thresholding.hard(cD,lambda_univ)
#   cD = modmax(cD)
#   cD = pywt.thresholding.hard(cD,lambda_univ)

    cD1m = modmax(cD1)
    cD1t = cD1m
    lambda_univ = np.std(cD1m) * math.sqrt(2.0*np.log2(len(cD1m)))
#   cD1t = pywt.thresholding.hard(cD1m,lambda_univ)
    cD1t = pywt.threshold(cD1m,lambda_univ,mode="hard")
#   lambda_univ = np.std(cD1) * math.sqrt(2.0*np.log2(len(cD1)))
#   print "lambda_univ: ", lambda_univ
#   cD1t = pywt.thresholding.hard(cD1,lambda_univ)

    cD2m = modmax(cD2)
    cD2t = cD2m
    lambda_univ2 = np.std(cD2m) * math.sqrt(2.0*np.log2(len(cD2m)))
    cD2t = pywt.threshold(cD2m,lambda_univ2,mode="hard")
#   lambda_univ2 = np.std(cD2) * math.sqrt(2.0*np.log2(len(cD2)))
#   print "lambda_univ2: ", lambda_univ2
#   cD2t = pywt.thresholding.hard(cD2,lambda_univ2)

    cD3m = modmax(cD3)
    cD3t = cD3m
    lambda_univ3 = np.std(cD3m) * math.sqrt(2.0*np.log2(len(cD3m)))
    cD3t = pywt.threshold(cD3m,lambda_univ3,mode="hard")
#   lambda_univ3 = np.std(cD3) * math.sqrt(2.0*np.log2(len(cD3)))
#   print "lambda_univ3: ", lambda_univ3
#   cD3t = pywt.thresholding.hard(cD3,lambda_univ3)

    cD4m = modmax(cD4)
    cD4t = cD4m
    lambda_univ4 = np.std(cD4m) * math.sqrt(2.0*np.log2(len(cD4m)))
    cD4t = pywt.threshold(cD4m,lambda_univ4,mode="hard")
#   lambda_univ4 = np.std(cD4) * math.sqrt(2.0*np.log2(len(cD4)))
#   print "lambda_univ4: ", lambda_univ4
#   cD4t = pywt.thresholding.hard(cD4,lambda_univ4)

    self.cD4t = cD4t
    self.cD3t = cD3t
    self.cD2t = cD2t
    self.cD1t = cD1t

    # note that I am artificially extending the signal by using the 2*i
    # index into the original data's timeline
#
#   for i in range(len(cA)):
#     self.pupildwt.append(Point(cA[i],cD[i],self.pupild[2*i].gettimestamp(),'None'))
#
    for i in range(len(cA2)):
      self.pupildwt.append(Point(cA2[i],cD2[i],self.pupild[4*i].gettimestamp(),'None'))

#   for i in range(len(cA3)):
#     self.pupildwt.append(Point(cA3[i],cD3[i],self.pupild[8*i].gettimestamp(),'None'))

#   for i in range(len(cA4)):
#     self.pupildwt.append(Point(cA4[i],cD4[i],self.pupild[16*i].gettimestamp(),'None'))

    outfile = open(fileName,'w')
#   for i in range(len(cD4)):
#     str = "%f %f\n" % (self.pupild[16*i].gettimestamp(), cD4[i])
    for i in range(len(cD2)):
      str = "%f %f\n" % (self.pupild[4*i].gettimestamp(), cD2[i])
      outfile.write(str)
    outfile.close()

  def pICA(self,fileName,w,h,herz,sfdegree,sfcutoff):

    # this function computes Sandra P. Marshall's Index of Cognitive Activity
    # using the DWT of the pupil diamter: pwdt (above) must be run first!

    ts = self.pupild[0].gettimestamp()
    te = self.pupild[-1].gettimestamp()
    tt = te - ts
#   print "ts: ", ts
#   print "te: ", te
#   print "total time: ", tt
    ctr = 0
#   for i in range(len(self.cD1t)):
#     if math.fabs(self.cD1t[i]) > 0:
#       ctr += 1
#   for i in range(len(self.cD2t)):
#     if math.fabs(self.cD2t[i]) > 0:
#       ctr += 1
#   for i in range(len(self.cD3t)):
#     if math.fabs(self.cD3t[i]) > 0:
#       ctr += 1
    for i in range(len(self.cD4t)):
      if math.fabs(self.cD4t[i]) > 0:
        ctr += 1
    self.ICA = float(ctr)/tt
    # 1/ICA gives good results, consistent with pupil diameter
    # but is meaningless...ICA appears to give results more consistent with
    # expected microsaccade rate
#   if self.ICA > 0:
#     self.ICA = 1.0 / self.ICA
#   print "total count: ", ctr
#   print "ICA: ", self.ICA

    outfile = open(fileName,'w')
    str = "%f\n" % (self.ICA)
    outfile.write(str)
    outfile.close()

  def pscat(self,pupils):

    # concatenate pupils into current scanpath

    for pt in pupils:
      err = pt.getStatus()

      # (smoothed and smoothed pfe) data we are concatenating
      self.pupils.append(Point(pt.at(0),pt.at(1),pt.gettimestamp(),err))

  def psmooth(self,fileName,w,h,herz,sfdegree,sfcutoff):

    # use Butterworth filter for smoothing
#   self.pupils = bw.applyBWFilter(self.pupild,sfdegree,herz,sfcutoff)

    nyq = herz/2.0

    # use Butterworth filter for smoothing (using scipy)
    # cf is the critical frequency
    # seemed ok for 500 Hz sampling
#   cf = (.008*herz)/nyq # way too smooth
#   cf = (.08*herz)/nyq
#   cf = (.16*herz)/nyq # more smooth
#   cf = (.25*herz)/nyq #
    cf = (.32*herz)/nyq # ok
#   cf = (.5*herz)/nyq
#   cf = herz/nyq # less smooth
    print("fs: %f,  sf: %f, cf: %f" % (herz,sfdegree,cf))
#   b, a = signal.butter(sfdegree,1.0/cf,'lowpass',False,'ba')
    b, a = signal.butter(sfdegree,cf,btype='lowpass',analog=False,output='ba',fs=herz)
#   b, a = signal.butter(sfdegree,1.0/cf,btype='low',analog=False,output='ba',fs=herz)
    cf = 4.0
    # seems to work with sosfilt (but not sosfiltfilt)
#   sos = signal.butter(sfdegree,1.0,btype='low',analog=False,output='sos',fs=herz)
    # for sosfiltfilt
    sos = signal.butter(sfdegree,cf,btype='low',analog=False,output='sos',fs=herz)
#   print("b: ",b)
#   print("a: ",a)
#   print("sos: ",sos)
#   print("len(sos): ",len(sos))
#   print("len(sos[0]): ",len(sos[0]))
#   print("len(sos[1]): ",len(sos[1]))
#   b, a = signal.butter(sfdegree,sfcutoff,'low',False,'ba')
#   sos = signal.butter(sfdegree,sfcutoff,'low',fs=herz,output='sos')

    d = []
    # construct 1D arrays (for raw data)
    for pt in self.pupild:
      d.append(float(pt.at(0)))
    sd = signal.filtfilt(b, a, np.asarray(d), padlen=len(d)-1)
#   sd = signal.filtfilt(b, a, d)
#   sd = signal.lfilter(b, a, d)
#   sd = signal.sosfilt(sos, d)
#   sd = signal.sosfilt(sos, np.asarray(d))
#   sd = signal.sosfiltfilt(sos, d, padtype=None, padlen=0)
#   sd = signal.sosfiltfilt(sos, d, padlen=len(d)-1)
#   sd = signal.sosfiltfilt(sos, np.asarray(d), padtype='constant', padlen=len(d)-1)
#   print("d: ",d)
#   print("sd: ",sd)
#   print("len(sd): ",len(sd))

    # delete d[], supposedly faster than 'del d[:]'
    d[:] = []
    # construct 1D arrays (for pfe data)
    for pt in self.pupild:
      d.append(pt.at(1))
    scd = signal.filtfilt(b, a, d, padlen=len(d)-1)
#   scd = signal.sosfiltfilt(sos, np.asarray(d), padlen=len(d)-1)

    i=0
    for pt in self.pupild:
      err = pt.getStatus()
      self.pupils.append(Point(sd[i],scd[i],pt.gettimestamp(),err))
      i += 1

    outfile = open(fileName,'w')
    for pt in self.pupils:
      d = pt.at(0)
      str = "%f %f\n" % (pt.gettimestamp(), d)
      outfile.write(str)
    outfile.close()

  def smooth(self,fileName,w,h,herz,sfdegree,sfcutoff,smooth):

    if smooth:
        # use Butterworth filter for smoothing
        self.smthpoints = bw.applyBWFilter(self.gazepoints,sfdegree,herz,sfcutoff)
    else:
        self.smthpoints = self.gazepoints
    outfile = open(fileName,'w')
    for pt in self.smthpoints:
      x = pt.at(0) * float(w)
      y = pt.at(1) * float(h)
      str = "%f %f %f\n" % (pt.gettimestamp(), x, y)
      outfile.write(str)
    outfile.close()

  def differentiate(self,fileName,w,h,screen,dist,herz,dfwidth,dfdegree,dfo):

    # differentiate smoothed points with Savitzky-Golay filter
    diffpoints = sg.applySGFilter(self.smthpoints,dfwidth,dfdegree,dfo)
    accelpoints = sg.applySGFilter(diffpoints, dfwidth, dfdegree, dfo)

    # sampling period in s
    period = float(1.0/float(herz))
    dt = period * float(dfwidth)
    # shouldn't this be 2 * dfwidth?  SG filter is 10 coefs long
#   dt = period * 2.0 * float(dfwidth)

    r = math.sqrt(float(w)*float(w) + float(h)*float(h))
    dpi = r/float(screen)

    D = float(dist)
#   fov = 2*math.degrees(math.atan(math.radians(screen/(2*D))))
    fov = 2*math.degrees(math.atan2(screen,2*D))
    fovx = 2*math.degrees(math.atan2(float(w)/dpi,2*D))
    fovy = 2*math.degrees(math.atan2(float(h)/dpi,2*D))

#   print "screen subtends %f (%f x %f) degrees" % \
#                   (float(fov),float(fovx),float(fovy))
#   print "screen aspect = %f" % float(float(w)/float(h))
#   print "sampling period = %f (ms)" % float(period * 1000.0)
#   print "filter window = %f (ms)" % float(dt * 1000.0)
#   print "dpi = %f" % dpi

    for pt in diffpoints:

      # distance covered in pixels (across diff filter window size)
      dx = pt.at(0) * float(w)
      dy = pt.at(1) * float(h)

      # assume D is in inches
#     degx = 2*math.degrees(math.atan(math.radians((dx/dpi)/(2*D))))
#     degy = 2*math.degrees(math.atan(math.radians((dy/dpi)/(2*D))))
      degx = 2*math.degrees(math.atan2((dx/dpi),(2*D)))
      degy = 2*math.degrees(math.atan2((dy/dpi),(2*D)))

      # degx, degy is degrees per filter window, div by dt to get per second
      velx = degx / dt
      vely = degy / dt

      self.magnitude.append(Point(degx,degy,pt.gettimestamp(),pt.getStatus()))
      self.velocity.append(Point(velx,vely,pt.gettimestamp(),pt.getStatus()))

#     print "pts: %d, smth: %d, diff: %d, vel: %d" % \
#           (len(self.gazepoints),len(self.smthpoints), \
#            len(diffpoints),len(self.velocity))

    for pt in accelpoints:
      dx = pt.at(0) * float(w)
      dy = pt.at(1) * float(h)

      degx = 2*math.degrees(math.atan2((dx/dpi),(2*D)))
      degy = 2*math.degrees(math.atan2((dy/dpi),(2*D)))

      # degx, degy is degrees per filter window, div by dt to get per second
      accelx = degx / dt
      accely = degy / dt

      self.acceleration.append(Point(accelx,accely,pt.gettimestamp(),pt.getStatus()))


    outfile = open(fileName,'w')
    for pt in self.velocity:
      # don't scale by w,h here, already did so above
      x = pt.at(0)
      y = pt.at(1)
      str = "%f %f %f\n" % (pt.gettimestamp(), x, y)
      outfile.write(str)
    outfile.close()

  def threshold(self,fileName,w,h,T,monitor,proximity):

    # state table
    # state         |   input    |  output
    # ------------------------------------
    # (1) fixation  | a < T (1)  |    1
    # (0) saccade   | a < T (1)  |    1
    # (1) fixation  | a > T (0)  |    0
    # (0) saccade   | a > T (0)  |    0
    #

    # fixation state enums
    fixation = 1
    saccade = 0

    # assuming starting state = saccade
    state = fixation

    st = 0.
    et = 0.
    tt = 0.
    ux = 0.
    uy = 0.
    k = 0

    # initialize list of microsaccade indeces (for each fixation)
    indeces = []

    # initialize list of saccade indeces (for each saccade)
    sindeces = []

    if(len(self.velocity) > 0):

      for i in range(len(self.smthpoints)):

        # (smoothed) data we are classifying
        x = self.smthpoints[i].at(0)
        y = self.smthpoints[i].at(1)

        # corresponding velocity
        vx = self.velocity[i].at(0)
        vy = self.velocity[i].at(1)

        # saccade amplitude
        amp = math.sqrt(vx*vx + vy*vy)

        if math.fabs(amp) > float(T)/2:

          # amplitude > T
          if state == fixation:
            # state transition from fixation to saccade (fixation ends)

            # get end time of current point, compute duration
            et = self.smthpoints[i].gettimestamp()
            tt = et - st
            # don't add fixation with -ve duration
            # (due to spurious timestamp in .raw file)
            if st > 0.0 and tt > 0.0:
              # need to clamp fixation centroid: it could have gone negative
              ux, uy = self.clamp(ux, uy)
              self.fixations.append(Fixation(ux,uy,st,tt))

              # append list of indeces (sequence of points within this fixation)
              if len(indeces) > 0:
                self.fixpoints.append(indeces)

            # saccade starts
            # reset list of osaccade indeces for this saccade
            sindeces = []

          else:
            # state transition from saccade to saccade (saccade continues)
#           pass

            # add in current index into sindeces
            sindeces.append(i)

          state = saccade

        else:

          # amplitude < T
          if state == fixation:
            # state transition from fixation to fixation (fixation continues)
            pass
          else:
            # state transition from saccade to fixation (fixation starts)

            # set start time
            st = self.smthpoints[i].gettimestamp()

            # reset running mean and counter
            ux = 0.
            uy = 0.
            k = 0

            # reset list of microsaccade indeces for this fixation
            indeces = []

            # saccade ends
            # append list of indeces (sequence of points within this saccade)
            if len(sindeces) > 0:
              self.sacpoints.append(sindeces)

          # compute running mean
          # from Brown, Robert Grover, "Introduction to Random Signal
          #   Analysis and Kalman Filtering", John Wiley & Sons, New York, NY
          #   1983 [p.182] TK5102.5 .B696
          ux = float(k)/float(k+1) * ux + 1.0/float(k+1) * x
          uy = float(k)/float(k+1) * uy + 1.0/float(k+1) * y
          k += 1

          # add in current index into indeces
          indeces.append(i)

          state = fixation

#     print "i = %d, st = %f, et = %f, tt = %f" % (i,st,et,tt)

    # compute mean fixation duration
    if(len(self.fixations) > 1):
      durations = [0.0]*len(self.fixations)
      for i in range(len(self.fixations)):
        durations[i] = self.fixations[i].getDuration()
      avedur = np.mean(durations)
      mindur = min(durations)
      maxdur = max(durations)
    elif(len(self.fixations) == 1):
      avedur = self.fixations[0].getDuration()
      mindur = avedur
      maxdur = avedur
    else:
      avedur = 0.0
      mindur = avedur
      maxdur = avedur

#   print "Found: %d fixations, %f sec mean duration, [%f,%f]" % \
#                 (len(self.fixations), avedur, mindur, maxdur)
#   print "Fixation point list length: %d" % len(self.fixpoints)
#   print "Fixation points (should be list of lists):", self.fixpoints

    # normalize
    for i in range(len(self.fixations)):
      self.fixations[i].normalizeDuration(mindur,maxdur)

    outfile = open(fileName,'w')
#    for pt in self.fixations:
#      x = pt.at(0) * float(w)
#      y = pt.at(1) * float(h)
#      str = "%f %f %f %f\n" % (pt.gettimestamp(), x, y, pt.getDuration())
##     str = "%f %f %f %f\n" % (pt.gettimestamp(), x, y, pt.getPercentDuration())
#      outfile.write(str)
    sacc_dur = 0.0
    for i in range(len(self.fixations)):
      x = self.fixations[i].at(0) * float(w)
      y = self.fixations[i].at(1) * float(h)
      t = self.fixations[i].gettimestamp()
      dur = self.fixations[i].getDuration()
      st = t   # fixation timestamp is its starttime (st)
      tt = dur # fixation duration is its endtime - startime (et - st)
      et = st + tt
      if i < len(self.fixations)-1:
        dx = x - (self.fixations[i+1].at(0) * float(w))
        dy = y - (self.fixations[i+1].at(1) * float(h))
        sacc_dur = self.fixations[i+1].gettimestamp() - et
        amp = math.sqrt(dx*dx + dy*dy)
      else:
        amp = 0.0
# what we used to do
#     str = "%f %f %f %f %f\n" % (t, x, y, dur, amp)
# now adding in saccade duration
      str = "%f %f %f %f %f %f\n" % (t, x, y, dur, amp, sacc_dur)
      outfile.write(str)

    outfile.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  def microsaccades(self,fileName,w,h,screen,dist,herz):

    outfile = open(fileName,'w')

    print("Found: %d fixations" % (len(self.fixations)))
#   print "Fixation point list length: %d" % len(self.fixpoints)
#   print "Fixation points (should be list of lists):", self.fixpoints

    # sampling period in s
    period = float(1.0/float(herz))

    # 1/6 is from Engbert and Kliegl (they use a 6-tap filter)
    # note: this is wrong since it assumes a constant sampling rate,
    #       which is unrealistic; see code below where we obtain dt
    #       from the data itself
    dt = period * 6.0

    r = math.sqrt(float(w)*float(w) + float(h)*float(h))
    dpi = r/float(screen)

    D = float(dist)
#   fov = 2*math.degrees(math.atan(math.radians(screen/(2*D))))
    fov = 2*math.degrees(math.atan2(screen,2*D))
    fovx = 2*math.degrees(math.atan2(float(w)/dpi,2*D))
    fovy = 2*math.degrees(math.atan2(float(h)/dpi,2*D))

    velocities_x = []
    velocities_y = []
    velocities_x2 = []
    velocities_y2 = []

    # record microsaccade rate per fixation
    self.msac_rates = [0.0]*len(self.fixpoints)

    # for each fixation (set of raw/smoothed gaze points)
    for k in range(len(self.fixpoints)):

      fx = self.fixations[k].at(0) * float(w)
      fy = self.fixations[k].at(1) * float(h)
      dx = 2*math.degrees(math.atan2((fx - float(w)/2)/dpi,2*D))
      dy = 2*math.degrees(math.atan2((fy - float(h)/2)/dpi,2*D))
      fdeg_away = math.sqrt(dx*dx + dy*dy)
#     print "oooooooooooo DEGREES AWAY FROM CENTER: ", fdeg_away

      # should be a list of lists
      indeces = self.fixpoints[k]
#     print "Fixation %d: " % (k)
#     print indeces

      # nuke vel lists and initialize to length of this fixation (no. of points)
      # (only do this if we don't want to save velocity info)
#     del velocities_x[:]
#     del velocities_y[:]
      velocities_x = [0.0]*len(indeces)
      velocities_y = [0.0]*len(indeces)
#     del velocities_x2[:]
#     del velocities_y2[:]
      velocities_x2 = [0.0]*len(indeces)
      velocities_y2 = [0.0]*len(indeces)

      # need at least length 5 for Engbert and Kleigl's algorithm
      if len(indeces) >= 5:

        # Step 1: compute velocities
        for i in range(len(indeces)):

#         print "index: %d" % (indeces[i])

          # compute velocity
          # v_i = \frac{x_{i+2} + x_{i+1} - x_{i-1} - x_{i-2}}{6 \Delta t}
          if 2 < i  and i < len(indeces) - 2:
            dx =  1.0 * self.smthpoints[indeces[i+2]].at(0) + \
                  1.0 * self.smthpoints[indeces[i+1]].at(0) + \
                 -1.0 * self.smthpoints[indeces[i-1]].at(0) + \
                 -1.0 * self.smthpoints[indeces[i-2]].at(0)

            dy =  1.0 * self.smthpoints[indeces[i+2]].at(1) + \
                  1.0 * self.smthpoints[indeces[i+1]].at(1) + \
                 -1.0 * self.smthpoints[indeces[i-1]].at(1) + \
                 -1.0 * self.smthpoints[indeces[i-2]].at(1)

            dt = self.smthpoints[indeces[i+2]].gettimestamp() - \
                 self.smthpoints[indeces[i-2]].gettimestamp()

            # distance covered in pixels (across diff filter window size)
            # normalize filter via scaling by 0.2?  Engbert and Kliegl don't
            dx = 1.0 * dx * float(w)
            dy = 1.0 * dy * float(h)

            # assume D is in inches
            degx = 2.0*math.degrees(math.atan2((dx/dpi),(2.0*D)))
            degy = 2.0*math.degrees(math.atan2((dy/dpi),(2.0*D)))

            # overwrite original magnitude with what we calculate above
            # the difference is largely in scale since the filters
            # used (6-tap vs. Savitzky-Golay) are different widths
            self.magnitude[indeces[i]].set(degx,degy)

            # degx, degy is deg per filter window, div by dt to get per second
            if dt < 0.0001:
              continue
            velx = degx / dt if dt > 0.0 else 0.0
            vely = degy / dt if dt > 0.0 else 0.0
#           print "velocity as per Engbert and Kliegl: (%f, %f)" % (velx,vely)

            velocities_x[i] = velx
            velocities_y[i] = vely
            velocities_x2[i] = velx*velx
            velocities_y2[i] = vely*vely

#     print "Velocities x: ", velocities_x
#     print "Velocities y: ", velocities_y

      # save the velocity data: same format as fixpoints (list of lists)
      self.fixpoints_vx.append(velocities_x)
      self.fixpoints_vy.append(velocities_y)
      self.fixpoints_vx2.append(velocities_x2)
      self.fixpoints_vy2.append(velocities_y2)

      # Step 2: medians of velocities squared and square median of velocities
#     median_vx = sorted(velocities_x)[len(velocities_x)/2]
#     median_vy = sorted(velocities_y)[len(velocities_y)/2]
#     median_vx2 = sorted(velocities_x2)[len(velocities_x2)/2]
#     median_vy2 = sorted(velocities_y2)[len(velocities_y2)/2]
      median_vx = np.median(velocities_x)
      median_vy = np.median(velocities_y)
      median_vx2 = np.median(velocities_x2)
      median_vy2 = np.median(velocities_y2)

#     print "Median velocity x: ", median_vx
#     print "Median velocity y: ", median_vy
#     print "Median velocity x2: ", median_vx2
#     print "Median velocity y2: ", median_vy2

      # compute median estimator
      # \sigma_{xy} = \langle v_{xy}^{2} \rangle - \langle v_{xy} \rangle^{2}
      sigma_x = median_vx2 - (median_vx*median_vx)
      sigma_y = median_vy2 - (median_vy*median_vy)

#     print "sigma_x: ", sigma_x
#     print "sigma_y: ", sigma_y

      if sigma_x < 0.00001:
        sigma_x = 0.0
      if sigma_y < 0.00001:
        sigma_y = 0.0

      # as per Engelbert 2006 paper: uses sqrt
      sigma_x = math.sqrt(sigma_x)
      sigma_y = math.sqrt(sigma_y)

      # Mergenthaler says as lambda increases, detected microsaccades decrease
      # compute detection threshold (with \lambda = 6)
      # \eta_{xy} = \lamba \sigma_{xy}
#     eta_x = 12.0 * sigma_x
#     eta_y = 12.0 * sigma_y
#     eta_x = 8.0 * sigma_x
#     eta_y = 8.0 * sigma_y
      eta_x = 6.0 * sigma_x
      eta_y = 6.0 * sigma_y
#     eta_x = 5.0 * sigma_x
#     eta_y = 5.0 * sigma_y
#     eta_x = 3.0 * sigma_x
#     eta_y = 3.0 * sigma_y

#     print "eta_x: ", eta_x
#     print "eta_y: ", eta_y

      # Step 3: according to Engbert and Kliegl, if right eye msac starts/ends
      #         at time r1, r2, then left msac start/end times l1, l2 should
      #         be l1 < r2 and l2 > r1
      #         ok, this makes sense, but we are operating on averaged
      #         cyclopean (x,y) data, so we can't really implement this
      #         unless we store left and right eye data

      # Step 4: iterate to find microsaccades
      ts = self.smthpoints[indeces[0]].gettimestamp()
      te = self.smthpoints[indeces[-1]].gettimestamp()
      tt = te - ts
#     print "ts: ", ts
#     print "te: ", te
#     print "total time: ", tt
      # BUG? above, tt shows up as 0.000 sometimes, which seems incorrect
      tt = self.fixations[k].getDuration()
      ctr = 0

      # loop through this set of fixation points and see how long the
      # microsaccade sequences are
      i = 0
      while i < len(indeces):

        # (smoothed) data we are classifying
        x = self.smthpoints[indeces[i]].at(0)
        y = self.smthpoints[indeces[i]].at(1)
        t = self.smthpoints[indeces[i]].gettimestamp()

        vx = velocities_x[i]
        vy = velocities_y[i]
#       vx2 = velocities_x2[i]
#       vy2 = velocities_y2[i]

        pvx = math.fabs(vx)
        pvy = math.fabs(vy)

        pvxr = math.pow(vx/eta_x,2) if eta_x > 0.000001 else 0.0
        pvyr = math.pow(vy/eta_y,2) if eta_y > 0.000001 else 0.0

        dx = self.magnitude[indeces[i]].at(0)
        dy = self.magnitude[indeces[i]].at(1)
        pmag = math.sqrt(dx*dx + dy*dy)
        minmag = pmag

        # compute running mean
        # from Brown, Robert Grover, "Introduction to Random Signal
        #   Analysis and Kalman Filtering", John Wiley & Sons, New York, NY
        #   1983 [p.182] TK5102.5 .B696
        umag = 0.0
        pmag_k = 0
        umag = float(pmag_k)/float(pmag_k+1) * umag + 1.0/float(pmag_k+1) * pmag
        pmag_k += 1

        pvel = math.sqrt(vx*vx + vy*vy)
        uvel = 0.0
        pvel_k = 0
        uvel = float(pvel_k)/float(pvel_k+1) * uvel + 1.0/float(pvel_k+1) * pvel
        pvel_k += 1

#       if pvx > eta_x or pvy > eta_y:
        if pvxr + pvyr > 1.0:

          # check consecutive lengths of points: need to find edges, not points
          j = i + 1

          while j < len(indeces):

            nvx = velocities_x[j]
            nvy = velocities_y[j]
#           nvx2 = velocities_x2[j]
#           nvy2 = velocities_y2[j]

            if math.fabs(nvx) > pvx:
              pvx = math.fabs(nvx)
            if math.fabs(nvy) > pvy:
              pvy = math.fabs(nvy)

            npvxr = math.pow(nvx/eta_x,2) if eta_x > 0.00001 else 0.0
            npvyr = math.pow(nvy/eta_y,2) if eta_y > 0.00001 else 0.0

            # keep going if we're over threshold, else break out of while loop
#           if math.fabs(nvx) > eta_x or math.fabs(nvy) > eta_y:
            if npvxr + npvyr > 1.0:
              j = j + 1

              dx = self.magnitude[indeces[j]].at(0)
              dy = self.magnitude[indeces[j]].at(1)
              nmag = math.sqrt(dx*dx + dy*dy)

              # mean magnitude
              # compute running mean
              # from Brown, Robert Grover, "Introduction to Random Signal
              #   Analysis and Kalman Filtering", Wiley & Sons, New York, NY
              #   1983 [p.182] TK5102.5 .B696
              umag = float(pmag_k)/float(pmag_k+1) * umag + \
                 1.0/float(pmag_k+1) * nmag
              pmag_k += 1

              nvel = math.sqrt(nvx*nvx + nvy*nvy)
              # mean velocity
              uvel = float(pvel_k)/float(pvel_k+1) * uvel + \
                 1.0/float(pvel_k+1) * nvel
              pvel_k += 1

              # peak magnitude
              if nmag > pmag:
                pmag = nmag
              if nmag < minmag:
                minmag = nmag

            else:
              break

          # in case we were at the last point?
          if j == len(indeces):
            j = j - 1

          # looking for runs with min. duration of 6 ms (2 ms sampling period)
#         if j - i > 2:
          if 2 < j - i:
#         if 2 < j - i and j - i < 6:

            # mean magnitude computed above
            mag = umag
#           mag = pmag
#           mag = minmag
#           print "mag min, avg, max: ", minmag, umag, pmag

            # get peak velocity from above
            velx = pvx
            vely = pvy

            # saccade amplitude
#           amp = math.sqrt(velx*velx + vely*vely)
#           amp = (math.fabs(velx) + math.fabs(vely))/2.0
#           amp = max(math.fabs(velx),math.fabs(vely))
            amp = max(velx,vely)

            dur = self.smthpoints[indeces[j]].gettimestamp() - \
                  self.smthpoints[indeces[i]].gettimestamp()
            # this is really MT = A/V_p
#           dur = umag / uvel
#           dur = umag / amp
#           dur = pmag / uvel

            # magnitude should be < 2 deg
#           if 1.0/60.0 < mag and mag < 2.0 and fdeg_away < 3.0:
            if 1.0/60.0 < mag and mag < 2.0:

              print("MICROSACCADE: vel (%f,%f), amp (%f), mag (%f), dur (%f)" \
                     % (velx, vely, amp, mag, dur))

#             print "mag less than 2 deg: %f" % (mag)
#             print "amp (deg/s): %f" % (amp)
              str = "%f %f %f %f %f %f\n" % (t, x, y, mag, amp, dur)
              outfile.write(str)

              ctr += 1

          # end of inner while; advance to next point in sequence past last one
#         i = j + 1
          # advance by 20 ms to not count overshoots, will differ with Hz rate
#         i = j + 10
          samples = 20.0/(period * 1000.0)
          i = int(j + samples)

        # if ith point wasn't a microsaccade, advance by one
        else:
          i = i + 1

      if tt > 0.0:
#       msac_rate = float(ctr)/tt
        msac_rate = float(ctr)
        self.msac_rates[k] = msac_rate
#       print "total micorsaccade count: ", ctr
#       print "microsaccade rate: ", msac_rate

    outfile.close()

#   print "fixpoints_vx:"
#   print self.fixpoints_vx

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  def microsaccade_rates(self,fileName,w,h,screen,dist,herz):
    outfile = open(fileName,'w')

    print("Found: %d fixations" % (len(self.fixations)))
    print("Found: %d fixation rate records" % (len(self.msac_rates)))

    # msac_rates[k] is now the number of microsaccades per fixation

    tt = 0.0
    nmsr = 0.0
    # compute total scanpath fixation duration
    for k in range(len(self.fixations)):
      tt += self.fixations[k].getDuration()
      nmsr += self.msac_rates[k]

    msrt = nmsr/tt if tt > 0.0 else 0.0

    for k in range(len(self.msac_rates)):
      x = self.fixations[k].at(0) * float(w)
      y = self.fixations[k].at(1) * float(h)
      t = self.fixations[k].gettimestamp()
      dur = self.fixations[k].getDuration()
      str = "%f %f %f %f %f %f\n" % (t, x, y, dur, self.msac_rates[k], msrt)
      outfile.write(str)

    outfile.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  def microsaccades_sg(self,fileName,w,h,screen,dist,herz):

    outfile = open(fileName,'w')

    print("Found: %d fixations" % (len(self.fixations)))
    print("Fixation point list length: %d" % len(self.fixpoints))

    # sampling period in s
    period = float(1.0/float(herz))

    # 1/6 is from Engbert and Kliegl (they use a 6-tap filter)
    dt = period * 6.0

    r = math.sqrt(float(w)*float(w) + float(h)*float(h))
    dpi = r/float(screen)

    D = float(dist)
#   fov = 2*math.degrees(math.atan(math.radians(screen/(2*D))))
    fov = 2*math.degrees(math.atan2(screen,2*D))
    fovx = 2*math.degrees(math.atan2(float(w)/dpi,2*D))
    fovy = 2*math.degrees(math.atan2(float(h)/dpi,2*D))

    velocities_x = []
    velocities_y = []
    velocities_x2 = []
    velocities_y2 = []

    # for each fixation (set of raw/smoothed gaze points)
    for k in range(len(self.fixpoints)):

      indeces = self.fixpoints[k]
      print("Fixation %d: " % (k))
#     print indeces

      # nuke vel lists and initialize to length of this fixation (no. of points)
      del velocities_x[:]
      del velocities_y[:]
      velocities_x = [0.0]*len(indeces)
      velocities_y = [0.0]*len(indeces)
      del velocities_x2[:]
      del velocities_y2[:]
      velocities_x2 = [0.0]*len(indeces)
      velocities_y2 = [0.0]*len(indeces)

      # Step 1: compute velocities
      for i in range(len(indeces)):

        print("index: %d" % (indeces[i]))

        # Step 1: obtain corresponding velocity (deg/s)
        velx = self.velocity[i].at(0)
        vely = self.velocity[i].at(1)

        print("velocity as per Savitzky-Golay: (%f, %f)" % (velx,vely))

        velocities_x[i] = velx
        velocities_y[i] = vely
        velocities_x2[i] = velx*velx
        velocities_y2[i] = vely*vely

      print("Velocities x: ", velocities_x)
      print("Velocities y: ", velocities_y)

      # Step 2: medians of velocities squared and square median of velocities
#     median_vx = sorted(velocities_x)[len(velocities_x)/2]
#     median_vy = sorted(velocities_y)[len(velocities_y)/2]
#     median_vx2 = sorted(velocities_x2)[len(velocities_x2)/2]
#     median_vy2 = sorted(velocities_y2)[len(velocities_y2)/2]
      median_vx = np.median(velocities_x)
      median_vy = np.median(velocities_y)
      median_vx2 = np.median(velocities_x2)
      median_vy2 = np.median(velocities_y2)

      print("Median velocity x: ", median_vx)
      print("Median velocity y: ", median_vy)
      print("Median velocity x2: ", median_vx2)
      print("Median velocity y2: ", median_vy2)

      # compute median estimator
      # \sigma_{xy} = \langle v_{xy}^{2} \rangle - \langle v_{xy} \rangle^{2}
      sigma_x = median_vx2 - (median_vx*median_vx)
      sigma_y = median_vy2 - (median_vy*median_vy)

      print("sigma_x: ", sigma_x)
      print("sigma_y: ", sigma_y)

      # compute detection threshold (with \lambda = 6)
      # \eta_{xy} = \lamba \sigma_{xy}
      eta_x = 6.0 * sigma_x
      eta_y = 6.0 * sigma_y

      print("eta_x: ", eta_x)
      print("eta_y: ", eta_y)

      # Step 3 & 4: skip step 3 (binocular check); iterate to find microsaccades
      for i in range(len(indeces)):

        # (smoothed) data we are classifying
        x = self.smthpoints[indeces[i]].at(0)
        y = self.smthpoints[indeces[i]].at(1)
        t = self.smthpoints[indeces[i]].gettimestamp()

        vx = velocities_x[i]
        vy = velocities_y[i]

        if vx > eta_x and vy > eta_y:

          # saccade amplitude
          amp = math.sqrt(vx*vx + vy*vy)

          # corresponding magnitude (deg)
          dx = self.magnitude[indeces[i]].at(0)
          dy = self.magnitude[indeces[i]].at(1)
          mag = math.sqrt(dx*dx + dy*dy)

          dur = self.smthpoints[indeces[j]].gettimestamp() - \
                self.smthpoints[indeces[i]].gettimestamp()

          print("MICROSACCADE: vel (%f, %f), amp (%f), mag (%f), dur (%f)" % \
                 (vx, vy, amp, mag, dur))

          # magnitude should be < 2 deg
          if mag < 2.0:
#           print "mag less than 2 deg: %f" % (mag)
#           print "amp (deg/s): %f" % (amp)
            str = "%f %f %f %f %f %f\n" % (t, x, y, mag, amp, dur)
            outfile.write(str)

    outfile.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  def saccades(self,fileName,w,h,screen,dist,herz):

    outfile = open(fileName,'w')

#   print "Found: %d saccades" % (len(self.sacpoints))
#   print "Saccade point list length: %d" % len(self.sacpoints)
#   print "Saccade points (should be list of lists):", self.sacpoints

    # sampling period in s
    period = float(1.0/float(herz))

    # 1/6 is from Engbert and Kliegl (they use a 6-tap filter)
    # note: this is wrong since it assumes a constant sampling rate,
    #       which is unrealistic; see code below where we obtain dt
    #       from the data itself
    dt = period * 6.0

    r = math.sqrt(float(w)*float(w) + float(h)*float(h))
    dpi = r/float(screen)

    D = float(dist)
#   fov = 2*math.degrees(math.atan(math.radians(screen/(2*D))))
    fov = 2*math.degrees(math.atan2(screen,2*D))
    fovx = 2*math.degrees(math.atan2(float(w)/dpi,2*D))
    fovy = 2*math.degrees(math.atan2(float(h)/dpi,2*D))

    velocities_x = []
    velocities_y = []
    velocities_x2 = []
    velocities_y2 = []

    magnitudes_x = []
    magnitudes_y = []

    # nuke vel lists and initialize to length of this fixation (no. of points)
    # (only do this if we don't want to save velocity info)
#   del velocities_x[:]
#   del velocities_y[:]
    velocities_x = [0.0]*len(self.smthpoints)
    velocities_y = [0.0]*len(self.smthpoints)

    velocities_x2 = [0.0]*len(self.smthpoints)
    velocities_y2 = [0.0]*len(self.smthpoints)

    veldiff_x = [0.0]*len(self.smthpoints)
    veldiff_y = [0.0]*len(self.smthpoints)

    magnitudes_x = [0.0]*len(self.smthpoints)
    magnitudes_y = [0.0]*len(self.smthpoints)

    # need at least length 5 for Engbert and Kleigl's algorithm
    if len(self.smthpoints) >= 5:

      # Step 1: compute velocities
      for i in range(len(self.smthpoints)):

#       print "index: %d" % (sindeces[i])

        # compute velocity
        # v_i = \frac{x_{i+2} + x_{i+1} - x_{i-1} - x_{i-2}}{6 \Delta t}
        if 2 < i  and i < len(self.smthpoints) - 2:
          dx =  1.0 * self.smthpoints[i+2].at(0) + \
                1.0 * self.smthpoints[i+1].at(0) + \
               -1.0 * self.smthpoints[i-1].at(0) + \
               -1.0 * self.smthpoints[i-2].at(0)

          dy =  1.0 * self.smthpoints[i+2].at(1) + \
                1.0 * self.smthpoints[i+1].at(1) + \
               -1.0 * self.smthpoints[i-1].at(1) + \
               -1.0 * self.smthpoints[i-2].at(1)

          dt = self.smthpoints[i+2].gettimestamp() - \
               self.smthpoints[i-2].gettimestamp()

#         print "dt = ", dt, " should be: ", period * 6.0
#         if dt < period * 3.0:
#           print "*******************"
#           print "dt = ", dt, " should be: ", period * 6.0

          # distance covered in pixels (across diff filter window size)
          # normalize filter via scaling by 0.2?  Engbert and Kliegl don't
          dx = 1.0 * dx * float(w)
          dy = 1.0 * dy * float(h)

          # assume D is in inches
          degx = 2.0*math.degrees(math.atan2((dx/dpi),(2.0*D)))
          degy = 2.0*math.degrees(math.atan2((dy/dpi),(2.0*D)))

          # overwrite original magnitude with what we calculate above
          # the difference is largely in scale since the filters
          # used (6-tap vs. Savitzky-Golay) are different widths
          magnitudes_x[i] = degx
          magnitudes_y[i] = degy
#         self.magnitude[sindeces[i]].set(degx,degy)

          # degx, degy is deg per filter window, div by dt to get per second
          if dt < 0.0001:
            continue
          velx = degx / dt
          vely = degy / dt
#         print "velocity as per Engbert and Kliegl: (%f, %f)" % (velx,vely)

          velocities_x[i] = velx
          velocities_y[i] = vely
          velocities_x2[i] = velx*velx
          velocities_y2[i] = vely*vely

#     print "Velocities x: ", velocities_x
#     print "Velocities y: ", velocities_y

      # Step 2: medians of velocities squared and square median of velocities
      # NOTE: using Engbert et al. 2016, where the median computations
      #       differ from other microsaccade detection papers
      median_vx = np.median(velocities_x)
      median_vy = np.median(velocities_y)
#     median_vx2 = np.median(velocities_x2)
#     median_vy2 = np.median(velocities_y2)

      for i in range(len(self.smthpoints)):
        veldiff_x[i] = math.pow(velocities_x[i] - median_vx, 2)
        veldiff_y[i] = math.pow(velocities_y[i] - median_vy, 2)

      median_vx2 = np.median(veldiff_x)
      median_vy2 = np.median(veldiff_y)

#     print "Median velocity x: ", median_vx
#     print "Median velocity y: ", median_vy
#     print "Median velocity x2: ", median_vx2
#     print "Median velocity y2: ", median_vy2

      # compute median estimator
      # \sigma_{xy} = \langle v_{xy}^{2} \rangle - \langle v_{xy} \rangle^{2}
#     sigma_x = median_vx2 - (median_vx*median_vx)
#     sigma_y = median_vy2 - (median_vy*median_vy)
      sigma_x = median_vx2
      sigma_y = median_vy2

      if sigma_x < 0.00001:
        sigma_x = 0.0
      if sigma_y < 0.00001:
        sigma_y = 0.0

      # as per Engelbert 2006: use sqrt
      sigma_x = math.sqrt(sigma_x)
      sigma_y = math.sqrt(sigma_y)

#     print "sigma_x: ", sigma_x
#     print "sigma_y: ", sigma_y

      # compute detection threshold (with \lambda = 6)
      eta_x = 8.0 * sigma_x
      eta_y = 8.0 * sigma_y

#     print "eta_x: ", eta_x
#     print "eta_y: ", eta_y

      # Step 3: according to Engbert and Kliegl, if right eye msac starts/ends
      #         at time r1, r2, then left msac start/end times l1, l2 should
      #         be l1 < r2 and l2 > r1
      #         ok, this makes sense, but we are operating on averaged
      #         cyclopean (x,y) data, so we can't really implement this
      #         unless we store left and right eye data

      # Step 4: iterate to find saccades
      # loop through this set of saccade points and see how long the
      # saccade sequences are
      i = 0
      while i < len(self.smthpoints):

        # (smoothed) data we are classifying
        x = self.smthpoints[i].at(0)
        y = self.smthpoints[i].at(1)
        t = self.smthpoints[i].gettimestamp()

        vx = velocities_x[i]
        vy = velocities_y[i]
#       vx2 = velocities_x2[i]
#       vy2 = velocities_y2[i]

        pvx = math.fabs(vx)
        pvy = math.fabs(vy)

        pvxr = math.pow(vx/eta_x,2) if eta_x > 0.000001 else 0.0
        pvyr = math.pow(vy/eta_y,2) if eta_y > 0.000001 else 0.0

#       dx = self.magnitude[i].at(0)
#       dy = self.magnitude[i].at(1)
        dx = magnitudes_x[i]
        dy = magnitudes_y[i]
        pmag = math.sqrt(dx*dx + dy*dy)
        minmag = pmag

        # compute running mean
        # from Brown, Robert Grover, "Introduction to Random Signal
        #   Analysis and Kalman Filtering", John Wiley & Sons, New York, NY
        #   1983 [p.182] TK5102.5 .B696
        umag = 0.0
        pmag_k = 0
        umag = float(pmag_k)/float(pmag_k+1) * umag + 1.0/float(pmag_k+1) * pmag
        pmag_k += 1

        pvel = math.sqrt(vx*vx + vy*vy)
        uvel = 0.0
        pvel_k = 0
        uvel = float(pvel_k)/float(pvel_k+1) * uvel + 1.0/float(pvel_k+1) * pvel
        pvel_k += 1

#       if pvx > eta_x or pvy > eta_y:
        if pvxr + pvyr > 1.0:

          # check consecutive lengths of points: need to find edges, not points
          j = i + 1

          while j < len(self.smthpoints):

            nvx = velocities_x[j]
            nvy = velocities_y[j]
#           nvx2 = velocities_x2[j]
#           nvy2 = velocities_y2[j]

            if math.fabs(nvx) > pvx:
              pvx = math.fabs(nvx)
            if math.fabs(nvy) > pvy:
              pvy = math.fabs(nvy)

            npvxr = math.pow(nvx/eta_x,2) if eta_x > 0.00001 else 0.0
            npvyr = math.pow(nvy/eta_y,2) if eta_y > 0.00001 else 0.0

            # keep going if we're over threshold, else break out of while loop
#           if math.fabs(nvx) > eta_x or math.fabs(nvy) > eta_y:
            if npvxr + npvyr > 1.0:
              j = j + 1

#             dx = self.magnitude[j].at(0)
#             dy = self.magnitude[j].at(1)
              dx = magnitudes_x[j]
              dy = magnitudes_y[j]
              nmag = math.sqrt(dx*dx + dy*dy)

              # mean magnitude
              # compute running mean
              # from Brown, Robert Grover, "Introduction to Random Signal
              #   Analysis and Kalman Filtering", Wiley & Sons, New York, NY
              #   1983 [p.182] TK5102.5 .B696
              umag = float(pmag_k)/float(pmag_k+1) * umag + \
                 1.0/float(pmag_k+1) * nmag
              pmag_k += 1

              nvel = math.sqrt(nvx*nvx + nvy*nvy)
              # mean velocity
              uvel = float(pvel_k)/float(pvel_k+1) * uvel + \
                 1.0/float(pvel_k+1) * nvel
              pvel_k += 1

              # peak magnitude
              if nmag > pmag:
                pmag = nmag
              if nmag < minmag:
                minmag = nmag

            else:
              break

          # in case we were at the last point?
          if j == len(self.smthpoints):
            j = j - 1

          # looking for runs with min. duration of 6 ms (2 ms sampling period)
#         if j - i > 2:
          if 2 < j - i:
#         if 2 < j - i and j - i < 6:

            # mean magnitude computed above
            mag = umag
#           mag = pmag
#           mag = minmag
#           print "mag min, avg, max: ", minmag, umag, pmag

            # get peak velocity from above
            velx = pvx
            vely = pvy

            # saccade amplitude
#           amp = math.sqrt(velx*velx + vely*vely)
#           amp = (math.fabs(velx) + math.fabs(vely))/2.0
#           amp = max(math.fabs(velx),math.fabs(vely))
#           if velx > 720:
#             amp = vely
#           elif vely > 720:
#             amp = velx
#           else:
#             amp = max(velx,vely)
            amp = max(velx,vely)

            dur = self.smthpoints[j].gettimestamp() - \
                  self.smthpoints[i].gettimestamp()
            # this is really MT = A/V_p
#           dur = mag / amp
#           dur = mag / uvel
#           dur = pmag / uvel

            # NOTE: amp is really peak velocity
            #       mag is really amplitude

#           if 1.0/60.0 < mag and mag < 2.0:
#           if 2.0 < mag:
            if 1.0/60.0 < mag:
#           if 1.0/60.0 < mag and 0.0001 < amp:
#           if 1.0/60.0 < mag and 10 < amp and amp < 720:
            # over the micro range, mag > 1.0 deg
#           if 1.0 < mag and 0.0001 < amp:
              print("SACCADE: vel (%f,%f), amp (%f), mag (%f), dur (%f)" \
                     % (velx, vely, amp, mag, dur))

              # (smoothed) data we are classifying
#             print "amp (deg/s): %f" % (amp)
              str = "%f %f %f %f %f %f %f %f %f\n" % \
                     (t, x, y, dt, dx, dy, mag, amp, dur)
              outfile.write(str)

          # end of inner while; advance to next point in sequence past last one
          i = j + 1

        # if ith point wasn't a microsaccade, advance by one
        else:
          i = i + 1

    outfile.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  def saccades_sg(self,fileName,w,h,screen,dist,herz):

    outfile = open(fileName,'w')

    print("Found: %d saccades" % (len(self.sacpoints)))
#   print "Saccade point list length: %d" % len(self.sacpoints)
#   print "Saccade points (should be list of lists):", self.sacpoints

    # sampling period in s
    period = float(1.0/float(herz))

    # 1/6 is from Engbert and Kliegl (they use a 6-tap filter)
    # note: this is wrong since it assumes a constant sampling rate,
    #       which is unrealistic; see code below where we obtain dt
    #       from the data itself
    dt = period * 6.0

    r = math.sqrt(float(w)*float(w) + float(h)*float(h))
    dpi = r/float(screen)

    D = float(dist)
#   fov = 2*math.degrees(math.atan(math.radians(screen/(2*D))))
    fov = 2*math.degrees(math.atan2(screen,2*D))
    fovx = 2*math.degrees(math.atan2(float(w)/dpi,2*D))
    fovy = 2*math.degrees(math.atan2(float(h)/dpi,2*D))

    velocities_x = []
    velocities_y = []

    # for each saccade (set of raw/smoothed gaze points)
    for k in range(len(self.sacpoints)):

      # should be a list of lists
      sindeces = self.sacpoints[k]
#     print "Saccde %d: " % (k)
#     print sindeces

      # nuke vel lists and initialize to length of this fixation (no. of points)
      # (only do this if we don't want to save velocity info)
#     del velocities_x[:]
#     del velocities_y[:]
      velocities_x = [0.0]*len(sindeces)
      velocities_y = [0.0]*len(sindeces)

      # Step 1: compute velocities
      for i in range(len(sindeces)):
        velocities_x[i] = self.velocity[sindeces[i]].at(0)
        velocities_y[i] = self.velocity[sindeces[i]].at(1)

#     print "Velocities x: ", velocities_x
#     print "Velocities y: ", velocities_y

      dx = 1.0 * self.smthpoints[sindeces[0]].at(0) + \
          -1.0 * self.smthpoints[sindeces[-1]].at(0)
      dy = 1.0 * self.smthpoints[sindeces[0]].at(1) + \
          -1.0 * self.smthpoints[sindeces[-1]].at(1)
      dt = 1.0 * self.smthpoints[sindeces[0]].gettimestamp() + \
          -1.0 * self.smthpoints[sindeces[-1]].gettimestamp()

      # distance covered in pixels (across diff filter window size)
      # normalize filter via scaling by 0.2?  Engbert and Kliegl don't
      dx = 1.0 * dx * float(w)
      dy = 1.0 * dy * float(h)

      # assume D is in inches
      degx = 2.0*math.degrees(math.atan2((dx/dpi),(2.0*D)))
      degy = 2.0*math.degrees(math.atan2((dy/dpi),(2.0*D)))

#     mag = math.sqrt(dx*dx + dy*dy)
      mag = math.sqrt(degx*degx + degy*degy)

      pvx = 0.0
      pvy = 0.0

      # Step 2: loop through this set of saccade points to get peak velocity
      for i in range(len(sindeces)):

        vx = velocities_x[i]
        vy = velocities_y[i]

        if math.fabs(vx) > pvx:
          pvx = math.fabs(vx)
        if math.fabs(vy) > pvy:
          pvy = math.fabs(vy)

      # get peak velocity from above
      velx = pvx
      vely = pvy

      # saccade amplitude
#     amp = math.sqrt(velx*velx + vely*vely)
#     amp = (math.fabs(velx) + math.fabs(vely))/2.0
#     amp = max(math.fabs(velx),math.fabs(vely))
      amp = max(velx,vely)

#     if 1.0/60.0 < mag and mag < 2.0:
#     if 2.0 < mag:
#     if 1.0/60.0 < mag:
      if 1.0/60.0 < mag and 0.0001 < amp:
        print("SACCADE w/ vel (%f,%f), amp (%f), mag (%f)" \
               % (velx, vely, amp, mag))

        # (smoothed) data we are classifying
#       print "amp (deg/s): %f" % (amp)
        str = "%f %f %f %f %f\n" % (dt, dx, dy, mag, amp)
        outfile.write(str)

    outfile.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  def amfoc(self,fileName,w,h):

    self.K = [0.0]*(len(self.fixations))

    # compute mean fixation duration
    if(len(self.fixations) > 1):
      durations = [0.0]*len(self.fixations)
      for i in range(len(self.fixations)):
        durations[i] = self.fixations[i].getDuration()
      avedur = np.mean(durations)
      stddur = np.std(durations)
      mindur = min(durations)
      maxdur = max(durations)

      # compute mean amplitudes
      amps = [0.0]*(len(self.fixations)-1)
      for i in range(len(self.fixations)-1):
        dx = self.fixations[i].at(0) - self.fixations[i+1].at(0)
        dy = self.fixations[i].at(1) - self.fixations[i+1].at(1)
        dist = math.sqrt(dx*dx + dy*dy)
        amps[i] = dist
      aveamp = np.mean(amps)
      stdamp = np.std(amps)

      # compute ambient/focal coefficients
      for i in range(len(self.fixations)-1):
        if(stdamp < 0.00001):
          self.K[i] = 0.0
        else:
          self.K[i] = (self.fixations[i].getDuration() - avedur)/stddur - \
                      (amps[i] - aveamp)/stdamp

      # set K for the last entry (should be same as penultimate K)
      self.K[len(self.fixations)-1] = self.K[len(self.fixations)-2]

    elif(len(self.fixations) == 1):
      avedur = self.fixations[0].getDuration()
      stddur = 1.0
      mindur = avedur
      maxdur = avedur
    else:
      avedur = 0.0
      stddur = 0.0
      mindur = avedur
      maxdur = avedur

    outfile = open(fileName,'w')
    for i in range(len(self.fixations)):
      str = "%f %f\n" % (self.fixations[i].gettimestamp(), self.K[i])
      outfile.write(str)
    outfile.close()

  def gridify(self,fileName,subj,cond,w,h,xtiles,ytiles):

    # AOIs: xtiles x ytiles overlaid on the image
    grid = Grid(xtiles,ytiles)

#   grid.fill(self.fixations,w,h)

    outfile = open(fileName,'w')
    str = "subj,cond,AOI,duration,order\n"
    outfile.write(str)

#   this prints out the number of fixations per AOI
#   for i in range(xtiles):
#     for j in range(ytiles):
#       str = "%c%c %d\n" % (chr(int(i+97)), chr(int(j+97)), grid.at(i,j))
#       outfile.write(str)
    ct = 1
    for fx in self.fixations:
      # ix,iy give the integer coords of the AOI
      # use float(w)-1 and float(h)-1 to avoid points on the last AOI boundary
      ix = (fx.at(0) * (float(w)-1))//((int(w)//xtiles))
#     iy = (fx.at(1) * (float(h)-1))//((int(h)//ytiles))
      # do the y-coord flip for rendering with (0,0) at bottom-left
      iy = ((1.0 - fx.at(1)) * (float(h)-1))//((int(h)//ytiles))
      # get character labels (1-26 set)
#     lx = "%c" % (chr(int(ix+97)))
#     ly = "%c" % (chr(int(iy+97)))
      # get character labels (1-52 set; note that range is [0,1] not [1,n])
      lx = "%c" % (chr(int(ix+97))) if ix <= 25 else (chr(int((ix-26)+65)))
      ly = "%c" % (chr(int(iy+97))) if iy <= 25 else (chr(int((iy-26)+65)))
      # output AOI label, fixation duration, count
      str = "%s,%s,%c%c,%f,%d\n" % (subj, cond, lx, ly, fx.getDuration(), ct)
      outfile.write(str)
      ct += 1
    outfile.close()

  def entropy(self,fileName,subj,cond,w,h):

    outfile = open(fileName,'w')
    str = "subj,cond,entropy\n"
    outfile.write(str)

    # diagonal
    d = math.sqrt(float(w)*float(w) + float(h)*float(h))
    sigma = 0.0

    # heatmap: a 32-bit floating point image, initially set to black
    lum = Image.new("F", (int(w),int(h)), 0.0)
    pix = lum.load()

    for fx in self.fixations:
      x = fx.at(0) * float(w)
      # do the y-coord flip for rendering with (0,0) at bottom-left
      y = (1.0 - fx.at(1)) * float(h)
      # hack: if there's only 1 fixation @ 0% duration: 1/6th of image
      sigma = fx.getPercentDuration() * d if fx.getPercentDuration() > 0 else d/6.0
 #    for i in range(int(h)):
 #      for j in range(int(w)):
      for i in range(int(y-2.0*sigma),int(y+2.0*sigma)):
        for j in range(int(x-2.0*sigma),int(x+2.0*sigma)):
          if( 0 <= i and i < int(h) and 0 <= j and j < int(w) ):
            sx = j - x
            sy = i - y
            heat = math.exp((sx*sx +sy*sy)/(-2.0*sigma*sigma))
            pix[j,i] = pix[j,i] + heat

    # get max value
    minlum, maxlum = lum.getextrema()

    # normalize
    lum = lum.point(lambda f: f * (1.0/maxlum) + 0)
    pix = lum.load()

    # compute lum entropy
    gray_freq = [0.0]*256;
    for i in range(int(h)):
      for j in range(int(w)):
        idx = int(pix[j,i] * 255.0)
 #      print "pix[%d,%d] = " % (j,i)
 #      print "%f, idx = %d" % (pix[j,i], idx)
        gray_freq[idx] = gray_freq[idx] + 1

    prob = 0.0
    entropy = 0.0
    for i in range(256):
      if(gray_freq[i] == 0):
        pass
      else:
        prob = float(gray_freq[i]/(float(w)*float(h)))
        if(prob > 0.0):
          entropy = entropy - (prob * np.log2(prob))

    # normalize entropy since max entropy is 2^8
    entropy = entropy / np.log2(256.0)

    # print entropy
    str = "%s,%s,%f\n" % (subj, cond, entropy)
    outfile.write(str)
    outfile.close()

  def dumpFixatedAOIs(self,fileName,w,h,aoilist):

    outfile = open(fileName,'w')
    # aoi_span: distance between this AOI and previous
    # aoi_order: same as what aoi_label used to be (numerical order of AOI)
    str = "timestamp,x,y,duration,prev_sacc_amplitude,aoi_span,aoi_label,aoi_order\n"
    outfile.write(str)
    i=0
    span=0
    prev_fixation=0
    prev_aoi=0
    for fx in self.fixations:
      x = fx.at(0) * float(w)
      y = fx.at(1) * float(h)
      inAOI = False
      for aoi in aoilist:
        if aoi.inside(x,y + aoi.getHeight()):
          inAOI = True
          break
        if inAOI:
          # compute the AOI "span": this is a numerical distance of
          # where we came from, e.g., if 0 it's a refixation, if it's < 0
          # we've regressed (gone left in the sentence) an dif it's > 0
          # then we've skipping ahead
          span = 0 if prev_aoi == 0 else int(aoi.getLocation()) - prev_aoi
          prev_aoi = int(aoi.getLocation())
          prev_sacc_amp = 0.0
          if prev_fixation > 0:
            sx = self.fixations[prev_fixation].at(0) * float(w)
            sy = self.fixations[prev_fixation].at(1) * float(h)
            dx = x - sx
            dy = y - sy
            prev_sacc_amp = math.sqrt(float(dx)*float(dx) + float(dy)*float(dy))
          str = "%f,%f,%f,%f,%f,%s,%s\n" % \
          (fx.gettimestamp(),x,y,fx.getDuration(),prev_sacc_amp,span,aoi.getAttributes_astring())
          outfile.write(str)
          prev_fixation=i
        i += 1
    outfile.close()

  def clamp(self, x, y):
    if x < 0:
      x = 0.0
    if y < 0:
      y = 0.0
    if x > 1:
      x = 1.0
    if y > 1:
      y = 1.0
    return x,y
