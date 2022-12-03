#!/usr/bin/env python3

# was using /sw/bin/python2.7 -tt
# filter.py
# this file is the driver for the python analysis script

# system includes
import sys
import os
import getopt
import glob
#import xml.etree.ElementTree as ET
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import tostring
import numpy as np
import math

# local includes
from scanpath import Scanpath
from aoi import AOI
from lagrange import Lagrange
from monitor import Monitor

import plotter
import point
from point import Point

from scipy import spatial 

def usage():
  print("Usage: python graph.py " \
        " --width=? --height=? --screen=? --dist=?" \
        " --xtiles=? --ytiles=?\n" \
        " --indir=? --imgdir=? --outdir=? --pltdir=?" \
        " --file=?\n"\
        " --hertz=?\n"\
        " --sfdegree=? --sfcutoff=?\n"\
        " --dfdegree=? --dfwidth=?\n"\
        " --vt=?\n" \
        " --baselineT=? --endT=?\n" \
        " --smooth=?\n"\
        " --proximity=?\n"\
        "   width, height: screen dimensions (in pixels)\n" \
        "   screen: screen diagonal dimension (in inches)\n" \
        "   dist: viewing distance (in inches)\n" \
        "   xtiles: number of AOI tiles (in x direction)\n" \
        "   ytiles: number of AOI tiles (in y direction)\n" \
        "   indir: a directory containing input files to be processed\n" \
        "   imgdir: a directory containing image files\n" \
        "   outdir: a directory containing output files\n" \
        "   pltdir: a directory containing plot files\n" \
        "   file: a single file to process\n" \
        "   image: a single image to use as background\n" \
        "   hertz: the sampling rate of the data\n" \
        "   sfdegree: butterworth filter smoothing degree \n" \
        "   sfcutoff: butterworth filter smoothing cutoff \n" \
        "   dfwidth: savitzky-golay filter width \n" \
        "   vt: min velocity for change to be a saccade\n"\
        "   baselineT: time for pupil diameter baseline\n"\
        "   endT: max time to average pupil diameter difference\n"\
        "   smooth: True enables butterworth smoothing\n"\
        "   proximity: fixation proximity check\n")

def parseAOI(aoidir,aoifile,aoilist):

  print("parsing: ", aoidir + '/' + aoifile)

  if os.path.isfile(aoidir + '/' + aoifile):
    pass
  else:
    return

  tree = ET.parse(aoidir + '/' + aoifile)
  #print("tree = %s" % tree)

  # should be something like Experiment, {}
  root = tree.getroot()
# print("root.tag = %s, root.attrib = %s" % (root.tag, root.attrib))

  # iterate through PAGEOBJECT objects, look for PTYPE=''6'
  print("PAGE INFO:")
  for obj in root.iter('MASTERPAGE'):
    px = obj.get('PAGEXPOS')
    py = obj.get('PAGEYPOS')
    pw = obj.get('PAGEWIDTH')
    ph = obj.get('PAGEHEIGHT')
    bl = obj.get('BORDERLEFT')
    br = obj.get('BORDERRIGHT')
    bt = obj.get('BORDERTOP')
    bb = obj.get('BORDERBOTTOM')
#   print("%s %s %s %s" % (px,py,pw,ph))
#   print("%s %s %s %s" % (bl,br,bt,bb))

  # iterate through PAGEOBJECT objects, look for PTYPE=''6'
  for obj in root.iter('PAGEOBJECT'):
    if obj.get('PTYPE') == '6':
      # x,y is the top-left corner
      x = obj.get('XPOS')
      y = obj.get('YPOS')
      w = obj.get('WIDTH')
      h = obj.get('HEIGHT')
# gone in v1.5.6
#     label = obj.get('ANNAME')
#     print("AOI: %s: %s %s %s %s" % (label,x,y,w,h))

# new in v1.5.6: attributes
#
# 'Name'     : user-specified
# 'Type'     : None, Boolean, Integer, Real Number, String
# 'Value'    : user-specified
      aoiattr = {}
      # <PageItemAttributes>
      pattr = obj.find('PageItemAttributes')
#     print("len(pattr): ", len(pattr))
#     print("pattr: ", pattr)
      for attr in pattr:
        # <ItemAttribute Name="????????" Type="Integer" Value="?" ... />
        name = attr.get('Name')
        type = attr.get('Type')
        value = attr.get('Value')
#       print("name: ",name)
#       print("type: ",type)
#       print("value: ",value)
        aoiattr[name] = value
      # </PageItemAttributes>

      # need to subtract the pagexpos and pageypos and need to add in height
      # to get bottom left corner
      fx = float(x) - float(px)
      fy = (float(y) - float(py)) + float(h)
      fw = float(w)
      fh = float(h)
      aoi = AOI(fx,fy,fw,fh)
      aoi.setAttributes(aoiattr)
      aoi.dump()
#     if aoidict.has_key(stimulus):
#       aoidict[stimulus].append(aoi)
#     else:
#       aoidict[stimulus] = [aoi]
      aoilist.append(aoi)

      del aoi

def main(argv):
# if not len(argv):
#   usage()

  try:
    opts, args = getopt.getopt(argv, '', \
                 ['width=','height=',\
                  'hertz=',\
                  'screen=','dist=',\
                  'xtiles=','ytiles=',\
                  'indir=','imgdir=','aoidir=','outdir=','pltdir=',\
                  'file=','aoifile=','image=', 'hertz=',\
                  'sfdegree=', 'sfcutoff=',\
                  'dfdegree=', 'dfwidth=',\
                  'baselineT=','endT=',\
                  'vt=', 'smooth=', 'proximity='])
  except getopt.GetoptError as err:
    usage()
    exit()

  # Enable/disable fixation proximity checking.
  proximity = False
  # Enable/disable butterworth smoothing.
  smooth = False
  # screen height in pixels
  width = 1920
  height = 1200
  # screen diagonal (in inches)
  screen = 22
  # viewing distance (in inches)
  dist = 21.65
  # sampling rate
  herz = 1200.0
  # smoothing (Butterworth) filter parameters: degree, cutoff
  sfdegree = 3
# sfcutoff = 0.5/herz # must be 0 < Wn < 1
  sfcutoff = 1.0/herz # must be 0 < Wn < 1
# sfcutoff = 1.0/herz # must be 0 < Wn < 1
# sfcutoff = 20.0/herz # must be 0 < Wn < 1
  # differentiation (SG) filter parameters: width, degree, order
  # 5, 3 for no smoothing, 3, 2 for smoothing
  if (smooth):
    dfwidth = 3
    dfdegree = 2
  else:
    dfwidth = 5
    dfdegree = 3
  dfo = 1
  # velocity threshold
# T = 5.0  # more fixations
# T = 7.0  # more fixations
# T = 10.0
# T = 20.0
# T = 25.0
# T = 30.0
# T = 35.0 # fewer fixations
  T = 36.0 # fewer fixations
# T = 40.0 # 
# T = 80.0 # last used
# T = 100.0 #
# T = 120.0 #
# T = 150.0 #
# T = 240.0 #

  baselineT = 2.0
  endT = 200.0

  file = None
  # initially don't use an image (just plain white bgnd for plots)
  image = None

  # set up AOI grid
  xtiles = 4
  ytiles = 3

  outdir = "./data"
  pltdir = "./plots"
  aoidir = "../../../../aois/"
  aoifile = "aoidefinition.sla"


  # Checked beforehand so that custom parameters could still be used...
  # Check if smooth is an option. We will set default parameters based on 
  # value. If others are provided via the command line, we will use them.
  try:
      arg = opts[[t[0] for t in opts].index('--smooth')][1]
      if arg.lower() == 'true':
          smooth = True
      elif arg.lower() == 'false':
          smooth = False
      else:
          print("Warning, invalid smooth value. Assuming default.")
      if (smooth):
        dfwidth = 3
        dfdegree = 2
      else:
        dfwidth = 5
        dfdegree = 3
  except Exception as e:
    print(e)
    sys.exit()
    pass

# try:
#     arg = opts[[t[0] for t in opts].index('--proximity')][1]
#     if arg.lower() == 'true':
#         proximity = True
#     elif arg.lower() == 'false':
#         proximity = False
#     else:
#         print("Warning, invalid proximity value. Assuming default.")
# except Exception as e:
#   print(e)
#   sys.exit()
#   pass

  dfwidth = 5
  dfdegree = 3

  for opt,arg in opts:
    opt = opt.lower()
    if(opt != '--file' and opt != '--image'):
      arg = arg.lower()

    if opt == '--indir':
      indir = arg
    elif opt == '--imgdir':
      imgdir = arg
    elif opt == '--outdir':
      outdir = arg
    elif opt == '--aoidir':
      aoidir = arg
    elif opt == '--aoifile':
      aoifile = arg
    elif opt == '--pltdir':
      pltdir = arg
    elif opt == '--width':
      width = arg
    elif opt == '--height':
      height = arg
    elif opt == '--screen':
      screen = float(arg)
    elif opt == '--dist':
      dist = float(arg)
    elif opt == '--xtiles':
      xtiles = int(arg)
    elif opt == '--ytiles':
      ytiles = int(arg)
    elif opt == '--file':
      file = arg
    elif opt == '--image':
      image = arg
    elif opt == '--hertz':
      herz = float(arg)
    elif opt == '--sfdegree':
      sfdegree = float(arg)
    elif opt == '--sfcutoff':
      sfcutoff = float(arg)
    elif opt == '--dfdegree':
      dfdegree = int(arg)
    elif opt == '--dfwidth':
      dfwidth = int(arg)
    elif opt == '--vt':
      T = float(arg)
    elif opt == '--baselinet':
      baselineT = float(arg)
    elif opt == '--endt':
      endT = float(arg)
      print("endT = ", endT)
    else:
      sys.argv[1:]

  basescanpath = None

  # get .raw input files to process
  if os.path.isdir(indir):
    files = glob.glob('%s/*.csv' % (indir))

  # if user specified --file="..." then we use that as the only one to process
  if file != None:
    file = indir + file
    if os.path.isfile(file):
      files = [file]
      print("overriding files with: ", files)

  # check to see if use specified image, if so, then don't compose filename
  # every time through loop; same image will be used (calibration image
  # presumably)...kind of a hack, but WTF...
  haveimage = False
# print("checking for imgdir: ", imgdir)
  if(image != None and os.path.isdir(imgdir)):
    image = os.path.join(imgdir,image)
#   print("checking for image: ", image)
    haveimage = True

  lagrange = Lagrange(int(width),int(height),screen,dist)
  monitor = Monitor(int(width),int(height),screen,dist)

  # process AOI file
  aoidict = {}
  aoilist = []

  parseAOI(aoidir,aoifile,aoilist)

  for file in files:

    # need to clear out matrices for each scanpath
    lagrange.reset()

    # don't process empty files
    if os.path.getsize(file) == 0:
      continue

#   base = os.path.basename(file)
    path, base = os.path.split(file)

    print("Processing: ", file, "[", base, "]")

    # split filename from extension
    filename, ext = os.path.splitext(base)

    # construct basescanpath that holds pupil diameter to be averaged later
#   base = os.path.basename(file)
#   path, ext = os.path.splitext(file)
    print("path, base, filename, ext: ", path, base, filename, ext)
    subj = filename.split('-')[0]
    exp_id = filename.split('-')[1]
    ses_id = filename.split('-')[2]
    marker = filename.split('-')[3]
    object = filename.split('-')[4]
    print("subj, exp_id, ses_id, marker, object: ", \
           subj, exp_id, ses_id, marker, object)

    # here would have to assemble file name 'shot' to be museum_signage_?

#   if len(aoilist):
#     aoilist.clear()
#   aoifile = aoidir + shot + ".sla"
#   parseAOI(aoidir,aoifile,aoilist)

    haveimage = False
    imagebase = 'white-1920x1080'

    # create filename of corresponding image
    if(haveimage == True):
      image = '{0}.jpg'.format(os.path.join("../../stimuli/",shot))
      print("Image: ", image)
    else:
      image = '{0}.png'.format(os.path.join("../../stimuli/",imagebase))
      print("Image: ", image, "[", imagebase, "]")

    scanpath = Scanpath()
    scanpath.parseFile(file,width,height,screen,dist,herz)

#   scanpath.pdwt("%s/%s-pdwt%s" % (outdir,filename,".dat"),\
#                     width,height,herz,sfdegree,sfcutoff)
    scanpath.pdwtlevel("%s/%s-pdwt%s" % (outdir,filename,".dat"),\
                    width,height,herz,sfdegree,sfcutoff,2)
    # scanpath.pdwtLH might override lof, hif
    lof, hif = 4, 1
    lof, hif = scanpath.pdwtLH("%s/%s-pdwtLH%s" % (outdir,filename,".dat"),\
                    width,height,herz,sfdegree,sfcutoff,lof,hif)

    scanpath.psmooth("%s/%s-pups%s" % (outdir,filename,".dat"),\
                    width,height,herz,sfdegree,sfcutoff)

#   dest = "%s/%s-pcpd%s" % (outdir,filename,".dat")
#   if not os.path.isfile(dest):
    scanpath.getpcpd("%s/%s-pcpd%s" % (outdir,filename,".dat"),\
                  width,height,baselineT,endT)

    scanpath.smooth("%s/%s-smth%s" % (outdir,filename,".dat"),\
                    width,height,herz,sfdegree,sfcutoff,smooth)
    scanpath.differentiate("%s/%s-diff%s" % (outdir,filename,".dat"),\
                            width,height,screen,dist,herz,dfwidth,dfdegree,dfo)
    scanpath.threshold("%s/%s-fxtn%s" % (outdir,filename,".dat"),\
                            width,height,T,monitor,proximity)
#   scanpath.microsaccades("%s/%s-msac%s" % (outdir,filename,".dat"),\
#                           width,height,screen,dist,herz)
    scanpath.amfoc("%s/%s-amfo%s" % (outdir,filename,".dat"),\
                            width,height)
#   scanpath.gridify("%s/%s-aois%s" % (outdir,filename,".csv"),\
#                           subj,cond,width,height,xtiles,ytiles)

    scanpath.dumpDAT("%s/%s%s" % (outdir,filename,".dat"),width,height)
#   scanpath.dumpXML("%s/%s%s" % (outdir,filename,".xml"),width,height)

#   plotter.renderPoints1D("%s/%s-%s" % (pltdir,filename,"gzpt"),\
#                          width,height,\
#                          scanpath.gazepoints,'x',\
#                          "Gaze point data")

#   plotter.renderPoints1D("%s/%s-%s" % (pltdir,filename,"gzpt"),\
#                          width,height,\
#                          scanpath.gazepoints,'y',\
#                          "Gaze point data")

    plotter.renderPupil1D("%s/%s-%s" % (pltdir,filename,"pupd"),\
                           width,height,\
                           scanpath.pupild,'d',\
                           "Raw pupil diameter")

    plotter.renderPupil1D("%s/%s-%s" % (pltdir,filename,"pups"),\
                           width,height,\
                           scanpath.pupils,'d',\
                           "Smoothed pupil diameter")

    plotter.renderPupil1D("%s/%s-%s" % (pltdir,filename,"pupc"),\
                           width,height,\
                           scanpath.pupils,'f',\
                           "Smoothed PFE pupil diameter")

    plotter.renderPupil1D("%s/%s-%s" % (pltdir,filename,"pcpd"),\
                           width,height,\
                           scanpath.pcpd,'p',\
                           "CPD")

    plotter.renderPupil1D("%s/%s-%s" % (pltdir,filename,"bpcpd"),\
                           width,height,\
                           scanpath.bpcpd,'p',\
                           "BCPD")

    plotter.renderPupil1D("%s/%s-%s" % (pltdir,filename,"cpcpd"),\
                           width,height,\
                           scanpath.pcpd,'c',\
                           "CCPD")

    plotter.renderPupil1D("%s/%s-%s" % (pltdir,filename,"cbpcpd"),\
                           width,height,\
                           scanpath.bpcpd,'c',\
                           "CBCPD")

#   plotter.renderPupilDWT("%s/%s-%s" % (pltdir,filename,"pDWT"),\
#                          width,height,\
#                          scanpath,'w',\
#                          "DWT wavelet coefficients")

#   plotter.renderPupilDWT("%s/%s-%s" % (pltdir,filename,"pDWT"),\
#                          width,height,\
#                          scanpath,'a',\
#                          "DWT wavelet approximation")

#   plotter.renderPupilDWT("%s/%s-%s" % (pltdir,filename,"pDWT"),\
#                          width,height,\
#                          scanpath,'t',\
#                          "DWT coefficients with hard thresholding")

    plotter.renderPupilDWTlevel("%s/%s-%s" % (pltdir,filename,"pDWT"),\
                           width,height,\
                           scanpath,'w',\
                           "DWT wavelet coefficients",2)

    plotter.renderPupilDWTlevel("%s/%s-%s" % (pltdir,filename,"pDWT"),\
                           width,height,\
                           scanpath,'a',\
                           "DWT wavelet approximation",2)

    plotter.renderPupilDWTLH("%s/%s-%s" % (pltdir,filename,"pDWTLH"),\
                           width,height,\
                           scanpath,'t',\
                           lof,hif,\
                           "DWT LH with hard thresholding")

    plotter.renderPoints2D("%s/%s-%s" % (pltdir,filename,"gzpt"),\
                           width,height,\
                           scanpath.gazepoints,\
                           "Gaze point data",\
                           image,True,\
                           xtiles,ytiles)

#   plotter.renderPoints1D("%s/%s-%s" % (pltdir,filename,"smth"),\
#                          width,height,\
#                          scanpath.smthpoints,'x',\
#                          "Smoothed gaze point data")

#   plotter.renderPoints1D("%s/%s-%s" % (pltdir,filename,"smth"),\
#                          width,height,\
#                          scanpath.smthpoints,'y',\
#                          "Smoothed gaze point data")

#   plotter.renderPoints2D("%s/%s-%s" % (pltdir,filename,"smth"),\
#                          width,height,\
#                          scanpath.smthpoints,\
#                          "Smoothed gaze point data",
#                          image,True,\
#                          xtiles,ytiles)

#   plotter.renderPoints1D("%s/%s-%s" % (pltdir,filename,"dfft"),\
#                         width,height,\
#                         scanpath.velocity,'x',\
#                         "Differentiated gaze point data",False)

#   plotter.renderPoints1D("%s/%s-%s" % (pltdir, filename, "accel"), \
#                          width, height, \
#                          scanpath.acceleration, 'x',\
#                          "Twice Differentiated gaze point data", False)
    
#    plotter.renderPoints1D("%s/%s-%s" % (pltdir,filename,"dfft"),\
#                          width,height,\
#                          scanpath.velocity,'y',\
#                          "Differentiated gaze point data",False)

#   plotter.renderPoints2D("%s/%s-%s" % (pltdir,filename,"dfft"),\
#                          width,height,\
#                          scanpath.velocity,\
#                          "Differentiated gaze point data",\
#                          None,False,\
#                          xtiles,ytiles)

    plotter.renderFixations("%s/%s-%s" % (pltdir,filename,"fxtn"),\
                            width,height,\
                            screen,dist,\
                            scanpath.fixations,\
                            "Fixations",\
                            image,\
                            lagrange,\
                            xtiles,ytiles)

#   plotter.renderKMicrosaccades("%s/%s-%s" % (pltdir,filename,"ksac"),\
#                          width,height,\
#                          screen,dist,herz,\
#                          scanpath.fixations,\
#                          scanpath.K,\
#                          scanpath.fixpoints,\
#                          scanpath.smthpoints,\
#                          scanpath.fixpoints_vx,\
#                          scanpath.fixpoints_vy,\
#                          scanpath.fixpoints_vx2,\
#                          scanpath.fixpoints_vy2,\
#                          scanpath.magnitude,\
#                          "Microsaccades",\
#                          image,True,\
#                          xtiles,ytiles)

#   plotter.renderMicrosaccades("%s/%s-%s" % (pltdir,filename,"msac"),\
#                          width,height,\
#                          screen,dist,herz,\
#                          scanpath.fixpoints,\
#                          scanpath.smthpoints,\
#                          scanpath.fixpoints_vx,\
#                          scanpath.fixpoints_vy,\
#                          scanpath.fixpoints_vx2,\
#                          scanpath.fixpoints_vy2,\
#                          scanpath.magnitude,\
#                          "Microsaccades",\
#                          image,True,\
#                          xtiles,ytiles)

    plotter.renderAOIs("%s/%s-%s" % (pltdir,filename,"aois"),\
                            width,height,\
                            aoilist,\
                            imagebase,\
                            "AOIs",\
                            image,\
                            xtiles,ytiles)

    plotter.renderAOIFixations("%s/%s-%s" % (pltdir,filename,"aoi-fxtn"),\
                            width,height,\
                            scanpath.fixations,\
                            aoilist,\
                            imagebase,\
                            "AOI Fixations",\
                            image,\
                            xtiles,ytiles)

    plotter.renderFixatedAOIs("%s/%s-%s" % (pltdir,filename,"fxtn-aoi"),\
                            width,height,\
                            scanpath.fixations,\
                            aoilist,\
                            imagebase,\
                            "AOI Fixations",\
                            image,\
                            xtiles,ytiles)

    plotter.renderAmfocFixations("%s/%s-%s" % (pltdir,filename,"affx"),\
                            width,height,\
                            scanpath.fixations,\
                            scanpath.K,\
                            "Ambient/Focal Fixations",\
                            image,\
                            xtiles,ytiles)

#   plotter.renderCalibFixations("%s/%s-%s" % (pltdir,filename,"clbf"),\
#                           width,height,\
#                           screen,dist,\
#                           scanpath.fixations,\
#                           "Mean Error and Fit of Fixated Calibration Locations",\
#                           image,\
#                           lagrange,\
#                           xtiles,ytiles)

#   plotter.renderHeatmap("%s/%s-%s" % (pltdir,filename,"heat"),\
#                           width,height,\
#                           scanpath.fixations,\
#                           "Heatmap",\
#                           image,\
#                           xtiles,ytiles)

    print(" ")
    del scanpath

if __name__ == "__main__":
  main(sys.argv[1:])
