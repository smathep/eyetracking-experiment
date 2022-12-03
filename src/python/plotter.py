#!/usr/bin/env python3
# Scanpath.py
# This class encapsulates the analysis program

# system includes
import sys
import os
import math
import re
import os.path
#import xml.etree.ElementTree as ET
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import tostring
import numpy as np
import xmlpp
from scipy import spatial 
# use agg non-gui backend
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.path import Path
import matplotlib.patches as patches
from PIL import Image
import pylab

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['pdf.fonttype'] = 42
# usetex to embed fonts
#plt.rcParams['text.usetex'] = True
plt.rcParams['text.usetex'] = False
plt.rcParams['axes.unicode_minus'] = False
plt.rc('font',family='sans-serif')
#plt.rc('text.latex',preamble='\usepackage{sfmath}')
#plt.rc('text.latex',preamble='\usepackage{cmbright}')

# local includes
from point import Point
from fixation import Fixation
from aoi import AOI

def renderPupilDWTLH(baseName,w,h,scanpath,coord,lof,hif,title,scale=True):
  # plot
  fileName = baseName + '-' + coord + '.pdf'
  # clear plot
  plt.clf()
  # axes dimensions
  ax = plt.axes()
# ax.set_ylim(-.5,.5)

# print("(len(scanpath.cD_LH)): ", len(scanpath.cD_LH))

  px = []
  for i in range(len(scanpath.cD_LH)):
    px.append(scanpath.pupild[(2**lof)*i].gettimestamp())
#   print("(2**lof)*i",(2**lof)*i)

# for pt in scanpath.pupildwt:
#   t = pt.gettimestamp()
#   px.append(t)

  if coord == 'w':
    py = scanpath.cD_LH
  elif coord == 't':
    py = scanpath.cD_LHt

  # title and labels
  plt.title(title,size="28",y=1.04)
  plt.ylabel("DWT LH",family="sans-serif",size="28")
  plt.xlabel("Time (s)",family="sans-serif",size="28")
  plt.tick_params(axis='both',which='major',labelsize="22")
  plt.tick_params(axis='both',which='minor',labelsize="22")
  plt.grid(True,'major',ls='solid',alpha=.1)

  # lines
# f, ax = plt.subplots(3,sharex=True,sharey=True)

  opt = {'antialiased':True,\
         'alpha':.8,\
         'color':"#3F003F",\
         'lw':1.0}
  if coord == 't':
    ax.vlines(px,[0],py,**opt)
  else:
    ax.plot(px,py,**opt)
  ax.set_title(title,size="26",y=1.04)
  str = "level %d" % (4)
  ax.set_ylabel(str,size="22")
# ax.set_ylim(-.3,.3)

  # margins
  plt.tight_layout()

  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderPupilDWTlevel(baseName,w,h,scanpath,coord,title,dwtlvl,scale=True):
  # plot
  fileName = baseName + '-' + coord + '.pdf'
  # clear plot
  plt.clf()
  # axes dimensions
  ax = plt.axes()
# ax.set_ylim(-.5,.5)

  px = []
  for i in range(len(scanpath.cA)):
    px.append(scanpath.pupild[(2**dwtlvl)*i].gettimestamp())

# for pt in scanpath.pupildwt:
#   t = pt.gettimestamp()
#   px.append(t)

  if coord == 'w':
    py = scanpath.cD
  elif coord == 'a':
    py = scanpath.cA
  elif coord == 't':
    py = scanpath.cDt

  # title and labels
  plt.title(title,size="28",y=1.04)
  plt.ylabel("DWT",family="sans-serif",size="28")
  plt.xlabel("Time (s)",family="sans-serif",size="28")
  plt.tick_params(axis='both',which='major',labelsize="22")
  plt.tick_params(axis='both',which='minor',labelsize="22")
  plt.grid(True,'major',ls='solid',alpha=.1)

  # lines
# f, ax = plt.subplots(3,sharex=True,sharey=True)

  opt = {'antialiased':True,\
         'alpha':.8,\
         'color':"#3F003F",\
         'lw':1.0}
  if coord == 't':
    ax.vlines(px,[0],py,**opt)
  else:
    ax.plot(px,py,**opt)
  ax.set_title(title,size="26",y=1.04)
  str = "level %d" % (dwtlvl)
  ax.set_ylabel(str,size="22")
# ax.set_ylim(-.3,.3)

  # margins
  plt.tight_layout()

  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderPupilDWT(baseName,w,h,scanpath,coord,title,scale=True):
  # plot
  fileName = baseName + '-' + coord + '.pdf'
  # clear plot
  plt.clf()
  # axes dimensions
  ax = plt.axes()
# ax.set_ylim(-.5,.5)

  px1 = []
  px2 = []
  px3 = []
  for i in range(len(scanpath.cA3)):
    px3.append(scanpath.pupild[8*i].gettimestamp())
  for i in range(len(scanpath.cA2)):
    px2.append(scanpath.pupild[4*i].gettimestamp())
  for i in range(len(scanpath.cA1)):
    px1.append(scanpath.pupild[2*i].gettimestamp())

# for pt in scanpath.pupildwt:
#   t = pt.gettimestamp()
#   px.append(t)

  if coord == 'w':
    py1 = scanpath.cD1
    py2 = scanpath.cD2
    py3 = scanpath.cD3
  elif coord == 'a':
    py1 = scanpath.cA1
    py2 = scanpath.cA2
    py3 = scanpath.cA3
  elif coord == 't':
    py1 = scanpath.cD1t
    py2 = scanpath.cD2t
    py3 = scanpath.cD3t

  # title and labels
  plt.title(title,size="28",y=1.04)
  plt.ylabel("DWT",family="sans-serif",size="28")
  plt.xlabel("Time (s)",family="sans-serif",size="28")
  plt.tick_params(axis='both',which='major',labelsize="22")
  plt.tick_params(axis='both',which='minor',labelsize="22")
  plt.grid(True,'major',ls='solid',alpha=.1)

  # lines
  f, ax = plt.subplots(3,sharex=True,sharey=True)

  opt = {'antialiased':True,\
         'alpha':.8,\
         'color':"#3F003F",\
         'lw':1.0}
  if coord == 't':
    ax[0].vlines(px1,[0],py1,**opt)
  else:
    ax[0].plot(px1,py1,**opt)
  ax[0].set_title(title,size="26",y=1.04)
  ax[0].set_ylabel("level 1",size="22")
# ax[0].set_ylim(-.3,.3)
  opt = {'antialiased':True,\
         'alpha':.8,\
         'color':"#003F3F",\
         'lw':1.0}
  if coord == 't':
    ax[1].vlines(px2,[0],py2,**opt)
  else:
    ax[1].plot(px2,py2,**opt)
  ax[1].set_ylabel("level 2",size="22")
# ax[1].set_ylim(-.3,.3)
  opt = {'antialiased':True,\
         'alpha':.8,\
         'color':"#3F3F00",\
         'lw':1.0}
  if coord == 't':
    ax[2].vlines(px3,[0],py3,**opt)
  else:
    ax[2].plot(px3,py3,**opt)
  ax[2].set_xlabel("Time (s)",size="22",y=1.04)
  ax[2].set_ylabel("level 3",size="22")
# ax[2].set_ylim(-.3,.3)

  # margins
  plt.tight_layout()

  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderPupil1D(baseName,w,h,points,coord,title,scale=True):
  # plot
  fileName = baseName + '-' + coord + '.pdf'
  # clear plot
  plt.clf()
  # axes dimensions
  ax = plt.axes()
# ax = plt.axes(aspect=1)
  if(coord == 'x' and scale):
    ax.set_ylim((0,int(w)))
  elif(coord == 'y' and scale):
    ax.set_ylim((0,int(h)))
  elif(coord == 'd' and scale):
#   ax.set_ylim((0,int(6)))
    pass
  elif(coord == 'p' and scale):
#   ax.set_ylim((-1,int(1)))
    pass
  # fill in data points
  px = []
  py = []
  for pt in points:
    t = pt.gettimestamp()
    if(coord == 'x'):
      x = pt.at(0) * float(w) if scale else pt.at(0)
    elif(coord == 'y'):
      x = pt.at(1) * float(h) if scale else pt.at(1)
    elif(coord == 'k'):
      x = pt.at(0)
    elif(coord == 'd'):
      x = pt.at(0)
    elif(coord == 'p'):
      x = pt.at(0)
    elif(coord == 'c'):
      x = pt.at(1)
    elif(coord == 'f'):
      x = pt.at(1)
    px.append(t)
    py.append(x)

  # lines
# opt = {'antialiased':True,\
#        'alpha':.6,\
#        'color':"#3F3F3F",\
#        'lw':1,\
#        'marker':"o",\
#        'markersize':0.5,\
#        'markeredgecolor': "#787878",\
#        'markeredgewidth':10}
# plt.plot(px,py,antialiased=True,alpha=.6,color="#787878")#, marker='o', markersize=2, linestyle='--')
  opt = {'antialiased':True,\
         'alpha':.8,\
         'color':"#3F3F3F",\
         'lw':2.0}
  plt.plot(px,py,**opt)
# line = plt.Line2D(px,py,**opt)
# ax.add_artist(line)
  if (coord == 'x'):
      labels = ['point{0}'.format(i) for i in range(len(px))]
      indexes = [i for i in range(len(px))]
  elif (coord == 'y'):
      labels = ['point{0}'.format(i) for i in range(len(py))]
      indexes = [i for i in range(len(py))]
  """
  for i, label, x, y in zip(indexes, labels, px, py):
      if (i >= 120 and i <= 125) or (i >= 240 and i <= 245):
          plt.annotate(
            label,
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad=0'))
  """
  # title and labels
  if(coord == 'x'):
    plt.title(title + ": $x$ coordinates",size="28",y=1.04)
    plt.ylabel("$x$-coordinate (pixels)",family="sans-serif",size="28")
  elif(coord == 'y'):
    plt.title(title + ": $y$ coordinates",size="28",y=1.04)
    plt.ylabel("$y$-coordinate (pixels)",family="sans-serif",size="28")
  elif(coord == 'k'):
    plt.title(title + ": $k$ coefficient",size="28",y=1.04)
    plt.ylabel("$k$-coefficient (standardized)",family="sans-serif",size="28")
  elif(coord == 'd'):
    plt.title(title,size="28",y=1.04)
    plt.ylabel("pupil diameter (mm)",family="sans-serif",size="28")
  elif(coord == 'f'):
    plt.title(title,size="28",y=1.04)
    plt.ylabel("PFE pupil diameter (mm)",family="sans-serif",size="28")
  elif(coord == 'p'):
    plt.title(title + ": change in pupil diameter",size="28",y=1.04)
    plt.ylabel("change in pupil diameter (mm)",family="sans-serif",size="28")
  elif(coord == 'c'):
    plt.title(title + ": change in PFE pupil diameter",size="28",y=1.04)
    plt.ylabel("change in PFE pupil diameter (mm)",family="sans-serif",size="28")
  plt.xlabel("Time (s)",family="sans-serif",size="28")
  plt.tick_params(axis='both',which='major',labelsize="22")
  plt.tick_params(axis='both',which='minor',labelsize="22")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
  plt.tight_layout()

  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderPoints1D(baseName,w,h,points,coord,title,scale=True):
  # plot
  fileName = baseName + '-' + coord + '.pdf'
  # clear plot
  plt.clf()
  # axes dimensions
  ax = plt.axes()
# ax = plt.axes(aspect=1)
  if(coord == 'x' and scale):
    ax.set_ylim((0,int(w)))
  elif(coord == 'y' and scale):
    ax.set_ylim((0,int(h)))
  # fill in data points
  px = []
  py = []
  for pt in points:
    t = pt.gettimestamp()
    if(coord == 'x'):
      x = pt.at(0) * float(w) if scale else pt.at(0)
    elif(coord == 'y'):
      x = pt.at(1) * float(h) if scale else pt.at(1)
    elif(coord == 'k'):
      x = pt.at(0)
    px.append(t)
    py.append(x)

  # lines
  opt = {'antialiased':True,\
         'alpha':.6,\
         'color':"#3F3F3F",\
         'lw':1,\
         'marker':"o",\
         'markersize':100,\
         'markeredgecolor': "#787878",\
         'markeredgewidth':10}
  plt.plot(px,py,antialiased=True,alpha=.6,color="#787878")#, marker='o', markersize=2, linestyle='--')
# line = plt.Line2D(px,py,**opt)
# ax.add_artist(line)
  if (coord == 'x'):
      labels = ['point{0}'.format(i) for i in range(len(px))]
      indexes = [i for i in range(len(px))]
  elif (coord == 'y'):
      labels = ['point{0}'.format(i) for i in range(len(py))]
      indexes = [i for i in range(len(py))]
  """
  for i, label, x, y in zip(indexes, labels, px, py):
      if (i >= 120 and i <= 125) or (i >= 240 and i <= 245):
          plt.annotate(
            label,
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad=0'))
  """
  # title and labels
  if(coord == 'x'):
    plt.title(title + ": $x$ coordinates")
    plt.ylabel("$x$-coordinate (pixels)",family="sans-serif")
  elif(coord == 'y'):
    plt.title(title + ": $y$ coordinates")
    plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  elif(coord == 'k'):
    plt.title(title + ": $k$ coefficient")
    plt.ylabel("$k$-coefficient (standardized)",family="sans-serif")
  plt.xlabel("Time (s)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
  plt.tight_layout()

  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderPoints2D(baseName,w,h,points,title,image=None,scale=False,xtiles=4,ytiles=3):
  # plot
  fileName = baseName + '.pdf'
  imagName = baseName + '.png'

  # clear plot
  plt.clf()

  # axes dimensions (scaled for gaze points, not scaled for velocity)
  if(scale):
    ax = plt.axes(aspect=1)
    ax.set_xlim((0,int(w)))
    ax.set_ylim((0,int(h)))

    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    # use +1 to get the last tick to show up
    ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
    ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
    plt.tick_params(labelsize="9")
  else:
    ax = plt.axes()

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("extent=",x0,x1,y0,y1)
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # fill in data points
  px = []
  py = []
  for pt in points:
    x = pt.at(0) * float(w) if scale else pt.at(0)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - pt.at(1)) * float(h) if scale else pt.at(1)
    px.append(x)
    py.append(y)

  if(not scale):
    ax.set_xlim(min(px),max(px))
    ax.set_ylim(min(py),max(py))

  # lines
  opt = {'antialiased':True,\
         'alpha':.3,\
         'color':"#3F3F3F",\
         'lw':1,\
         'marker':"o",\
         'markersize':2,\
         'markeredgecolor':"#787878",\
         'markeredgewidth':1}
  line = plt.Line2D(px,py,**opt)
  ax.add_artist(line)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)

  plt.title('')
  plt.axis('off')
  plt.savefig(imagName,transparent=False)

  plt.close('all')

def renderKMicrosaccades(baseName,w,h,screen,dist,herz,fixations,K,fixpoints,smthpoints,fixpoints_vx,fixpoints_vy,fixpoints_vx2,fixpoints_vy2,magnitude,title,image=None,scale=False,xtiles=4,ytiles=3):
  # plot
  fileName = baseName + '.pdf'
  imagName = baseName + '.png'

  # sampling period in s
  period = float(1.0/float(herz))

  # clear plot
  plt.clf()

  # axes dimensions (scaled for gaze points, not scaled for velocity)
  if(scale):
    ax = plt.axes(aspect=1)
    ax.set_xlim((0,int(w)))
    ax.set_ylim((0,int(h)))

    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    # use +1 to get the last tick to show up
    ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
    ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
    plt.tick_params(labelsize="9")
  else:
    ax = plt.axes()

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("extent=",x0,x1,y0,y1)
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # calculate K
  if len(K) > 0:
    Kmin = min(K)
    Kmax = max(K)
    Krange = Kmax - Kmin
    print("K min, max, range = ", Kmin, Kmax, Krange)
    nK = [0.0]*len(K)
      # normalize K
    for i in range(len(K)):
      if Krange > 0.0:
        nK[i] = (K[i] - Kmin)/(Krange)
      else:
        nK[i] = (K[i] - Kmin)

  # circles (fixations)
  i=0
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    r = fx.getPercentDuration() * 100.0
    circ = plt.Circle((x,y),radius=r,fc=plt.cm.Oranges(nK[i]),ec='#393939',alpha=.6)
#   circ = plt.Circle((x,y),radius=r,fc=plt.cm.Blues(nK[i]),ec='#393939',alpha=.6)
#   circ = plt.Circle((x,y),radius=r,fc=plt.cm.Greys(nK[i]),ec='#393939',alpha=.4)
    ax.add_patch(circ)
    i += 1

  # fill in data points (basic gazepoint rendering)
  px = []
  py = []
  for k in range(len(fixpoints)):

    # should be list of lists
    indeces = fixpoints[k]

    for i in range(len(indeces)):
      pt = smthpoints[indeces[i]]
      x = pt.at(0) * float(w) if scale else pt.at(0)
      # do the y-coord flip for rendering with (0,0) at bottom-left
      y = (1.0 - pt.at(1)) * float(h) if scale else pt.at(1)
      px.append(x)
      py.append(y)

  # lines
  opt = {'antialiased':True,\
         'alpha':.3,\
#        'color':"#3F3F3F",\
#        'color':"#737373",\
         'color':"#a63603",\
         'lw':1,\
         'marker':"o",\
         'markersize':2,\
         'markeredgecolor':"#787878",\
         'markeredgewidth':1}
  line = plt.Line2D(px,py,**opt)
  ax.add_artist(line)

  # find microsaccades
  for k in range(len(fixpoints)):
    mx = []
    my = []

    # should be list of lists
    indeces = fixpoints[k]

    velocities_x = fixpoints_vx[k]
    velocities_y = fixpoints_vy[k]
    velocities_x2 = fixpoints_vx2[k]
    velocities_y2 = fixpoints_vy2[k]

#   median_vx = sorted(velocities_x)[len(velocities_x)/2]
#   median_vy = sorted(velocities_y)[len(velocities_y)/2]
#   median_vx2 = sorted(velocities_x2)[len(velocities_x2)/2]
#   median_vy2 = sorted(velocities_y2)[len(velocities_y2)/2]
    median_vx = np.median(velocities_x)
    median_vy = np.median(velocities_y)
    median_vx2 = np.median(velocities_x2)
    median_vy2 = np.median(velocities_y2)

    sigma_x = median_vx2 - (median_vx*median_vx)
    sigma_y = median_vy2 - (median_vy*median_vy)

    if sigma_x < 0.00001:
      sigma_x = 0.0
    if sigma_y < 0.00001:
      sigma_y = 0.0

    # as per Engelbert 2006 paper: uses sqrt
    sigma_x = math.sqrt(sigma_x)
    sigma_y = math.sqrt(sigma_y)

    # compute detection threshold (with \lambda = 6)
    # \eta_{xy} = \lamba \sigma_{xy}
    eta_x = 6.0 * sigma_x
    eta_y = 6.0 * sigma_y

    # loop through this set of fixation points and see how long the
    # microsaccade sequences are
    i = 0
    while i < len(indeces):

      # (smoothed) data we are classifying
      x = smthpoints[indeces[i]].at(0)
      y = smthpoints[indeces[i]].at(1)
      t = smthpoints[indeces[i]].gettimestamp()

      vx = velocities_x[i]
      vy = velocities_y[i]
#     vx2 = velocities_x2[i]
#     vy2 = velocities_y2[i]

      pvx = math.fabs(vx)
      pvy = math.fabs(vy)

      pvxr = math.pow(vx/eta_x,2) if eta_x > 0.000001 else 0.0
      pvyr = math.pow(vy/eta_y,2) if eta_y > 0.000001 else 0.0

      dx = magnitude[indeces[i]].at(0)
      dy = magnitude[indeces[i]].at(1)
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

      if pvxr + pvyr > 1.0:

        # check consecutive lengths of points: need to find edges, not points
        j = i + 1

        while j < len(indeces):

          nvx = velocities_x[j]
          nvy = velocities_y[j]
#         nvx2 = velocities_x2[j]
#         nvy2 = velocities_y2[j]

          if math.fabs(nvx) > pvx:
            pvx = math.fabs(nvx)
          if math.fabs(nvy) > pvy:
            pvy = math.fabs(nvy)

          npvxr = math.pow(nvx/eta_x,2) if eta_x > 0.00001 else 0.0
          npvyr = math.pow(nvy/eta_y,2) if eta_y > 0.00001 else 0.0

          dx = magnitude[indeces[j]].at(0)
          dy = magnitude[indeces[j]].at(1)
          nmag = math.sqrt(dx*dx + dy*dy)

          # mean magnitude
          # compute running mean
          # from Brown, Robert Grover, "Introduction to Random Signal
          #   Analysis and Kalman Filtering", John Wiley & Sons, New York, NY
          #   1983 [p.182] TK5102.5 .B696
          umag = float(pmag_k)/float(pmag_k+1) * umag + \
             1.0/float(pmag_k+1) * nmag
          pmag_k += 1

          # peak magnitude
          if nmag > pmag:
            pmag = nmag
          if nmag < minmag:
            minmag = nmag

          # keep going if we're over threshold, else break out of while loop
#         if math.fabs(nvx) > eta_x or math.fabs(nvy) > eta_y:
          if npvxr + npvyr > 1.0:
            j = j + 1
          else:
            break

        # in case we were at the last point?
        if j == len(indeces):
          j = j - 1

        # looking for runs with min. duration of 6 ms (2 ms sampling period)
#       if j - i > 2:
        if 2 < j - i:
#       if 2 < j - i and j - i < 6:

          # mean magnitude computed above
          mag = umag
#         mag = pmag
#         mag = minmag
#         print("mag min, avg, max: ", minmag, umag, pmag)

          # get peak velocity from above
          velx = pvx
          vely = pvy

          # saccade amplitude
#         amp = math.sqrt(velx*velx + vely*vely)
#         amp = (math.fabs(velx) + math.fabs(vely))/2.0
#         amp = max(math.fabs(velx),math.fabs(vely))
          amp = max(velx,vely)

          # magnitude should be < 2 deg
          if 1.0/60.0 < mag and mag < 2.0:

            pt = smthpoints[indeces[i]]
            x = pt.at(0) * float(w) if scale else pt.at(0)
            # do the y-coord flip for rendering with (0,0) at bottom-left
            y = (1.0 - pt.at(1)) * float(h) if scale else pt.at(1)
            mx.append(x)
            my.append(y)

            pt = smthpoints[indeces[j]]
            x = pt.at(0) * float(w) if scale else pt.at(0)
            # do the y-coord flip for rendering with (0,0) at bottom-left
            y = (1.0 - pt.at(1)) * float(h) if scale else pt.at(1)
            mx.append(x)
            my.append(y)

        # end of inner while; advance to next point in sequence past last one
#       i = j + 1
        # advance by 20 ms to not count overshoots, will differ with Hz rate
#       i = j + 10
        samples = 20.0/(period * 1000.0)
        i = int(j + samples)

      # if ith point wasn't a microsaccade, advance by one
      else:
        i = i + 1

    # lines
    opt = {'antialiased':True,\
#          'alpha':.6,\
           'alpha':.7,\
#          'color':"#333333",\
#          'color':"#252525",\
           'color':"#ffeda0",\
#          'color':plt.cm.Greys(nK[k]),\
#          'color':plt.cm.Blues(nK[k]),\
#          'color':plt.cm.Oranges(nK[k]),\
#          'lw':6,\
           'lw':1,\
           'marker':"o",\
#          'markersize':6,\
           'markersize':2,\
           'markeredgecolor':"#ffeda0",\
           'markeredgewidth':1}
    line = plt.Line2D(mx,my,**opt)
    ax.add_artist(line)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)

  plt.title('')
  plt.axis('off')
  plt.savefig(imagName,transparent=False)
  plt.close('all')

def renderMicrosaccades(baseName,w,h,screen,dist,herz,fixpoints,smthpoints,fixpoints_vx,fixpoints_vy,fixpoints_vx2,fixpoints_vy2,magnitude,title,image=None,scale=False,xtiles=4,ytiles=3):
  # plot
  fileName = baseName + '.pdf'

  # sampling period in s
  period = float(1.0/float(herz))

  # clear plot
  plt.clf()

  # axes dimensions (scaled for gaze points, not scaled for velocity)
  if(scale):
    ax = plt.axes(aspect=1)
    ax.set_xlim((0,int(w)))
    ax.set_ylim((0,int(h)))

    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    # use +1 to get the last tick to show up
    ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
    ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
    plt.tick_params(labelsize="9")
  else:
    ax = plt.axes()

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("extent=",x0,x1,y0,y1)
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # fill in data points (basic gazepoint rendering)
  px = []
  py = []
  for k in range(len(fixpoints)):

    # should be list of lists
    indeces = fixpoints[k]

    for i in range(len(indeces)):
      pt = smthpoints[indeces[i]]
      x = pt.at(0) * float(w) if scale else pt.at(0)
      # do the y-coord flip for rendering with (0,0) at bottom-left
      y = (1.0 - pt.at(1)) * float(h) if scale else pt.at(1)
      px.append(x)
      py.append(y)

  # lines
  opt = {'antialiased':True,\
         'alpha':.3,\
         'color':"#3F3F3F",\
         'lw':1,\
         'marker':"o",\
         'markersize':2,\
         'markeredgecolor':"#787878",\
         'markeredgewidth':1}
  line = plt.Line2D(px,py,**opt)
  ax.add_artist(line)

  # find microsaccades
  for k in range(len(fixpoints)):
    mx = []
    my = []

    # should be list of lists
    indeces = fixpoints[k]

    velocities_x = fixpoints_vx[k]
    velocities_y = fixpoints_vy[k]
    velocities_x2 = fixpoints_vx2[k]
    velocities_y2 = fixpoints_vy2[k]

#   median_vx = sorted(velocities_x)[len(velocities_x)/2]
#   median_vy = sorted(velocities_y)[len(velocities_y)/2]
#   median_vx2 = sorted(velocities_x2)[len(velocities_x2)/2]
#   median_vy2 = sorted(velocities_y2)[len(velocities_y2)/2]
    median_vx = np.median(velocities_x)
    median_vy = np.median(velocities_y)
    median_vx2 = np.median(velocities_x2)
    median_vy2 = np.median(velocities_y2)

    sigma_x = median_vx2 - (median_vx*median_vx)
    sigma_y = median_vy2 - (median_vy*median_vy)

    if sigma_x < 0.00001:
      sigma_x = 0.0
    if sigma_y < 0.00001:
      sigma_y = 0.0

    # as per Engelbert 2006 paper: uses sqrt
    sigma_x = math.sqrt(sigma_x)
    sigma_y = math.sqrt(sigma_y)

    # compute detection threshold (with \lambda = 6)
    # \eta_{xy} = \lamba \sigma_{xy}
    eta_x = 6.0 * sigma_x
    eta_y = 6.0 * sigma_y

    # loop through this set of fixation points and see how long the
    # microsaccade sequences are
    i = 0
    while i < len(indeces):

      # (smoothed) data we are classifying
      x = smthpoints[indeces[i]].at(0)
      y = smthpoints[indeces[i]].at(1)
      t = smthpoints[indeces[i]].gettimestamp()

      vx = velocities_x[i]
      vy = velocities_y[i]
#     vx2 = velocities_x2[i]
#     vy2 = velocities_y2[i]

      pvx = math.fabs(vx)
      pvy = math.fabs(vy)

      pvxr = math.pow(vx/eta_x,2) if eta_x > 0.000001 else 0.0
      pvyr = math.pow(vy/eta_y,2) if eta_y > 0.000001 else 0.0

      dx = magnitude[indeces[i]].at(0)
      dy = magnitude[indeces[i]].at(1)
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

      if pvxr + pvyr > 1.0:

        # check consecutive lengths of points: need to find edges, not points
        j = i + 1

        while j < len(indeces):

          nvx = velocities_x[j]
          nvy = velocities_y[j]
#         nvx2 = velocities_x2[j]
#         nvy2 = velocities_y2[j]

          if math.fabs(nvx) > pvx:
            pvx = math.fabs(nvx)
          if math.fabs(nvy) > pvy:
            pvy = math.fabs(nvy)

          npvxr = math.pow(nvx/eta_x,2) if eta_x > 0.00001 else 0.0
          npvyr = math.pow(nvy/eta_y,2) if eta_y > 0.00001 else 0.0

          dx = magnitude[indeces[j]].at(0)
          dy = magnitude[indeces[j]].at(1)
          nmag = math.sqrt(dx*dx + dy*dy)

          # mean magnitude
          # compute running mean
          # from Brown, Robert Grover, "Introduction to Random Signal
          #   Analysis and Kalman Filtering", John Wiley & Sons, New York, NY
          #   1983 [p.182] TK5102.5 .B696
          umag = float(pmag_k)/float(pmag_k+1) * umag + \
             1.0/float(pmag_k+1) * nmag
          pmag_k += 1

          # peak magnitude
          if nmag > pmag:
            pmag = nmag
          if nmag < minmag:
            minmag = nmag

          # keep going if we're over threshold, else break out of while loop
#         if math.fabs(nvx) > eta_x or math.fabs(nvy) > eta_y:
          if npvxr + npvyr > 1.0:
            j = j + 1
          else:
            break

        # in case we were at the last point?
        if j == len(indeces):
          j = j - 1

        # looking for runs with min. duration of 6 ms (2 ms sampling period)
#       if j - i > 2:
        if 2 < j - i:
#       if 2 < j - i and j - i < 6:

          # mean magnitude computed above
          mag = umag
#         mag = pmag
#         mag = minmag
#         print("mag min, avg, max: ", minmag, umag, pmag)

          # get peak velocity from above
          velx = pvx
          vely = pvy

          # saccade amplitude
#         amp = math.sqrt(velx*velx + vely*vely)
#         amp = (math.fabs(velx) + math.fabs(vely))/2.0
#         amp = max(math.fabs(velx),math.fabs(vely))
          amp = max(velx,vely)

          # magnitude should be < 2 deg
          if 1.0/60.0 < mag and mag < 2.0:

            pt = smthpoints[indeces[i]]
            x = pt.at(0) * float(w) if scale else pt.at(0)
            # do the y-coord flip for rendering with (0,0) at bottom-left
            y = (1.0 - pt.at(1)) * float(h) if scale else pt.at(1)
            mx.append(x)
            my.append(y)

            pt = smthpoints[indeces[j]]
            x = pt.at(0) * float(w) if scale else pt.at(0)
            # do the y-coord flip for rendering with (0,0) at bottom-left
            y = (1.0 - pt.at(1)) * float(h) if scale else pt.at(1)
            mx.append(x)
            my.append(y)

        # end of inner while; advance to next point in sequence past last one
#       i = j + 1
        # advance by 20 ms to not count overshoots, will differ with Hz rate
#       i = j + 10
        samples = 20.0/(period * 1000.0)
        i = int(j + samples)

      # if ith point wasn't a microsaccade, advance by one
      else:
        i = i + 1

    # lines
    opt = {'antialiased':True,\
           'alpha':.7,\
#          'color':"#333333",\
           'color':"#ffeda0",\
           'lw':1,\
           'marker':"o",\
           'markersize':2,\
           'markeredgecolor':"#ffeda0",\
           'markeredgewidth':1}
    line = plt.Line2D(mx,my,**opt)
    ax.add_artist(line)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderFixations(baseName,w,h,screen,viewdist,fixations,title,image=None,lagrange=None,xtiles=4,ytiles=3):
  # plot
  fileName = baseName + '.pdf'

  # clear plot
  plt.clf()

  # axes dimensions
  ax = plt.axes(aspect=1)
  ax.set_xlim((0,int(w)))
  ax.set_ylim((0,int(h)))

  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  # use +1 to get the last tick to show up
  ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
  ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
  plt.tick_params(labelsize="9")

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  cx = []
  cy = []
  # calibration dots
  for i in range(lagrange.n):
    x = lagrange.S[i,0] * float(w)
    # don't do the y-coord flip for cal dots rendering with (0,0) at bottom-left
    # no I think we do need the y-flip
    y = (1.0 - lagrange.S[i,1]) * float(h)
    cx.append(x)
    cy.append(y)

  # init kd-tree with calibration points
  # (each fixation point will then query for nearest neighbor)
  kdtree = spatial.KDTree(list(zip(cx,cy)))

  # for each fixation, find nearest neighbor and distance to it (in pixels)
  X_err_i = [0]*(lagrange.n)
  X_err = [0.0]*(lagrange.n)
  distances = []
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)

    # see: http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
    # get two nearest neigbhors: the first will be itself, with distance 0
    # return (lists of) nearest neighbor distances and indeces in tree data

    # use fixation point only if it is inside image borders
    # (image may be smaller than screen w,h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):

      # get nearest neighbor
      nndist,nnidxs = kdtree.query(np.array([[x,y]]),1)
#     print("nndist: ", nndist)
#     print("nnidxs: ", nnidxs)
#     print("dist to nearest neighbor (in pixels): ", nndist[0])

      # get index to kdtree element (same index as lagrange calib points)
      i = nnidxs[0]
      # get distance to closest kdtree element (calibration point)
      dist = float(nndist[0])
  
      # add in to list of distances
      distances.append(dist)

      # compute running mean for each calibration point as we encounter it
      X_err[i] = float(X_err_i[i])/float(X_err_i[i]+1) * X_err[i] + \
                               1.0/float(X_err_i[i]+1) * dist
      # compute centroid (running mean) of error at each calibration point
      lagrange.X_bar[i,0] = float(X_err_i[i])/float(X_err_i[i]+1) * \
                                              lagrange.X_bar[i,0] + \
                               1.0/float(X_err_i[i]+1) * x
      lagrange.X_bar[i,1] = float(X_err_i[i])/float(X_err_i[i]+1) * \
                                              lagrange.X_bar[i,1] + \
                               1.0/float(X_err_i[i]+1) * y

      # increment how many times we've been close to this calib point
      X_err_i[i] += 1

# print("X_err: ", X_err)

  if len(fixations) > 0:
    # compute observed mean distance between each feature (fixation) and
    # its nearest neihgbor (calibration point)
    avedist = np.mean(distances)
    stddist = np.std(distances)
    print("mean dist b/ween any fixation and clb: %f (pixels)" % (avedist))

    r = math.sqrt(float(w)*float(w) + float(h)*float(h))
    dpi = r/float(screen)

    D = float(viewdist)

    fov = 2*math.degrees(math.atan2(screen,2*D))
    fovx = 2*math.degrees(math.atan2(float(w)/dpi,2*D))
    fovy = 2*math.degrees(math.atan2(float(h)/dpi,2*D))

    avedist_deg = 2*math.degrees(math.atan2((avedist/dpi),(2*D)))
    stddist_deg = 2*math.degrees(math.atan2((stddist/dpi),(2*D)))

    print("mean dist b/ween any fixation and clb: %f (degrees)" % (avedist_deg))
  
    strout = "view distance: %5.2f (inches), screen: %3.0f (inches), %5.2f$^{\circ}$ (visual angle), dpi: %5.2f" % \
             (D,screen,fov,dpi)
    ax.text(10,int(h)-40,strout,fontsize=10)
    strout = "mean error (accuracy): %5.2f$^{\circ}$ (degrees visual angle), standard deviation (precision): %5.2f$^{\circ}$" % \
             (avedist_deg, stddist_deg)
    ax.text(10,10,strout,fontsize=10)

  # fill in data points
  px = []
  py = []
  i=0
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):
      px.append(x)
      py.append(y)
      if i>0:
        sx = fixations[i-1].at(0) * float(w)
        # do the y-coord flip for rendering with (0,0) at bottom-left
        sy = (1.0 - fixations[i-1].at(1)) * float(h)
        if image is None or (x0 < sx and sx < x1 and y0 < sy and sy < y1):
          dx = x - sx
          dy = y - sy
#         arrow = plt.Arrow(sx,sy,dx,dy,\
#                           width=35,\
#                           alpha=.2,\
#                           fill=False,\
#                           fc="none",\
#                           ec="#101010")
#         ax.add_patch(arrow)
          opt = {'shape':"full",\
                 'fill':False,\
                 'fc':"none",\
                 'ec':"#101010",\
                 'width':1,\
                 'head_width':15,\
                 'head_length':45,\
                 'length_includes_head':True}
          if(abs(dx) > 0.0 and abs(dy) > 0.0):
            plt.arrow(sx,sy,dx,dy,alpha=.3,**opt)
    i += 1
  # old way, straight plot
# plt.plot(px,py,antialiased=True,alpha=.6,color="#756BB1")
  # or using Artist
  line = plt.Line2D(px,py,antialiased=True,alpha=.3,color="#494949",lw=1)
  ax.add_artist(line)

  # circles
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):
      r = fx.getPercentDuration() * 100.0
      circ = plt.Circle((x,y),radius=r,fc='#BFBFBF',ec='#393939',alpha=.6)
      ax.add_patch(circ)
  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderAmfocFixations(baseName,w,h,fixations,K,title,image=None,xtiles=4,ytiles=3):
  # plot
  fileName = baseName + '.pdf'

  # clear plot
  plt.clf()

  # axes dimensions
  ax = plt.axes(aspect=1)
  ax.set_xlim((0,int(w)))
  ax.set_ylim((0,int(h)))

  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  # use +1 to get the last tick to show up
  ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
  ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
  plt.tick_params(labelsize="9")

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # fill in data points
  px = []
  py = []
  i=0
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    px.append(x)
    py.append(y)
    if i>0:
      sx = fixations[i-1].at(0) * float(w)
      # do the y-coord flip for rendering with (0,0) at bottom-left
      sy = (1.0 - fixations[i-1].at(1)) * float(h)
      dx = x - sx
      dy = y - sy
#     arrow = plt.Arrow(sx,sy,dx,dy,\
#                       width=35,\
#                       alpha=.2,\
#                       fill=False,\
#                       fc="none",\
#                       ec="#101010")
#     ax.add_patch(arrow)
      opt = {'shape':"full",\
             'fill':False,\
             'fc':"none",\
             'ec':"#101010",\
             'width':1,\
             'head_width':15,\
             'head_length':45,\
             'length_includes_head':True}
      if(abs(dx) > 0.0 and abs(dy) > 0.0):
        plt.arrow(sx,sy,dx,dy,alpha=.3,**opt)
    i += 1
  # old way, straight plot
# plt.plot(px,py,antialiased=True,alpha=.6,color="#756BB1")
  # or using Artist
  line = plt.Line2D(px,py,antialiased=True,alpha=.3,color="#494949",lw=1)
  ax.add_artist(line)

  if len(K) > 0:
    Kmin = min(K)
    Kmax = max(K)
    Krange = Kmax - Kmin
    print("K min, max, range = ", Kmin, Kmax, Krange)
    nK = [0.0]*len(K)
      # normalize K
    for i in range(len(K)):
      if Krange > 0.0:
        nK[i] = (K[i] - Kmin)/(Krange)
      else:
        nK[i] = (K[i] - Kmin)

  # circles
  i=0
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    r = fx.getPercentDuration() * 100.0
#   if K[i] > 0.3:
#     # focal
#     circ = plt.Circle((x,y),radius=r,fc='#fc8d59',ec='#393939',alpha=.6)
#   elif K[i] < -0.3:
#     # ambient
#     circ = plt.Circle((x,y),radius=r,fc='#91bfdb',ec='#393939',alpha=.6)
#   else:
#     # neither ambient nor focal, use grey color
#     circ = plt.Circle((x,y),radius=r,fc='#BFBFBF',ec='#393939',alpha=.6)
#   circ = plt.Circle((x,y),radius=r,fc=plt.cm.Oranges(nK[i]),ec='#393939',alpha=.6)
    circ = plt.Circle((x,y),radius=r,fc=plt.cm.Blues(nK[i]),ec='#393939',alpha=.6)
    ax.add_patch(circ)
    i += 1
  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderCalibFixations(baseName,w,h,screen,viewdist,fixations,title,image=None,lagrange=None,xtiles=4,ytiles=3):

  # diagonal
  d = math.sqrt(float(w)*float(w) + float(h)*float(h))

  # plot
  fileName = baseName + '.pdf'

  # clear plot
  plt.clf()

  # axes dimensions
  ax = plt.axes(aspect=1)
  ax.set_xlim((0,int(w)))
  ax.set_ylim((0,int(h)))

  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  # use +1 to get the last tick to show up
  ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
  ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
  plt.tick_params(labelsize="9")

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  cx = []
  cy = []
  # calibration dots
  for i in range(lagrange.n):
    x = lagrange.S[i,0] * float(w)
    # don't do the y-coord flip for cal dots rendering with (0,0) at bottom-left
    # no I think we do need the y-flip
    y = (1.0 - lagrange.S[i,1]) * float(h)
#   print("cb ",i,": (x,y): (", x, ",", y, ")")
    r = 1.0/256.0 * d
    # circle for calibrtation dot
#   circ = plt.Circle((x,y),radius=r,fc='#bababa',ec='#393939',alpha=.6)
#   ax.add_patch(circ)
    # square box for calibration dot
#   path = Path([(x-r,y-r),(x+r,y-r),(x+r,y+r),(x-r,y+r),(x-r,y-r)],\
#               [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY])
#   patch = patches.PathPatch(path,fc='#bababa',ec='#393939',alpha=.6)
#   ax.add_patch(patch)
    # an x for centroid
    path = Path([(x-r,y-r),(x+r,y+r)],[Path.MOVETO,Path.LINETO])
    patch = patches.PathPatch(path,fc='#4d4d4d',ec='#4d4d4d',alpha=.8)
    ax.add_patch(patch)
    path = Path([(x-r,y+r),(x+r,y-r)],[Path.MOVETO,Path.LINETO])
    patch = patches.PathPatch(path,fc='#4d4d4d',ec='#4d4d4d',alpha=.8)
    ax.add_patch(patch)
    # don't append the last calibration point: it is the same as the first
    cx.append(x)
    cy.append(y)

  # init kd-tree with calibration points
  # (each fixation point will then query for nearest neighbor)
  kdtree = spatial.KDTree(list(zip(cx,cy)))

  # for each fixation, find nearest neighbor and distance to it (in pixels)
  X_err_i = [0]*(lagrange.n)
  X_err = [0.0]*(lagrange.n)
  distances = []
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)

    # see: http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
    # get two nearest neigbhors: the first will be itself, with distance 0
    # return (lists of) nearest neighbor distances and indeces in tree data
#   nndist,nnidxs = kdtree.query(np.array([[x,y]]),2)
#   print("query: (x,y): (", x, ",", y, ")")
#   print("nndist: ", nndist)
#   print("nnidxs: ", nnidxs)
#   print("dist to nearest neighbor (in pixels): ", nndist[0][0])
#   print("kdtree.data: ", kdtree.data)
#   print("kdtree.data[nnidxs[0][0]]: ", kdtree.data[nnidxs[0][0]])
#   print("kdtree.data[nnidxs[0][0]][0]: ", kdtree.data[nnidxs[0][0]][0])
#   print("kdtree.data[nnidxs[0][0]][1]: ", kdtree.data[nnidxs[0][0]][1])
#
#   # get index to kdtree element (same index as lagrange calib points)
#   i = nnidxs[0][0]
#   # get distance to closest kdtree element (calibration point)
#   dist = float(nndist[0][0])

    # use fixation point only if it is inside image borders
    # (image may be smaller than screen w,h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):
      # get nearest neighbor
      nndist,nnidxs = kdtree.query(np.array([[x,y]]),1)
#     print("nndist: ", nndist)
#     print("nnidxs: ", nnidxs)
#     print("dist to nearest neighbor (in pixels): ", nndist[0])

      # get index to kdtree element (same index as lagrange calib points)
      i = nnidxs[0]
      # get distance to closest kdtree element (calibration point)
      dist = float(nndist[0])

      # add in to list of distances
      distances.append(dist)

      # compute running mean for each calibration point as we encounter it
      X_err[i] = float(X_err_i[i])/float(X_err_i[i]+1) * X_err[i] + \
                               1.0/float(X_err_i[i]+1) * dist
      # compute centroid (running mean) of error at each calibration point
      lagrange.X_bar[i,0] = float(X_err_i[i])/float(X_err_i[i]+1) * \
                                              lagrange.X_bar[i,0] + \
                               1.0/float(X_err_i[i]+1) * x
      lagrange.X_bar[i,1] = float(X_err_i[i])/float(X_err_i[i]+1) * \
                                              lagrange.X_bar[i,1] + \
                               1.0/float(X_err_i[i]+1) * y
      lagrange.X_bar_n[i] += 1

      # increment how many times we've been close to this calib point
      X_err_i[i] += 1

# print("X_err: ", X_err)

  # render error centroid
  for i in range(lagrange.n):
    sx = lagrange.S[i,0] * float(w)
    sy = (1.0 - lagrange.S[i,1])* float(h)
    # check for no points at calib point, if none, set error centroid to
    # calib point coords (0 error)
    if lagrange.X_bar_n[i] == 0:
      lagrange.X_bar[i,0] = sx
      lagrange.X_bar[i,1] = sy
    x = lagrange.X_bar[i,0]
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = lagrange.X_bar[i,1]
#   print("er ",i,": (x,y): (", x, ",", y, ")")
    r = 1.0/256.0 * d
    # circle for centroid, not good...
#   circ = plt.Circle((x,y),radius=r,fc='#fdae61',ec='#393939',alpha=.6)
#   ax.add_patch(circ)
    # square box for centroid
    path = Path([(x-r,y-r),(x+r,y-r),(x+r,y+r),(x-r,y+r),(x-r,y-r)],\
                [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY])
    patch = patches.PathPatch(path,fc='#fddbc7',ec='#393939',alpha=.6)
    ax.add_patch(patch)
    # an x for centroid
#   path = Path([(x-r,y-r),(x+r,y+r)],[Path.MOVETO,Path.LINETO])
#   patch = patches.PathPatch(path,fc='#fdae61',ec='#393939',alpha=.6)
#   ax.add_patch(patch)
#   path = Path([(x-r,y+r),(x+r,y-r)],[Path.MOVETO,Path.LINETO])
#   patch = patches.PathPatch(path,fc='#fdae61',ec='#393939',alpha=.6)
#   ax.add_patch(patch)
    # line from centroid calibration point to fixation centroid
    # from: http://matplotlib.org/users/path_tutorial.html
    path = Path([(sx,sy),(x,y)],[Path.MOVETO,Path.LINETO])
    patch = patches.PathPatch(path,fc='#878787',ec='#878787',lw=1,alpha=.8)
    ax.add_patch(patch)

  if len(fixations) > 0:
    # compute observed mean distance between each feature (fixation) and
    # its nearest neihgbor (calibration point)
    avedist = np.mean(distances)
    stddist = np.std(distances)
    print("mean dist b/ween any fixation and clb: %f (pixels)" % (avedist))

    r = math.sqrt(float(w)*float(w) + float(h)*float(h))
    dpi = r/float(screen)

    D = float(viewdist)

    fov = 2*math.degrees(math.atan2(screen,2*D))
    fovx = 2*math.degrees(math.atan2(float(w)/dpi,2*D))
    fovy = 2*math.degrees(math.atan2(float(h)/dpi,2*D))

    avedist_deg = 2*math.degrees(math.atan2((avedist/dpi),(2*D)))
    stddist_deg = 2*math.degrees(math.atan2((stddist/dpi),(2*D)))

    print("mean dist b/ween any fixation and clb: %f (degrees)" % (avedist_deg))
  
    strout = "view distance: %5.2f (inches), screen: %3.0f (inches), %5.2f$^{\circ}$ (visual angle), dpi: %5.2f" % \
             (D,screen,fov,dpi)
    ax.text(10,int(h)-40,strout,fontsize=10)
    strout = "mean error (accuracy): %5.2f$^{\circ}$ (degrees visual angle), standard deviation (precision): %5.2f$^{\circ}$" % \
             (avedist_deg, stddist_deg)
    ax.text(10,10,strout,fontsize=10)

# don't bother; for n-back this won't work since we only have observed points
#               at 1 "calibration point"
# lagrange.solve(w,h)

  # draw arrowheads
  px = []
  py = []
  i=0
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):
      px.append(x)
      py.append(y)
      if i>0:
        sx = fixations[i-1].at(0) * float(w)
        # do the y-coord flip for rendering with (0,0) at bottom-left
        sy = (1.0 - fixations[i-1].at(1)) * float(h)
        if image is None or (x0 < sx and sx < x1 and y0 < sy and sy < y1):
          dx = x - sx
          dy = y - sy
#         arrow = plt.Arrow(sx,sy,dx,dy,\
#                           width=35,\
#                           alpha=.2,\
#                           fill=False,\
#                           fc="none",\
#                           ec="#101010")
#         ax.add_patch(arrow)
          opt = {'shape':"full",\
                 'fill':False,\
                 'fc':"none",\
                 'ec':"#101010",\
                 'width':1,\
                 'head_width':15,\
                 'head_length':45,\
                 'length_includes_head':True}
          if(abs(dx) > 0.0 and abs(dy) > 0.0):
            plt.arrow(sx,sy,dx,dy,alpha=.2,**opt)
      i += 1
    # old way, straight plot
# plt.plot(px,py,antialiased=True,alpha=.6,color="#756BB1")
  # or using Artist
  line = plt.Line2D(px,py,antialiased=True,alpha=.3,color="#494949",lw=1)
  ax.add_artist(line)

  # circles (corrected): won't work in this (n-back) case because
  # we only have samples at 1 "calibration point"
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):
      r = fx.getPercentDuration() * 100.0
      circ = plt.Circle((x,y),radius=r,fc='#fddbc7',ec='#393939',alpha=.2)
      ax.add_patch(circ)
#     (fitx, fity) = lagrange.transform(x,y)
#     circ = plt.Circle((fitx,fity),radius=r,fc='#d6604d',ec='#393939',alpha=.4)
#     ax.add_patch(circ)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderAOIs(baseName,w,h,aoilist,key,title,image=None,xtiles=4,ytiles=3):
  # plot
  fileName = baseName + '.pdf'

  # clear plot
  plt.clf()

  # axes dimensions
  ax = plt.axes(aspect=1)
  ax.set_xlim((0,int(w)))
  ax.set_ylim((0,int(h)))

  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  # use +1 to get the last tick to show up
  ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
  ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
  plt.tick_params(labelsize="9")

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # AOIs
  for aoi in aoilist:
    rect = plt.Rectangle(xy=aoi.getXY(h),\
                         width=aoi.getWidth(),\
                         height=aoi.getHeight(),\
                         fc='#BFBFBF',ec='#393939',alpha=.4)
    ax.add_patch(rect)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderAOIFixations(baseName,w,h,fixations,aoilist,key,title,image=None,xtiles=4,ytiles=3):
  # plot
  fileName = baseName + '.pdf'

  # clear plot
  plt.clf()

  # axes dimensions
  ax = plt.axes(aspect=1)
  ax.set_xlim((0,int(w)))
  ax.set_ylim((0,int(h)))

  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  # use +1 to get the last tick to show up
  ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
  ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
  plt.tick_params(labelsize="9")

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # fill in data points
  px = []
  py = []
  i=0
  for fx in fixations:
    x = fx.at(0) * float(w)
    y = fx.at(1) * float(h)
#   this code omits fixations that are outside of AOIs...trouble with it
#   is that it messes up the fixation scanpath sequence
#   inAOI = False
#   for aoi in aoidict[key]:
#     if aoi.inside(x,y + aoi.getHeight()):
#       inAOI = True
#       break
    inAOI = True
    if inAOI:
      # do the y-coord flip for rendering with (0,0) at bottom-left
      y = (1.0 - fx.at(1)) * float(h)
      px.append(x)
      py.append(y)
      if i>0:
        sx = fixations[i-1].at(0) * float(w)
        # do the y-coord flip for rendering with (0,0) at bottom-left
        sy = (1.0 - fixations[i-1].at(1)) * float(h)
        dx = x - sx
        dy = y - sy
#       arrow = plt.Arrow(sx,sy,dx,dy,\
#                         width=35,\
#                         alpha=.2,\
#                         fill=False,\
#                         fc="none",\
#                         ec="#101010")
#       ax.add_patch(arrow)
        opt = {'shape':"full",\
               'fill':False,\
               'fc':"none",\
               'ec':"#101010",\
               'width':1,\
               'head_width':15,\
               'head_length':45,\
               'length_includes_head':True}
        if(abs(dx) > 0.0 and abs(dy) > 0.0):
          plt.arrow(sx,sy,dx,dy,alpha=.3,**opt)
      i += 1
  # old way, straight plot
# plt.plot(px,py,antialiased=True,alpha=.6,color="#756BB1")
  # or using Artist
  line = plt.Line2D(px,py,antialiased=True,alpha=.3,color="#494949",lw=1)
  ax.add_artist(line)

  # circles
  for fx in fixations:
    x = fx.at(0) * float(w)
    y = fx.at(1) * float(h)
#   this code omits fixations that are outside of AOIs...trouble with it
#   is that it messes up the fixation scanpath sequence
#   inAOI = False
#   for aoi in aoidict[key]:
#     if aoi.inside(x,y + aoi.getHeight()):
#       inAOI = True
#       break
    inAOI = True
    if inAOI:
      # do the y-coord flip for rendering with (0,0) at bottom-left
      y = (1.0 - fx.at(1)) * float(h)
      r = fx.getPercentDuration() * 100.0
#     circ = plt.Circle((x,y),radius=r,fc='#BFBFBF',ec='#393939',alpha=.6)
      circ = plt.Circle((x,y),radius=r,fc='#fdae61',ec='#393939',alpha=.6)
      ax.add_patch(circ)

  # AOIs
# for aoi in aoidict[key]:
  for aoi in aoilist:
    inside = False
    for fx in fixations:
      x = fx.at(0) * float(w)
      # don't do the y-flip here (same ref. frame as AOIs)
      y = fx.at(1) * float(h)
      # add in AOI height---AOIs apear to be off-by-one line height
      if aoi.inside(x,y + aoi.getHeight()):
        inside = True
        break
    if inside:
      rect = plt.Rectangle(xy=aoi.getXY(h),\
                           width=aoi.getWidth(),\
                           height=aoi.getHeight(),\
                           fc='#d7191c',ec='#393939',alpha=.4)
    else:
      rect = plt.Rectangle(xy=aoi.getXY(h),\
                           width=aoi.getWidth(),\
                           height=aoi.getHeight(),\
                           fc='#abdda4',ec='#393939',alpha=.4)
    ax.add_patch(rect)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderFixatedAOIs(baseName,w,h,fixations,aoilist,key,title,image=None,xtiles=4,ytiles=3):
  # plot
  fileName = baseName + '.pdf'

  # clear plot
  plt.clf()

  # axes dimensions
  ax = plt.axes(aspect=1)
  ax.set_xlim((0,int(w)))
  ax.set_ylim((0,int(h)))

  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  # use +1 to get the last tick to show up
  ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
  ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
  plt.tick_params(labelsize="9")

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # fill in data points
  px = []
  py = []
  i=0
  prev_fixation=0
  for fx in fixations:
    x = fx.at(0) * float(w)
    y = fx.at(1) * float(h)
    inAOI = False
    for aoi in aoilist:
      if aoi.inside(x,y + aoi.getHeight()):
        inAOI = True
        break
    if inAOI:
#     rect = plt.Rectangle(xy=aoi.getXY(h),\
#                          width=aoi.getWidth(),\
#                          height=aoi.getHeight(),\
#                          fc='#d7191c',ec='#393939',alpha=.4)
      rect = plt.Rectangle(xy=aoi.getXY(h),\
                           width=aoi.getWidth(),\
                           height=aoi.getHeight(),\
                           fc='#BFBFBF',ec='#393939',alpha=.2)
      ax.add_patch(rect)
      # do the y-coord flip for rendering with (0,0) at bottom-left
      y = (1.0 - fx.at(1)) * float(h)
      r = fx.getPercentDuration() * 100.0
#     circ = plt.Circle((x,y),radius=r,fc='#BFBFBF',ec='#393939',alpha=.6)
      circ = plt.Circle((x,y),radius=r,fc='#fdae61',ec='#393939',alpha=.6)
      ax.add_patch(circ)
      px.append(x)
      py.append(y)
      if prev_fixation>0:
        sx = fixations[prev_fixation].at(0) * float(w)
        # do the y-coord flip for rendering with (0,0) at bottom-left
        sy = (1.0 - fixations[prev_fixation].at(1)) * float(h)
        dx = x - sx
        dy = y - sy
        opt = {'shape':"full",\
               'fill':False,\
               'fc':"none",\
               'ec':"#101010",\
               'width':1,\
               'head_width':15,\
               'head_length':45,\
               'length_includes_head':True}
        if(abs(dx) > 0.0 and abs(dy) > 0.0):
          plt.arrow(sx,sy,dx,dy,alpha=.3,**opt)
      prev_fixation=i
    i += 1
  # or using Artist
  line = plt.Line2D(px,py,antialiased=True,alpha=.3,color="#494949",lw=1)
  ax.add_artist(line)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.close('all')

def renderHeatmap(baseName,w,h,fixations,title,image=None,xtiles=4,ytiles=3):
  # see also: http://www.pmavridis.com/misc/heatmaps/

  # diagonal
  d = math.sqrt(float(w)*float(w) + float(h)*float(h))
  sigma = 0.0

  # plot
  fileName = baseName + '.pdf'
  fileName_png = baseName + '.png'

  # clear plot
  plt.clf()

  # axes dimensions
  ax = plt.axes(aspect=1)
  ax.set_xlim((0,int(w)))
  ax.set_ylim((0,int(h)))

  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  # use +1 to get the last tick to show up
  ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
  ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
  plt.tick_params(labelsize="9")

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),alpha=0.5,cmap=pylab.gray(),
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,cmap=pylab.gray(),
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # heatmap: a 32-bit floating point image, initially set to black
  lum = Image.new("F", (int(w),int(h)), 0.0)
  pix = lum.load()
  print("processing ", len(fixations), " fixations...")
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    # hack: if there's only 1 fixation @ 0% duration: 1/6th of image
    sigma = fx.getPercentDuration() * d/6.0 if fx.getPercentDuration() > 0 else d/6.0
#   print("percent duration = ", fx.getPercentDuration(), "d = ", d)
#   print("sigma = ", sigma)
    # lum.putdata() might be faster!!
#   for i in range(int(h)):
#     for j in range(int(w)):
    for i in range(int(y-2.0*sigma),int(y+2.0*sigma)):
      for j in range(int(x-2.0*sigma),int(x+2.0*sigma)):
        if( 0 <= i and i < int(h) and 0 <= j and j < int(w) ):
          sx = j - x
          sy = i - y
          heat = math.exp((sx*sx + sy*sy)/(-2.0*sigma*sigma))
          pix[j,i] = pix[j,i] + heat
  print("done.")

  # get max value
  minlum, maxlum = lum.getextrema()
  print("minlum, maxlum = ", (minlum, maxlum))

  # normalize
# for i in range(int(h)):
#   for j in range(int(w)):
#     pix[j,i] = pix[j,i] / maxlum
  if(abs(maxlum) < 0.00001):
    maxlum = 1.0
  lum = lum.point(lambda f: f * (1.0/maxlum) + 0)
  print("done normalizing")

  # convert to grayscale
# out = lum.point(lambda f: f * 255.0 + 0,"L")
  out = lum.point(lambda f: f * 255.0 + 0)
  print("done converting")

  # plot
  # imshow can handle lum image, default colormap is "jet", e.g., heatmap
# plt.imshow(np.asarray(lum),cmap="gray")
  plt.imshow(np.asarray(lum),cmap="jet",alpha=0.5)
# plt.imshow(np.asarray(lum),cmap="rainbow",alpha=0.5)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.savefig(fileName_png,transparent=True)
  plt.close('all')

def renderNbackFixations(baseName,w,h,screen,viewdist,fixations,title,image=None,monitor=None,ttype=None,xtiles=4,ytiles=3):

  # diagonal
  d = math.sqrt(float(w)*float(w) + float(h)*float(h))

  # plot
  fileName = baseName + '.pdf'

  # clear plot
  plt.clf()

  # axes dimensions
  ax = plt.axes(aspect=1)
  ax.set_xlim((0,int(w)))
  ax.set_ylim((0,int(h)))

  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  # use +1 to get the last tick to show up
  ax.set_xticks(np.arange(xmin,xmax+1,xmax/xtiles))
  ax.set_yticks(np.arange(ymin,ymax+1,ymax/ytiles))
  plt.tick_params(labelsize="9")

  # set background image
  # from: http://matplotlib.org/users/image_tutorial.html
  if image is not None:
#   img = mpimg.imread(image)
#   plt.imshow(img)
    # using PIL, see: http://effbot.org/imagingbook/image.htm
#   img = Image.open(image)
    img = Image.open(image).rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    (imw, imh) = img.size
    x0 = (int(w)-imw)/2
    y0 = (int(h)-imh)/2
    x1 = x0+imw
    y1 = y0+imh
#   print("Image mode: ", img.mode)
    if img.mode == "L":
#     plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),cmap=pylab.gray(),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))
    else:
#     plt.imshow(np.asarray(img),alpha=0.5,
#                origin='None',aspect='auto',extent=(x0,x1,y0,y1))
      plt.imshow(np.asarray(img),alpha=0.5,
                 origin='lower',aspect='auto',extent=(x0,x1,y0,y1))

  # calibration dots
  for i in range(monitor.n):
    x = monitor.S[i,0] * float(w)
    # don't do the y-coord flip for cal dots rendering with (0,0) at bottom-left
    # no I think we do need the y-flip
    y = (1.0 - monitor.S[i,1]) * float(h)
#   print("cb ",i,": (x,y): (", x, ",", y, ")")
    r = monitor.deg2pix(3.0)
    # circle for target dot
    if i == monitor.ndict[ttype]: 
      circ = plt.Circle((x,y),radius=r,fc='#d6604d',ec='#393939',alpha=.4)
      ax.add_patch(circ)
    # an x for centroid
    r = 1.0/256.0 * d
    path = Path([(x-r,y-r),(x+r,y+r)],[Path.MOVETO,Path.LINETO])
    patch = patches.PathPatch(path,fc='#4d4d4d',ec='#4d4d4d',alpha=.8)
    ax.add_patch(patch)
    path = Path([(x-r,y+r),(x+r,y-r)],[Path.MOVETO,Path.LINETO])
    patch = patches.PathPatch(path,fc='#4d4d4d',ec='#4d4d4d',alpha=.8)
    ax.add_patch(patch)

  # for each fixation, find nearest neighbor and distance to it (in pixels)
  X_err_i = 0
  X_err = 0.0
  distances = []
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)

    # use fixation point only if it is inside image borders
    # (image may be smaller than screen w,h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):

      # get intended monitor point
#     if ttype == 'center':
#       sx = monitor.S[0,0] * float(w)
#       sy = (1.0 - monitor.S[0,1])* float(h)
#     elif ttype == 'bottom_right':
#       sx = monitor.S[1,0] * float(w)
#       sy = (1.0 - monitor.S[1,1])* float(h)
#     elif ttype == 'bottom_left':
#       sx = monitor.S[2,0] * float(w)
#       sy = (1.0 - monitor.S[2,1])* float(h)
#     elif ttype == 'top_left':
#       sx = monitor.S[3,0] * float(w)
#       sy = (1.0 - monitor.S[3,1])* float(h)
#     elif ttype == 'top_right':
#       sx = monitor.S[4,0] * float(w)
#       sy = (1.0 - monitor.S[4,1])* float(h)
#     else:
#       print("WARNING: UNKNOWN NBACK POINT")

      sx = monitor.S[monitor.ndict[ttype],0] * float(w)
      sy = (1.0 - monitor.S[monitor.ndict[ttype],1])* float(h)

      dx = x - sx
      dy = y - sy
      dist = math.sqrt(dx*dx + dy*dy)

      # add in to list of distances
      distances.append(dist)

      # compute running mean for each calibration point as we encounter it
      X_err = float(X_err_i)/float(X_err_i+1) * X_err + \
                            1.0/float(X_err_i+1) * dist

      # increment how many times we've been close to this monitor point
      X_err_i += 1

# print("X_err: ", X_err)

  if len(fixations) > 0:
    # compute observed mean distance between each feature (fixation) and
    # its nearest neihgbor (calibration point)
    avedist = np.mean(distances)
    stddist = np.std(distances)
    print("mean dist b/ween fixations and monitor points: %f (pixels)" % (avedist))

    r = math.sqrt(float(w)*float(w) + float(h)*float(h))
    dpi = r/float(screen)

    D = float(viewdist)

    fov = 2*math.degrees(math.atan2(screen,2*D))
    fovx = 2*math.degrees(math.atan2(float(w)/dpi,2*D))
    fovy = 2*math.degrees(math.atan2(float(h)/dpi,2*D))

    avedist_deg = 2*math.degrees(math.atan2((avedist/dpi),(2*D)))
    stddist_deg = 2*math.degrees(math.atan2((stddist/dpi),(2*D)))

    print("mean dist b/ween fixations and monitor points: %f (degrees)" % (avedist_deg))
  
    strout = "view distance: %5.2f (inches), screen: %3.0f (inches), %5.2f$^{\circ}$ (visual angle), dpi: %5.2f" % \
             (D,screen,fov,dpi)
    ax.text(10,int(h)-40,strout,fontsize=10)
    strout = "mean error (accuracy): %5.2f$^{\circ}$ (degrees visual angle), standard deviation (precision): %5.2f$^{\circ}$" % \
             (avedist_deg, stddist_deg)
    ax.text(10,10,strout,fontsize=10)

  # draw arrowheads
  px = []
  py = []
  i=0
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):
      px.append(x)
      py.append(y)
      if i>0:
        sx = fixations[i-1].at(0) * float(w)
        # do the y-coord flip for rendering with (0,0) at bottom-left
        sy = (1.0 - fixations[i-1].at(1)) * float(h)
        if image is None or (x0 < sx and sx < x1 and y0 < sy and sy < y1):
          dx = x - sx
          dy = y - sy
#         arrow = plt.Arrow(sx,sy,dx,dy,\
#                           width=35,\
#                           alpha=.2,\
#                           fill=False,\
#                           fc="none",\
#                           ec="#101010")
#         ax.add_patch(arrow)
          opt = {'shape':"full",\
                 'fill':False,\
                 'fc':"none",\
                 'ec':"#101010",\
                 'width':1,\
                 'head_width':15,\
                 'head_length':45,\
                 'length_includes_head':True}
          if(abs(dx) > 0.0 and abs(dy) > 0.0):
            plt.arrow(sx,sy,dx,dy,alpha=.2,**opt)
      i += 1
    # old way, straight plot
# plt.plot(px,py,antialiased=True,alpha=.6,color="#756BB1")
  # or using Artist
  line = plt.Line2D(px,py,antialiased=True,alpha=.3,color="#494949",lw=1)
  ax.add_artist(line)

  # circles (corrected): won't work in this (n-back) case because
  # we only have samples at 1 "calibration point"
  for fx in fixations:
    x = fx.at(0) * float(w)
    # do the y-coord flip for rendering with (0,0) at bottom-left
    y = (1.0 - fx.at(1)) * float(h)
    if image is None or (x0 < x and x < x1 and y0 < y and y < y1):
      r = fx.getPercentDuration() * 100.0
      circ = plt.Circle((x,y),radius=r,fc='#fddbc7',ec='#393939',alpha=.2)
      ax.add_patch(circ)

  # title and labels
  plt.title(title)
  plt.ylabel("$y$-coordinate (pixels)",family="sans-serif")
  plt.xlabel("$x$-coordinate (pixels)",family="sans-serif")
  plt.grid(True,'major',ls='solid',alpha=.1)
  # margins
# plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.1)
  plt.tight_layout()

# fig = fileName[:-4] + ".pdf"
  plt.savefig(fileName,transparent=True)
  plt.close('all')

