#!/usr/bin/env python3

import platform,sys,os,getopt,glob
import locale
import numpy as np
import csv
import math
from statistics import mode

def is_odd(a):
  return bool(a - ((a>>1)<<1))

def isfloat(value):
  try:
#   float(value)
    float(locale.atof(value))
    return True
  except:
    return False

def usage():
  print("Usage: python csv2raw.py " \
        " --indir=? -outdir=?\n" \
        " --width=? -height=?\n" \
        " --group=?\n" \
        " --sheet=?\n" \
        " --file=?\n" \
        "   indir: a directory containing input files to be processed\n" \
        "   outdir: a directory containing output files\n" \
        "   width: horiz. gaze coord extents\n" \
        "   height: vert. gaze coord extents\n" \
        "   group: which data set\n" \
        "   sheet: which sheet in Excel spreadsheet\n" \
        "   file: a single file to process\n")

def readidx(infile):

  infile = infile.replace('\\','/')

  # get infile base
# base = infile[:infile.rfind('.')]
  base = os.path.basename(infile)
  print("Processing: ", infile, "[", base, "]")

  # strip out extension
  filename, ext = os.path.splitext(base)
  print('Processing: ' + filename)

  idxdict = []

  csvfile = open(infile,newline='')
  
  # using csv DictReader
  rows = csv.DictReader(csvfile)

  # here we read through the .csv file and copy into a separate list of dicts
  # this seems redundant since why not just process through the .csv file
  # directly?
  # we can only do this if the .csv file stays open...I want to pass the
  # dict around
  for row in rows:
#   print(row)
#   print(row['file'],row['shot'],row['cond'],row['subj'],row['ordr'])

    elem = {}
    elem['file'] = row['file']
    elem['shot'] = row['shot']
    elem['cond'] = row['cond']
    elem['subj'] = row['subj']
    elem['ordr'] = row['ordr']
    idxdict.append(elem)

  csvfile.close()

  return idxdict

def dump(data):

# print(data)
  print("Number of entries: ",len(data))

  for row in data:
#   print(row)
    print(row['file'],row['shot'],row['cond'],row['subj'],row['ordr'])

# convert asc file to raw
def csv2raw(idxdict,indir,outdir,width,height,dist):

  # idxdict is the index.csv file which has file names, shot, cond, subj info
  for row in idxdict:

    file = row['file']
    shot = row['shot']
    cond = row['cond']
    subj = row['subj']
    ordr = row['ordr']

    # skip 'Marta', just a test file
    if subj == 'Marta':
      continue

    # get input file name (file with data from VR run)
    infile = indir + '/' + file + '_AOI.csv'

    print("infile: ", infile)

    outfile = None

    print("subj, shot, cond, ordr: ", \
           subj, shot, cond, ordr)

    # set up new data file name
    oname = "%s/%s-%s-%s-%s.raw" % \
            (outdir, subj, shot, cond, ordr)

    print("oname: ", oname)

    if os.path.isfile(infile):

      # open output file
      outfile = open(oname,'w+')

      # write out header
#     header = "x,y,d,t"
#     outfile.write(header + '\n')

      csvfile = open(infile,newline='')
      lines = csv.DictReader(csvfile)

      # print header
#     print(lines.fieldnames)

      for line in lines:

        t = eval(line['timestamp'])
        x = eval(line['x'])
        y = eval(line['y'])
        fxd = eval(line['fixation_duration'])
        amp = eval(line['amplitude'])
        sxd = eval(line['sacc_duration'])
        vel = eval(line['velocity'])
        d_l = eval(line['pupil_d_left'])
        d_r = eval(line['pupil_d_right'])
        grp = eval(line['aoi_group'])
        aoi = eval(line['aoi'])
        wrd = eval(line['word'])

        # pupil diameter in mm
        d = (d_l + d_r)/2.0

        if float(x) > 0.0 and float(y) > 0.0:
          strout = "%f %f %f %f" % \
                    (x,y,d,t)
          outfile.write(strout + '\n')

      outfile.close()

def main(argv):
# if not len(argv):
#   usage()

  try:
    opts, args = getopt.getopt(argv, '', \
                 ['indir=','outdir=',\
                  'idxfile=',\
                  'width=','height=','dist='])
  except getopt.GetoptError as err:
    usage()
    exit()

  idxfile = None

  for opt,arg in opts:
    opt = opt.lower()
    if opt != '--indir':
      arg = arg.lower()

    if opt == '--idxfile':
      idxfile = arg
    elif opt == '--indir':
      indir = arg
    elif opt == '--outdir':
      outdir = arg
    elif opt == '--width':
      width = float(arg)
    elif opt == '--height':
      height = float(arg)
    elif opt == '--dist':
      # convert distance to cm and then to mm
      dist = float(arg) * 2.54 * 10.0
    else:
      sys.argv[1:]

  if idxfile is not None:
    idxdict = readidx(idxfile)
    dump(idxdict)
    csv2raw(idxdict,indir,outdir,width,height,dist)

#######################
# on Windows, see bottom of this web page:
# http://stackoverflow.com/questions/19709026/how-can-i-list-all-available-windows-locales-in-python-console

locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' )
#locale.setlocale( locale.LC_ALL, 'de_DE' )
#if platform.system() == 'Windows':
#  locale.setlocale( locale.LC_ALL, 'German' )
#else:
#  locale.setlocale( locale.LC_ALL, 'de_DE' )

if __name__ == "__main__":
  main(sys.argv[1:])
