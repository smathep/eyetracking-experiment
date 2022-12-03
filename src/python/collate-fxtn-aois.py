#!/usr/bin/env python3

import sys, os, math

def catCSVFile(infile,df,ct):
  try:
    f = open(infile,'r')
  except IOError:
    print("Can't open file: " + infile)
    return

  path, base = os.path.split(infile)

  print("Processing: ", infile, "[", base, "]")

  # split filename from extension
  filename, ext = os.path.splitext(base)

  print("path, base, filename, ext: ", path, base, filename, ext)

  # extract stimulus name and subj id
  # filename now has the form 'date-rest_of_it', extract just the second part
  subj = filename.split('-')[0]
  exp_id = filename.split('-')[1]
  ses_id = filename.split('-')[2]
  marker = filename.split('-')[3]
  object = filename.split('-')[4]
  print("subj, exp_id, ses_id, marker, object: ", \
         subj, exp_id, ses_id, marker, object)

  # read lines, throwing away first one (header)
# linelist = f.readlines()
# linelist = f.read().split('\r')
  linelist = f.read().splitlines()
  header = linelist[0].split(',')
  linelist = linelist[1:]

  # timestamp,x,y,duration,prev_sacc_amplitude,aoi_label
  for idx, label in enumerate(header):
    if label.strip() == "timestamp":
      TIMESTAMP = idx
    if label.strip() == "x":
      X = idx
    if label.strip() == "y":
      Y = idx
    if label.strip() == "duration":
      DURATION = idx
    if label.strip() == "prev_sacc_amplitude":
      PREV_SACC_AMPLITUDE = idx
    if label.strip() == "aoi_span":
      AOI_SPAN = idx
    if label.strip() == "aoi_label":
      AOI_LABEL = idx
    if label.strip() == "aoi_order":
      AOI_ORDER = idx

  for line in linelist:
    entry = line.split(',')

    # get line elements
    timestamp = entry[TIMESTAMP]
    x  = entry[X]
    y  = entry[Y]
    duration  = entry[DURATION]
    prev_sacc_amplitude  = entry[PREV_SACC_AMPLITUDE]
    aoi_span  = entry[AOI_SPAN]
    aoi_label  = entry[AOI_LABEL]
    aoi_order  = entry[AOI_ORDER]

    str = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % ( \
                         subj, \
                         exp_id, \
                         ses_id, \
                         marker, \
                         object, \
                         timestamp,\
                         x,y,\
                         duration,\
                         prev_sacc_amplitude,\
                         aoi_span,\
                         aoi_label,\
                         aoi_order,\
                         ct)
    print(str, file=df)
    ct += 1

  return ct

###############################################################################

# clear out output file
df = open("fxtn-aois.csv",'w')
print("subj,exp_id,ses_id,marker,object,timestamp,x,y,duration,prev_sacc_amplitude,aoi_span,aoi_label,aoi_order,order", file=df)

dir = './data/'

# find all files in dir with .csv extension
lst = [a for a in os.listdir(dir) if a.endswith('-fxtn-aoi.csv')]

lineno = 1

for item in lst:

  file = dir + item
  print('Processing ', file)

  # cat csv files into one
  lineno = catCSVFile(file,df,lineno)

df.close()
