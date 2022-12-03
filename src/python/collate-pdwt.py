#!/usr/bin/env python3

import sys, os, math

def catDATFile(infile,df,ct):
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
# header = linelist[0].split(',')
# linelist = linelist[1:]

  k = 0
  up = 0.0

  for line in linelist:
    entry = line.split(' ')

    # get line elements
    pdwt = float(entry[1])

    # compute running mean
    # from Brown, Robert Grover, "Introduction to Random Signal
    #   Analysis and Kalman Filtering", John Wiley & Sons, New York, NY
    #   1983 [p.182] TK5102.5 .B696
    up = float(k)/float(k+1) * up + 1.0/float(k+1) * pdwt
    k += 1

    str = "%s,%s,%s,%s,%s,%s" % ( \
                         subj, \
                         exp_id, \
                         ses_id, \
                         marker, \
                         object, \
                         up)
  print(str, file=df)
  ct += 1

  return ct

###############################################################################

# clear out output file
df = open("pdwt.csv",'w')
print("subj,exp_id,ses_id,marker,object,pdwt", file=df)

dir = './data/'

# find all files in dir with .csv extension
lst = [a for a in os.listdir(dir) if a.endswith('-pdwt.dat')]

lineno = 1

for item in lst:

  if "VALIDATION" in item:
    continue

  file = dir + item
  print('Processing ', file)

  # cat csv files into one
  lineno = catDATFile(file,df,lineno)

df.close()
