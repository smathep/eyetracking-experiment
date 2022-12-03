#!/usr/bin/env python3

import sys,os,getopt,glob
import locale
import numpy as np
import math
from tqdm import tqdm
import csv
import h5py
import tables

def usage():
  print("Usage: python hdf52raw.py " \
        " --indir=? -outdir=?\n" \
        " --width=? -height=?\n" \
        " --file=?\n" \
        "   indir: a directory containing input files to be processed\n" \
        "   outdir: a directory containing output files\n" \
        "   width: horiz. gaze coord extents\n" \
        "   height: vert. gaze coord extents\n" \
        "   file: a single file to process\n")

# convert csv file to raw
def hdf52raw(infile,outdir,width,height,dist):

  try:
    h5file = tables.open_file(infile,mode='r')
  except IOError:
    print("Can't open file with tables.open_file: " + infile)
    return

  infile = infile.replace('\\','/')

  # parse the subject name from the file name
# base = infile[:infile.rfind('.')]
  base = os.path.basename(infile)
  dname = os.path.dirname(infile)
  print("Processing: ", infile, "[", base, "]", "<", dname, ">")

  # strip out extension
  filename, ext = os.path.splitext(base)
  print('Processing: ' + filename)

  # subj id from filename
  subj = filename.split('_')[0]
  print('subj: ', subj)

  ################ open csv file

  # csvfile (expected)
  incsvfile = dname + '/' + filename + '.csv'

  if os.path.isfile(incsvfile):
    print('csv file: ' + incsvfile + ' exists')
    # utf-8-sig encoding removes leading '\ufeff'
    csvfile = open(incsvfile,newline='',encoding='utf-8-sig')
    # get number of rows, skipping header
#   nrows = sum(1 for _ in csvfile) - 1
#   print("nrows: ",nrows)
    csvfile.seek(0)
    csvrows = csv.DictReader(csvfile,delimiter=',')
    csvlist = list(csvrows)
    # print header
#   print("header: ",csvrows.fieldnames)
  else:
    print('csv file: ' + incsvfile + ' does NOT exist')

  ############### init outfile

  outfile = None

  ################ open HDF5 file

  # processing hdf5 file using pytables
# print("h5file: ", h5file)

  # find session_data_data...
  metatable = h5file.root.data_collection.session_meta_data
# print("metatable: ", metatable)

  # first use class_table_mapping to get at paths to various datasets
  pathtable = h5file.root.class_table_mapping
# print("pathtable: ",pathtable)
# for rec in pathtable.iterrows():
#   print("rec['table_path']: ", rec['table_path'].decode('utf-8'))
# bespath = [rec['table_path'] for rec in \
#           pathtable.where("""(class_name == 'BinocularEyeSampleEvent')""")]
# print("bespath: ", bespath)

  # from https://www.pytables.org/usersguide/tutorials.html
  # rows in class_table_mapping table give us path to other data objects

  # find gaze data...
  # in theory there may be more than one such node but not likely in practice
# for record in pathtable.where('(class_name == b"BinocularEyeSampleEvent")'):
  for record in pathtable.where('(class_name == b"GazepointSampleEvent")'):
    # pull out row from class_table_mapping to GazepointSampleEvent node 
    bespath = record['table_path'].decode('utf-8')
    # get table
    bestable = h5file.get_node(bespath)
#   print("bestable: ", bestable)

  # this is if we manually add in code to issue messages in PsychoPy
  # find message data...
  # in theory there may be more than one such node but not likely in practice
  for record in pathtable.where('(class_name == b"MessageEvent")'):
    # pull out row from class_table_mapping to MessageEvent node 
    mespath = record['table_path'].decode('utf-8')
    # get table
    mestable = h5file.get_node(mespath)
#   print("mestable: ", mestable)

  # init time vars
  st = 0.0
  et = math.inf

  # init condition vars
  marker = None
  object = None

#  # alternative method: direct SQL-like query
#  for row in mestable.where('(category == "trial") & \
#                                 (text == "start")' \
#                           ):
#    st = row['time']
#    exp_id = row['experiment_id']
#    ses_id = row['session_id']
#
#  for row in mestable.where('(category == "trial") & \
#                                 (text == "stop")' \
#                           ):
#    et = row['time']
#    exp_id = row['experiment_id']
#    ses_id = row['session_id']

  # this is if we manually add in code to issue messages in PsychoPy
  # process each pair of rows in message table: they contain start/end times
  for msg in mestable.iterrows():
    # get record info: first of pair of lines should be
    # category:'trial', text:'start'
    if msg['category'].decode('utf-8') == 'trial' and \
       msg['text'].decode('utf-8') == 'start':

      # pick up start time
      st = msg['time']

      # pick up experiment_id and session_id
      exp_id = msg['experiment_id']
      ses_id = msg['session_id']

    # these are additional things we sent here from PsychoPy (conditions)
    if msg['category'].decode('utf-8') == 'marker':
      marker = msg['text'].decode('utf-8')
      marker = marker.replace("'","")

    if msg['category'].decode('utf-8') == 'object':
      object = msg['text'].decode('utf-8')
      object = object.replace("'","")

    # get record info: second of pair of lines should be 'end'
    if msg['category'].decode('utf-8') == 'trial' and \
       msg['text'].decode('utf-8') == 'stop':
      et = msg['time']

      # if we got 'trial':stop then we have st, et, and condition vars
      # dump out eye movement data
 
      print("exp_id: ", exp_id)
      print("ses_id: ", ses_id)
      print("marker: ", marker)
      print("object: ", object)
      print("st: ", st)
      print("et: ", et)

      outfilename = "%s-%s-%s-%s-%s.raw" % (subj,exp_id,ses_id,marker,object)
      outfile = open(outdir + outfilename,'w+')
      print('Outfile: %s' % (outfilename))
      # write header (not for legacy .raw files)
#     strout = "x,y,d,t"
#     outfile.write(strout + '\n')

#     for row in bestable.where('(experiment_id == exp_id) & \
#                                (session_id == ses_id) & \
#                                (st <= time)  &  (time <= et)' \
      for row in bestable.where('(st <= time)  &  (time <= et)' \
                                 ):
#       print("%d %d %f %f %f %f %f %f %f" % (\
#             row['experiment_id'], \
#             row['session_id'], \
#             row['time'], \
#             row['left_gaze_x'], \
#             row['left_gaze_y'], \
#             row['left_pupil_measure1'], \
#             row['right_gaze_x'], \
#             row['right_gaze_y'], \
#             row['right_pupil_measure1']))

# see:
# https://www.psychopy.org/api/iohub/device/eyetracker_interface/GazePoint_Implementation_Notes.html#psychopy.iohub.devices.eyetracker.BinocularEyeSampleEvent
        t = row['time']
# these are the eye coords of the pupils in the camera images, don't use
#       x_l = row['left_raw_x']
#       y_l = row['left_raw_y']
#       x_r = row['right_raw_x']
#       y_r = row['right_raw_y']
# gaze position in psychopy window in the window coord unit being used
# (e.g., norm, height, pix, etc.)
        x_l = row['left_gaze_x']
        y_l = row['left_gaze_y']
        x_r = row['right_gaze_x']
        y_r = row['right_gaze_y']
        d_l = row['left_pupil_measure2']
        d_r = row['right_pupil_measure2']
        v = row['status']

        # height unit mapping [0,w,0,h] -> [-.5*(w/h),.5*(w/h),-.5,.5]
        # x' = ( (x - 0)/w -.5 )*(w/h)
        # y' = ( (y - 0)/h -.5 )*(h/h)

        # convert height units to normalized coords
        # [-.5*(h/w), .5*(h/w)] -> [0,1]
        x_l = (x_l * height/width + 0.5)
        x_r = (x_r * height/width + 0.5)

        # [-.5, .5] -> [0,1]
        y_l = (y_l * 1.0 + 0.5)
        y_r = (y_r * 1.0 + 0.5)

        # y-flip
        y_l = 1.0 - y_l
        y_r = 1.0 - y_r

        x = (x_l + x_r)/2.0
        y = (y_l + y_r)/2.0
        d = (d_l + d_r)/2.0

#       if(v < 6 and \
#          x < 1 and y < 1 and \
#          x >= 0 and y >= 0):
        if(v < 6):
          strout = "%f %f %f %f" % (x,y,d,t)
#         strout = "%f %f %f %f" % (x,1-y,d,t)
          outfile.write(strout + '\n')

      if outfile is not None:
        outfile.close()

  h5file.close()

  return

  """ this is my first attempt using h5py...ended up being too low-level
      above attempt using pytables was much easier
  """

  try:
    f = h5py.File(infile,'r')
  except IOError:
    print("Can't open file with h5py: " + infile)
    return

  # processing hdf5 file using h5py
  print("keys:")
  print(f.keys())

  if 'class_table_mapping' in f.keys():
    ctm = f['class_table_mapping']
    print("ctm: ", ctm)
#   tpath = f['class_table_mapping'].value # deprecated
    tpath = f['class_table_mapping'][()]
    print("tpath: ", tpath)

# df = f['data_collection/session_meta_data']
# print(df)

  # attributes
  print("attributes:")
  for item in f.attrs.keys():
    print(item + ":", f.attrs[item])

  # meta data
  print("data_collection/session_meta_data:")
  md = f['data_collection/session_meta_data']
  for i in range(len(md)):
#   print("%d\t %g" % (i, md[i]))
    print(md[i])

  # eye tracker
  print("data_collection/events/eyetracker:")
  print(f['data_collection/events/eyetracker'])
  for g in f['data_collection/events/eyetracker']:
    print(g)
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['experiment_id'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['session_id'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['time'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['status'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['left_raw_x'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['left_raw_y'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['left_pupil_measure1'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['right_raw_x'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['right_raw_y'])
  print(f['data_collection/events/eyetracker/BinocularEyeSampleEvent']['right_pupil_measure1'])

  # gaze data -- but how to query it?
  gd = f['data_collection/events/eyetracker/BinocularEyeSampleEvent']
  print(gd)
  print("len(gd): ",len(gd))
  print("gd.shape: ",gd.shape)

# for i in range(len(gd)):
#   print("%d\t %g" % (i, gd[i]))
#   print(gd['time'][i])
#   print(gd['status'][i])
#   print(gd['left_raw_x'][i])
#   print(gd['left_raw_y'][i])
#   print(gd['left_pupil_measure1'][i])
#   print(gd['right_raw_x'][i])
#   print(gd['right_raw_y'][i])
#   print(gd['right_pupil_measure1'][i])
#   print('')

def main(argv):
# if not len(argv):
#   usage()

  try:
    opts, args = getopt.getopt(argv, '', \
                 ['indir=','outdir=','file=',\
                  'width=','height=','dist=',\
                  'hertz=','sfdegree=','sfcutoff='])
  except getopt.GetoptError as err:
    usage()
    exit()

  file = None
  files = []
  indir = './'
  outdir = './'

  for opt,arg in opts:
    opt = opt.lower()
    if(opt != '--file' and opt != '--indir'):
      arg = arg.lower()

    if opt == '--indir':
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

  # get .hdf5 input files to process
  if os.path.isdir(indir):
#   files = glob.glob('%s/*.hdf5' % (indir))
    for top, dirs, listing in os.walk(indir):

#     print("top: %s" % (top))
#     print("dirs: %s" % (dirs))
#     print("listing: %s" % (listing))

      for file in listing:
        path = os.path.join(top,file)
        # get second-to-last path entry
        dname = path.split('/')[-2]
        base = os.path.basename(path)

#       print("path: %s" % (path))
#       print("dname: %s" % (dname))
#       print("base: %s" % (base))

        if base.endswith('.hdf5'):
          basename = base[:base.rfind('.hdf5')]
#         print path, top, dname, base, "[",basename,"]","[",file,"]"
          files.extend([path])

  # if user specified --file="..." then we use that as the only one to process
  if(file != None and os.path.isfile(file)):
    files = [file]

  # set up progress bar
# pbar = tqdm(total=len(files))

  for file in files:
    hdf52raw(file,outdir,width,height,dist)
    # update progress bar
#   pbar.update(1)


#######################
#locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' )

if __name__ == "__main__":
  main(sys.argv[1:])
