#!/usr/bin/env python3

# system includes
import sys
import os
#import xml.etree.ElementTree as ET
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import tostring
import numpy as np

# local includes
from aoi import AOI

def parseAOI(aoifile,df):

# aoidict = {}
  aoilist = []

  if(os.path.isfile(aoifile)):

    print("parsing: ", aoifile)
    path, base = os.path.split(aoifile)
    print("Processing: ", aoifile, "[", base, "]")
    filename, ext = os.path.splitext(base)
    print("path, base, filename, ext: ", path, base, filename, ext)
    stim = filename.strip()
    print("stim: ", stim)

    tree = ET.parse(aoifile)
    #print "tree = %s" % tree

    # should be something like Experiment, {}
    root = tree.getroot()
  # print "root.tag = %s, root.attrib = %s" % (root.tag, root.attrib)

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
      print("%s %s %s %s" % (px,py,pw,ph))
      print("%s %s %s %s" % (bl,br,bt,bb))

    print("IMAGE INFO:")
    for obj in root.iter('PAGEOBJECT'):
      if obj.get('PTYPE') == '2':
        x = obj.get('XPOS')
        y = obj.get('YPOS')
        w = obj.get('WIDTH')
        h = obj.get('HEIGHT')
        print("%s %s %s %s" % (x,y,w,h))

    # iterate through PAGEOBJECT objects, look for PTYPE=''6'
    for obj in root.iter('PAGEOBJECT'):
      if obj.get('PTYPE') == '6':

#       print("obj tag, attr")
        x = obj.get('XPOS')
        y = obj.get('YPOS')
        w = obj.get('WIDTH')
        h = obj.get('HEIGHT')
# gone in v1.5.6
#       label = obj.get('ANNAME')
# new in v1.5.6: attributes
#
# location   : integer (which AOI in sequence)
# key        : boolean (key word or not)
# syllables  : integer (how many syllables)
# frequency  : string  (not sure what this is)
# length     : integer (length of word)
# function   : boolean (is it a function or not?)
# length/key : string  ('long' or 'short' only on key AOIs)
        aoiattr = {}
#       aoiattr['label'] = 'None'
        # <PageItemAttributes>
        pattr = obj.find('PageItemAttributes')
#       print("len(pattr): ", len(pattr))
        for attr in pattr:
          # <ItemAttribute Name="location" Type="Integer" Value="1" ... />
          name = attr.get('Name')
          type = attr.get('Type')
          value = attr.get('Value')
#         print("name: ",name)
#         print("type: ",type)
#         print("value: ",value)
          # this is the key idea: set up dict from XML file as per above
          aoiattr[name] = value
        # </PageItemAttributes>

        x = str(float(x) - float(px))
        y = str(float(y) - float(py))
#       print("%s %s %s %s" % (x,y,w,h))

#       for (0,0) at bottom
#       aoi = AOI(x,y,w,h-y)
#       for (0,0) at top
        aoi = AOI(x,y,w,h)
        aoi.setAttributes(aoiattr)
#       aoi.dump()
#       if aoidict.has_key(stimulus):
#         aoidict[stimulus].append(aoi)
#       else:
#         aoidict[stimulus] = [aoi]
        aoilist.append(aoi)

        del aoi

# print aoidict
# for key in aoidict:
#   print key
#   print "number of AOIs: ",len(aoidict[key])
#   for aoi in aoidict[key]:
#     aoi.dump()
  print("number of AOIs: ",len(aoilist))
  for i, aoi in enumerate(aoilist):
    print("================= AOI %d:" % (i+1))
    # see aoi.y for how to get at key, value pairs
    aoi.dump()
    attrstr = aoi.getAttributes_astring()
    print(stim + ',' + attrstr,file=df)

def main(argv):

# aoi = AOI(179,764.0,211,108)
# aoi.dump()
# print("in" if aoi.inside(180.0,780) else "out")
# print("in" if aoi.inside(0.0,0.0) else "out")

  # process AOI file
  aoidir = "../../../../aois/1920x1080/"
# aoifile = aoidir + "1.sla"

  file = None
  files = []

  # get .csv input files to process
  if os.path.isdir(aoidir):
#   files = glob.glob('%s/*.sla' % (aoidir))
    for top, dirs, listing in os.walk(aoidir):
      for file in listing:
        path = os.path.join(top,file)
        dname = path.split('/')[2]
        base = os.path.basename(path)
        if base.endswith('.sla'):
          basename = base[:base.rfind('.sla')]
#         print(path, top, dname, base, "[",basename,"]","[",file,"]")
          files.extend([path])

  # if user specified --file="..." then we use that as the only one to process
  if(file != None and os.path.isfile(file)):
    files = [file]

  df = open("aois.csv",'w')
  print("stim,aoi",file=df)

  for aoifile in files:
    print("aoifile = ", aoifile)
    parseAOI(aoifile,df)

  df.close()

if __name__ == "__main__":
  main(sys.argv[1:])
