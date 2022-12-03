#!/usr/bin/env python3

# system includes
import sys
import os
#import xml.etree.ElementTree as ET
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import tostring

import lxml.etree
import lxml.builder

import numpy as np

def readxml():

  tree = ET.parse('output.xml')
  root = tree.getroot()

  # read from .xml file how many fields file name should have
  fields = root.find('schema')
  # len(fields) tells us how many
  print("schema has %d fields" % (len(fields)))
  for field in fields:
    name = field.get('name')
    posn = field.get('posn')
    print("name: ",name,"posn: ",posn)

def dumpxml():

  # we start with this to come up with an XML file for as many fields as
  # intended file name should have, e.g.,
  #
  # block_subj_school_shot_datetime
  #
  root = ET.Element("root")
  schema = ET.SubElement(root,"schema")
  ET.SubElement(schema,"field", name="block", posn='0')
  ET.SubElement(schema,"field", name="subj", posn='1')
  ET.SubElement(schema,"field", name="shot", posn='2')
  ET.SubElement(schema,"field", name="school", posn='3')
  ET.SubElement(schema,"field", name="datetime", posn='4')
  tree = ET.ElementTree(root)
  tree.write('output.xml')

# tree = lxml.builder.ElementMaker()
# root = tree.root
# doc = tree.doc
# f1 = tree.field1 
# f2 = tree.field2 
# f3 = tree.field3 
# f4 = tree.field4 
# f5 = tree.field5 
# the_doc = root(
#             doc(
#               f1('',name='blok'),
#               f2('',name='subj'),
#               f3('',name='skul'),
#               f4('',name='stim'),
#               f5('',name='date')
#             )
#           )
# print(lxml.etree.tostring(the_doc,pretty_print=True))

def main(argv):

  dumpxml()
  readxml()

if __name__ == "__main__":
  main(sys.argv[1:])
