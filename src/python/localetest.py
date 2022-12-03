#!/usr/bin/env python3

import sys,os,getopt,glob
import locale
import numpy as np
import math
#######################
#locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' )
locale.setlocale( locale.LC_ALL, 'de_DE' )
#locale.setlocale( locale.LC_ALL, '' )

def main(argv):

  local = locale.localeconv()
  print("LC_NUMERIC: ", local['decimal_point'])

  print("%f" % (locale.atof("513,30")))

if __name__ == "__main__":
  main(sys.argv[1:])
