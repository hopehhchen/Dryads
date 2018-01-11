import sys
import os
from collections import *

#
import numpy as np
import scipy

#
from astropy.io import fits
import astropy.wcs as wcs
import astropy.units as u
import astropy.constants as c
import FITS_tools as fits_tools

#
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from matplotlib import rcParams

#
import pandas as pd

##
from Dryads import Dryads
import styles


##### template from Droplets
#from Droplets import *
#from stat_tools import *
#from constants import *
#import styles


#### Data from the GAS DR1 and from Herschel ####
# data folder (within the current directory)
#direcData = os.getcwd()+'/data/'
