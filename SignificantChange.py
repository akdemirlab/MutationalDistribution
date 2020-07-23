import os,sys
from pylab import *
from math  import sqrt, isnan, floor, ceil, pi
from numpy import log2, log10, array, max
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker  import MultipleLocator
from matplotlib.patches import Polygon, Rectangle, Circle
from scipy.signal import argrelextrema
from scipy import ndimage
from matplotlib import pyplot as plt
import numpy as np
import argparse
import bisect
import logging
import pybedtools
import glob
from shutil import copyfile
import scipy.stats
from scipy import stats
from matplotlib.ticker import FormatStrFormatter
from scipy.signal import savgol_filter
from numpy import diff
from scipy.signal import find_peaks

def derivative(mutation='',output=''):	

	chrs = {}
	data = []
	doms = {}
	
	resolution = 25000
	
	boun = {}
	prev = 100
	fone=open(mutation,'r')
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if tags[0] not in chrs.keys():
			chrs[tags[0]]=[]
			chrs[tags[0]].append(float(tags[3]))
		else: 
			chrs[tags[0]].append(float(tags[3]))
		prev = float(tags[3])

	writer = open(output+'-secondDerivative.bed','w')
	
	start = 1
	end = 22
			
	for item in range(start,end):
		
		name="chr"+str(item)
		
		data = chrs["chr"+str(item)]
		yhat = savgol_filter(data, 7, 3)
		peaks, _ = find_peaks(yhat, height=6)

 #		for xx in range(0,len(yhat)):
 #  			writer.write("chr"+str(item)+'\t'+str(xx*resolution)+'\t'+str(xx*resolution)+'\t'+str(yhat[xx])+'\n') 		

  		for xx in range(0,len(peaks)):
  			writer.write("chr"+str(item)+'\t'+str(peaks[xx]*resolution)+'\t'+str(peaks[xx]*resolution)+'\t'+str(data[peaks[xx]])+'\n') 		
 
if __name__=='__main__':
	
	parser = argparse.ArgumentParser(usage='SignificantChange.py -m Mutation.bedGraph -o output',add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
	
	group = parser.add_argument_group("Required Parameters")
	group.add_argument('-m','--mutation',default='', help='',metavar='',required=True)
	group.add_argument('-o', '--output',default='',metavar='',required=True)
	
	group1 = parser.add_argument_group("Optional Parameters")
	group1.add_argument('-h', '--help', action="help")
# 	group1.add_argument('-r', '--region',default='',metavar='',help='')
# 	group1.add_argument('-d', '--domain',default='',metavar='',help='')
# 	group1.add_argument('-f', '--folder',default='',metavar='',help='')
	
	args = vars(parser.parse_args())
	
	derivative(**args)