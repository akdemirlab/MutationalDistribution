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
import seaborn as sns


def read_bedGraph(filename,resolution): # add stopping after certain chromosome passed
	
	'''
    reads bedGraph files for various file type plottings

    parameters:

    filename: file name. format could be either "chr\tstart\tend" or "chr\tstart\tend\tvalue..."
    resolution: bin size for the matrix
	
	returns:
	x_scores = location along the given chromosome - start sites
	x_scores2 = location along the given chromosome - end sites
	y_scores = signal scores for the assay
	colors = allow for colors option
    '''
	
	try:
		fone=open(filename,'r')
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	x_scores={}
	x_scores2={}
	y_scores={} 

	average = 0.0
	count = 1
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if line[0]!='t':
			if tags[0] not in x_scores.keys() and tags[0] != 'chrY':
				x_scores[tags[0]]=[]
				x_scores2[tags[0]]=[]
				y_scores[tags[0]]=[]
				x_scores[tags[0]].append(float(tags[1])/resolution)
				x_scores2[tags[0]].append(float(tags[2])/resolution)
				if len(tags) > 3: 
					if tags[3]=='.': y_scores[tags[0]].append(0.0)
					else: y_scores[tags[0]].append(float(tags[3]))
			elif tags[0] != 'chrY':
				x_scores[tags[0]].append(float(tags[1])/resolution)
				x_scores2[tags[0]].append(float(tags[2])/resolution)
				if len(tags) > 3: 
					if tags[3]=='.': y_scores[tags[0]].append(0.0)
					else: y_scores[tags[0]].append(float(tags[3]))
			
			average += (float(tags[2])-float(tags[1]))/resolution
			count +=1
	average = int(round(average/count,0))	
	
	return x_scores,x_scores2,y_scores,average


def read_bed(filename,resolution): # add stopping after certain chromosome passed
	
	'''
    reads bedGraph files for various file type plottings

    parameters:

    filename: file name. format could be either "chr\tstart\tend" or "chr\tstart\tend\tvalue..."
    resolution: bin size for the matrix
	
	returns:
	x_scores = location along the given chromosome - start sites
	x_scores2 = location along the given chromosome - end sites
	y_scores = signal scores for the assay
	colors = allow for colors option
    '''
	
	try:
		fone=open(filename,'r')
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	x_scores={}
	x_scores2={} 

	average = 0.0
	ave = []
	count = 1
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if line[0]=='#': continue
		if line[0]!='t':
			if tags[0] not in x_scores.keys() and tags[0] != 'chrY':
				x_scores[tags[0]]=[]
				x_scores2[tags[0]]=[]
				x_scores[tags[0]].append(float(tags[1])/resolution)
				x_scores2[tags[0]].append(float(tags[2])/resolution)
			elif tags[0] != 'chrY':
				x_scores[tags[0]].append(float(tags[1])/resolution)
				x_scores2[tags[0]].append(float(tags[2])/resolution)
			
			average+=(float(tags[2])-float(tags[1]))/resolution
			count +=1
	average = int(round(average/count,0))	
	
	return x_scores,x_scores2,average
	
def where(start,end,arr):
    """Find where the start location and end location indexes in an array"""
    
    astart = bisect.bisect_left(arr, start)
    aend = bisect.bisect_right(arr[start:], end) + start
        
    return astart, aend


def profiler(file1,file2): ## rechange to 20/41 from 10/21
	resolution = 25000
	x_comps,x_comps2,y_comps,_ = read_bedGraph(file1,resolution)
	
	x2_comps,x2_comps2,average = read_bed(file2,resolution)
	pseudocount = 0
	regions = []

	for key in x2_comps.keys():
		for item in range(0,len(x2_comps[key])):
			midpoint = int(round(x2_comps[key][item],0))+int(round((x2_comps2[key][item]-x2_comps[key][item])/2,0))
			if midpoint-20 > 0 : ystart,yend = where(midpoint-20,midpoint+20,x_comps[key])
			else: print 'Discarding a region with midpoint', midpoint #ystart,yend = where(0,midpoint+20,x_comps[key])
			
 			if len(y_comps[key][ystart:yend]) < 41:
 				continue
 				print 'Short region in', file2
 				if len(regions)==0 : regions = np.lib.pad(y_comps[key][ystart:yend], (0,21-len(y_comps[key][ystart:yend])), 'constant', constant_values=(0))
 				else: regions = np.vstack([regions,np.lib.pad(y_comps[key][ystart:yend], (0,21-len(y_comps[key][ystart:yend])), 'constant', constant_values=(0))])
 			else:	
 				if len(regions)==0 : regions = y_comps[key][ystart:yend]#;writer.write(str(y_comps[key][ystart:yend]));writer.write('\n')
 				else: regions = np.vstack([regions,np.array(y_comps[key][ystart:yend])])#;writer.write(str(y_comps[key][ystart:yend]));writer.write('\n')

	final = np.mean(regions,axis=0)
	return final,average

def ind_profiler(file1,file2,write): ## rechange to 20/41 from 10/21
	
	resolution = 25000
	x_comps,x_comps2,y_comps,_ = read_bedGraph(file1,resolution)
	
	x2_comps,x2_comps2,average = read_bed(file2,resolution)
	pseudocount = 0
	regions = []
	upstream = []
	downstream = []
	if write=='yes':
		writer = open(file1+file2+".txt",'w')
	for key in x2_comps.keys():
		for item in range(0,len(x2_comps[key])):
			midpoint = int(round(x2_comps[key][item],0))+int(round((x2_comps2[key][item]-x2_comps[key][item])/2,0))
			if midpoint-20 > 0 : ystart,yend = where(midpoint-20,midpoint+20,x_comps[key])
			else: print 'Discarding a region with midpoint', midpoint #ystart,yend = where(0,midpoint+20,x_comps[key])
			
			#print len(y_comps[key][ystart:yend])
			if write=='yes':
				median = np.median(np.array(y_comps[key][ystart:yend]))
				divided = np.divide(np.array(y_comps[key][ystart:yend]),float(median))
				logged = np.log2(divided)
				print median
				logged.tofile(writer,sep='\t')
				writer.write('\n')
 			if len(y_comps[key][ystart:yend]) < 41:
 				print 'Short region in', file2
 				
 			else:
 				upstream.append(np.average(y_comps[key][ystart+4:ystart+18]))
 				downstream.append(np.average(y_comps[key][ystart+22:ystart+36]))
 	
	return upstream,downstream

def mut_profiler():
	

	fig, ax = plt.subplots(1, 1)

	mutations = 'Pancancer-Cumulative-MutBurden.bedGraph'
	
	boundaries = 'ActiveToInactive-Boundaries.bed'
	

	p_average3,l_average3 = profiler(mutations,boundaries)

	x = np.linspace(0, len(p_average3), 1)
	
	ax.plot(p_average3,linewidth = 4, color='#008080',label='Mutations')
	ax.set_ylabel('Mut Load / 25 Kb',fontsize=20)
	ax.locator_params(axis='y',tight=False, nbins=5)
	ax.locator_params(axis='x',tight=False, nbins=3)
	ticks= ax.get_xticks().tolist()
	ticks = ['-500Kb','Boundary','+500Kb']
	ax.set_xticklabels(ticks,fontsize=20)
	ax.axvspan(0,20,facecolor='#fc9929',alpha=0.45)
	ax.axvspan(20,40,facecolor='#cbcbcb',alpha=0.45)
	ax.set_xlim(0,40)
	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	
	plt.savefig('MutationDistribution.png',dpi=200)
	
mut_profiler()	
