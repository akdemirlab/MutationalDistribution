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
#import seaborn as sns
from matplotlib import colors as clrs


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


def profiler(file1,file2):
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


def mut_profiler(mutations='',boundaries='',output='',name=''):
	
	fig, ax = plt.subplots(1, 1)	

	p_average3,l_average3 = profiler(mutations,boundaries)

	x = np.linspace(0, len(p_average3), 1)
	
	ax.plot(p_average3,linewidth = 4, color='#008080',label='Mutations')
	ax.set_ylabel('Mut Load / 25 Kb',fontsize=20)
	ax.locator_params(axis='y',tight=False, nbins=5)
	ax.locator_params(axis='x',tight=False, nbins=3)
	ticks= ax.get_xticks().tolist()
	ticks = ['-500Kb','Boundary','+500Kb']
	ax.set_xticklabels(ticks,fontsize=20)
	#ax.axvspan(0,20,facecolor='#fc9929',alpha=0.45)
	#ax.axvspan(20,40,facecolor='#cbcbcb',alpha=0.45)
	ax.set_xlim(0,40)
	
	ax.set_title(name)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	
	plt.savefig(output+'.png',dpi=200)
	
	
def dom_profiler(mutations='',domain='',output='',name=''):
	
	
	mutations = pybedtools.BedTool(mutations).sort()
	doms = pybedtools.BedTool(domain).sort()
	overlap = doms.map(mutations,c=4,o='sum')
	
	fig, ax = plt.subplots(1, 1)

	domains = {}
	data = []
	
	for tags in overlap:
		if tags[3] not in domains.keys():
			domains[tags[3]]=[]
			domains[tags[3]].append(float(tags[5])*25000/(int(tags[2])-int(tags[1])))
		else: 
			domains[tags[3]].append(float(tags[5])*25000/(int(tags[2])-int(tags[1])))

	data = [domains['0'],domains['1'],domains['2'],domains['3'],domains['4']]

	colors=["#8a91d0","#cbcbcb","#4da6ff","#fc9929","#ff4444"]

	for i in [1,2,3,4,5]:
		y = data[i-1]
		x = np.random.normal(i, 0.04, len(y))
		ax.plot(x, y, colors[i-1], alpha=0.4,mec='k', ms=7, marker="o", linestyle="None")

	 
	ax.set_facecolor("white") 
	#if counter != len(files)-1: plt.setp(ax.get_xticklabels(), visible=False)
	#ax.set_title(name,fontsize=10)
	#ax.get_yaxis().set_label_coords(-0.2,0.5)
	ax.locator_params(axis='y',tight=False, nbins=4)
	
	print scipy.stats.ranksums(domains['1'],domains['3'])

	box = ax.boxplot(data, notch=False, patch_artist=True,widths=(0.6,0.6,0.6,0.6,0.6),showfliers=False)
	for patch, color in zip(box['boxes'], colors):
		patch.set_facecolor('none');patch.set_edgecolor('none');
	for median in box['medians']: median.set(color='black', linewidth=2)
	#for median in box['fliers']: median.set(color='gray', marker='o')
	plt.setp(box['whiskers'], color='black')	
	
	doms = ['Heterochromatin','Inactive','Repressed','Active','Active-2']
	ax.set_xticks(range(1,6,1))
 	ticks= ax.get_xticks().tolist()
	for ditem in range(0,len(ticks)): ticks[ditem]=doms[ditem]	
 	ax.set_xticklabels(ticks, fontsize=8, fontweight='bold')

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.set_title(name)
	
	plt.savefig(output+'.png',dpi=200)	


def Plotter(mutation='',output='',region='',domain='',folder='',name=''):

	if len(mutation)>0 and len(domain)==0 and len(folder)==0: mut_profiler(mutation,region,output,name)
	if len(mutation)>0 and len(domain)>0: dom_profiler(mutation,domain,output,name)
	if len(folder)>0: 
		
		tog,l_average3 = profiler(mutation,region)
		
		files = glob.glob(folder+"/*.bedGraph")
		data = []
	
		for item in files:
			p_average3,l_average3 = profiler(item,region)
 		 	if len(data)==0 : data = p_average3[:]
	 		else: data = np.vstack([data,np.array(p_average3[:])])
	 			
		x = np.linspace(0, len(data[0]), 1)
		norm_data = (data - np.mean(data, axis=1)[:, np.newaxis]) / np.ptp(data, axis=1)[:, np.newaxis]
		scmap = clrs.ListedColormap(['#f662ff', '#55adff'])
		bounds=[0,1]
	
		fig, (ax,ax1) = plt.subplots(2,1,figsize=(8,16))
		fig.subplots_adjust(hspace=0.25,wspace=0.25)
		
		ax.plot(tog,linewidth = 4, color='#008080',label='Mutations')
		ax.set_ylabel('Mut Load / 25 Kb',fontsize=20)
		ax.locator_params(axis='y',tight=False, nbins=5)
		ax.locator_params(axis='x',tight=False, nbins=3)
		ticks= ax.get_xticks().tolist()
		ticks = ['-500Kb','Boundary','+500Kb']
		ax.set_xticklabels(ticks,fontsize=20)
		ax.set_xlim([0,40])
		
		ax.set_title(name)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		
		with np.errstate(divide='ignore'): img = ax1.imshow(norm_data,interpolation='nearest',aspect='auto',cmap='binary')
		ax1.set_xlim([0,40])
		ax1.xaxis.set_ticks([])
		ax1.yaxis.set_ticks([])
		ax1.set_frame_on(False)
		ax1.set_ylabel(name,fontsize=15)
		ax1.get_yaxis().set_label_coords(-0.075,0.5)
		plt.setp(ax1.get_xticklabels(), visible=False)
		ax1.axvspan(19, 21, facecolor='#5781A4', alpha=0.30, linestyle='dashed')
		
		plt.savefig(output+'.png',dpi=200)	
	
if __name__=='__main__':
	
	parser = argparse.ArgumentParser(usage='MutationAggregate.py -m Mutation.bedGraph -r Boundaries.bed -o output -n title',add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
	
	group = parser.add_argument_group("Required Parameters")
	group.add_argument('-m','--mutation',default='', help='',metavar='',required=True)
	group.add_argument('-o', '--output',default='',metavar='',required=True)
	group.add_argument('-n', '--name',default='',metavar='',required=True)
	
	group1 = parser.add_argument_group("Optional Parameters")
	group1.add_argument('-h', '--help', action="help")
	group1.add_argument('-r', '--region',default='',metavar='',help='')
	group1.add_argument('-d', '--domain',default='',metavar='',help='')
	group1.add_argument('-f', '--folder',default='',metavar='',help='')
	
	args = vars(parser.parse_args())
	
	if len(args['domain'])==0 and len(args['region'])==0 and len(args['folder'])==0:
		print >>sys.stderr, 'Upps!! Please activate one of the modules with -r, -d or -r & -h'
		raise SystemExit
	if len(args['domain'])>0 and len(args['folder'])>0:
		print >>sys.stderr, 'Upps!! Please either provide a domain file or a folder.'
		raise SystemExit

	
	Plotter(**args)