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
from scipy.signal import find_peaks
from collections import OrderedDict 
from scipy import signal
import scipy.sparse as sps
import pickle

def read_bedGraph(filename,resolution,chromosome): # add stopping after certain chromosome passed
	
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
	
	x_scores=[]
	x_scores2=[]
	y_scores=[] 
	colors=[]
	texts=[]
	
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if tags[0]==chromosome:
			x_scores.append(float(tags[1])/resolution)
			x_scores2.append(float(tags[2])/resolution)
			if len(tags) > 3:
				y_scores.append(float(tags[3]))
				if len(tags) > 4:
					hex = '#%02x%02x%02x' % (int(tags[4].split(',')[0]), int(tags[4].split(',')[1]), int(tags[4].split(',')[2]))
					colors.append(hex)
					if len(tags) > 5:
						texts.append(tags[5])
				
	if len(y_scores) !=0 and len(y_scores)!=len(x_scores):
		print >>sys.stderr, 'BedGraph('+filename+') has some missing values'
		raise SystemExit
	if len(x_scores)==0 or len(x_scores2)==0:
		print >>sys.stderr, 'BedGraph('+filename+') has some missing values'
		raise SystemExit
	# color and text controls
	return x_scores,x_scores2,y_scores,colors,texts
	

def read_genes(filename,resolution,chromosome,start,end):
		
	try:
		fone=open(filename,'r')
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	start = resolution * start
	end = resolution * end
	
	minDist = 8000
	
	genes = {}
	row_list = []
	row_genes = {}
	current_start=0;current_end=0;prev_end=0;
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if tags[0]==chromosome:
			if int(tags[1]) >= start and int(tags[2]) <= end and tags[1]+'-'+tags[2] not in genes.keys():
				
				if len(row_list)==0:
					current_start = int(tags[1])
					current_end = int(tags[2])
					prev_start = int(tags[1])
					genes[tags[1]+'-'+tags[2]]=[]
					genes[tags[1]+'-'+tags[2]].append(1)
					genes[tags[1]+'-'+tags[2]].append(tags[3])
					row_list.append(current_end)
					row_genes[1]=[]
					row_genes[1].append(current_end)
				else:
					if prev_start > int(tags[1]):
						print prev_end, int(tags[1])
						print >>sys.stderr, 'Gene File ('+filename+') is not sorted.'
						raise SystemExit
					else:
						current_end = int(tags[2])
						current_start = int(tags[1])
						execute=0
						genes[tags[1]+'-'+tags[2]]=[]
						for item in range(0,len(row_list)):
							if current_start > row_list[item]+minDist: 
								row_list[item]=current_end
								execute=1
								genes[tags[1]+'-'+tags[2]].append(item+1)
								if item+1 not in row_genes.keys(): row_genes[item+1]=[]
								row_genes[item+1].append(current_end)
								break
						if execute == 0:
							genes[tags[1]+'-'+tags[2]].append(len(row_list)+1)
							row_list.append(current_end)
							if len(row_list) not in row_genes.keys(): row_genes[len(row_list)]=[]
							row_genes[len(row_list)].append(current_end)
							
						genes[tags[1]+'-'+tags[2]].append(tags[3])
						if len(tags)>5:
							genes[tags[1]+'-'+tags[2]].append(tags[4])
							genes[tags[1]+'-'+tags[2]].append(tags[5])
							genes[tags[1]+'-'+tags[2]].append(tags[6])
							#if len(tags)>7:
							#	hex = '#%02x%02x%02x' % (int(tags[7].split(',')[0]), int(tags[7].split(',')[1]), int(tags[7].split(',')[2]))
							#	genes[tags[1]+'-'+tags[2]].append(hex)
							#	if len(tags)>8:
							#		genes[tags[1]+'-'+tags[2]].append(float(tags[8]))
								
					prev_start = current_start
							
	if len(genes.keys()) ==0:
		print >>sys.stderr, 'Gene File ('+filename+') has some missing values'
		raise SystemExit
	
	return genes,len(row_list)+1,row_genes


def compare(typeA='',typeB='',output=''):
	
	first = {}
	fone=open(typeA,'r')
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if line[0]!='t':
			if tags[0] not in first.keys():
				first[tags[0]]=[]
				first[tags[0]].append(float(tags[3]))
			else:
				first[tags[0]].append(float(tags[3]))
	second = {}
	fone=open(typeB,'r')
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if line[0]!='t':
			if tags[0] not in second.keys():
				second[tags[0]]=[]
				second[tags[0]].append(float(tags[3]))
			else:
				second[tags[0]].append(float(tags[3]))
	
	sigs = {}
	sigsO = {}
	called = {}
	resolution = 25000	
	firstCalls = []
		
	blacklist = {1:[500,4500,5000,5500],2:[3500],3:[3500],4:[],5:[],6:[2000],7:[2000,2500],8:[1500],9:[1500,2000,2500],10:[1500],11:[2000],12:[1000],13:[0,500],14:[0,500],15:[0,500],16:[1000,1500],17:[500],18:[],19:[1000],20:[1000],21:[0,500]}
	
	upR = []
	downR = []
	sel = []
	
	log2s = []
	rands = []
	
	for chrI in range(1,22):
		chromosome = 'chr'+str(chrI)
		start = 0
		end = 1000
		sigs[chrI]=[]
		sigsO[chrI]=[]
		called[chrI]=[]
		
		while end < len(first[chromosome]):
			length = end-start
			print chromosome, start,end
			if start in blacklist[chrI]: 
				start = end; 
				if end+500 > len(first[chromosome]):
					end = len(first[chromosome])
				else: end += 500
				continue
			
			y1 = np.array(first[chromosome][start:end])
			y2 = np.array(second[chromosome][start:end])
			n = len(y1)
			
			y1 = np.divide(y1*100,np.sum(y1))
			y2 = np.divide(y2*100,np.sum(y2))
			z = np.subtract(y1,y2)
			
			x = savgol_filter(z, 11, 3)
			
			decreasing  = x < 0
			increasing = x > 0

			changes = np.zeros_like(x)

			change_increasing = np.logical_and(increasing[1:] , decreasing[:-1])
			change_decreasing = np.logical_and(decreasing[1:] , increasing[:-1])

			changes[0] = (1 * increasing[0]) + (-1 * decreasing[0])
			changes[1:][change_increasing] = 1
			changes[1:][change_decreasing] = -1

			regions = []
			regions1 = np.where(x[:-1] * x[1:] < 0 )[0] +1


			zx = savgol_filter(z, 41, 3,deriv=1)
			peaks, _ = find_peaks(zx, height=0.015)
			npeaks, _ = find_peaks(np.array(zx)*-1, height=0.015)
			
			peaks = np.append(peaks,npeaks,axis=0)
			xpeaks = sorted(np.append(peaks,npeaks,axis=0))
			
			
			rand = []
			for xx in range(0,len(np.array(xpeaks))):
				occured = 1
				idx = (np.abs(np.array(regions1)-xpeaks[xx])).argmin()
				if abs(regions1[idx]-xpeaks[xx]) < 5: 
					
					if (abs(np.sum(x[xpeaks[xx]:xpeaks[xx]+5]))) not in upR : upR.append((abs(np.sum(x[xpeaks[xx]:xpeaks[xx]+5]))))
					if (abs(np.sum(x[xpeaks[xx]-5:xpeaks[xx]]))) not in downR : downR.append((abs(np.sum(x[xpeaks[xx]-5:xpeaks[xx]]))))
					if (abs(np.sum(x[xpeaks[xx]-5:xpeaks[xx]]))) not in rands : rands.append((abs(np.sum(x[xpeaks[xx]-115:xpeaks[xx]-110]))))
					
					if (abs(np.sum(x[xpeaks[xx]:xpeaks[xx]+5]))) > (abs(np.sum(x[xpeaks[xx]-5:xpeaks[xx]]))): sel.append((abs(np.sum(x[xpeaks[xx]:xpeaks[xx]+5]))))
					else: sel.append((abs(np.sum(x[xpeaks[xx]-5:xpeaks[xx]]))))
					
					#print abs(np.sum(x[xpeaks[xx]-5:xpeaks[xx]])), abs(np.sum(x[xpeaks[xx]:xpeaks[xx]+5]))
					
					if abs(np.sum(x[xpeaks[xx]-5:xpeaks[xx]])) + abs(np.sum(x[xpeaks[xx]:xpeaks[xx]+5])) > 0.5:
						regions.append(xpeaks[xx])
						sigs[chrI].append((xpeaks[xx]+start)*25000);rand.append(xpeaks[xx])
						firstCalls.append(chromosome+'\t'+str((xpeaks[xx]+start)*25000)+'\t'+str((xpeaks[xx]+start)*25000))

								
			start = end
			if end+500 > len(first[chromosome]):
				end = len(first[chromosome])
			else: end += 500
		sigsO[chrI] = list(OrderedDict.fromkeys(sorted(sigs[chrI])))
		
	print sigsO

### Plotting	
##	 			if not os.path.exists('Comparisons/'+flabel+'-vs-'+slabel+'/'): os.makedirs('Comparisons/'+flabel+'-vs-'+slabel+'/')
# 
# # 			if len(rand)>0:
# # 				fmatrix=read_sparseHiCdata('IndvMutChange/'+fhic+'.'+chromosome+'.txt',chromosome,'IndvMutChange/'+chromosome+'-bins.txt',start,end)
# # 				smatrix=read_sparseHiCdata('IndvMutChange/'+shic+'.'+chromosome+'.txt',chromosome,'IndvMutChange/'+chromosome+'-bins.txt',start,end)
# # 				
# # 				for xz in range(0,len(rand)):
# # 					matrix1 = np.nan_to_num(np.matrix(fmatrix[rand[xz]-6:rand[xz],rand[xz]-6:rand[xz]], dtype=np.float64))
# # 					matrix2 = np.nan_to_num(np.matrix(smatrix[rand[xz]-6:rand[xz],rand[xz]-6:rand[xz]], dtype=np.float64))
# # 					with np.errstate(divide='ignore',invalid='ignore'): upmatrix = np.absolute(np.nan_to_num(np.matrix(log2(matrix1/(matrix2*1.0)), dtype=np.float64)))
# # 					matrix1 = np.nan_to_num(np.matrix(fmatrix[rand[xz]:rand[xz]+6,rand[xz]:rand[xz]+6], dtype=np.float64))
# # 					matrix2 = np.nan_to_num(np.matrix(smatrix[rand[xz]:rand[xz]+6,rand[xz]:rand[xz]+6], dtype=np.float64))
# # 					with np.errstate(divide='ignore',invalid='ignore'): downmatrix = np.absolute(np.nan_to_num(np.matrix(log2(matrix1/(matrix2*1.0)), dtype=np.float64)))
# # 					
# # 					
# # 					matrix1 = np.nan_to_num(np.matrix(fmatrix[rand[xz]+120:rand[xz]+126,rand[xz]+120:rand[xz]+126], dtype=np.float64))
# # 					matrix2 = np.nan_to_num(np.matrix(smatrix[rand[xz]+120:rand[xz]+126,rand[xz]+120:rand[xz]+126], dtype=np.float64))
# # 					with np.errstate(divide='ignore',invalid='ignore'): random = np.absolute(np.nan_to_num(np.matrix(log2(matrix1/(matrix2*1.0)), dtype=np.float64)))
# # 					
# # 					log2s.append(np.nansum(upmatrix)); log2s.append(np.nansum(downmatrix)); rands.append(np.nansum(random))
# 			
# # 			if len(rand)>0:
# # 				fig, (ax_orig, ax_noise, ax_diff, ax_corr, ax_genes) = plt.subplots(5, 1)
# # 				fig.subplots_adjust(hspace=1,wspace=1)
# # 	
# # 				ax_orig.plot(y1,color='gray')
# # 				ax_orig.set_title(flabel)
# # 				ax_orig.set_ylim(0,1)
# # 				ax_orig.spines['right'].set_visible(False)
# # 				ax_orig.spines['top'].set_visible(False)
# # 				ax_orig.xaxis.set_ticks_position('bottom')
# # 				ax_orig.yaxis.set_ticks_position('left')
# # 				
# # 				for ritem in range(0,len(cdoms)):
# # 					if cdoms[ritem][0] >= start and cdoms[ritem][1] <= end:
# # 						#print cdoms[ritem][0]-start,cdoms[ritem][1]-cdoms[ritem][0],cdoms[ritem][2]
# # 						rect = Rectangle((cdoms[ritem][0]-start,0.9), (cdoms[ritem][1]-cdoms[ritem][0]), 0.05, color=cdoms[ritem][2])
# # 						ax_orig.add_patch(rect)
# # 				ax_orig.set_xlim(0,length)
# # 				
# # 				ax_noise.plot(y2,color='black')
# # 				ax_noise.set_title(slabel)
# # 				ax_noise.set_ylim(0,1)
# # 				ax_noise.set_xlim(0,length)
# # 				ax_noise.spines['right'].set_visible(False)
# # 				ax_noise.spines['top'].set_visible(False)
# # 				ax_noise.xaxis.set_ticks_position('bottom')
# # 				ax_noise.yaxis.set_ticks_position('left')
# # 	
# # 				ax_diff.plot(x,color='white')
# # 				ax_diff.set_title('Diff')
# # 				ax_diff.axhline(0,color='black')
# # 				with np.errstate(all='ignore'):ax_diff.fill_between(np.arange(0,length),z, where=z>0, color='gray', interpolate=True)
# # 				with np.errstate(all='ignore'):ax_diff.fill_between(np.arange(0,length),z, where=z<0, color='black', interpolate=True)
# # 				ax_diff.set_xlim(0,length)
# # 				ax_diff.spines['right'].set_visible(False)
# # 				ax_diff.spines['top'].set_visible(False)
# # 				ax_diff.xaxis.set_ticks_position('bottom')
# # 				ax_diff.yaxis.set_ticks_position('left')
# # 	
# # 				ax_corr.set_title('Derv')
# # 				ax_corr.plot(zx,color='black')
# # 				ax_corr.plot(rand, zx[rand], "x")
# # 				#ax_corr.plot(sigsO[20], zx[sigsO[20]], "x")
# # 				ax_corr.axhline(0,color='black')
# # 				ax_corr.set_xlim(0,length)
# # 				ax_corr.spines['right'].set_visible(False)
# # 				ax_corr.spines['top'].set_visible(False)
# # 				ax_corr.xaxis.set_ticks_position('bottom')
# # 				ax_corr.yaxis.set_ticks_position('left')
# # 
# # 				for xz in range(0,len(regions)):
# # 					ax_corr.axvspan(regions[xz], regions[xz]+2, facecolor='g', alpha=0.10, linestyle='dashed')
# # 
# # 				for xz in range(0,len(rand)):
# # 					ax_orig.axvspan(rand[xz], rand[xz]+2, facecolor='g', alpha=0.10, linestyle='dashed')
# # 					ax_noise.axvspan(rand[xz], rand[xz]+2, facecolor='g', alpha=0.10, linestyle='dashed')
# # 					ax_diff.axvspan(rand[xz], rand[xz]+2, facecolor='g', alpha=0.10, linestyle='dashed')
# # 
# # 				
# # 				if exp==0: ax_genes.set_ylabel('Genes')
# # 				ax_genes.get_yaxis().set_label_coords(-0.125,0.5)
# # 				genes,trackCount,nearest = read_genes('genes.sorted.bed',25000,chromosome,start,end)
# # 				plength = (end-start)*float(resolution)/1000000
# # 			
# # 				gcolor = '#3C3C8C';icolor='#0C0C78'
# # 						
# # 				for item in genes.keys():
# # 								
# # 					gstart = float(item.split('-')[0])/resolution-start
# # 					gend = float(item.split('-')[1])/resolution-start
# # 					gtrack = genes[item][0]
# # 					
# # 					if len(genes[item])>5: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=genes[item][5])
# # 					else: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=gcolor)
# # 					ax_genes.add_patch(rect)
# # 					
# # 				ax_genes.set_xlim(0,length)
# # 				ax_genes.set_ylim(0,trackCount+1)
# # 				ax_genes.spines['right'].set_visible(False)
# # 				ax_genes.spines['top'].set_visible(False)
# # 				ax_genes.xaxis.set_ticks_position('bottom')
# # 				ax_genes.yaxis.set_ticks_position('left')					
# # 				
# # 				ticks= ax_orig.get_xticks().tolist()
# # 				for item in range(0,len(ticks)): ticks[item]=round((ticks[item]+start)*25000/1000000,3)
# # 				ax_orig.set_xticklabels(ticks)
# # 				ax_noise.set_xticklabels(ticks)
# # 				ax_diff.set_xticklabels(ticks)
# # 				ax_corr.set_xticklabels(ticks)
# # 				ax_genes.set_xticklabels(ticks)
# # 	
# # 				ax_genes.set_xlabel('Chromosome %s (Mb)' % (chromosome.replace("chr","")))
# # 				plt.savefig('Comparisons/'+flabel+'-vs-'+slabel+'/'+chromosome+'-'+str(start)+'-'+str(end)+'.png',dpi=200)

if __name__=='__main__':
	
	parser = argparse.ArgumentParser(usage='MutLoadCompare.py -a CancerA.bedGraph -b CancerB.bedGraph -o output',add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
	
	group = parser.add_argument_group("Required Parameters")
	group.add_argument('-a','--typeA',default='', help='',metavar='',required=True)
	group.add_argument('-b','--typeB',default='', help='',metavar='',required=True)
	group.add_argument('-o', '--output',default='',metavar='',required=True)
	
	group1 = parser.add_argument_group("Optional Parameters")
	group1.add_argument('-h', '--help', action="help")
# 	group1.add_argument('-r', '--region',default='',metavar='',help='')
# 	group1.add_argument('-d', '--domain',default='',metavar='',help='')
# 	group1.add_argument('-f', '--folder',default='',metavar='',help='')
	
	args = vars(parser.parse_args())
	
	compare(**args)