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
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import scipy.signal as signal
from scipy.signal import argrelextrema
import seaborn as sns
import math
Polynomial = np.polynomial.Polynomial
from scipy.interpolate import UnivariateSpline



def SigCont(sigs='',output='',domain='',folder=''):

	doms = pybedtools.BedTool(domain).sort()
	
	files = glob.glob(folder+"/*.bedGraph")
	data = []
	
	for item in files:
		mutations = pybedtools.BedTool(item).sort()
		overlap = doms.map(mutations,c=4,o='sum')
	
	fsig = collections.OrderedDict()
	scounts = {}
	
	reader = open(sigs,'r')	
	for line in reader.readlines():
		tabs=line.strip().split('\t')
		if tabs[0] != 'Cancer Types' and tabs[1] in IDs.keys():
			scounts[IDs[tabs[1]][0]]={}
			counter = 0
			ratio = 0 
			for x in range(3,len(tabs)): counter+=float(tabs[x])
			if IDs[tabs[1]][0] in slopes.keys(): ratio = slopes[IDs[tabs[1]][0]]
			for x in range(3,len(tabs)): scounts[IDs[tabs[1]][0]][fsig.keys()[x-3]]=float(tabs[x])*100/counter;writer.write(str(float(tabs[x])*100/counter)+'\t')
		elif tabs[0] == 'Cancer Types':
			for x in range(3,len(tabs)):fsig[tabs[x]]=[];writer.write(tabs[x]+'\t')
		
	slope = []

	for item in scounts.keys():
		if item in slopes.keys(): slope.append(slopes[item])
		else: print item

			
	for item in scounts.keys():
		for x in fsig.keys():
			if item in slopes.keys():
				fsig[x].append(scounts[item][x])


	for item in range(0,49):
		if item in removed: continue
		print fsig.keys()[item],item
		sns.set_style("white")
		fig, (ax3) = plt.subplots(1, 1,figsize=(6,8))
		sig = fsig[fsig.keys()[item]]
		
		fslope = []
		ffsig = []
		for xx in range(0,len(sig)):
			fslope.append(slope[xx])
			ffsig.append(sig[xx])
				
		if len(ffsig)<3: continue
			

		cmin, cmax = min(fslope), max(fslope)
		pfit, stats = Polynomial.fit(fslope, np.log2(np.add(ffsig,1)), 1, full=True, window=(cmin, cmax), domain=(cmin, cmax))
		coefs = np.polyfit(fslope,np.log2(np.add(ffsig,1)),1)

		sns.despine(right=True)
		wwriter.write(fsig.keys()[item]+'\t'+str(coefs[0])+'\t'+str((np.log2(ffsig[fslope.index(max(fslope))])- np.log2(ffsig[fslope.index(min(fslope))]))/(max(fslope)-min(fslope)))+'\t'+str((np.log2(max(ffsig)) - np.log2(min(ffsig))) /	(fslope[ffsig.index(max(ffsig))]-fslope[ffsig.index(min(ffsig))]))+'\t'+str(np.median(fslope))+'\t'+str(np.average(fslope))+'\t'+str(100*scipy.stats.pearsonr(fslope, np.log2(np.add(ffsig,1)))[0])+'\t'+str(math.log10(scipy.stats.pearsonr(fslope, np.log2(np.add(ffsig,1)))[1]))+'\n')
		
		predict = np.poly1d(coefs)
		
		ax3.plot(fslope, pfit(np.array(fslope)), color='k')
		ax3.scatter(fslope, np.log2(np.add(ffsig,1)), c='#404040',s=80, alpha=0.6,label='Mutations')	
		ax3.scatter([fslope[ffsig.index(max(ffsig))], fslope[ffsig.index(min(ffsig))]], [np.log2(max(ffsig)), np.log2(min(ffsig))], c='red',s=80, alpha=0.6,label='Mutations')
		ax3.scatter([max(fslope), min(fslope)], [np.log2(ffsig[fslope.index(max(fslope))]), np.log2(ffsig[fslope.index(min(fslope))])], c='blue',s=80, alpha=0.6,label='Mutations')

		ax3.spines['right'].set_visible(False)
		ax3.spines['top'].set_visible(False)
		ax3.xaxis.set_ticks_position('bottom')
		ax3.yaxis.set_ticks_position('left')
		res = scipy.stats.theilslopes(slope, np.log2(np.add(sig,1)), 0.90)
		ax3.plot(slope, res[1] + res[0] * np.array(slope), color='black',linestyle=':',linewidth=1)
		at = AnchoredText('Pearson Co: '+str(round(scipy.stats.pearsonr(slope, sig)[0],3))+'\nP-value: '+str(scipy.stats.pearsonr(slope, sig)[1]),prop=dict(size=12), frameon=False,loc=1)
		ax3.add_artist(at)
		ax3.set_ylim(-0.1,7.0)
		ax3.set_ylabel('log2(% Mutation Load)',fontsize=18)
		ax3.set_xlabel('Active/Inactive Mut Ratio',fontsize=18)
		ax3.set_title(fsig.keys()[item],fontsize=30)
		plt.savefig(output+'-'+fsig.keys()[item]+'_vs_Slope.png',dpi=200)


if __name__=='__main__':
	
	parser = argparse.ArgumentParser(usage='MutationDistribution.py -f folderName -d domains.txt -s signatureCont.txt -o output',add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
	
	group = parser.add_argument_group("Required Parameters")
	group.add_argument('-o','--output',default='', help='',metavar='',required=True)
	group.add_argument('-s','--sigs',default='', help='',metavar='',required=True)
	group.add_argument('-d', '--domain',default='',metavar='',required=True)
	group.add_argument('-f', '--folder',default='',metavar='',required=True)
	
	group1 = parser.add_argument_group("Optional Parameters")
	group1.add_argument('-h', '--help', action="help")
	
	args = vars(parser.parse_args())
	
	SigCont(**args)