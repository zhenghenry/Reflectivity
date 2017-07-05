#!/usr/bin/python
from pylab import *
import pylab as P
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from numpy import genfromtxt

import sys

np.seterr(divide='ignore', invalid='ignore')

# #==============================================================
# #===> define some colors
# #==============================================================
RED = "#880000"; BLU = "#0088cc"; BLK = "#000000"; GRY = "#555555"


# #==============================================================
# #===> define some parameters for figure
# #==============================================================
plt.rc("font", size=12)
rcParams['figure.figsize'] = 8, 5
#plt.rcParams['font.family']='M+ 2c'


# #==============================================================
# #===> create the graphics
# #==============================================================
gs = gridspec.GridSpec(1, 1)
ax1 = plt.subplot(gs[0, 0])


# #==============================================================
# #===> some axis properties
# #==============================================================
ax1.set_xlabel("q$_\mathregular{z}$ (Å$^\mathregular{-1}$)")
ax1.set_ylabel("reflectivity")


#==============================================================
#===> Argparse to use optional command line parameters
#==============================================================
import argparse

parser = argparse.ArgumentParser(description='reduce reflectivity data.')

#parser.add_argument('paramFile',    metavar='PARAMFILE',  type=str, help='reflectivity parameter file from D8')
parser.add_argument('inputFiles',   metavar='INPUTFILE', type=str, nargs='+',   help='all Data files (refl and background)')
parser.add_argument('-b',           dest='beam',          type=float, default=0.011, help='beam height')
parser.add_argument('-s',           dest='size',          type=float, default=10,  help='sample size')
parser.add_argument('-i',           dest='I0',            type=float, default=1.76e7,   help='incident inentsity') # for liquid this is 1.3e7
parser.add_argument('-l',           dest='lam',        type=float, default=0.56356,   help='wavelength')
parser.add_argument('-o',           dest='outputFile',    type=str, default="output.xrr",   help='output File')
parser.add_argument('-q', dest='quiet', action='store_false', help='quiet output')

args = parser.parse_args()


# #==============================================================
# #===> Assign variables
# #==============================================================
ss = args.beam
size = args.size
i0 = args.I0
lam = args.lam




# #==============================================================
# #===> Set the title of output-graph
# #==============================================================
ax1.set_title(sys.argv[1] + " - " + sys.argv[2])


# #==============================================================
# #===> Absorbers
# #==============================================================

foil_thickness = 46;  # thickness of the one foil in micron
attlength = 43.6;     # 1/e attenuation length  in micron
thc=[1.0,1.01600,1.02343,1.05248,1.03896,1.05287,0.98626,1.00244,1.00671,1.0226,1.01913,1.01977,1.01052,1.00231,1.00142,1.00546,1.,1.,1.,1.] # thc array of the filters thickness corre



# #==============================================================
# #===> load the data!!!
# #==============================================================
# individual files for intensity from unspeced data
chi, mon, curratt, det = np.loadtxt(args.inputFiles[0], unpack=True, usecols = (0,4,9,-1))

def chi_data(paramsFile):
	chi = np.loadtxt(paramsFile[0], unpack = True, usecols = 0)
	return chi

def mon_data(paramsFile):
	mon = np.loadtxt(paramsFile[0], unpack = True, usecols = 4)
	return mon

def curratt_data(paramsFile):
	curratt = np.loadtxt(paramsFile[0], unpack = True, usecols = 9)
	return curratt

def det_data(paramsFile):
	det = np.loadtxt(paramsFile[0], unpack = True, usecols = -1)
	ii = det
	for i in range(len(det)):
		ii[i] = det[i] * exp(curratt_data(paramsFile)[i] * foil_thickness/attlength * thc[int(curratt_data(paramsFile)[i])])/mon_data(paramsFile)[i]
	return ii

# correct for "wrong" absorbers
ii = det
for i in range(len(det)):
	ii[i] = det[i] * exp(curratt[i]*foil_thickness/attlength*thc[int(curratt[i])])/mon[i]

# individual files for background from unspeced data
bg_mon, bg_curratt, bg_det = np.loadtxt(args.inputFiles[1], unpack=True, usecols = (4,9,-1))

# make array for bg
bg2 = bg_det

# correct for "wrong" absorbers
for i in range(len(bg_det)):
	bg2[i] = bg_det[i] * exp(bg_curratt[i]*foil_thickness/attlength*thc[int(bg_curratt[i])])/bg_mon[i]

# make array of length where background was not measured, and fill with zeros
bg1 = np.zeros(len(ii)-len(bg2))

# calculate qz
qz = 4 * 3.1415 / lam * np.sin (chi  /180 * 3.1415)

# concatenate array of zeros and background measured
bg = np.concatenate((bg1, bg2), axis=1)
bg_corr = ii-bg

# plot intensity, background, and intensity - background
ax1.semilogy(qz, ii,marker="o" ,linestyle='none', markerfacecolor='m', markeredgecolor='m', markersize=3, markeredgewidth=0.8, label='int = det x exp(curratt * foil_thickness / attlength * thc)  / mon')
ax1.semilogy(qz, bg,marker="o" ,linestyle='none', markerfacecolor=RED, markeredgecolor=RED, markersize=3, markeredgewidth=0.8, label='bg = det x exp(curratt * foil_thickness / attlength * thc)  / mon')
ax1.semilogy(qz, bg_corr,marker="o" ,linestyle='none', markerfacecolor='g', markeredgecolor='g', markersize=3, markeredgewidth=0.8, label='bg_corr = int -bg')


# #==============================================================
# #===> Correct I0
# #==============================================================
i0_corr = bg_corr/i0

# plot I0 corrected
ax1.semilogy(qz, i0_corr,marker="o" ,linestyle='none', markerfacecolor='r', markeredgecolor='r', markersize=3, markeredgewidth=0.8, label='I0 corrected')


# #==============================================================
# #===> Correct footprint
# #==============================================================
fp_corr = i0_corr

x=ss*4*3.1415/lam/qz
fp_corr[x>=float(size)]=(fp_corr*(x/float(size)))[x>=float(size)]

xrr = fp_corr

# #==============================================================
# #===> And plot!
# #==============================================================
ax1.semilogy(qz, xrr,marker="o" ,linestyle='none', markerfacecolor='k', markeredgecolor='k', markersize=3, markeredgewidth=0.8, label='footprint corrected = final xrr --> is saved')

ax1.axvline(x=amax(qz[x>=float(size)]), color='k')

# make a legend
lg = ax1.legend(loc=1, numpoints = 1, fontsize=8)
lg.draw_frame(False)

# #==============================================================
# #===> Save it!
# #==============================================================
output = np.column_stack((qz,xrr))
np.savetxt("%s.xrr"%(args.outputFile),output)

#gs.tight_layout(fig, w_pad=0.1, h_pad=0.1)
savefig("%s.png"%(args.outputFile), bbox_inches='tight')
if args.quiet: show()