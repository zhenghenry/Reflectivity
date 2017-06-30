#!/usr/bin/python

import warnings
warnings.simplefilter("ignore", RuntimeWarning) 


from pylab import *
import sys
import matplotlib.gridspec as gridspec

import pyRefFit
import fitfunc_layers as fitfunc



data   = pyRefFit.loadRefl(sys.argv[1])
params = pyRefFit.loadParams(sys.argv[2])

#test = fitfunc.loadParams(sys.argv[2])
#print(test)



fitResult, fitFOM = pyRefFit.doFit(data, fitfunc.reflectivity, params, qmin=00.11, qmax=0.48, verbose=True)

fitfunc.saveParams(sys.argv[3], fitResult, footer='%.5f'% fitFOM)




############################## Load data to be plotted
dataX, dataY = data

############################## Range for density profile
r1 = np.arange(-30,300,0.01)

############################## unpack parameters
values, lower, upper = params




############################## Make subfigures
rcParams['figure.figsize'] = 8, 6
plt.rc("font", size=12)
plt.rcParams['font.family']='M+ 2c'
fig=figure()

suptitle(sys.argv[3])

gs = gridspec.GridSpec(2, 2, width_ratios=[6,2])
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[1, 0])
ax3 = plt.subplot(gs[:, 1])

############################## Make plot for reflectivity
ax1.xaxis.set_ticks(np.arange(0, 1.2, 0.2))
ax1.set_ylim(1e-8,2)
ax1.set_xlabel("q$_\mathregular{z}$ (\u00c5\u207b\u00b9)")
ax1.set_ylabel("R/R$_\mathregular{F}$")
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position("top")

ax1.semilogy(dataX, dataY  /      abs( ( dataX - sqrt(dataX**2 - (0.0375*sqrt(values[3]-values[2]))**2 ) ) / ( dataX + sqrt(dataX**2 - (0.0375*sqrt(values[3]-values[2]))**2 ) ) ) **2    , marker="o", linestyle="none", markerfacecolor='None', markeredgecolor='k', markersize=5.5, markeredgewidth=2)
#ax1.semilogy(dataX, dataY/((0.0375*sqrt(values[12]-values[1]))**4/16/dataX**4), marker="o", linestyle="none", markerfacecolor='None', markeredgecolor='k', markersize=5.5, markeredgewidth=2)

#p1 = arange(0.001,0.8,0.001)
#ax1.plot(p1, fitfunc(p1, values)/((0.0375*sqrt(values[12]-values[1]))**4/16/p1**4), color='b', linewidth=2)

t1 = arange(0.001,1.2,0.001)
ax1.semilogy(t1, fitfunc.reflectivity(t1, fitResult[0])/((0.0375*sqrt(values[3]-values[2]))**4/16/t1**4), color='r', linewidth=2)

z1 = arange(0.001,1.2,0.001)
ax1.semilogy(z1, fitfunc.reflectivity(z1, params[0])/((0.0375*sqrt(values[3]-values[2]))**4/16/t1**4), color='b', linewidth=1)

output = np.column_stack((t1, fitfunc.reflectivity(t1, fitResult[0])))
np.savetxt("%s.r"%(sys.argv[3]),output)

#output = np.column_stack((z1, fitfunc.reflectivity(t1, params[0])))
#np.savetxt("%s.r.initial"%(sys.argv[3]),output)

############################## Make plot for density
ax2.set_xlabel("z (\u00c5)")
ax2.set_ylabel("\u03C1 (e/\u00c5\u00b3)")
ax2.xaxis.set_ticks(np.arange(-100, 300, 10))
ax2.yaxis.set_ticks(np.arange(0, 0.8, 0.2))
ax2.set_ylim(0,1.3)

ax2.plot(r1,fitfunc.density(r1, fitResult[0]), color='k', linewidth=2)
ax2.plot(r1,fitfunc.density(r1, fitResult[0], roughness=False), color='k', linewidth=2, ls="dashed")

output = np.column_stack((r1,fitfunc.density(r1, fitResult[0])))
np.savetxt("%s.ed"%(sys.argv[3]),output)

output = np.column_stack((r1,fitfunc.density(r1, fitResult[0], roughness=False)))
np.savetxt("%s.ned"%(sys.argv[3]),output)

############################## Make plot for values

tmp0=np.around(fitResult[0],decimals=3)
tmp1=list(tmp0)
tmp1=map(str,tmp1)
results="\n".join(tmp1)

ax3.set_ylim(0,10)
ax3.set_xlim(0,10)
ax3.annotate(sys.argv[1], xy=(1,9.5))
ax3.annotate(results, xy=(1,1), fontsize=9)
ax3.annotate("chi2 = %0.2f" %(fitFOM), xy=(1,0.3))
ax3.set_xticks([]) 
ax3.set_yticks([]) 

gs.tight_layout(fig, w_pad=0.1, h_pad=0.1)
savefig("%s.pdf"%(sys.argv[3]), bbox_inches='tight')

############################## Show plots for rho(z) and reflectivity
show()


