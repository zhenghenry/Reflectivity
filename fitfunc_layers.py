import numpy as np
import scipy.special

import pyRefFit


def loadParams ( paramsFile ):

	params = np.loadtxt(paramsFile, unpack=True)

	paramsDict = {}

	for n,k in enumerate(('I0', 'bkg', 'rho0', 'rhoInf', 'sigInf')):
		paramsDict[k] =  pyRefFit.parameter(params[0, n], params[1, n], params[2, n], False)

	for l in range(5, len(params[0,:]), 3):
		for i,k in enumerate(('r', 'd', 's')):
			paramsDict['%s%d'%(k,((l-5)/3)+1)] = pyRefFit.parameter(params[0, l + i], params[1, l + i], params[2, l + i], False)

	print(paramsDict['rhoInf'].v, paramsDict['rhoInf'].l, paramsDict['rhoInf'].u, paramsDict['rhoInf'].f)

	return paramsDict



	




def saveParams ( paramsFile, params, header='# values \tlower \tupper', footer=''):

	output = np.column_stack(params)
	np.savetxt(paramsFile, output, header=header, footer=footer)



def reflectivity(q, params, isSLD=False):

	I0, bkg, rho0, rhoInf, sigInf, *params = params

	#if len(params) % 3 != 0: raise RuntimeError('Please make sure you have 5 + (n * 3) parameters in your set!'

	layers = [(rho0, 0, 0)]
	for n in range(0, len(params), 3):
		layers.append(tuple(params[n:n+3]))
	layers.append((rhoInf, 0, sigInf))

	sldConversion = 1
	if isSLD: sldConversion /= 28.179


	r = 0j
	depth =0

	for n in range(0, len(layers)-1):
		
		depth += layers[n][1] 
		r += ( layers[n][0] - layers[n+1][0] ) * sldConversion * np.exp( 1j*q*depth ) * np.exp( -q**2 * layers[n+1][2]**2 / 2 )

	r /= rhoInf - rho0

	r = I0 * abs(r)**2 * (0.0375*np.sqrt(rhoInf-rho0))**4 / 16/ q**4

	return r



def density(z, params, roughness=True):

	I0, bkg, rho0, rhoInf, sigInf, *params = params

	#if len(params) % 3 != 0: raise RuntimeError('Please make sure you have 5 + (n * 3) parameters in your set!'

	layers = [(rho0, 0, 0)]
	for n in range(0, len(params), 3):
		layers.append(tuple(params[n:n+3]))
	layers.append((rhoInf, 0, sigInf))

	density = layers[0][0]
	depth = 0


	for n in range(0, len(layers)-1):

		depth += layers[n][1]

		if roughness:
			density += ( layers[n+1][0] - layers[n][0] ) * ( 1 + scipy.special.erf((z-depth)/np.sqrt(2)/layers[n+1][2]) ) / 2 
		else:
			density += ( layers[n+1][0] - layers[n][0] ) * ( 1 + scipy.special.erf((z-depth)/np.sqrt(2)/0.0001) ) / 2 

	return density
