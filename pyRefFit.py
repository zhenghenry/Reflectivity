from pylab import *
from scipy import optimize
from scipy.special import *
import cmath


class parameter:

	def __init__(self, value, lower, upper, fixed):

		self.v = value
		self.l = lower
		self.u = upper
		self.f = fixed

	def __repr__(self):

		return '<value: %0.3f bounds: %0.3f / %0.3f>' % (self.v, self.l, self.u)






fitData     = None
fitFunction = None 
fitResult   = None



def loadRefl( reflFile ):

	return np.loadtxt(reflFile, unpack=True, usecols=(0,1))




def loadParams ( paramsFile ):

	return loadtxt(paramsFile, unpack=True)




def doFit( data, function, params, qmin=0.1, qmax=0.7, verbose=False ):

	"""
	Perform fit using logarithmic figure of merit using custom fit function:

	
	Parameters
	----------
	data: 2-D array (2,n)
	    Reflectivity data to be fitted; use numpy.array as generated by loadtxt
	function: callable
	    Function describing the reflectivity curve. Must be in the form of f(x, *args), where
	    'x' is the argument in the form of a 1-D array and 'args' is a tupel of parameters
	     used to evaluate the function.
	params: 2-D array (3,n)
	    Array that contains the initial guess for the parameters of 'function' (first column),
	    the lower and the upper boundaries (second and third column).
	qmin: float, optional
	    Minimum q-Value for evaluating the fit against the data.
	qmax: float, optional
	    Maximum q-Value for evaluationg the fit against the data.

	"""
	
	global fitData, fitFunction, fitResult

	fitData     = data
	fitFunction = function

	values, lower, upper = params
	
	# restrict q-range to qmin <= q <= qmax
	dataX, dataY = data
	dataX, dataY = dataX[dataX >= qmin], dataY[dataX >= qmin]
	dataX, dataY = dataX[dataX <= qmax], dataY[dataX <= qmax]
	data = dataX, dataY


		# equalize dimensions of parameter space to intervals of [0; 1]
	startParams = np.nan_to_num((values - lower) / (values - upper))
	bounds = np.column_stack((np.zeros(lower.shape), np.ones(upper.shape)))

	# perform fit
	minimizer_kwargs = dict(method="L-BFGS-B", bounds=np.column_stack((lower, upper)), args=(function, data))
	result = optimize.basinhopping(FOM, values, niter=100, minimizer_kwargs=minimizer_kwargs, niter_success=10)

	fitResult   = result

	if verbose: print(result)

	return (result.x, lower, upper), result.fun





def FOM( params, *args ):

	refFunction, data = args
	dataX, dataY = data

	chi2 = (  (dataY - refFunction(dataX, params))  /  dataY  )**2
	#chi2 = (np.log(dataY) - np.log(refFunction(dataX, params)))**2
	#chi2 = (dataY - refFunction(dataX, params))**2

	return np.sum(chi2)
