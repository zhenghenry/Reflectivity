import numpy as np
import scipy.special
from astropy.convolution import Gaussian1DKernel, convolve


# missing loadParams function
def saveParams(paramsFile, params, header='# values \tlower \tupper', footer=''):
    output = np.column_stack(params)
    np.savetxt(paramsFile, output, header=header, footer=footer)


def reflectivity(q, params, isSLD=False):
    # make sure we have a numpy array
    params = np.array(params)

    # introduce constrained parameters
    params[params < 0] = params[-params[params < 0].astype(int)]

    # extract parameters
    I0, bkg = params[0:2]
    rho = np.concatenate((np.ones(1) * params[2], params[4::3]))
    sig = params[3::3]
    d = np.concatenate((np.zeros(1), params[5::3]))

    numInterfaces = sig.shape[0]

    # create mesh of qs; one row per layer (including ambient/substrate), one column per q-value
    qs = np.sqrt(np.subtract(*np.meshgrid((q ** 2), (8 * 2 * np.pi * 2.82e-5 * rho))))

    # create empty arrays of the right shape for r, p; one row per interface, one coulmn per q-value
    r = np.zeros((numInterfaces, q.shape[0]), dtype=np.complex64)
    p = np.copy(r)

    # calculate reflective indexes for each interface, phase terms for each layer
    for ii in np.arange(0, numInterfaces):
        r[ii] = (qs[ii] - qs[ii + 1]) / (qs[ii] + qs[ii + 1]) * np.exp(-sig[ii] ** 2 * qs[ii] ** 2 / 2)
        p[ii] = np.exp(1j * qs[ii] * d[ii])

        # recursively build the reflective index of the entire system from the bottom up
        rr = r[numInterfaces - 1]
    for ii in np.arange(0, numInterfaces - 1)[::-1]:
        rr = (r[ii] + rr * p[ii + 1]) / (1 + r[ii] * rr * p[ii + 1])

    xrr = I0 * np.abs(rr) ** 2

    np.place(xrr, np.isnan(rr), 1)

    sigG = 0.0058 / 2.355
    g = Gaussian1DKernel(sigG / (q[1] - q[0]))
    xrr_conv = convolve(xrr, g, boundary='extend')

    # xrr = I0 * np.abs(xr)**2


    return xrr_conv


def density(z, params, roughness=True):
    # make sure we have a numpy array
    params = np.array(params)

    # introduce constrained parameters
    params[params < 0] = params[-params[params < 0].astype(int)]
    # params[params<0] = -params[params<0]

    # extract parameters
    I0, bkg = params[0:2]

    rho = np.concatenate((np.ones(1) * params[2], params[4::3]))

    sig = params[3::3]
    if not roughness:
        sig = 0 * sig + 0.0001

    d = np.concatenate((np.zeros(1), params[5::3]))
    Z = np.cumsum(d)

    numInterfaces = sig.shape[0]

    # sum up density contributions each layer
    density = z * 0 + rho[0]
    for ii in np.arange(0, numInterfaces):
        density += (rho[ii + 1] - rho[ii]) * (1 + scipy.special.erf((z - Z[ii]) / np.sqrt(2) / sig[ii])) / 2

    return density
