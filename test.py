#!/usr/bin/python
from pylab import *
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt


def theta_data(paramsFile, filler=',', param="Theta"):
    name = np.array(np.genfromtxt(paramsFile, comments="!!!", dtype=None, delimiter=',', max_rows=1, names=True))
    name = name.dtype.names
    for i in name:
        if i == param:
            col = name.index(i)
    theta = np.loadtxt(paramsFile, delimiter=filler, skiprows=1, unpack=True, usecols=col)
    return theta


def monitor_data(paramsFile, filler=',', param="Monitor"):
    name = np.array(np.genfromtxt(paramsFile, comments="!!!", dtype=None, delimiter=',', max_rows=1, names=True))
    name = name.dtype.names
    for i in name:
        if i == param:
            col = name.index(i)
    monitor = np.loadtxt(paramsFile, delimiter=filler, skiprows=1, unpack=True, usecols=col)
    return monitor


def foil_data(paramsFile, filler=',', param="Foils"):
    name = np.array(np.genfromtxt(paramsFile, comments="!!!", dtype=None, delimiter=',', max_rows=1, names=True))
    name = name.dtype.names
    for i in name:
        if i == param:
            col = name.index(i)
    foil = np.loadtxt(paramsFile, delimiter=filler, skiprows=1, unpack=True, usecols=col)
    new_foil = []
    for i in foil:
        i = np.array([int(d) for d in str(int(i))])
        n = 0
        b = len(i)
        while n < 4 - b:
            i = np.insert(i, 0, 0)
            n += 1
        new_foil.append(i)
    return new_foil


def detector_data(paramsFile, bkgFile = None, filler=',', param="Detector", param_2="Foils", param_3 = "Detector", param_4="Monitor", sep_bkg = False, abs=[exp(1.444), exp(2.74), exp(5.316), exp(10.622)]):
    name = np.array(np.genfromtxt(paramsFile, comments="!!!", dtype=None, delimiter=filler, max_rows=1, names=True))
    name = name.dtype.names
    for i in name:
        if i == param:
            col = name.index(i)
    detector = np.loadtxt(paramsFile, delimiter=filler, skiprows=1, unpack=True, usecols=col)
    new_detector = detector
    i = 0
    while i < len(detector):
        abs_corrected = np.array(abs) * foil_data(paramsFile, param=param_2)[i]
        abs_corrected = abs_corrected[abs_corrected != 0]
        new_detector[i] = detector[i] * np.prod(abs_corrected) / monitor_data(paramsFile, param=param_4)[i]
        i += 1

    if sep_bkg == True:
        bkg_name = np.array(np.genfromtxt(bkgFile, comments="!!!", dtype=None, delimiter=filler, max_rows=1, names=True))
        bkg_name = bkg_name.dtype.names
        for i in bkg_name:
            if i == param_3:
                bkg_col = bkg_name.index(i)
        bkg_detector = np.loadtxt(bkgFile, delimiter=filler, skiprows=1, unpack=True, usecols=bkg_col)
        new_bkg_detector = bkg_detector
        i = 0
        while i < len(bkg_detector):
            abs_corrected = np.array(abs) * foil_data(bkgFile, param=param_2)[i]
            abs_corrected = abs_corrected[abs_corrected != 0]
            new_bkg_detector[i] = bkg_detector[i] * np.prod(abs_corrected)/monitor_data(bkgFile, param=param_4)[i]
            i += 1
        final_bkg_detector = np.concatenate((np.zeros(len(new_detector)-len(new_bkg_detector)), new_bkg_detector), axis=0)
        final_detector = new_detector - final_bkg_detector


    if sep_bkg == False:
        bkg_name = np.array(np.genfromtxt(paramsFile), comments="!!!", dtype=None, delimiter=filler, max_rows=1, names=True)
        bkg_name = bkg_name.dtype.names
        for i in bkg_name:
            if i == param_3:
                bkg_col = bkg_name.index(i)
        bkg_detector = np.loadtxt(paramsFile, delimiter=filler, skiprows=1, unpack=True, usecols=bkg_col)
        new_bkg_detector = bkg_detector
        i = 0
        while i < len(bkg_detector):
            abs_corrected = np.array(abs) * foil_data(bkgFile, param=param_2)[i]
            abs_corrected = abs_corrected[abs_corrected != 0]
            new_bkg_detector[i] = bkg_detector[i] * np.prod(abs_corrected) / monitor_data(bkgFile, param=param_4)[i]
            print (new_bkg_detector[i])
            i += 1
        final_bkg_detector = np.concatenate((np.zeros(len(new_detector)-len(new_bkg_detector)), new_bkg_detector), axis=0)
        final_detector = new_detector - final_bkg_detector
    return final_detector, new_detector, final_bkg_detector

def plotData(final_detector, new_detector, final_bkg_detector, theta, i0 = 1.76e7, param = "I0", lam = 1.078, beam = 0.06, size = 3.14, save = True, filename = 'outputFile'):
    final_detector = np.array(final_detector)
    new_detector = np.array(new_detector)
    final_bkg_detector = np.array(final_bkg_detector)
    theta = np.array(theta)
    qz = 4 * np.pi / lam * np.sin(theta / 180 * np.pi)
    i0_corrected = final_detector/i0
    fp_corrected = i0_corrected

    x = beam * 4 * np.pi / lam / qz
    fp_corrected[x >= float(size)] = (fp_corrected * (x / float(size)))[x >= float(size)]

    xrr = fp_corrected
    np.seterr(divide='ignore', invalid='ignore')

    # #==============================================================
    # #===> define some colors
    # #==============================================================
    RED = "#880000";
    BLU = "#0088cc";
    BLK = "#000000";
    GRY = "#555555"

    # #==============================================================
    # #===> define some parameters for figure
    # #==============================================================
    plt.rc("font", size=12)
    rcParams['figure.figsize'] = 8, 5
    # plt.rcParams['font.family']='M+ 2c'


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
    ax1.semilogy(qz, new_detector, marker="o", markerfacecolor='m', markeredgecolor='m', markersize=3,
                 markeredgewidth=0.8, label='int = det x exp(curratt * foil_thickness / attlength * thc)  / mon')
    ax1.semilogy(qz, final_bkg_detector, marker="o", markerfacecolor=RED, markeredgecolor=RED, markersize=3,
                 markeredgewidth=0.8, label='bg = det x exp(curratt * foil_thickness / attlength * thc)  / mon')
    ax1.semilogy(qz, final_detector, marker="o", markerfacecolor='g', markeredgecolor='g', markersize=3,
                 markeredgewidth=0.8, label='bg_corr = int -bg')
    ax1.semilogy(qz, xrr, marker="o", markerfacecolor='k', markeredgecolor='k', markersize=3,
                 markeredgewidth=0.8, label='footprint corrected = final xrr --> is saved')

    ax1.axvline(x=amax(qz[x >= float(size)]), color='k')

    # make a legend
    lg = ax1.legend(loc=1, numpoints=1, fontsize=8)
    lg.draw_frame(False)

    # #==============================================================
    # #===> Save it!
    # #==============================================================
    if save == True:
        output = np.column_stack((qz, xrr))
        np.savetxt("%s.xrr" % (filename), output)

        # gs.tight_layout(fig, w_pad=0.1, h_pad=0.1)
        savefig("%s.png" % (filename), bbox_inches='tight')
    show()