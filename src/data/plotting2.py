# Plotting program for plotting every file that corresponds to the shortcut
# on the variable filepath, it search for the directory where the program is

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
import glob
# For latex stuff
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathpazo}',
        )
        #  \usepackage{foo-name} `...')
#matplotlib.verbose.level = 'debug-annoying'

filepath = "GaussSpectrakmax10000_a_*_k0_1_sigma_*.txt"
filepath2 = "powerSpectra_a_1_n_initial_*.txt"
txt = glob.glob(filepath)
txt2 = glob.glob(filepath2)
print(txt2)
i=0
# 0: gaussian, 1: dark spectrum, 2: peak graph
determine_program = 2

if determine_program ==0:
    for textfile in txt:
        plt.figure(i)
        data = np.loadtxt(textfile, skiprows = 8)
        k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables

        maxBorn = np.max(BornApprox)
        indexMax = np.where(BornApprox == maxBorn)
        indexMax = indexMax[0][0]
        kMax = k[indexMax]

        get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]
        array_index_ = get_indexes("_",textfile)
        t_index = get_indexes("t",textfile)
        a = textfile[(array_index_[1]+1):array_index_[2] ]
        k0 = textfile[(array_index_[3]+1):array_index_[4] ]
        sigma =textfile[(array_index_[5]+1):(t_index[1]-1) ]
        a = float(a)
        k0 = float(k0)
        sigma = float(sigma)

        lin_k = np.log(k[(indexMax +5):(indexMax + 20)])
        linBornApprox =  BornApprox[(indexMax +5):(indexMax + 20)]
        linBornApprox = np.log(linBornApprox)
        slope = stats.linregress(lin_k,linBornApprox) # it returns an array with various values, 0 is the slope
        lin_kCurly = np.log(k[100:118]) # with this you start with k = 10
        linCurlyP = np.log(curlyP[100:118])
        slopeCurly =  stats.linregress(lin_kCurly,linCurlyP)

        print('%.2f & %.4f & %f $ \pm $ %f & %f $ \pm $ %f & %f \\\\  ' %(a, sigma, slope[0],
            slope[4], slopeCurly[0], slopeCurly[4],kMax))
        plt.plot(k, BornApprox, k, curlyP)
        plt.title(r'Gaussian power spectrum, $ k_0 = %f \, , \ \sigma = %f $ $ [h \, \mathrm{Mpc}^{-1} ] $ $ a = %f$ ' %(k0,sigma, a) )
        plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend((r'Born approximated $ \bar{\mathcal{P}} $',
            r'Free non-linearly evolved $ \mathcal{P} $'), loc='lower left')
        plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
        i=i+1
if determine_program ==1:
    for textfile in txt2:
        plt.figure(i)
        data = np.loadtxt(textfile, skiprows = 8)
        k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables

        maxBorn = np.max(BornApprox)
        indexMax = np.where(BornApprox == maxBorn)
        indexMax = indexMax[0][0]
        kMax = k[indexMax]
        # to make the slope computations not go out of bounds if the max is on
        # a too big k
        if indexMax > 106:
            indexMax = 106

        get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]
        array_index_ = get_indexes("_",textfile)
        t_index = get_indexes("t",textfile)
        a = textfile[(array_index_[1]+1):array_index_[2] ]
        n = textfile[(array_index_[4]+1):(t_index[2] -1)]
        a = float(a)
        n = float(n)

        lin_k = np.log(k[(indexMax +5):(indexMax + 20)])
        linBornApprox =  BornApprox[(indexMax +5):(indexMax + 20)]
        linBornApprox = np.log(linBornApprox)
        slope = stats.linregress(lin_k,linBornApprox) # it returns an array with various values, 0 is the slope
        lin_kCurly = np.log(k[100:118]) # with this you start with k = 10
        linCurlyP = np.log(curlyP[100:118])
        slopeCurly =  stats.linregress(lin_kCurly,linCurlyP)

        print('%.2f & %f $ \pm $ %f & %f $ \pm $ %f & %f \\\\  ' %(a, slope[0],
            slope[4], slopeCurly[0], slopeCurly[4],kMax))
        plt.plot(k, BornApprox, k, curlyP)
        plt.title(r'Gaussian power spectrum, $ a = %f$, $ n = %.1f $ ' %(a, n) )
        plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend((r'Born approximated $ \bar{\mathcal{P}} $',
            r'Free non-linearly evolved $ \mathcal{P} $'), loc='lower left')
        plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
        i=i+1

if determine_program == 2:
    kplot = []
    n_plot = []
    get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]
    for textfile in txt2:
        data = np.loadtxt(textfile, skiprows = 8)
        k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables

        maxBorn = np.max(BornApprox)
        indexMax = np.where(BornApprox == maxBorn)
        indexMax = indexMax[0][0]
        kMax = k[indexMax]

        array_index_ = get_indexes("_",textfile)
        t_index = get_indexes("t",textfile)
        n = textfile[(array_index_[4]+1):(t_index[2] -1)]
        if len(n) == 1:
            n=int(n)
            determine = True
        elif n[1]=='.':
            n = float(n)
            determine = False
        elif n[1] =='k':
            n= n[0:1]
            n = int(n)
            determine = False
        elif len(n) >2 and n[2] =='k':
            n= n[0:2]
            n = int(n)
            determine = False
        else:
            n=int(n)
            determine = True
        if determine == True:
            n_plot.append(n)
            kplot.append(kMax)
    plt.xlabel('$ n $')
    plt.ylabel('$ k_\mathrm{max} $')
    print(kplot)
    print(n_plot)
    plt.yscale('log')
    plt.scatter(n_plot,kplot)

plt.show()


