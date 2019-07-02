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
txt = glob.glob(filepath)
i=0
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

    #  print('a = %f , k0 = %f , sigma = %f , Born peak = %f' %( a, k0 , sigma , kMax))
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

#  plt.show()


# To change font size
# plt.rcParams.update({'font.size': 14})
