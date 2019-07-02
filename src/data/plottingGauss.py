import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
# For latex stuff
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathpazo}',
        )
        #  \usepackage{foo-name} `...')
#matplotlib.verbose.level = 'debug-annoying'

plt.rcParams.update({'font.size': 14})

# 1: slope computations, 2: comparison with linearly evolved
# 3: classic
determine_program = 3

filename = sys.argv[1]
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

if determine_program == 1:
    plt.figure(1)
    maxBorn = np.max(BornApprox)
    indexMax = np.where(BornApprox == maxBorn)
    indexMax = indexMax[0][0]
    kMax = k[indexMax]

    lin_k = np.log(k[(indexMax +5):(indexMax + 20)])
    linBornApprox =  BornApprox[(indexMax +5):(indexMax + 20)]
    linBornApprox = np.log(linBornApprox)
    slope = stats.linregress(lin_k,linBornApprox) # it returns an array with various values, 0 is the slope
    lin_kCurly = np.log(k[100:118]) # with this you start with k = 10
    linCurlyP = np.log(curlyP[100:118])
    slopeCurly =  stats.linregress(lin_kCurly,linCurlyP)

    print('a = %f , k0 = %f , sigma = %f , Born peak = %f' %( a, k0 , sigma , kMax))
    print('slope Born = %f , slope Curly = %f' %(slope[0],slopeCurly[0]) )
    print(' %f & %f $ \pm $ %f & %f $ \pm $ %f & %f \\\\ \n ' %(sigma, slope[0],
        slope[4], slopeCurly[0], slopeCurly[4],kMax))
    plt.plot(k, BornApprox, k, curlyP)
    plt.title(r'Gaussian power spectrum, $ k_0 = %f \, , \ \sigma = %f $ $ [h \, \mathrm{Mpc}^{-1} ] $ $ a = %f$ ' %(k0,sigma, a) )
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend((r'Born approximated $ \bar{\mathcal{P}} $',
        r'Free non-linearly evolved $ \mathcal{P} $'), loc='lower left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')

if determine_program == 2:
    plt.figure(2)
    plt.plot(k, linPower, k, curlyP, k, BornApprox)
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    # to set the lower limit
    plt.gca().set_ylim([1e-20, 1e10])
    plt.gca().set_xlim(right=100)
    plt.legend((r'Linearly evolved $ P_\delta^\mathrm{lin} $',r'Free non-linearly evolved $  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )

    # Just check maximum of linPower (should be very near k0)
    maxlinear= np.max(linPower)
    indexMax = np.where(linPower == maxlinear)
    indexMax = indexMax[0][0]
    kMax = k[indexMax]
    print(kMax)

if determine_program == 3:
    plt.figure(2)
    plt.plot(k, curlyP, k, BornApprox)
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    # to set the lower limit
    #  plt.gca().set_ylim([1e-20, 1e10])
    #  plt.gca().set_xlim(right=1000)
    plt.legend((r'Free non-linearly evolved $  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='upper left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )

plt.show()


# To change font size
# plt.rcParams.update({'font.size': 14})
