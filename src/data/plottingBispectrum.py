# Program for bispectrum plots

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
# For latex stuff
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathpazo}',
        )

# 1 = bispectrum, 2 = bispectrum fixed damping, 3 = trial curly P_ij
# 4 = bispectrum gaussian
determine_program = 1

filename = sys.argv[1]
#  filename2 = sys.argv[2]
data = np.loadtxt(filename, skiprows = 1)
#  data2 = np.loadtxt(filename2, skiprows = 1)

if determine_program ==1 or determine_program ==4:
    mu, linear_Eulerian, bispectrum, QN_denominator, bispectrum = data.T  # trick: columns to variables
    linear_Eulerian = linear_Eulerian*2
    #  mu2, linear_Eulerian2,bispectrum2,QN_denominator2, bispectrum2 = data2.T

if determine_program ==2:
    mu, damping_corrected, P31_3, QN_denominator, bispectrum = data.T  # trick: columns to variables

if determine_program ==3:
    k, P1, P2, P3, P4 = data.T  # trick: columns to variables

get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]

t_index = get_indexes("t",filename)
array_index_ = get_indexes("_",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
if determine_program == 1 or determine_program ==2:
    ratio= filename[(array_index_[6]+1):array_index_[7] ]
    n_initial = filename[(array_index_[4]+1):array_index_[5] ]
    k1 = filename[(array_index_[8]+1):(t_index[3]-1) ]
    ratio = float(ratio)
    k1 = float(k1)
    n_initial = float(n_initial)
if determine_program == 3:
    k2= filename[(array_index_[6]+1):array_index_[7] ]
    mu= filename[(array_index_[8]+1):(t_index[2]-1) ]
    n_initial = filename[(array_index_[4]+1):array_index_[5] ]
    k2 = float(k2)
    mu = float(mu)
    n_initial = float(n_initial)
if determine_program ==4:
    sigma= filename[(array_index_[3]+1):array_index_[4] ]
    k0= filename[(array_index_[5]+1):array_index_[6] ]
    ratio = filename[(array_index_[7]+1):(array_index_[8]) ]
    k1 = filename[(array_index_[9]+1):(t_index[2]-1) ]
    ratio = float(ratio)
    k1 = float(k1)
    sigma = float(sigma)
    k0 = float(k0)

a = float(a)

if determine_program == 1:
    plt.figure(1)
    plt.plot(mu,bispectrum/QN_denominator) #  mu2, bispectrum2/QN_denominator2)
    plt.title(r'Bispectrum, $ a = %.2f \, ,$ $ n = %.0f , \ k_1 = %.3f \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(a, n_initial, k1, ratio))
    plt.xlabel('$\mu $ ')
    plt.ylabel('Bispectrum $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.yscale('log')
    
    plt.figure(2)
    theta = np.arccos(mu)/np.pi
    plt.plot(theta,bispectrum/QN_denominator, theta, linear_Eulerian)
    plt.title(r'Bispectrum, $ a = %.2f \, ,$ $ n = %.0f , \ k_1 = %.3f \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(a, n_initial, k1, ratio))
    plt.xlabel(r'$\theta / \pi $ ')
    plt.ylabel('Bispectrum $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.legend(('KFT','Eulerian'))

if determine_program == 2:
    plt.figure(1)
    plt.plot(mu,damping_corrected/QN_denominator)
    plt.title(r'$ Q_\mathrm{D} $, $ a = %.2f \, ,$ $ n = %.0f , \ k_1 = %.3f \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(a, n_initial, k1, ratio))
    plt.xlabel('$\mu $ ')
    plt.ylabel('Bispectrum $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    
    plt.figure(2)
    plt.plot(mu,P31_3)
    plt.title(r'$ \mathcal{P} $ example, $ a = %.2f \, ,$ $ n = %.0f , \ k_1 = %.3f \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(a, n_initial, k1, ratio))
    plt.xlabel('$\mu $ ')
    plt.ylabel(r'$ \mathcal{P} \ [h^{-3}\,\mathrm{Mpc}^3 ]$')

if determine_program == 3:
    plt.figure(1)
    plt.plot(k,P1, k,P2,k,P3,k,P4)
    plt.xlabel('$ k $ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(('P1','P2','P3','P4'))

if determine_program == 4:
    plt.figure(1)
    plt.plot(mu,bispectrum/QN_denominator) #  mu2, bispectrum2/QN_denominator2)
    plt.title(r'Bispectrum, $ a = %.2f \, ,$ $ \sigma = %.3f , \ k_0 = %.3f , \ k_1 = %.3f , \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(a, sigma, k0, k1, ratio))
    plt.xlabel('$\mu $ ')
    plt.ylabel('Bispectrum $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.yscale('log')
plt.show()
