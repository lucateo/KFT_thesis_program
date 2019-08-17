import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
# For latex stuff
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathpazo}',
        )

# 1 = bispectrum, 2 = bispectrum fixed damping, 3 = trial curly P_ij 
determine_program = 3

filename = sys.argv[1]
data = np.loadtxt(filename, skiprows = 1)

if determine_program ==1:
    mu, first_factor, second_factor, QN_denominator, bispectrum = data.T  # trick: columns to variables

if determine_program ==2:
    mu, damping_corrected, P31_3, QN_denominator, bispectrum = data.T  # trick: columns to variables

if determine_program ==3:
    k, P1, P2, P3, P4 = data.T  # trick: columns to variables

get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]

t_index = get_indexes("t",filename)
array_index_ = get_indexes("_",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
n_initial = filename[(array_index_[4]+1):array_index_[5] ]
if determine_program == 1 or determine_program ==2:
    ratio= filename[(array_index_[6]+1):array_index_[7] ]
    k1 = filename[(array_index_[8]+1):(t_index[3]-1) ]
    ratio = float(ratio)
    k1 = float(k1)
if determine_program == 3:
    k2= filename[(array_index_[6]+1):array_index_[7] ]
    mu= filename[(array_index_[8]+1):(t_index[2]-1) ]
    k2 = float(k2)
    mu = float(mu)

a = float(a)
n_initial = float(n_initial)

if determine_program == 1:
    plt.figure(1)
    plt.plot(mu,bispectrum/QN_denominator)
    plt.title(r'Bispectrum, $ a = %.2f \, ,$ $ n = %.0f , \ k_1 = %.3f \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(a, n_initial, k1, ratio))
    plt.xlabel('$\mu $ ')
    plt.ylabel('Bispectrum $[h^{-3}\,\mathrm{Mpc}^3 ]$')

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
plt.show()
