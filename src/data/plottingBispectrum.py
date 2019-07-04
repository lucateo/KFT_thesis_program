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

#  plt.rcParams.update({'font.size': 14})

# 1: slope computations, 2: comparison with linearly evolved
# 3: classic
determine_program = 1

filename = sys.argv[1]
data = np.loadtxt(filename, skiprows = 7)
k, meanField, curlyP, FixedBornApprox = data.T  # trick: columns to variables
get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]

t_index = get_indexes("t",filename)
array_index_ = get_indexes("_",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
n_initial = filename[(array_index_[4]+1):array_index_[5] ]
kprime = filename[(array_index_[6]+1):array_index_[7] ]
mu = filename[(array_index_[8]+1):(t_index[2]-1) ]

a = float(a)
n_initial = float(n_initial)
kprime = float(kprime)
mu = float(mu)

if determine_program == 1:
    plt.figure(1)
    plt.plot(k, np.exp(meanField+FixedBornApprox)* curlyP , k, curlyP)
    plt.title(r'Bispectrum, $ a = %.2f \, ,$ $ n = %.0f , \ \mu = %.3f , \ k_2 = %.3f $ $ [h \, \mathrm{Mpc}^{-1} ] $ ' %(a, n_initial, mu, kprime))
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend((r'Born approximated bispectrum',
        r'Free non-linearly evolved bispectrum'), loc='lower left')
    plt.ylabel('Bispectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.show()
