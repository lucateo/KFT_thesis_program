import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
# For latex stuff
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathpazo}',
        )

get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]

# n=4
filename = "powerSpectra_a_1_n_initial_4.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(1)
plt.plot(k, linPower, k, curlyP, k, BornApprox)
plt.title(r'Dark matter power spectrum, $ n $ = 4 ')
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Linearly evolved $ P_\delta^\mathrm{lin} $',r'Free non-linearly evolved $ \exp(-Q_\mathrm{D}) \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.savefig("/home/luca/thesis_master/img_tikz/power_a_1_n_4_full.pdf")

# n=17
filename = "powerSpectra_a_1_n_initial_17.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(2)
plt.plot(k, curlyP, k, BornApprox)
plt.title(r'Dark matter power spectrum, $ n $ = 17 ')
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $\exp(-Q_\mathrm{D})  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.savefig("/home/luca/thesis_master/img_tikz/power_a_1_n_17.pdf")

# n=10 with SI
filename = "powerSpectra_a_1_n_initial_10kmax_10000.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
fig = plt.figure(3, figsize=(10,5))
ax1 = fig.add_subplot(121) #subplot(nrows, ncols, index)
ax2 = fig.add_subplot(122) #subplot(nrows, ncols, index)
ax1.plot(k, BornApprox, k, curlyP)
ax2.plot(k, meanField)#, k, linPower, k, curlyP, k, BornApprox)
ax1.set_xlabel('$ k \,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
ax2.set_xlabel('$ k \,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
ax1.set_ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax1.legend((r'Born approximated $ \bar{\mathcal{P}} $',r'Free non-linearly evolved $\exp(-Q_\mathrm{D}) \mathcal{P} $'), loc='lower left')
ax1.set_title(r'$\exp(-Q_\mathrm{D}) \mathcal{P} $ and $ \bar{\mathcal{P}} $')
ax2.set_title('$ S_\mathrm{I} $')
plt.suptitle(r'Dark matter power spectrum, $ n $ = 10 ')
plt.savefig("/home/luca/thesis_master/img_tikz/power_a_1_n_10_fullkmax_10000_SI.pdf")

# n=3.2
#  filename = "powerSpectrakmax10000_a_1_n_initial_3.2.txt"
#  data = np.loadtxt(filename, skiprows = 8)
#  k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
#  plt.figure(4)
#  plt.plot(k, linPower, k, curlyP, k, BornApprox)
#  plt.title(r'Dark matter power spectrum, $ n $ = 3.2 ')
#  plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
#  plt.xscale('log')
#  plt.yscale('log')
#  plt.legend((r'Linearly evolved $ P_\delta^\mathrm{lin} $',r'Free non-linearly evolved $ \exp(-Q_\mathrm{D}) \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
#  plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
#  plt.savefig("/home/luca/thesis_master/img_tikz/power_a_1_n_3.2fullkmax_10000.pdf")

# Gauss plots
filename = "GaussSpectra_a_0.029907_k0_0.01_sigma_4.64159.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(5)

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

plt.rcParams.update({'font.size': 14})
plt.plot(k, curlyP, k, BornApprox)
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $ \exp(-Q_\mathrm{D}) \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
plt.savefig("/home/luca/thesis_master/img_tikz/gauss_a_0.03_k0_0.01_sigma_4.64.pdf")



filename = "GaussSpectra_a_1_k0_0.01_sigma_4.64159.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(6)

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

plt.plot(k, curlyP, k, BornApprox)
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $ \exp(-Q_\mathrm{D}) \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='upper left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
plt.savefig("/home/luca/thesis_master/img_tikz/gauss_a_1_k0_0.01_sigma_4.64.pdf")


# n=6, a=5
filename = "powerSpectra_a_5_n_initial_6.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(7)
plt.plot(k, curlyP, k, BornApprox)
plt.title(r'Dark matter power spectrum at $ a=5 $, $ n $ = 6 ')
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $\exp(-Q_\mathrm{D})  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='upper left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.savefig("/home/luca/thesis_master/img_tikz/power_a_5_n_6.pdf")

# n=6, a=1
filename = "powerSpectra_a_1_n_initial_6.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(8)
plt.plot(k, curlyP, k, BornApprox)
plt.title(r'Dark matter power spectrum at $ a=1 $, $ n $ = 6 ')
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $ \exp(-Q_\mathrm{D}) \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.savefig("/home/luca/thesis_master/img_tikz/power_a_1_n_6.pdf")


# Gauss, a=0.01
filename = "GaussSpectrakmax10000_a_0.01_k0_0.01_sigma_0.1.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(10)

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

plt.plot(k, curlyP, k, BornApprox)
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $\exp(-Q_\mathrm{D})  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
plt.savefig("/home/luca/thesis_master/img_tikz/gausskmax10000_a_0.01_k0_0.1_sigma_0.1.pdf")

# gauss, a=0.1
filename = "GaussSpectrakmax10000_a_0.1_k0_0.01_sigma_0.1.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(11)

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

plt.plot(k, curlyP, k, BornApprox)
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $\exp(-Q_\mathrm{D})  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
plt.savefig("/home/luca/thesis_master/img_tikz/gausskmax10000_a_0.1_k0_0.01_sigma_0.1.pdf")


# Gauss a=1, sigma = 0.1
filename = "GaussSpectrakmax10000_a_1_k0_0.01_sigma_0.1.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(12)

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

plt.plot(k, curlyP, k, BornApprox)
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $\exp(-Q_\mathrm{D})  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
plt.savefig("/home/luca/thesis_master/img_tikz/gausskmax10000_a_1_k0_0.01_sigma_0.1.pdf")


# a = 0.01
filename = "GaussSpectrakmax10000_a_0.01_k0_1_sigma_0.5.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(13)

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

plt.plot(k, curlyP, k, BornApprox)
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $ \exp(-Q_\mathrm{D}) \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
plt.savefig("/home/luca/thesis_master/img_tikz/gausskmax10000_a_0.01_k0_1_sigma_0.5.pdf")

# a = 0.1, sigma = 0.5
filename = "GaussSpectrakmax10000_a_0.1_k0_1_sigma_0.5.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(14)

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

plt.plot(k, curlyP, k, BornApprox)
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $\exp(-Q_\mathrm{D})  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
plt.savefig("/home/luca/thesis_master/img_tikz/gausskmax10000_a_0.1_k0_1_sigma_0.5.pdf")

# a = 1, sigma = 0.5
filename = "GaussSpectrakmax10000_a_1_k0_1_sigma_0.5.txt"
data = np.loadtxt(filename, skiprows = 8)
k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
plt.figure(15)

array_index_ = get_indexes("_",filename)
t_index = get_indexes("t",filename)
a = filename[(array_index_[1]+1):array_index_[2] ]
k0 = filename[(array_index_[3]+1):array_index_[4] ]
sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
a = float(a)
k0 = float(k0)
sigma = float(sigma)

plt.plot(k, curlyP, k, BornApprox)
plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'Free non-linearly evolved $ \exp(-Q_\mathrm{D}) \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='upper left')
plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
plt.savefig("/home/luca/thesis_master/img_tikz/gausskmax10000_a_1_k0_1_sigma_0.5.pdf")

plt.show()
