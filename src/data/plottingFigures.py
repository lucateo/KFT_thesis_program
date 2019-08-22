import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
# For latex stuff
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathpazo}',
        )

get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]
# 0 = dark, 1 = gauss, 2= second dark, 3 = second gaussian, 4 = bispectrum gauss,
# 5 = bispectrum dark
determine_program = 4

if determine_program ==0:
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
    plt.legend((r'Linearly evolved $ P_\delta^\mathrm{lin} $',r'Free non-linearly evolved $ \mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
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
    plt.legend((r'Free non-linearly evolved $\mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
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
    ax1.legend((r'Born approximated $ \bar{\mathcal{P}} $',r'Free non-linearly evolved $\mathcal{P}_\mathrm{d} $'), loc='lower left')
    ax1.set_title(r'$\mathcal{P}_\mathrm{d} $ and $ \bar{\mathcal{P}} $')
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

if determine_program ==1:
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
    plt.legend((r'Free non-linearly evolved $ \mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
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
    plt.legend((r'Free non-linearly evolved $ \mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='upper left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
    plt.savefig("/home/luca/thesis_master/img_tikz/gauss_a_1_k0_0.01_sigma_4.64.pdf")

if determine_program ==2:
    # n=6, a=5
    plt.rcParams.update({'font.size': 14})
    filename = "powerSpectra_a_5_n_initial_6.txt"
    data = np.loadtxt(filename, skiprows = 8)
    k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
    plt.figure(7)
    plt.plot(k, curlyP, k, BornApprox)
    plt.title(r'Dark matter power spectrum at $ a=5 $, $ n $ = 6 ')
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend((r'Free non-linearly evolved $\mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='upper left')
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
    plt.legend((r'Free non-linearly evolved $ \mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.savefig("/home/luca/thesis_master/img_tikz/power_a_1_n_6.pdf")

if determine_program==3:
    # Gauss, a=0.01
    filename = "GaussSpectrakmax10000_a_0.01_k0_0.01_sigma_0.1.txt"
    plt.rcParams.update({'font.size': 14})
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
    plt.legend((r'Free non-linearly evolved $\mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
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
    plt.legend((r'Free non-linearly evolved $\mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
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
    plt.legend((r'Free non-linearly evolved $\mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
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
    plt.legend((r'Free non-linearly evolved $ \mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
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
    plt.legend((r'Free non-linearly evolved $\mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
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
    plt.legend((r'Free non-linearly evolved $ \mathcal{P}_\mathrm{d} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='upper left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.title(r'Gaussian power spectrum, $ a = %f  \\   k_0 = %f \, , \ \sigma = %f \  [h \, \mathrm{Mpc}^{-1} ] \ , $ ' %(a, k0,sigma,) )
    plt.savefig("/home/luca/thesis_master/img_tikz/gausskmax10000_a_1_k0_1_sigma_0.5.pdf")

if determine_program==4:
    filename1 = "powerSpectraHigherOrder_a_0.1_sigma_0.0464159_k0_0.01_ratio_2_k1_2.5.txt"
    filename2 = "powerSpectraHigherOrder_a_0.1_sigma_0.0464159_k0_0.1_ratio_2_k1_2.5.txt"
    filename3 = "powerSpectraHigherOrder_a_0.1_sigma_0.599484_k0_0.01_ratio_2_k1_2.5.txt"
    filename4 = "powerSpectraHigherOrder_a_0.1_sigma_0.599484_k0_0.1_ratio_2_k1_2.5.txt"
    filename5 = "powerSpectraHigherOrder_a_0.05_sigma_0.0464159_k0_0.01_ratio_2_k1_2.5.txt"
    filename6 = "powerSpectraHigherOrder_a_0.05_sigma_0.0464159_k0_0.1_ratio_2_k1_2.5.txt"
    filename7 = "powerSpectraHigherOrder_a_0.05_sigma_0.599484_k0_0.01_ratio_2_k1_2.5.txt"
    filename8 = "powerSpectraHigherOrder_a_0.05_sigma_0.599484_k0_0.1_ratio_2_k1_2.5.txt"
    
    t_index = get_indexes("t",filename1)
    array_index_ = get_indexes("_",filename1)

    data1 = np.loadtxt(filename1, skiprows = 1)
    mu1, linear_Eulerian1, bispectrum1, QN_denominator1, bispectrum1 = data1.T  
    data2 = np.loadtxt(filename2, skiprows = 1)
    mu2, linear_Eulerian2, bispectrum2, QN_denominator2, bispectrum2 = data2.T  
    data3 = np.loadtxt(filename3, skiprows = 1)
    mu3, linear_Eulerian3, bispectrum3, QN_denominator3, bispectrum3 = data3.T  
    data4 = np.loadtxt(filename4, skiprows = 1)
    mu4, linear_Eulerian4, bispectrum4, QN_denominator4, bispectrum4 = data4.T
    data5 = np.loadtxt(filename5, skiprows = 1)
    mu5, linear_Eulerian5, bispectrum5, QN_denominator5, bispectrum5 = data5.T  
    data6 = np.loadtxt(filename6, skiprows = 1)
    mu6, linear_Eulerian6, bispectrum6, QN_denominator6, bispectrum6 = data6.T  
    data7 = np.loadtxt(filename7, skiprows = 1)
    mu7, linear_Eulerian7, bispectrum7, QN_denominator7, bispectrum7 = data7.T  
    data8 = np.loadtxt(filename8, skiprows = 1)
    mu8, linear_Eulerian8, bispectrum8, QN_denominator8, bispectrum8 = data8.T

    a = filename1[(array_index_[1]+1):array_index_[2] ]
    ratio = filename1[(array_index_[7]+1):(array_index_[8]) ]
    k1 = filename1[(array_index_[9]+1):(t_index[2]-1) ]
    ratio = float(ratio)
    k1 = float(k1)
    a=float(a)
    theta1 = np.arccos(mu1)/np.pi
    theta2 = np.arccos(mu2)/np.pi
    theta3 = np.arccos(mu3)/np.pi
    theta4 = np.arccos(mu4)/np.pi
    theta5 = np.arccos(mu5)/np.pi
    theta6 = np.arccos(mu6)/np.pi
    theta7 = np.arccos(mu7)/np.pi
    theta8 = np.arccos(mu8)/np.pi

    plt.rcParams.update({'font.size': 14})
    
    plt.figure(1)
    plt.plot(theta1,bispectrum1/QN_denominator1, theta2, bispectrum2/QN_denominator2) 
    plt.title(r'$ Q_3 $, $ a = %.2f \, , r = %.3f \\ \sigma =0.046 , \ k_1 = %.3f  \ [h \, \mathrm{Mpc}^{-1} ] $  ' %(a, ratio, k1 ))
    plt.xlabel(r'$\theta/\pi $ ')
    plt.ylabel(r'$ Q_3 $')
    plt.legend((r'$k_0 = 0.01 \, h \, \mathrm{Mpc}^{-1}$',r'$k_0 = 0.1 \, h \, \mathrm{Mpc}^{-1}$'))
    plt.savefig("/home/luca/thesis_master/img_tikz/bispectrumGauss_a_0.1_sigma_0.046_k0_0.01.pdf")
    
    plt.figure(2)
    plt.plot(theta5,bispectrum5/QN_denominator5, theta6,bispectrum6/QN_denominator6) 
    plt.title(r'$ Q_3 $, $ a = 0.05 \, ,$ $  r = %.3f \\ \sigma =0.046 , \ k_1 = %.3f , \ [h \, \mathrm{Mpc}^{-1} ] $  ' %(ratio, k1))
    plt.xlabel(r'$\theta/\pi $ ')
    plt.ylabel(r'$ Q_3 $')
    plt.legend((r'$k_0 = 0.01 \, h \, \mathrm{Mpc}^{-1}$',r'$k_0 = 0.1 \, h \, \mathrm{Mpc}^{-1}$'))
    plt.savefig("/home/luca/thesis_master/img_tikz/bispectrumGauss_a_0.05_sigma_0.046_k0_0.1.pdf")
    
    plt.figure(3)
    plt.plot(theta3,bispectrum3/QN_denominator3, theta4,bispectrum4/QN_denominator4) 
    plt.title(r'$ Q_3 $, $ a = %.2f \, ,$ $  r = %.3f \\  \sigma =0.599 , \ k_1 = %.3f , \ [h \, \mathrm{Mpc}^{-1} ] $  ' %(a, ratio, k1))
    plt.xlabel(r'$\theta/\pi $ ')
    plt.ylabel(r'$ Q_3 $')
    plt.legend((r'$k_0 = 0.01 \, h \, \mathrm{Mpc}^{-1}$',r'$k_0 = 0.1 \, h \, \mathrm{Mpc}^{-1}$'))
    plt.savefig("/home/luca/thesis_master/img_tikz/bispectrumGauss_a_0.1_sigma_0.599_k0_0.01.pdf")
    
    plt.figure(4)
    plt.plot(theta7,bispectrum7/QN_denominator7, theta8,bispectrum8/QN_denominator8) 
    plt.title(r'$ Q_3 $, $ a = 0.05 \, ,$ $  r = %.3f \\ \sigma =0.599 , k_1 = %.3f , \ [h \, \mathrm{Mpc}^{-1} ] $  ' %(ratio, k1))
    plt.xlabel(r'$\theta/\pi $')
    plt.ylabel(r'$ Q_3 $')
    plt.legend((r'$k_0 = 0.01 \, h \, \mathrm{Mpc}^{-1}$',r'$k_0 = 0.1 \, h \, \mathrm{Mpc}^{-1}$'))
    plt.savefig("/home/luca/thesis_master/img_tikz/bispectrumGauss_a_0.05_sigma_0.599_k0_0.1.pdf")

if determine_program ==5:
    filename11 = "powerSpectraHigherOrder_a_0.1_n_initial_11_ratio_2_k1_2.5.txt"
    filename14 = "powerSpectraHigherOrder_a_0.1_n_initial_14_ratio_2_k1_2.5.txt"
    filename20 = "powerSpectraHigherOrder_a_0.1_n_initial_20_ratio_2_k1_2.5.txt"
    filename4 = "powerSpectraHigherOrder_a_0.1_n_initial_4_ratio_2_k1_2.5.txt"
    filename5 = "powerSpectraHigherOrder_a_0.1_n_initial_5_ratio_2_k1_2.5.txt"
    filename8 = "powerSpectraHigherOrder_a_0.1_n_initial_8_ratio_2_k1_2.5.txt"
    
    filename = "powerSpectraHigherOrder_a_1_n_initial_4_ratio_2_k1_2.5.txt"
    filename7 = "powerSpectraHigherOrder_a_1_n_initial_7_ratio_2_k1_2.5.txt"
    
    t_index = get_indexes("t",filename11)
    array_index_ = get_indexes("_",filename11)

    plt.rcParams.update({'font.size': 14})
    data11 = np.loadtxt(filename11, skiprows = 1)
    mu11, linear_Eulerian11, bispectrum11, QN_denominator11, bispectrum11 = data11.T  
    data14 = np.loadtxt(filename14, skiprows = 1)
    mu14, linear_Eulerian14, bispectrum14, QN_denominator14, bispectrum14 = data14.T  
    data4 = np.loadtxt(filename4, skiprows = 1)
    mu4, linear_Eulerian4, bispectrum4, QN_denominator4, bispectrum4 = data4.T  
    data5 = np.loadtxt(filename5, skiprows = 1)
    mu5, linear_Eulerian5, bispectrum5, QN_denominator5, bispectrum5 = data5.T  
    data8 = np.loadtxt(filename8, skiprows = 1)
    mu8, linear_Eulerian8, bispectrum8, QN_denominator8, bispectrum8 = data8.T  
    data20 = np.loadtxt(filename20, skiprows = 1)
    mu20, linear_Eulerian20, bispectrum20, QN_denominator20, bispectrum20 = data20.T

    data = np.loadtxt(filename, skiprows = 1)
    mu, linear_Eulerian, bispectrum, QN_denominator, bispectrum = data.T
    data7 = np.loadtxt(filename7, skiprows = 1)
    mu7, linear_Eulerian7, bispectrum7, QN_denominator7, bispectrum7 = data.T


    a = filename11[(array_index_[1]+1):array_index_[2] ]
    ratio= filename11[(array_index_[6]+1):array_index_[7] ]
    k1= filename11[(array_index_[8]+1):(t_index[3]-1) ]
    ratio = float(ratio)
    k1 = float(k1)
    a=float(a)
    theta11 = np.arccos(mu11)/np.pi
    theta20 = np.arccos(mu20)/np.pi
    theta14 = np.arccos(mu14)/np.pi
    theta4 = np.arccos(mu4)/np.pi
    theta5 = np.arccos(mu5)/np.pi
    theta8 = np.arccos(mu8)/np.pi
    
    theta = np.arccos(mu)/np.pi
    theta7 = np.arccos(mu7)/np.pi
    
    plt.figure(1)
    plt.plot(theta4,bispectrum4/QN_denominator4, theta5, bispectrum5/QN_denominator5) 
    plt.title(r'$ Q_3 $, $ a = %.2f \, ,$ $ k_1 = %.3f , \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(a, k1, ratio))
    plt.xlabel(r'$\theta/\pi $ ')
    plt.ylabel(r'$ Q_3 $')
    plt.yscale('log')
    plt.legend((r'$ n $ = 4', r'$ n=5 $'))
    plt.savefig("/home/luca/thesis_master/img_tikz/bispectrum_a_0.1_n_4_5_ratio_2_k1_2.5.pdf")

    plt.figure(2)
    plt.plot(theta8,bispectrum8/QN_denominator8, theta11, bispectrum11/QN_denominator11,
            theta20, bispectrum20/QN_denominator20) 
    plt.title(r'$ Q_3 $ , $ a = %.2f \, ,$ $ k_1 = %.3f , \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(a, k1, ratio))
    plt.xlabel(r'$\theta/\pi $ ')
    plt.ylabel(r'$ Q_3 $')
    plt.yscale('log')
    plt.legend((r'$ n $ = 8', r'$ n=11 $', r'$ n = 20$'))
    plt.savefig("/home/luca/thesis_master/img_tikz/bispectrum_a_0.1_n_8_11_14_20_ratio_2_k1_2.5.pdf")

    plt.figure(3)
    plt.plot(theta,bispectrum/QN_denominator) 
    plt.title(r'$ Q_3 $ , $ a = 1.0 \, ,$ $ k_1 = %.3f , \ [h \, \mathrm{Mpc}^{-1} ], \ r = %.3f $' %(k1, ratio))
    plt.xlabel(r'$\theta/\pi $ ')
    plt.ylabel(r'$ Q_3 $')
    #  plt.yscale('log')
    plt.savefig("/home/luca/thesis_master/img_tikz/bispectrum_a_1_n_4_ratio_2_k1_2.5.pdf")

plt.show()
