import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats

# The name of the text file you want to plot
filename = sys.argv[1]

# dark matter case, this reads the first two letters
if filename[0:2] == "po":
    #Loading data
    data = np.loadtxt(filename, skiprows = 8)
    k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
    #Plotting
    plt.figure(1)
    plt.plot(k, meanField)#, k, linPower, k, curlyP, k, BornApprox)
    #  plt.legend(('$S_\mathrm{I}$', '$ D^2_+ P_\delta [h^{-3}\,\mathrm{Mpc}^3 ] $','$ \mathcal{P} [h^{-3}\,\mathrm{Mpc}^3 ] $','$ \bar{\mathcal{P}} [h^{-3}\,\mathrm{Mpc}^3 ] $'), loc='upper left')
   # plt.plot(k, BornApprox)
    plt.xlabel('$k \, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('$S_\mathrm{I} $')
    n = 4 # Change this to put the proper index in the title
    plt.title(f'Initial condition index $n = {n} $ ')

    plt.figure(2)
    plt.plot(k, BornApprox)
    #  plt.plot( k, linPower, k, curlyP, k, BornApprox)
    plt.xlabel('$k \, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Born approx case')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    #  plt.legend((r'$ D^2_+ P_\delta $','$ \mathcal{P} $','$ \overline{\mathcal{P}}  $'), loc='upper left')

    plt.figure(3)
    plt.plot(k, (BornApprox / linPower) -1)
    plt.xlabel('$k \, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Comparison linear vs Born')
    plt.ylabel('$ \overline{\mathcal{P}}/D^2_+ P_\delta  -1  $')

    # I want to determine the slope on large scales
    lin_k = np.log(k[50:-1])
    linBornApprox =  BornApprox[50:-1]
    # I remove the values that corresponds to too low values (giving error when
    # you do the log
    for x in linBornApprox:
        if x < 1e-7: # Change it to remove the probably garbage data at the end of the curve
            index = np.where(linBornApprox == x)
            index = index[0][0]
            linBornApprox = np.delete(linBornApprox, index)
            lin_k = np.delete(lin_k, index)
    linBornApprox = np.log(linBornApprox)
    slope = stats.linregress(lin_k,linBornApprox) # it returns an array with various values, 0 is the slope
    print(slope[0])
    plt.title("slope = %f " %slope[0])

    plt.figure(4)
    plt.plot(k, linPower)
    plt.title('Linear Power')
    plt.xlabel('$k/k_0$ ')
    plt.xscale('log')
    plt.yscale('log')

    plt.figure(5)
    plt.plot(k, curlyP)
    plt.title(r'$\mathcal{P}$ spectrum')
    plt.xlabel('$k/k_0$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

# Gaussian case
if filename[0:2] == "Ga":
    #Loading data
    data = np.loadtxt(filename, skiprows = 8)
    k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
    #Plotting
    plt.figure(1)
    plt.plot(k, meanField)
    plt.xlabel('$k \, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('$S_\mathrm{I} $')
    norm = 4
    plt.title(f'Normalization factor = {norm} ')

    plt.figure(2)
    plt.plot(k, BornApprox)
    plt.xlabel('$k/k_0 $ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Born approx case')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    plt.show()

    plt.figure(3)
    plt.plot(k, (BornApprox / linPower) -1)
    plt.xlabel('$k/k_0$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('$\overline{\mathcal{P}}/D^2_+ P_\delta  -1$ ')

    plt.figure(4)
    plt.plot(k, linPower)
    plt.title('Linear Power')
    plt.xlabel('$k/k_0$ ')
    plt.xscale('log')
    plt.yscale('log')

    plt.figure(5)
    plt.plot(k, curlyP)
    plt.title(r'$\mathcal{P}$ spectrum')
    plt.xlabel('$k/k_0$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()


if filename[0:4] == "Test":
    #Loading data
    data = np.loadtxt(filename)
    k, P = data.T  # trick: columns to variables

    # I want to determine the slope on large scales
    lin_k = np.log(k[50:-1])
    linP =  np.log(P[50:-1])
    # I remove the values that corresponds to too low values (giving error when
    # you do the log
    slope = stats.linregress(lin_k,linP) # it returns an array with various values, 0 is the slope
    print(slope[0])
    #Plotting
    plt.plot(k, P)
    plt.xlabel('$k$')
    plt.xscale('log')
    plt.yscale('log')
    plt.title("slope = %f " %slope[0])
    plt.show()



#curlyP plotting
if filename[0:7] == "curlyP.":
    #Loading data
    data = np.loadtxt(filename)
    q, B_1, B_2 = data.T  # trick: columns to variables

    #Plotting
    plt.plot(q, B_1, B_2)
    plt.legend(('$B_1$','$B_2$'), loc='upper right')
    plt.xlabel('$q$')
    plt.xscale('log')
    plt.title("plotting")
    plt.show()

#curlyP plotting
if filename[0:7] == "curlyP2":
    #Loading data
    data = np.loadtxt(filename)
    q, integral_mu = data.T  # trick: columns to variables

    #Plotting
    plt.plot(q, integral_mu)
    #plt.legend(('$B_1$','$B_2$'), loc='upper right')
    plt.xlabel('$q$')
    plt.xscale('log')
    plt.title("plotting")
    plt.show()

if filename[0:2] == "Q_":
    #Loading data
    data = np.loadtxt(filename)
    k, Q = data.T  # trick: columns to variables

    #Plotting
    plt.plot(k, Q)
    #plt.legend(('$B_1$','$B_2$'), loc='upper right')
    plt.xlabel('$k$')
    plt.xscale('log')
    plt.ylabel('$Q_\mathrm{D}$')
    plt.yscale('log')
    plt.title("plotting")
    plt.show()

if filename[0:2] == "Al":
    #Loading data
    data = np.loadtxt(filename, skiprows = 2)
    k, Q, S = data.T  # trick: columns to variables
    first_line = np.loadtxt(open(filename, 'rt').readlines(1))
    a, n_in, A, k_0 = first_line.T
    print(a,n_in,A,k_0)
    #Plotting
    plt.plot(k, Q, S)
    plt.legend(('$Q_\mathrm{D}$','$S_\mathrm{I}$'), loc='upper right')
    plt.xlabel('$k$')
    plt.xscale('log')
    #plt.ylabel('$Q_\mathrm{D}$')
    #plt.yscale('log')
    plt.title(f'a = {a}, n_in = {n_in}, A = {A}, $k_0$ = {k_0} ')
    plt.show()

