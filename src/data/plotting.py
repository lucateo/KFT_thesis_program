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

# To change font size
# plt.rcParams.update({'font.size': 14})
# The name of the text file you want to plot
filename = sys.argv[1]

# dark matter case, this reads the first two letters
if filename[0:2] == "po" or filename[0:4]=="luca":
    #Loading data
    data = np.loadtxt(filename, skiprows = 8)
    k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
    #Plotting
    fig = plt.figure(1, figsize=(10,5))
    ax1 = fig.add_subplot(121) #subplot(nrows, ncols, index)
    ax2 = fig.add_subplot(122) #subplot(nrows, ncols, index)
    ax1.plot(k, BornApprox, k, curlyP)
    ax2.plot(k, meanField)#, k, linPower, k, curlyP, k, BornApprox)
    ax1.set_xlabel('$ k \,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax1.legend((r'Born approximated $ \bar{\mathcal{P}} $',r'Free non-linearly evolved $ \mathcal{P} $'), loc='lower left')
    ax1.set_title(r'$ \mathcal{P} $ and $ \bar{\mathcal{P}} $')
    ax2.set_title('$ S_\mathrm{I} $')
    #  plt.title(f'Mean interacting action, initial condition index $n = {n} $ ')

    plt.figure(2)
    plt.plot(k, BornApprox)
    plt.xlabel('$k \,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    # I want to determine the slope on large scales
    lin_k = np.log(k[110:122]) # with this you start with k = 10
    linBornApprox =  BornApprox[110:122]
    # I remove the values that corresponds to too low values (giving error when
    # you do the log
    #  for x in linBornApprox:
    #      if x < 1e-7: # Change it to remove the probably garbage data at the end of the curve
    #          index = np.where(linBornApprox == x)
    #          index = index[0][0]
    #          linBornApprox = np.delete(linBornApprox, index)
    #          lin_k = np.delete(lin_k, index)
    linBornApprox = np.log(linBornApprox)
    slope = stats.linregress(lin_k,linBornApprox) # it returns an array with various values, 0 is the slope
    linCurlyP = np.log(curlyP[110:122])
    slopeCurly =  stats.linregress(lin_k,linCurlyP)

    plt.figure(3)
    plt.plot(k, (BornApprox / linPower) -1)
    plt.xlabel('$k \,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Comparison linear vs Born')
    plt.ylabel(r'$ \bar{\mathcal{P}}/P_\delta^\mathrm{lin}  -1  $')

    plt.figure(4)
    plt.plot(k, BornApprox, k, curlyP)
    maxBorn = np.max(BornApprox)
    indexMax = np.where(BornApprox == maxBorn)
    indexMax = indexMax[0][0]
    kMax = k[indexMax]

    get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]
    array_index_ = get_indexes("_",filename)
    point_index = get_indexes(".",filename)
    n = filename[(array_index_[4]+1):point_index[-1] ]
    if len(n) == 1:
        n=int(n)
        determine = False
    elif n[1]=='.':
        n = float(n)
        determine = True # determines if n is float or int
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
        determine = False
    print(n)

    print('Born peak = %f' %kMax)
    plt.title(r'Dark matter power spectrum, $ n $ = %s' %n)
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend((r'Born approximated $ \bar{\mathcal{P}} $',r'Free non-linearly evolved $ \mathcal{P} $'), loc='lower left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    print("Slope Born = %f , error = %f"%(slope[0],slope[4] ))
    print("Slope curlyP = %f, error = %f" %(slopeCurly[0],slopeCurly[4] ) )

    plt.figure(5)
    plt.plot(k, linPower, k, curlyP, k, BornApprox)
    plt.title(r'Dark matter power spectrum, $ n $ = %s' %n)
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend((r'Linearly evolved $ P_\delta^\mathrm{lin} $',r'Free non-linearly evolved $  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    if determine == True:
        print("%.1f & %f $ \pm $ %f & %f $ \pm $ %f & %f \\\\ " %(n,slope[0],slope[4], slopeCurly[0],slopeCurly[4], kMax))
    else:
        print("%i & %f $ \pm $ %f & %f $ \pm $ %f & %f \\\\ " %(n,slope[0],slope[4], slopeCurly[0],slopeCurly[4], kMax))
    plt.show()

# Gaussian case
if filename[0:2] == "Ga":
    #Loading data
    data = np.loadtxt(filename, skiprows = 8)
    k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
    #Plotting
    fig = plt.figure(1, figsize=(10,5))
    ax1 = fig.add_subplot(121) #subplot(nrows, ncols, index)
    ax2 = fig.add_subplot(122) #subplot(nrows, ncols, index)
    ax1.plot(k, BornApprox, k, curlyP)
    ax2.plot(k, meanField)#, k, linPower, k, curlyP, k, BornApprox)
    ax1.set_xlabel('$ k \,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax1.legend((r'Born approximated $ \bar{\mathcal{P}} $',r'Free non-linearly evolved $ \mathcal{P} $'), loc='lower left')
    ax1.set_title(r'$ \mathcal{P} $ and $ \bar{\mathcal{P}} $')
    ax2.set_title('$ S_\mathrm{I} $')

    plt.figure(2)
    plt.plot(k, BornApprox)
    plt.xlabel('$k \,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    # I want to determine the slope on large scales
    lin_k = np.log(k[100:118]) # with this you start with k = 10
    linBornApprox =  BornApprox[100:118]
    linBornApprox = np.log(linBornApprox)
    slope = stats.linregress(lin_k,linBornApprox) # it returns an array with various values, 0 is the slope
    linCurlyP = np.log(curlyP[100:118])
    slopeCurly =  stats.linregress(lin_k,linCurlyP)
    plt.title('Born')

    plt.figure(3)
    plt.plot(k, (BornApprox / linPower) -1)
    plt.xlabel('$k \,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Comparison linear vs Born')
    plt.ylabel(r'$ \bar{\mathcal{P}}/P_\delta^\mathrm{lin}  -1  $')

    plt.figure(4)
    plt.plot(k, BornApprox, k, curlyP)
    maxBorn = np.max(BornApprox)
    indexMax = np.where(BornApprox == maxBorn)
    indexMax = indexMax[0][0]
    kMax = k[indexMax]

    #  index_ = np.where(filename == 'a')
    get_indexes = lambda filename, x: [i for (y, i) in zip(x, range(len(x))) if filename == y]
    array_index_ = get_indexes("_",filename)
    t_index = get_indexes("t",filename)
    a = filename[(array_index_[1]+1):array_index_[2] ]
    k0 = filename[(array_index_[3]+1):array_index_[4] ]
    sigma =filename[(array_index_[5]+1):(t_index[1]-1) ]
    a = float(a)
    k0 = float(k0)
    sigma = float(sigma)
    print('a = %f' %a)
    print('k0 = %f' %k0)
    print('sigma = %f'%sigma)

    print('Born peak = %f' %kMax)
    print(' %f & %f $ \pm $ %f & %f $ \pm $ %f & %f \\\\ ' %(sigma, slope[0], slope[4], slopeCurly[0], slopeCurly[4],kMax))
    plt.title(r'Gaussian power spectrum, $ k_0 = %f \, , \ \sigma = %f $ $ [h \, \mathrm{Mpc}^{-1} ] $ ' %(k0,sigma) )
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend((r'Born approximated $ \bar{\mathcal{P}} $',r'Free non-linearly evolved $ \mathcal{P} $'), loc='lower left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')
    print("Slope Born = %f , error = %f"%(slope[0],slope[4] ))
    print("Slope curlyP = %f, error = %f" %(slopeCurly[0],slopeCurly[4] ) )

    plt.figure(5)
    plt.plot(k, linPower, k, curlyP, k, BornApprox)
    plt.title(r'Gaussian power spectrum')
    plt.xlabel('$k\,\, [h \, \mathrm{Mpc}^{-1} ]$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend((r'Linearly evolved $ P_\delta^\mathrm{lin} $',r'Free non-linearly evolved $  \mathcal{P} $', r'Born approximated $ \bar{\mathcal{P}} $'), loc='lower left')
    plt.ylabel('Power spectra $[h^{-3}\,\mathrm{Mpc}^3 ]$')

    plt.show()


if filename[0:4] == "Test" or filename[0:4]=="test":
    #Loading data
    data = np.loadtxt(filename)
    k, P = data.T  # trick: columns to variables

    # I want to determine the slope on large scales
    lin_k = np.log(k[100:-1])
    linP =  np.log(P[100:-1])
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

# Try with mainTest of libKFT
if filename[0:4] == "Bart":
    #Loading data
    data = np.loadtxt(filename, skiprows = 8)
    k, meanField, linPower, curlyP, BornApprox = data.T  # trick: columns to variables
    #Plotting
    plt.figure(1)
    plt.plot(k, linPower, k, curlyP, k, BornApprox)
    plt.legend((r'$ P_\delta^\mathrm{lin} $',r'$ \mathcal{P} $', r'$ \bar{\mathcal{P}} $'), loc='lower left')
    plt.xlabel(r'$k$ ')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'$\bar{\mathcal{P}} $')
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

