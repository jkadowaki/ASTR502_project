import numpy as np
import matplotlib.pyplot as plt


def main():
    sep, xi, error = np.loadtxt('autocorr.dat', skiprows=1, usecols=[0,3,4], unpack=True)

    plt.errorbar(sep, xi, yerr=error, fmt='.')
    plt.axhline(0.0, color='k')
    plt.xscale('log')

    plt.ylabel(r'$\xi(r)$')
    plt.xlabel(r'$r\;\;\mathrm{arcmin}$')

    plt.savefig('correlationfunction.png')

    return


if __name__ == '__main__':
    main()