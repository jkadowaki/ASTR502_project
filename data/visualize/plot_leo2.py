import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

def main():
    t = Table.read('../leo2.fits')

    plt.plot(t['ra'], t['dec'], '.k', alpha=0.2, markeredgewidth=0)
    plt.xlabel('ra')
    plt.ylabel('dec')

    plt.savefig('leo2.png')








if __name__ == "__main__":
    main()
