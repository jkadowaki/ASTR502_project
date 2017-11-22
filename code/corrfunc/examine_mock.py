import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table

def main():
    mock = Table().read('../../data/mock_fields/test_field1.fits')
    ra = np.array(mock['ra'])
    dec = np.array(mock['dec'])

    np.savetxt('data.dat', np.transpose([ra, dec]))






    return


if __name__ == '__main__':
    main()