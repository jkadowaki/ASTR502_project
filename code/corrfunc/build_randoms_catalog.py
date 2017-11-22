import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def main():
    nobj, good_obj = 1000000, 0
    radius = 0.5

    random_catalog = []
    while good_obj < nobj:
        # Draw random from [-0.5, 0.5)
        ra, dec = np.random.random(2)-0.5

        if ra**2 + dec**2 < radius**2:
            random_catalog.append([ra, dec])
            good_obj +=1

    random_catalog = np.asarray(random_catalog)


    plt.plot(random_catalog[:,0], random_catalog[:,1], '.k')
    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.5, 0.5)

    plt.savefig('randoms_catalog.png')
    np.savetxt('randoms.dat', random_catalog)


if __name__ == '__main__':
    main()