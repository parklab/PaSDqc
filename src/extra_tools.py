import pandas as pd
import numpy as np
import astropy.stats
import statsmodels.nonparametric
import matplotlib.pyplot as plt
import seaborn as sns
import re
import pathlib2

def chroms_from_build(build):
    chroms = {'grch37': [i for i in range(1, 23)] + ['X', 'Y'],
    }

    try:
        return chroms[build]
    except KeyError:
        raise ValueError("Oops, I don't recognize the build {}".format(build))

def symmetric_KL(psd_1, psd_2):
    """ Calculate the approximate symmetric KL divergence between to power spectra
    """
    KL = (psd_1 / psd_2 + psd_2 / psd_1 - 2).sum()

    return KL

def plot_KL_div_by_chrom(j):
    """ Plot KL divergence by chrom for grch37 sample
    """
    f = plt.figure()
    ax = f.add_subplot(111)

    mu = j.median()
    mad = np.median(np.abs(j - mu))

    chroms = j.index.tolist()
    chroms[-2:] = ['23', '24']
    chroms = [int(c) for c in chroms]
    chroms = pd.Series(chroms)

    # mask = (j <= (mu+mad))
    # print(mask)
    # print(len(chroms[(j <= (mu+mad)).tolist()]))
    # print(len(j[j <= (mu+mad)]))

    # return chroms

    ax.plot(chroms[(j <= (mu+mad)).tolist()], j[j <= (mu+mad)], 'o', label='Consistent')
    ax.plot(chroms[(j > (mu+mad)).tolist()], j[j > (mu+mad)], 'o', label='Abberant', color='red')
    ax.plot(sorted(chroms), [mu+mad for i in range(len(chroms))], '--', label='MAD cutoff')
    ax.legend(loc='upper left')
    ax.set_xlabel('chromosome')
    ax.set_ylabel('Symmetric KL Divergence')

    chroms = sorted(chroms)
    # chroms = chroms.tolist()
    ax.xaxis.set_ticks(chroms)
    chroms[-2:] = ['X', 'Y']
    ax.xaxis.set_ticklabels(chroms)

    return f

def PSD_to_ACF(freq, psd, lags):
    """ Convert PSD to ACF using the Wiener-Khinchin theorem.
        Approximates the integral using the right rectangle rule
    """
    freq_sym = np.append(-freq[::-1], freq) 
    psd_sym = np.append(psd[::-1], psd)

    steps = freq_sym[1:] - freq_sym[:-1]
    height = psd_sym[1:]

    nd = np.tile(freq_sym[1:], (len(lags), 1)).T

    acf = np.cos(-2*np.pi*nd*lags)*(height*steps)[:, np.newaxis]
    acf = acf.sum(axis=0)

    return acf

def PSD_sym_KL(u, v):
    """ Calculate the symmetric KL divergence between two spectral densities
        j = sum( u / v + v / u - 2)
    """
    j = (u / v + v / u - 2).sum()

def hclust(nd, method='ward'):
    """ Perform heirarchical clustering of PSDs using the symmetric KL divergence
    """
    dist = hc.distance.pdist(nd, PSD_sym_KL)
    link = hc.linkage(dist, method)

