# PSDTools.py - classes for MDAqc Power Spectral Density calculations
#
# v 0.0.7
# rev 2017-03-19 (MS: summarize chrom KL divergences)
# Notes:

import pandas as pd
import numpy as np
import astropy.stats
import statsmodels.nonparametric
import matplotlib.pyplot as plt
import seaborn as sns
import re
import pathlib2
import scipy.cluster.hierarchy as hc

from . import PSDTools

def chroms_from_build(build):
    """ Get list of chromosomes from a particular genome build

        Args:
            build           str

        Returns:
            chrom_list      list
    """
    chroms = {'grch37': [i for i in range(1, 23)] + ['X', 'Y'],
    }

    try:
        return chroms[build]
    except KeyError:
        raise ValueError("Oops, I don't recognize the build {}".format(build))

def summarize_KL_div_by_chrom(j_list, sample_list):
    chrom_pass = []
    chrom_warn = []
    chrom_fail = []

    for j, sample in zip(j_list, sample_list):
        mu = j.median()
        mad = np.median(np.abs(j - mu))
        sd_1 = mu + mad
        sd_2 = mu + 2*mad

        chrom_pass.append(j[j <= sd_1].index.tolist())
        chrom_warn.append(j[(j > sd_1) & (j <= sd_2)].index.tolist())
        chrom_fail.append(j[j > sd_2].index.tolist())

    cols = ['Chrom: pass', 'Chrom: warn', 'Chrom: fail']

    df = pd.DataFrame.from_items(zip(cols, [chrom_pass, chrom_warn, chrom_fail]))
    df.index = sample_list

    return df

def plot_KL_div_by_chrom(j):
    """ Plot KL divergence by chrom for grch37 sample

        Args:
            j       pandas series       KL divergence of each chromosome

        Returns:
            f       matplotlib figure
    """
    f = plt.figure()
    ax = f.add_subplot(111)

    mu = j.median()
    mad = np.median(np.abs(j - mu))
    sd_1 = mu + mad
    sd_2 = mu + 2 * mad

    chroms = j.index.tolist()
    chroms[-2:] = ['23', '24']
    chroms = [int(c) for c in chroms]
    chroms = pd.Series(chroms)

    ax.plot(chroms[(j <= (sd_1)).tolist()], j[j <= (sd_1)], 'o', label='Pass')
    ax.plot(chroms[((j > sd_1) & (j <= sd_2)).tolist()], j[(j > sd_1) & (j <= sd_2)], 'o', label='Warn', color='orange')
    ax.plot(chroms[(j > (sd_2)).tolist()], j[j > (sd_2)], 'o', label='Fail', color='red')
    ax.plot(sorted(chroms), [sd_1 for i in range(len(chroms))], '--', label='one std')
    ax.plot(sorted(chroms), [sd_2 for i in range(len(chroms))], '--', label='two std')
    ax.legend(loc='upper left')
    ax.set_xlabel('chromosome')
    ax.set_ylabel('Symmetric KL Divergence')

    chroms = sorted(chroms)
    ax.xaxis.set_ticks(chroms)
    chroms[-2:] = ['X', 'Y']
    ax.xaxis.set_ticklabels(chroms)

    return f

def PSD_to_ACF(freq, psd, lags):
    """ Convert PSD to ACF using the Wiener-Khinchin theorem.
        Approximates the integral using the right rectangle rule

        Args:
            freq        1d array        array of frequencies
            psd         1d array        PSD values at specified frequencies
            lags        1d array        array of lags at which to calcualte ACF

        Returns:
            acf         1d array        ACF at specified lags
    """
    freq_sym = np.append(-freq[::-1], freq) 
    psd_sym = np.append(psd[::-1], psd)

    steps = freq_sym[1:] - freq_sym[:-1]
    height = psd_sym[1:]

    nd = np.tile(freq_sym[1:], (len(lags), 1)).T

    acf = np.cos(-2*np.pi*nd*lags)*(height*steps)[:, np.newaxis]
    acf = acf.sum(axis=0)

    return acf

def mk_ndarray(dir_in):
    """ Load PSD dataframes from a directory, calculate each sample's average
        and form into an ND array for clustering.

        Args:
            dir_in      str         input directory

        Returns:
            freq        1d array    frequencies of PSDs
            nd          nd array    array of PSDs
            sample_list list        names of samples
    """
    p = pathlib2.Path(dir_in)
    file_list = sorted(p.glob("*chroms.spec"))
    sample_list = [f.name.split('.chroms.spec')[0] for f in file_list]

    df_list = [pd.read_table(str(f), index_col=0) for f in file_list]
    psd_list = [PSDTools.SamplePSD(df, s) for df, s in zip(df_list, sample_list)]
    avg_list = [psd.avg_PSD() for psd in psd_list]

    nd = np.array(avg_list)
    freq = np.array(df_list[0].index.tolist())

    return freq, nd, sample_list

def PSD_sym_KL(u, v):
    """ Calculate the symmetric KL divergence between two spectral densities
        j = sum( u / v + v / u - 2)

        Args:
            u       1d array        PSD array
            v       1d array        PSD array
    """
    n = len(v)
    j = (1 / n) * (u / v + v / u - 2).sum()

    return j

def hclust(nd, method='ward'):
    """ Perform heirarchical clustering of PSDs using the symmetric KL divergence

        Args:
            nd          nd array    array of PSDs

        Kwargs
            method      str         clustering method

        Returns:
            link        nd array    array of linkage relationships
    """
    dist = hc.distance.pdist(nd, PSD_sym_KL)
    link = hc.linkage(dist, method)

    return link

def mk_categorical_spectra(freq, nd, labels):
    """ Construct categorical spectra from an nd array of average spectra and a list of labels
        labels should be of the form ['good', 'good', 'bad', ...]

        Args:
            nd      nd array    array of PSDs
            labels  list        list of labels, one for each PSD

        Returns:
            --      dataframe   dataframe of categorical spectra with one column for each label
    """
    labels_uniq = set(labels)

    d = dict.fromkeys(labels_uniq)

    for label in labels_uniq:
        mask = np.array([label == l for l in labels])
        d[label] = np.median(nd[mask, :], axis=0)

    return pd.DataFrame(d, index=freq)

def classify_samples(nd, sample_list, cat_spec):
    """ Classify new samples as belonging to a catergory based on distance to categorical spectra

        Args:
            nd          nd array        array of PSDs
            sample_list list            sample names
            cat_spec    dataframe       categorical spectra

        Returns:
            --          dataframe       sample  category    prob cat 1  prob cat 2 ... prob cat n
    """
    tmp = []

    for key in cat_spec:
        tmp.append([PSD_sym_KL(psd, cat_spec[key]) for psd in nd])

    KL = np.array(tmp).T

    # This is a confusing formula
    # amounts to: for a given sample, prob of belonging to class k is:
    #   (1 / KL_k) / sum_k(KL_i) = sum_k\i(KL_i) / sum_k(KL_i)
    prob = ((1 / KL).T / (1 / KL).sum(axis=1)).T

    row_masks = np.array([row == row.max() for row in prob])
    cats = [cat_spec.columns[mask][0] for mask in row_masks]

    items = [('label', cats)] + [('P({})'.format(lab), p) for lab, p in zip(cat_spec.columns, prob.T)]
    df = pd.DataFrame.from_items(items)
    df.index = sample_list

    return df
