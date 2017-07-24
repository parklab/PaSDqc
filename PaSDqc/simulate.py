# simulate.py - methods for simulating amplification profiles
#
# v 1.0.8
# rev 2017-07-24 (MS: created)
# Notes:

import numpy as np
import pandas as pd
import scipy.stats

from . import PSDTools

def simulate_logis_profile(start, end, mu, sigma, depth=1):
    """ Simulate linear amplification using a logistic distribution

        Start position of an amplicon is randomly chosen from a uniform distribution
        over the [start, end]
    """

    chrom_len = end - start
    read_depth = np.zeros(chrom_len, dtype=int)
    read_total = 0.

    logis = scipy.stats.logistic(loc=mu, scale=sigma)
    unif = scipy.stats.randint(low=start, high=end)

    while read_total / chrom_len < 1:
        # Amplicon length
        amp_log_len = logis.rvs(1)[0]
        while amp_log_len < 0:
            amp_log_len = logis.rvs(1)[0]

        amp_len = np.int(10**logis.rvs(1)[0])

        # Amplicon start
        amp_start = unif.rvs(1)[0]

        while amp_start + amp_len > end:
            amp_start = unif.rvs(1)[0]

        amp_end = amp_start + amp_len

        # Assign values
        read_depth[amp_start:amp_end+1] += 1
        read_total += amp_len

    return read_depth

def simulate_erf_profile(start, end, mu, sigma, depth=1):
    """ Simulate linear amplification using a normal (erf) distribution

        Start position of an amplicon is randomly chosen from a uniform distribution
        over the [start, end]
    """

    chrom_len = end - start
    read_depth = np.zeros(chrom_len, dtype=int)
    read_total = 0.

    norm = scipy.stats.norm(loc=mu, scale=sigma)
    unif = scipy.stats.randint(low=start, high=end)

    while read_total / chrom_len < 1:
        # Amplicon length
        amp_log_len = norm.rvs(1)[0]
        while amp_log_len < 0:
            amp_log_len = norm.rvs(1)[0]

        amp_len = np.int(10**norm.rvs(1)[0])

        # Amplicon start
        amp_start = unif.rvs(1)[0]

        while amp_start + amp_len > end:
            amp_start = unif.rvs(1)[0]

        amp_end = amp_start + amp_len

        # Assign values
        read_depth[amp_start:amp_end+1] += 1
        read_total += amp_len

    return read_depth

def simulate_gamma_profile(start, end, alpha, beta, shift=3, depth=1):
    """ Simulate linear amplification using a gamma distribution

        Start position of an amplicon is randomly chosen from a uniform distribution
        over the [start, end]
    """

    chrom_len = end - start
    read_depth = np.zeros(chrom_len, dtype=int)
    read_total = 0.

    gamma = scipy.stats.gamma(a=alpha, scale=1./beta)
    unif = scipy.stats.randint(low=start, high=end)

    while read_total / chrom_len < 1:
        # Amplicon length
        amp_log_len = gamma.rvs(1)[0]
        while amp_log_len < 0:
            amp_log_len = gamma.rvs(1)[0]

        amp_len = np.int(10**(gamma.rvs(1)[0]+shift))

        # Amplicon start
        amp_start = unif.rvs(1)[0]

        while amp_start + amp_len > end:
            amp_start = unif.rvs(1)[0]

        amp_end = amp_start + amp_len

        # Assign values
        read_depth[amp_start:amp_end+1] += 1
        read_total += amp_len

    return read_depth

def restrict_to_uniq_pos(read_depth, start, end, pos):
    pos_cut = pos[(pos > start) & (pos < end)]
    pos_shift = pos_cut - end

    return read_depth[pos_shift]
    
def simulate(f_sim, start, end, n_obs, *args, **kwargs):
    """ Peform a complete amplification simulation using a read depth density estimate
        Also performs PSD estimation

        Args:
            f_sim: simulation function (see above fns)
            start: start position of chromosome arm
            end:   end position of chromosome arm
            n_obs: number of observations to simulate
            *args: other arguments to pass to f_sim

        Kwargs:
            **kwargs: kwargs to pass to f_sim

        Returns:
            read depth data frame
            psd estimate
    """
    print("Simulating read depth")
    rd = f_sim(start, end+100000, *args, **kwargs)
    rd_s = pd.DataFrame(rd[:end], columns['depth'])
    rd_uniq = rd_s[rd_s.depth > 0].sample(n=n_obs, replace=False).sort_index()
    rd_uniq['pos'] = rd_uniq.index

    print("Performing PSD analysis")
    df = pd.DataFrame({'chrom': 3, 'pos': rd_uniq.index.tolist(), 'depth': rd_uniq.depth.as_matrix()})
    freq = np.linspace(1e-6, 5e-3, 8000)
    pwr, count = PSDTools.ChromPSD.PSD_LS_manual(df, freq, l_seg=1e6)

    return df, pwr / count
