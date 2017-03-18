# PSDTools.py - classes for MDAqc Power Spectral Density calculations
#
# v 0.0.6
# rev 2017-03-18 (MS: Added documentation)
# Notes:

import pandas as pd
import numpy as np
import astropy.stats
import statsmodels.nonparametric
import matplotlib.pyplot as plt
import seaborn as sns
import re
import pathlib2

class ChromPSD(object):
    """ Lombe-Scargle PSD estimation for a single Chromosome
        
        Normalizes depth to be copy ratio neutral
        Splits by p- and q-arm to avoid centromere problems
        Performs a (modified) Welch estimation procedure to smooth the periodogram

        Args:
            f       file path to chromosome coverage file

        Attributes:
            chrom       chromosome number
            df          data frame containing position and depth information for the chromosome
    """
    def __init__(self, f):
        self.df = pd.read_table(f, header=None, names=['chrom', 'pos', 'depth'], index_col=None)
        self.chrom = str(self.df.chrom.iloc[0])
        self._norm_depth()

    def _norm_depth(self):
        """ Normalize depth by mean of all observations
        """
        self.df.depth = self.df.depth / self.df.depth.mean()

    @staticmethod
    def _welch_seg_bounds(pos, l_seg, p_overlap):
        """ Define boundaries of segments for Welch PSD based on length and desired overlap

            Args:
                pos         series      chromosome positions of observations
                l_seq       int         length of segments
                p_overlap   float       percentage by which segments should overlap

            Returns:
                starts      1d array    array of segment start positions
                ends        1d array    array of segment end positions
        """
        step = l_seg - p_overlap * l_seg
        starts = np.arange(pos.iloc[0], pos.iloc[-1], step)
        ends = np.arange(pos.iloc[0]+l_seg, pos.iloc[-1], step)
        ends[-1] = pos.iloc[-1]

        return starts, ends

    @staticmethod
    def _welch_PSD(freq, pwr, df, starts, ends, count=1, verbose=False):
        """ Perform Welch spectral density estimation
            NOTE: requires an initial freq / pwr vector

            Args:
                freq        1d array        frequencies at which to calculate PSD
                pwr         1d array        initialized pwr array
                df          dataframe       read depth at each position
                starts      1d array        Welch segment start positions
                ends        1d array        Welch segment end positions

            Kwargs:
                counts      int             number of PSDs already included in pwr
                verbose     bool            calculate verbosely

            Returns:
                freq        1d array        frequencies at which PSD calculates
        """
        for start, end in zip(starts, ends):
            if verbose:
                print(start, end)

            mask = (df.pos >= start) & (df.pos < end) 
            df_tmp = df[mask]

            if(len(df_tmp) > 0):
                count += 1
                pwr += astropy.stats.LombScargle(df_tmp.pos, df_tmp.depth).power(freq, normalization='psd')

        return freq, pwr, count

    @staticmethod
    def PSD_LS_auto(df, l_seg=1e7, p_overlap=0.5, verbose=False):
        """ Lomb-Scargle PSD estimate at automatically determined frequencies

            Args:
                df          dataframe   read depth at each position

            Kwargs:
                l_seg       int         length of Welch segments
                p_overlap   float       amount Welch segments should overlap
                verbose     bool        calculate verbosely

            Returns:
                freq_agg    1d array    frequencies at which PSD calculate
                pwr_agg     1d array    power of PSD at each frequency
                count       int         number of Welch PSDs included in pwr_agg
        """

        # Create boundaries for Welch PSD segments
        starts, ends = ChromPSD._welch_seg_bounds(df.pos, l_seg, p_overlap)

        # Initialize PSD with first segment
        df_red = df[df.pos < ends[0]]
        freq_agg, pwr_agg = astropy.stats.LombScargle(df_red.pos, df_red.depth).autopower(normalization='psd')

        # Welch PSD
        freq_agg, pwr_agg, count = chromPSD._welch_PSD(freq_agg, pwr_agg, df, starts[1:], ends[1:], verbose=verbose) 

        return freq_agg, pwr_agg, count

    @staticmethod
    def PSD_LS_manual(df, freq, l_seg=1e7, p_overlap=0.5, verbose=False, inplace=True):
        """ Lomb-Scargle PSD estimate at manually specified frequencies

            Args:
                df          dataframe   read depth at each position
                freq        1d array    frequencies at which to calculate PSD

            Kwargs:
                l_seg       int         length of Welch segments
                p_overlap   float       amount Welch segments should overlap
                verbose     bool        calculate verbosely

            Returns:
                freq_agg    1d array    frequencies at which PSD calculate
                pwr_agg     1d array    power of PSD at each frequency
                count       int         number of Welch PSDs included in pwr_agg
        """
        # Create boundaries for Welch PSD segments
        starts, ends = ChromPSD._welch_seg_bounds(df.pos, l_seg, p_overlap)

        # Initialize PSD with first segment
        df_red = df[df.pos < ends[0]]
        pwr_agg = astropy.stats.LombScargle(df_red.pos, df_red.depth).power(freq, normalization='psd')

        # Welch PSD
        freq_agg, pwr_agg, count = ChromPSD._welch_PSD(freq, pwr_agg, df, starts[1:], ends[1:], verbose=verbose) 

        return pwr_agg, count

    def PSD_LS_chrom(self, f_cent, freq, l_seg=1e6, p_overlap=0.5, verbose=False):
        """ Centromere aware PSD estimation

            Args:
                f_cent      str         path to file of centromere locations
                freq        1d array    frequencies at which to calcualte PSD

            Kwargs:
                l_seg       int         length of Welch segments
                p_overlap   float       amount Welch segments should overlap
                verbose     bool        calculate verbosely

            Returns:
                freq_agg    1d array    frequencies at which PSD calculate
                pwr_agg     1d array    power of PSD at each frequency
                count       int         number of Welch PSDs included in pwr_agg
        """
        # load centromere file
        # split df by pos < start of cent and > end of cent
        # peform PSD analysis on two arms
        # merge
        df_cent = pd.read_table(f_cent, names=['start', 'end'], index_col=0)
        cent_start = df_cent.loc[self.chrom, :].start
        cent_end = df_cent.loc[self.chrom, :].end

        df_p = self.df[self.df.pos < cent_start]
        df_q = self.df[self.df.pos > cent_end]

        # count = 0

        print("{}.p arm".format(self.chrom))
        if len(df_p > 1e7):
            pwr_p, count_p = self.PSD_LS_manual(df_p, freq, l_seg, p_overlap, verbose=verbose)
        else:
            pwr_p = 0
            count_p = 0

        print("{}.q arm".format(self.chrom))
        if len(df_q > 1e7):
            pwr_q, count_q = self.PSD_LS_manual(df_q, freq, l_seg, p_overlap, verbose=verbose)
        else:
            pwr_q = 0
            count_q = 0

        self.freq = freq
        self.pwr = (pwr_p + pwr_q) / (count_p + count_q)

class SamplePSD(object):
    """ Lombe-Scargle PSD estimation for an entire sample (all chromosomes)
        
        Normalizes depth to be copy ratio neutral
        Splits by p- and q-arm to avoid centromere problems
        Performs a (modified) Welch estimation procedure to smooth the periodogram

        Args:
            f       file path to chromosome coverage file

        Attributes:
            name        name of sample
            df          data frame containing position and depth information for the chromosome
    """

    name = None
    df = None

    def __init__(self, df, name):
        self.name = name
        self.df = df

    @classmethod
    def build_from_dir(cls, d_path, sample=None):
        """ Build SamplePSD object from a directory of .cov files

            Args:
                d_path      str     path to directory of .cov files

            Kwargs:
                sample      str     sample name
        """
        pattern = "*.cov"
        if sample:
            pattern = sample + pattern

        p = pathlib2.Path(d_path)
        file_list = sorted(p.glob(pattern))
        name = cls.name_from_file(file_list[0])

        df = cls._build_dataframe(file_list)

        # chrom_list = cls.chroms_from_files(file_list, build)

        return cls(df, name)

    @staticmethod
    def chroms_from_files(file_list, build='grch37'):
        """ Determine chromosomes from file list
        """
        chrom_list = [re.search('chr[0-9XY]{1,}', str(f)) for f in file_list]

        if build.startswith('grch'):
            chrom_list = [re.sub('chr', '', s.group()) for s in chrom_list]

        return chrom_list

    @staticmethod
    def name_from_file(f):
        """ Get sample name from a file path

            Args:
                f       pathlib Path

            Returns
                name    str             name of sample
        """
        name = str(f.name).split('.')[0]

        return name

    @staticmethod
    def _build_dataframe(f_list):
        """ Build dataframe of a sample's PSDs -- one column per chromosome

            Args:
                f_list      list        pathlib paths to .cov files

            Returns:
                df          dataframe   dataframe of sample PSDs
        """
        # create ChromPSD objects
        psd_list = [ChromPSD(str(f)) for f in f_list]

        # perform PSD estimation
        freq = np.linspace(1e-6, 5e-3, 8000)
        [psd.PSD_LS_chrom("db/grch37.centromeres.bed", freq=freq) for psd in psd_list]

        # Assemble into dataframe
        items = [('freq', freq)]
        items += [(psd.chrom, psd.pwr) for psd in psd_list]

        df = pd.DataFrame.from_items(items, )
        df = df.set_index('freq')

        return df

    # def avg_PSD(self, f_chrom_sizes):
    def avg_PSD(self):
        """ Calculate median of chromosomal PSDs
        """
        # chrom_list = self.df.columns
        # chrom_sizes = pd.read_table(f_chrom_sizes, header=None, index_col=0, names=['size'])
        # weights are size of chromosome divided by total size of genome
        # sizes = np.array([chrom_sizes.loc[c].as_matrix()[0] for c in chrom_list])
        # weights = pd.Series(sizes.astype(float) / sizes.sum(), index=chrom_list)
        # print(weights)

        # multiply each chrom column by its weight and sum to get avg
        # df = self.df.iloc[:, 1:].copy() * weights
        avg = self.df.median(axis=1)

        return avg

    def KL_div_by_chrom(self):
        """ Calculate each chromosome's PSD KL divergence from the sample's average PSD
        """
        avg = self.avg_PSD()

        kl_for = self.df.multiply(1 / avg, 'index')
        kl_back = (1 / self.df).multiply(avg, 'index')

        j = (kl_for + kl_back - 2).sum(axis='index')
        j *= 1 / len(avg)

        return j

    def save(self, f_out):
        """ Save the sample's PSD dataframe

            Args:
                f_out   str     name of output file
        """
        self.df.to_csv(f_out, sep="\t", header=True, index=True)
