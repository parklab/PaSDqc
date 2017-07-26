# amplicon.py - classes for fitting amplicon distributions from PSDs
#
# v 1.0.10
# rev 2017-07-26 (MS: more robust curve fitting)
# Notes:

import scipy.optimize
import scipy.special
import scipy.stats
import numpy as np
import pandas as pd

from . import extra_tools

def func_logis(x, inter, asym, xmid, scal):
    """ Return logistic curve values for a given set of point
    """
    return inter + asym / (1 + np.exp(-(x - xmid) / scal))

class PSDLogis(object):
    # psd = {}
    # freq = {}
    # inter = 0
    # asym = 0
    # xmid = 0
    # scal = 0
    # var = 0

    def __init__(self, freq, psd):
        self.freq = {}
        self.psd = {}

        self.freq['raw'] = freq
        self.psd['raw'] = psd

        # self.freq['scale'] = 10 * np.log10(freq)
        # self.freq['cut'] = self.freq['scale'][freq < 1e-3]

    def fit_logistic(self, bulk="db/bulk_1x.spec"):
        """ Fit a logistic curve with intercept to power spectral density
            The logistic is fit to the "flipped" psd where the x-axis is in BP not 1/BP
        """
        psd = self.psd['raw']
        freq = self.freq['raw']

        # Load bulk spectrum
        psd_bulk = pd.Series.from_csv(bulk, index_col=0, header=None, sep="\t").as_matrix()
        # psd_bulk = np.loadtxt(bulk)

        # Normalize psd by bulk spectrum
        # psd_norm = 10 * np.log10(psd / psd_bulk)
        psd_norm = 10 * np.log10(psd / psd_bulk)
        psd_cut = psd_norm[freq < 1e-3]

        freq_cut = freq[freq < 1e-3]
        freq_scale = -np.log10(freq_cut)

        # Fit the curve
        try:
            popt, pcov = scipy.optimize.curve_fit(func_logis, freq_scale, psd_cut)

            # Make sure the scale parameter is positive
            # May be negative for bulk samples
            if popt[-1] > 0:
                self.success=True
            else:
                self.success=False

        # Sometimes curve fitting fails for bulk samples
        except RuntimeError:
            popt = [np.nan, np.nan, np.nan, np.nan]
            self.success=False

        self.psd['logis'] = psd_cut
        self.freq['logis'] = freq_scale

        self.popt = popt
        self.inter, self.asym, self.xmid, self.scal = popt

    def predict_vals(self, point=None):
        """ Return predictions for values of the logistic curve
        """
        if point:
            if isinstance(point, str):
                vals = getattr(self, point)
            else:
                 vals = point

        else:
            vals = self.freq['logis']

        return func_logis(vals, *self.popt)

    def amplicon_range(self):
        """ Calculate mean and 95% upper and lower bounds on amplicon sizes using logistic fit
        """
        # var = self.scal**2 * np.pi**2 / 3
        if self.success:
            print(self.xmid, self.scal)
            x = scipy.stats.logistic.rvs(loc=self.xmid, scale=self.scal, size=100000)
            x_scale = 10**x
            self.mean = np.mean(x_scale)
            self.median = np.median(x_scale)
            self.lower_95 = np.percentile(x_scale, 5)
            self.upper_95 = np.percentile(x_scale, 95)

        else:
            print('pass')
            self.mean = 0
            self.median = 0
            self.lower_95 = 0
            self.upper_95 = 0

        return [self.median, self.mean, self.lower_95, self.upper_95]

        # self.upper_q = self.xmid + self.scal * np.log(0.95 / 0.05)
        # self.lower_q = self.xmid + self.scal * np.log(0.05 / 0.95)

        # self.mean = 10**(self.xmid)
        # self.upper_95 = 10**(self.upper_q)
        # self.lower_95 = 10**(self.lower_q)

        # return [self.mean, self.lower_95, self.upper_95]

    def logistic_dist(self, point=None):
        """ Return the logistic distribution of the amplicon sizes
            (in log10 coordinates)
        """
        if point:
            freq = point

        else:
            freq = -np.log10(self.freq['raw'][self.freq['raw'] < 3e-3])

        if self.success:
            pdf = scipy.stats.logistic.pdf(freq, loc=self.xmid, scale=self.scal)

        else:
            pdf = np.zeros(len(freq))

        self.freq['dist'] = freq

        return pdf

class AmplDist(object):
    def __init__(self, freq, psd):
        self.freq = {}
        self.psd = {}
        self.popt = {}

        self.freq['raw'] = freq
        self.psd['raw'] = psd

    @staticmethod
    def func_logis(x, inter, asym, xmid, scal):
        """ Return logistic curve values for a given set of point
        """
        return inter + asym / (1 + np.exp(-(x - xmid) / scal))

    @staticmethod
    def func_erf(x, inter, asym, mu, sigma):
        """ Return erf curve values for a given set of point
            The erf function is fit s.t. mu and sigma have a Gaussian interpretation
        """
        return inter + asym * scipy.special.erf((x-mu) / (np.sqrt(2) * sigma))

    @staticmethod
    def func_gamma(x, inter, asym, alpha, beta):
        """ Return logistic curve values for a given set of point
        """
        return inter + asym * scipy.special.gammainc(alpha, beta*x)

    def fit_curve(self, method='erf', bulk="bulk_1x.smooth3.spec", shift=0):
        """ Fit a curve of specified type to power spectral density
            The logistic is fit to the "flipped" psd where the x-axis is in BP not 1/BP
        """
        psd = self.psd['raw']
        freq = self.freq['raw']

        f_fit = getattr(self, "func_{}".format(method))

        # Load bulk spectrum
        f_bulk = extra_tools.get_data_file(bulk)
        psd_bulk = pd.Series.from_csv(f_bulk, index_col=0, header=None, sep="\t").as_matrix()
        # psd_bulk = np.loadtxt(bulk)

        # Normalize psd by bulk spectrum and shift
        psd_norm = 10 * np.log10(psd / psd_bulk)
        psd_cut = psd_norm[freq < 1e-3]

        freq_cut = freq[freq < 1e-3]
        period = -np.log10(freq_cut) - shift

        # Fit the curve
        try:
            popt, pcov = scipy.optimize.curve_fit(f_fit, period, psd_cut, method='trf')

            # Make sure the scale parameter is positive
            if popt[-1] < 0:
                popt[-1] = np.abs(popt[-1])

            self.success=True

        # Sometimes curve fitting fails for bulk samples
        except RuntimeError:
            popt = [np.nan, np.nan, np.nan, np.nan]
            self.success=False

        self.psd[method] = psd_cut
        self.freq[method] = period

        popt2 = np.append(popt, shift)
        self.popt[method] = popt2
        # self.inter, self.asym, self.xmid, self.scal = popt

    def predict_vals(self, method='erf', point=None):
        """ Return predictions for values of the logistic curve
        """
        fit_fnc = getattr(self, "func_{}".format(method))

        if point:
            if isinstance(point, str):
                vals = getattr(self, point)
            else:
                 vals = point

        else:
            vals = self.freq[method]

        return fit_fnc(vals, *self.popt[method])

    def amplicon_range(self, method='erf'):
        """ Calculate mean and 95% upper and lower bounds on amplicon sizes using logistic fit
        """

        if self.success:
            popt = self.popt[method]
            shift = popt[-1]

            if method == 'erf':
                dist = scipy.stats.norm(loc=popt[2], scale=popt[3])
            elif method == 'logis':
                dist = scipy.stats.logistic(loc=popt[2], scale=popt[3])
            elif method == 'gamma':
                dist = scipy.stats.gamma(a=popt[2], scale=1/popt[3])

            x = dist.rvs(size=100000)
            x_scale = 10**(x + shift)
            self.mean = np.mean(x_scale)
            self.median = np.median(x_scale)
            self.lower_95 = np.percentile(x_scale, 5)
            self.upper_95 = np.percentile(x_scale, 95)

        else:
            print('pass')
            self.mean = 0
            self.median = 0
            self.lower_95 = 0
            self.upper_95 = 0

        return [self.median, self.mean, self.lower_95, self.upper_95]

    def amplicon_dist(self, method='erf'):
        """ Calculate the distribution of amplicon sizes based on a curve fit
        """
        if self.success:
            popt = self.popt[method]

            if method == 'erf':
                dist = scipy.stats.norm(loc=popt[2], scale=popt[3])

            elif method == 'logis':
                dist = scipy.stats.logistic(loc=popt[2], scale=popt[3])
            elif method == 'gamma':
                dist = scipy.stats.gamma(a=popt[2], scale=1/popt[3])

            lower = dist.ppf(0.001)
            upper = dist.ppf(0.999)
            vals = np.linspace(lower, upper, 100)
            pdf = dist.pdf(vals)

        else:
            pdf = np.zeros(100)
            vals = np.linspace(2, 6, 100)

        self.freq['dist'] = 10**vals

        return pdf

    @staticmethod
    def param_names(method):
        """ Return parameter names for a curve fit method
        """
        d = {'erf': ['intercept', 'scale', 'mu', 'sigma', 'shift'],
             'logis': ['intercept', 'asymptote', 'xmid', 'scale', 'shift'],
             'gamma': ['intercept', 'scale', 'alpha', 'beta', 'shift'],
            }

        return d[method]

