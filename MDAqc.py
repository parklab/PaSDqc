#!/usr/bin/env python

# MDAqc.py - top level executible for MDAqc package
#
# v 0.0.12
# rev 2017-03-23 (MS: Add gold standard spectra to dendrogram plot)
# Notes:

import os
import argparse
import pathlib2
import pandas as pd
import numpy as np
import multiprocessing as mp
import plotly.offline as py

from src import mappable_positions
from src import PSDTools
from src import extra_tools
from src import plotly_tools
from src import report_writer

def extract(args):
    """ Extract coverage at uniquely mappable positions in a bam file
    """
    # Make out dir if it does not exist
    out_dir = os.path.join(args.dir, 'cov')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get bam paths from file if necessary
    if args.file_list:
        with open(args.file_list, 'r') as f:
            args.in_file = [line.rstrip('\n') for line in f]

        f.close()

    # Search a dir for bams if necessary
    elif args.dir_in:
        p = pathlib.Path(args.dir_in)
        file_list = sorted(p.glob("*.bam"))
        args.in_file = [str(f) for f in file_list]

    chrom_list = extra_tools.chroms_from_build(args.build)

    pool = mp.Pool(processes=args.n_procs)
    for f_in in args.in_file:
        fn_args = (f_in, out_dir, chrom_list)
        kwargs = {'build': args.build, 'map_qual': args.map_qual, 'clean': args.clean}
        pool.apply_async(mappable_positions.extract_coverage, args=fn_args, kwds=kwargs)

    pool.close()
    pool.join()

def PSD(args):
    """ Calculate power spectral densities from a directory of coverage files
    """
    dir_in = pathlib2.Path(args.dir)
    dir_search = dir_in / 'cov'

    pattern = "*.cov"
    if args.pattern:
        pattern = '*' + args.pattern + pattern

    f_list = sorted(dir_search.glob(pattern))
    samples = [f.name.split('.chr')[0] for f in f_list]
    samples = set(samples)
    # print(samples)

    pool = mp.Pool(processes=args.n_procs)

    for sample in samples:
        # _build_and_save(dir_in, sample=sample)
        fname = sample + ".chroms.spec"
        fout = dir_in / fname

        # _build_and_save(dir_search, sample, args.clean, fout)
        p = pool.apply_async(_build_and_save, (dir_search, sample, args.clean, fout))
        # p.get()
        # pool.apply_async(_build_and_save, args=(dir_in), kwds={'sample': sample})

    pool.close()
    pool.join()

def _build_and_save(dir_in, sample, clean, fout):
    """ Wrapper to build and save SamplePSD objects
    """

    print(fout)
    psd = PSDTools.SamplePSD.build_from_dir(str(dir_in), sample=sample, clean=clean)
    psd.save(str(fout))

def report(args):
    p = pathlib2.Path(args.dir)
    file_list = sorted(p.glob('*.chroms.spec'))

    # Chrom warnings
    sample_list = [f.name.split('.chroms.spec')[0] for f in file_list]
    df_list = [pd.read_table(str(f), index_col=0) for f in file_list]
    psd_list = [PSDTools.SamplePSD(df, name=s) for df, s in zip(df_list, sample_list)]
    j_list = [psd.KL_div_by_chrom() for psd in psd_list]
    df_chrom = extra_tools.summarize_KL_div_by_chrom(j_list, sample_list)

    # Categorization
    avg_list = [psd.avg_PSD() for psd in psd_list]
    nd = np.array(avg_list)
    freq = np.array(df_list[0].index.tolist())
    cat_spec = pd.read_table(args.cat_spec, index_col=0)
    df_clust = extra_tools.classify_samples(nd, sample_list, cat_spec)

    # ACF
    lags = np.array([0, 1e2, 2e2, 3e2, 4e2, 5e2, 7e2, 8e2, 1e3, 5e3, 1e4, 2e4, 3e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 1e6])
    nd_acf = np.array([extra_tools.PSD_to_ACF(freq, avg, lags) for avg in nd]) 
    df_var = pd.DataFrame(nd_acf[:, 0], columns=['variance'], index=sample_list)

    # Cluster plot
    dend = plotly_tools.dendrogram(nd, sample_list, args.cat_spec)
    div_dend = py.plot(dend, output_type='div')

    # PSD plot
    psd = plotly_tools.PSD_plot(freq, nd, sample_list)
    div_psd = py.plot(psd, output_type='div')

    # ACF plot
    acf = plotly_tools.ACF_plot(lags, nd_acf, sample_list)
    div_acf = py.plot(acf, output_type='div')

    # Chrom plot
    chrom = plotly_tools.chrom_KL_plot(j_list, sample_list)
    div_chrom = py.plot(chrom, output_type='div')

    # Report generation
    df = df_var.join(df_clust).join(df_chrom)
    fout_name = args.rep_name + '.html'
    fout = p / fout_name
    
    report_writer.writer(df, div_dend, div_psd, div_acf, div_chrom, str(fout))

def QC(args):
    extract(args)
    PSD(args)
    report(args)

def parse_args():
    parser = {}
    parser['argparse'] = argparse.ArgumentParser(description='MDAqc: tools for MDA quality control checks')
    parser['subparse'] = parser['argparse'].add_subparsers()

    parser['extract'] = parser['subparse'].add_parser('extract', help='Extract depth at uniquely mappable positions from bam file(s)')
    parser['extract'].add_argument('-i', '--in_file', default=None, nargs='+', help='input bam file(s)')
    parser['extract'].add_argument('-f', '--file_list', default=None, help='file of bam paths to use instead of -i')
    parser['extract'].add_argument('-d', '--dir_in', default=None, help='dir to search for bam files (instead of -i or -f)')
    parser['extract'].add_argument('-o', '--out_dir', default='.', help='output directory')
    parser['extract'].add_argument('-n', '--n_procs', default=1, type=int, help='number of cores to use')
    parser['extract'].add_argument('-b', '--build', default='grch37', help='genome build')
    parser['extract'].add_argument('-q', '--map_qual', default=30, help='mapping quality threshold')
    parser['extract'].add_argument('--noclean', dest='clean', action='store_false', default=True, help='Do not delete pos / cov files automatically')
    parser['extract'].set_defaults(func=extract)

    parser['PSD'] = parser['subparse'].add_parser('PSD', help='Estimate power spectral densities from coverage files')
    parser['PSD'].add_argument('-d', '--dir', default=None, help='directory in which to search for coverage files')
    parser['PSD'].add_argument('-n', '--n_procs', default=1, type=int, help='number of cores to use')
    parser['PSD'].add_argument('-p', '--pattern', default=None, help='pattern to match when finding samples')
    parser['PSD'].add_argument('--noclean', dest='clean', action='store_false', default=True, help='Do not delete pos / cov files automatically')
    parser['PSD'].set_defaults(func=PSD)

    parser['report'] = parser['subparse'].add_parser('report', help='Generate HTML report of MDA sample quality')
    parser['report'].add_argument('-d', '--dir', required=True,  default=None, help='directory in which to search for sample PSD files')
    parser['report'].add_argument('-c', '--cat_spec', default='db/categorical_spectra_1x.txt', help='path to file of categorical spectra')
    parser['report'].add_argument('-o', '--rep_name', default='MDAqc_report', help='file name for html report')
    parser['report'].set_defaults(func=report)

    parser['QC'] = parser['subparse'].add_parser('QC', help='Run all MDAqc steps, starting with BAM files and ending with an html report')
    parser['QC'].add_argument('-i', '--in_file', default=None, nargs='+', help='input bam file(s)')
    parser['QC'].add_argument('-f', '--file_list', default=None, help='file of bam paths to use instead of -i or -d')
    parser['QC'].add_argument('-d', '--dir_in', default=None, help='dir to search for bam files (instead of -i or -f)')
    parser['QC'].add_argument('-o', '--dir', default='data/', help='output directory')
    parser['QC'].add_argument('-n', '--n_procs', default=1, type=int, help='number of cores to use')
    parser['QC'].add_argument('-b', '--build', default='grch37', choices=['grch37, hg19'], help='genome build')
    parser['QC'].add_argument('-q', '--map_qual', default=30, type=int, help='mapping quality threshold')
    parser['QC'].add_argument('-c', '--cat_spec', default='db/categorical_spectra_1x.txt', help='path to file of categorical spectra')
    parser['QC'].add_argument('-r', '--rep_name', default='MDAqc_report', help='file name for html report')
    parser['QC'].add_argument('-p', '--pattern', default=None, help='pattern to match when finding samples')
    parser['QC'].add_argument('--noclean', dest='clean', action='store_false', default=True, help='Do not delete pos / cov files automatically')
    parser['QC'].set_defaults(func=QC)

    return parser['argparse'].parse_args()

if __name__ == "__main__":
    args = parse_args()
    args.func(args)
