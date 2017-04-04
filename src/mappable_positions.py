# mappable_positions.py - methods to extract coverage from bam files
#
# v 0.0.8
# rev 2017-03-21 (MS: removal of pos files after cov extraction)
# Notes:

import pathlib2
import pandas as pd
import numpy as np
# import pybedtools as pybt
import subprocess
import os
import argparse

def exec_cmd(cmd, verbose=True):
    if verbose:
        print("Executing: {}".format(cmd))

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()

def map_to_bed(f_map, dir_out, build='grch37'):
    """ Split a file listing unique regions across genome into bed files per chrom
        assumes columns are chrom   start   end ...
        columns beyond the third are ignored
    """

    df = pd.read_table(f_map, header=None, index_col=None)
    df = df.rename(columns=dict(zip([0, 1, 2], ['chrom', 'start', 'end'])))

    # chrm_list = df.chrom.unique()
    for name, group in df.groupby('chrom'):
        print('processing: {}'.format(name))
        name = to_hg19_format(name)

        if build.startswith('grch'):
            group = group.replace('chr', '', regex=True)

        elif build.starswith('hg'):
            chrom = group.chrom.astype(str)

            if not chrom.iloc[0].startswith('chr'):
                group.chrom = chrom.replace('^', 'chr', regex=True)

        f_name = "{}.{}.map.bed".format(build, name)
        f_out = os.path.join(dir_out, f_name)

        group.to_csv(f_out, columns=['chrom', 'start', 'end'], sep="\t", header=None, index=None)

def uniq_pos_by_chrm(bam_file, chrom, out_file, map_qual=30, build='hg19'):
    """ Extract the uniquely mappable positions which occur at the start of a read
        in a chromosome specific fashion
    """
    map_file = lookup_map_file(build, chrom)

    kwargs = {'bam_file': bam_file, 
              'chrom': chrom, 
              'out_file': out_file, 
              'map_qual': map_qual,
              'map_file': str(map_file)
    }

    cmd = "samtools view -q {map_qual} -L {map_file} {bam_file} {chrom} | cut -f 3,4 | uniq > {out_file}".format(**kwargs) 

    exec_cmd(cmd)
    # exec_cmd(cmd, verbose=True)

def depth_at_pos(bam_file, pos_file, chrom, out_file, map_qual=30, clean=True):
    """ Extract depth at positions specified in a file.
    """
    kwargs = {'bam_file': bam_file, 
              'pos_file': pos_file,
              'chrom': chrom, 
              'out_file': out_file, 
              'map_qual': map_qual,
    }

    cmd = "samtools view -uh -q {map_qual} {bam_file} {chrom} | samtools depth -b {pos_file} - > {out_file}".format(**kwargs)

    if clean:
        cmd += "; rm {}".format(pos_file)

    exec_cmd(cmd)
    # exec_cmd(cmd, verbose=True)

def extract_coverage(bam_file, out_dir, chrom_list, build='hg19', map_qual=30, clean=True):
    """ Efficiently extract coverage at uniquely mappable positions 
        occuring at the starts of reads in a chromosome specific manner

        Proceeds in two steps
          1) Extract uniq map'ble positions from bam file
          2) Determine coverage at extracted positions
    """
    # print('something')
    p_bam = pathlib2.Path(bam_file)
    p_out = pathlib2.Path(out_dir)

    # Should make this parallel!!!
    for chrom in chrom_list:
        print("Processing: {}".format(chrom))
        chrom_lab = to_hg19_format(chrom)
        pos_file = p_out / p_bam.with_suffix('.{}.map.pos'.format(chrom_lab)).name
        cov_file = p_out / p_bam.with_suffix('.{}.map.pos.cov'.format(chrom_lab)).name

        uniq_pos_by_chrm(bam_file, chrom, pos_file, map_qual, build)
        depth_at_pos(bam_file, pos_file, chrom, cov_file, map_qual, clean)

def to_hg19_format(chrom):
    try:
        if not chrom.startswith('chr'):
            chrom = 'chr{}'.format(chrom)

    except AttributeError:
        chrom = 'chr{}'.format(chrom)

    return chrom

def to_grch_format(chrom):
    try:
        if chrom.startswith('chr'):
            chrom = chrom.split('chr')[-1]

    except AttributeError:
        pass
        
    return chrom

def lookup_map_file(build, chrom):

    chrom = to_hg19_format(chrom)

    path = pathlib2.Path("db/")
    map_file = sorted(path.glob('{}.{}.map.bed'.format(build, chrom)))[0]

    if not map_file:
        raise IOError("Uh oh, there's no map file matching {} and {}".format(build, chrom))

    return map_file

def parse_args():
    parser = argparse.ArgumentParser(description='Extract coverage at uniquely mappable positions from bam file(s)')
    parser.add_argument('-i', '--in_file', default=None, nargs='+', help='input bam file(s)')
    parser.add_argument('-f', '--file_list', default=None, help='file of bam paths to use instead of -i')
    parser.add_argument('-o', '--out_dir', default='.', help='output directory')
    parser.add_argument('-n', '--n_procs', default=1, help='number of cores to use')
    parser.add_argument('-b', '--build', default='grch37', help='genome build')
    parser.add_argument('-q', '--map_qual', default=30, help='mapping quality threshold')

    args = parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
