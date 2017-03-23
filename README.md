# MDAqc
Comprehensive evaluation of single cell MDA library quality.

## Purpose
Anyone who works with single cell mulitple displacement amplified (MDA) libraries knows that this technology is a fickle beast. Sometimes the amplification works great. Sometimes it works terribly. A central challenge is thus to identify high quality libraries before spending a bunch of time / money on deep sequencing and variant analysis.

MDAqc provides a simple command line tool to rapidly evaluate very low coverage (0.1 - 1x) sequencing of MDA libraries for quality using power spectral density estimation techniques. MDAqc provides:
+ A simple "good" / "bad" label for each library evaluated
+ Sample clustering based on library quality similarity
+ Identification of individual chromosomes displaying particularly poor amplification
+ An html report full of fun, interactive plots that you can explore to your hearts content.

## Requirements
0. linux / macOS
1. python 3.5+
    + numpy
    + scipy
    + pandas
    + matplotlib
    + seaborn
    + plotly
    + astropy
2. samtools

## Installation
1. clone this repository `git clone https://github.com/parklab/MDAqc.git`
2. Install the necessary python packages
    + e.g. `pip install numpy` or `conda install numpy`
3. Install [samtools](http://www.htslib.org/download/) and make sure its on your PATH.

## Usage
0. Preprocessing bam files: MDAqc requires indexed BAM files aligned to either the GrCh37 or hg19 reference genomes as input. To ensure consistency across samples, we recommend sequencing (or downsampling) all libraries to either 1x or 0.1x prior to running MDAqc.

1. Quick usage: `python MDAqc.py QC -d <dir/of/bams> -o <my/output/dir> -c db/categorical_spectra_[XX].txt`
   * `db/categorical_spectra_[XX].txt` is the gold-standard spectra corresponding to the depth of your samples (see FAQ for more)

2. Long usage: `python MDAqc.py QC [OPTIONS]`

   ```
   Options:
   -h        display help message and exit
   -i        list of bam files to analyze [bam1 [bam2 ...]]
   -f        file containing paths to bam files to analyze (instead of -i, -d)
   -d        directory to search for bam files to analyze (instead of -i, -f)
   -o        output directory (will be created if it does not exist)
   -c        categorical spectra file
                NOTE: if using generic spectra, be sure to specify the correct depth (see FAQ).
   -n        number of threads to use
   -b        Build {grch37, hg19}
   -q        Mapping quality for read extraction
   -r        name of html report
   ```

3. Even longer useage `python MDAqc.py -h`

## FAQ
1. *How do I interpret the results in the html report?*
   * Summary table:
      + __Variance__: the higher this number, the worse the overall quality of the library
      + __label__: inferred quality of the library based on gold-standard performance
      + __P(good) / P(bad)__: the probability of label assignment
      + __Chrom: pass__: list of chromosomes with consistent quality
      + __Chrom: warn__: list of chromosomes with slightly outlying quality
      + __Chrom: fail__: list of chromosomes with signficantly outlying quality
   * Sample Clustering: hopefully this is self-explanatory
   * Chromosome outlier plots: chromosome quality by sample
      + The higher the y-value, the worse the quality. Warn chromosomes are one standard deviation above the average. Fail chromosomes are two standard deviations above.
      + NOTE: due to a short coming in the plotly API, all points for all samples are plotted at startup. Selecting a particular sample from the dropdown menu will display a sane plot.
   * Sample Periodograms:
      + Please see our paper for details on interpretting this plot.
   * Sample Autocorrelation:
      + The most interesting feature of this plot is where it approaches zero. This provides an estimate of the maximum size of the amplicons in a library. 

2. *How does MDAqc categorize a sample as good or bad?*

   MDAqc compares the behavior of the new sample to previously calculated gold-standard spectra. MDAqc is distributed with generic gold standard-spectra for depths corresponding to 1x (`db/categorical_spectra_1x.txt`) and 0.1x (`db/categorical_spectra_0.1x.txt`). If you are using these generic spectra, please be sure to specify the correct file corresponding to the depth of your samples using the `-c` option.

3. *Can I define my own gold-standard spectra?*

   Yes! MDAqc is also a python library. It includes a method to generate new gold standard spectra using previously analyzed samples. Within the MDAqc dir, start python
   ```
   >>> from src import extra_tools
   >>> freq, nd, sample_list = extra_tools.mk_ndarray("/path/to/analysis/dir")
   >>> cat_spec = extra_tools.mk_categorical_spectra(freq, nd, labels)
   ```
   where `labels` is a list of the form `['good', 'good', 'bad', ...]` corresponding to user-assigned classifications to each sample in `sample_list`. Then save the `cat_spec` to a file. Tell MDAqc to use this new file using the `-c` option the next time you run the program.

4. *I deleted the html report. Do I have to rerun the whole pipeline?*

   No. The steps of MDAqc can be run separately. To regenerate the report, simply run `python MDAqc.py report -d /path/to/analysis/dir -c db/categorical_spectra_XX.txt`. Regenerating a report generally takes less than 1 minute.

   To see all functions of MDAqc, run `python MDAqc.py -h`.

5. *Do you have some example data I can try out on?*

   Yes. The source code is distributed with an `examples` directory. For size reasons, we do not distribute BAM files, but we do provide example power spectral densities generated by MDAqc. You can compile these into an example report by running
   ```
   python MDAqc.py report -d examples/example_01 -c db/categorical_spectra_0.1x.txt
   ```
   This will recreate the html report in the directory `examples/example_01`.

6. *Doesn't plotly upload its plots to the interwebs? Can anyone see my MDAqc plots??*

   Rest easy, none of your plots are uploaded to the internet. We use the plotly javascript API to embed the plots directly into the html report. This means that a report is entirely self-contained, so you can send it around to collaborators without fear of plots not loading. It also means that the reports can be fairly large (several MBs), as the data required to generate the plots must also be embedded in the html file.

5. *I want to use MDAqc to analyze the quality of samples aligned to a different genome build. Is this possible*

   Yes, but you'll need to a bed file of uniquely mappable positions for the genome. Contact us for help with this.

## Support
__Maxwell Sherman__: maxwell\_sherman {at} hms.harvard.edu

__Peter Park__: peter\_park {at} hms.harvard.edu
