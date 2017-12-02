# PaSD-qc
PaSD-qc ("Power Spectral Density-qc", pronounced "passed-qc") provides comprehensive evaluation of the properties and qualities of single cell whole-genome sequencing (scWGS) data.

## Purpose
Anyone who works with scWGS libraries knows that the technology is a fickle beast. Sometimes the whole-genome amplification works great. Sometimes it works terribly. Biases abound, and these biases can affect the accuracy of variant calls. A central challenge is thus to determine the quality of amplification prior to deep sequencing and to measure the biases prior to variant calling.

PaSD-qc provides a simple command line tool to rapidly evaluate the amplification properties and quality of scWGS sequencing using a custom, highly accurate power spectral density estimation technique. PaSD-qc provides:
+ A robust measure of read depth variance for each sample
+ An accurate estimate of the autocovariance patterns each sample
+ An estimate of the full distribution of amplicon sizes in each sample
+ Identification of copy number altered and poorly amplified chromosomes 
+ A comparison of samples based on quality
+ A simple "good" / "bad" label for each library evaluated
+ An html report full of fun, interactive plots that you can explore to your hearts content.

You can read more about the powers of PaSD-qc in [our paper at Nucleic Acids Research](https://doi.org/10.1093/nar/gkx1195) and by exploring the `examples` directory of this repository.

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
**PaSDqc is currently in beta. A stable release will be available soon.**
1. clone this repository `git clone https://github.com/parklab/PaSDqc.git`
2. cd into the newly created `PaSDqc` directory and run `make`
    * This will create the python distribution and install it as a site-package
    * It will also automatically (try to) install any missing python dependencies
3. Install [samtools](http://www.htslib.org/download/) and make sure its on your PATH.

## Usage
0. Preprocessing bam files: PaSD-qc requires indexed BAM files aligned to either the GrCh37 or hg19 reference genomes as input. To ensure consistency across samples, we recommend sequencing (or downsampling) all libraries to the same depth prior to running PaSD-qc. Depth recommendations:
    * 0.1X - basic evaluation of very low-coverage data to identify high quality samples for deep sequencing
    * 1X - General purpose quality evaluation. Best for identifying potential copy-altered chromosomes.
    * 5X - Accurate amplicon distribution evaluation for calibrating correlation-based variant calling methods. 

1. Quick usage: `PaSDqc QC -d <dir/of/bams> -o <my/output/dir> -c db/categorical_spectra_[XX].txt`
   * `db/categorical_spectra_[XX].txt` is the gold-standard spectra corresponding to the depth of your samples (see FAQ for more)

2. Long usage: `PaSDqc QC [OPTIONS]`

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

3. Even longer useage `PaSDqc -h`
    ```
    positional arguments:
        extract         Extract depth at uniquely mappable positions from bam file(s)
        PSD             Estimate power spectral densities from coverage files
        report          Generate HTML report of single-cell sample quality
        QC              Run all PaSDqc steps, starting with BAM files and ending with an html report
    ```
    PaSD-qc is comprised of three steps: 1) read depth extraction; 2) power spectral density estimation; 3) report generation. These steps can be run independently by the user or all at together by using the `QC` method.

## FAQ
0. *Why should I used PaSD-qc?*

    PaSD-qc offers unparalleled insight into the properties of single-cell whole-genome sequencing. While other QC approaches generally provide information only about read depth dispersion (often at a fixed-genomic scale), PaSD-qc provides comprehensive quality control. For example, PaSD-qc can infer the full distribution of amplicon sizes in a sample. This provides a principled method to determine what bandwidth to use when employing a single-cell specific caller which uses allele balance of nearby loci to perform accurate variant calling.

1. *How do I interpret the results in the html report?*
   * Summary table:
      + __Variance__: the higher this number, the worse the overall quality of the library
      + __label__: inferred quality of the library based on gold-standard performance
      + __P(good) / P(bad)__: the probability of label assignment
      + __Amplicon median__: estimated median amplicon size
      + __mean size__: estimated mean amplicon size
      + __lower size__: 5% lower bound on amplicon size range
      + __upper size__: 95% upper bound on amplicon size range 
      + __Chrom: pass__: list of chromosomes with consistent quality
      + __Chrom: gain?__: list of chromosomes with potential copy gains
      + __Chrom: loss?__: list of chromosomes with potential copy losses
      + __Chrom: fail__: list of chromosomes with signs of poor amplification quality
   * Sample Clustering: hopefully this is self-explanatory
   * Chromosome classification plots: chromosome quality by sample
      + The y-value is the KL-divergence of the chromosome from the sample average. Colors indicate pass (grey), possible gain (green), possible loss (blue), and poor quality (red).
   * Amplicon size distributions
      + Estimated distributions of amplicon sizes (in log coordinates)
   * Sample Periodograms:
      + Please see our paper for details on interpretting this plot.
   * Sample Autocorrelation:
      + The most interesting feature of this plot is where it approaches zero. This provides an estimate of the maximum size of the amplicons in a library.
      
2. *Where can I find more data on my samples not included in the HTML report (e.g. chromosomal amplicon properties)?*

    In the output directory, PaSD-qc creates a sub-directory `data`. With this directory you will find
    + <name>.table.txt: a tab-deliminated version of the table in the html report
    + <name>.fit.txt: a table of the parameters fit when estimating the amplicon distribution for each sample
    + chrom: directory containing a file for each sample summarizing that sample's chromosome properties
        + <name>_<sample_name>_chromProps.txt: each row is a chromosome; the columns include: KL-divergence, sub-amplicon variance, supra-amplicon variance, min_pos, median amplicon size, mean amplicon size, upper amplicon size, lower amplicon size, and classification (pass, possible gain, possible loss, fail).            

3. *How do I know if my sample is high or low quality?*

    The answer to this question depends on the application, so unfortunately no single answer exists. In general, however, the more closely a single-cell PSD resembles a bulk PSD, the higher the quality. That is, total variance should be small; sub-amplicon variance should approach that of bulk sequencing; and the lower the supra-amplicon variance, the better. Additionally, the smaller the variance of the amplicon size distribution, generally the higher the quality of the sample. Most importantly, the PSD of all samples should be highly similar. Any sample with a noticeably different PSD should be discarded as low quality.
    
   Additionally, PaSD-qc can compare the behavior of new samples to previously calculated gold-standard spectra. PaSD-qc is distributed with generic gold standard-spectra for depths corresponding to 1x (`db/categorical_spectra_1x.txt`) and 0.1x (`db/categorical_spectra_0.1x.txt`). However, the shape of high-quality PSDs is highly dependent on amplification kit and protocol, so we recommend that you generate your own gold-standard spectra based on your own samples which you believe to be high quality. Instructions for how to do this can be found in `example/08_example_gold_standard_spectra`.
   
   If you need help interpretting the results of your analysis, please feel free to email us (contact info at bottom of page).
   
4. *Can I use PaSD-qc to evaluate whether sub-chromosomal regions are poorly amplified?*

    Yes. See `examples/07_example_region_blacklist` for instructions.
    
5. *What should I supply as a bulk sample with the `-u` flag?*

    The `-u` flag is optional. By default, PaSD-qc will use a generic bulk PSD to normalize your single-cell samples. This will work fine for general purpose evaluation of samples. However, if you performed bulk sequencing at the same time as you sequenced the single cells, you may get more accurate estimates by using your own bulk sample. To do this, simply run PaSD-qc on the bulk sample and then supply the path to the resulting PSD file using the `-u` flag when profiling the single cells. E.g 
    ```
    PaSDqc QC -u /path/to/your/<bulk_name>.chroms.spec ...
    ```

6. *I deleted the html report. Do I have to rerun the whole pipeline?*

   No. The steps of PaSD-qc can be run separately. To regenerate the report, simply run `./PaSDqc.py report -d /path/to/analysis/dir -c categorical_spectra_XX.txt`. Regenerating a report generally takes less than 1 minute.

   To see all functions of PaSD-qc, run `./PaSD-qc.py -h`.

7. *Do you have some example data I can try out on?*

   Yes. The source code is distributed with an `examples` directory, which contains comprehensive instructions on basic and advanced usage of PaSD-qc. The tutorials are distributed as Jupyter notebooks, which can run once PaSDqc has been installed. To check that your installation is functional, you can re-compile the report in `01_example_report` navigating to the cloned repository and running
   ```
   PaSD-qc.py report -d examples/01_example_report -c categorical_spectra_0.1x.txt
   ```
   This will recreate the html report in the directory `examples/01_examples_report`. If you want to start with BAMs, the samples used in our [paper](https://doi.org/10.1093/nar/gkx1195) are available on SRA. See the Availability section of the paper for accession numbers.

8. *Doesn't plotly upload its plots to the interwebs? Can anyone see my PaSD-qc plots??*

   Rest easy, none of your plots are uploaded to the internet. We use the plotly javascript API to embed the plots directly into the html report. This means that a report is entirely self-contained, so you can send it around to collaborators without fear of plots not loading. It also means that the reports can be fairly large (several MBs), as the data required to generate the plots must also be embedded in the html file.

9. *Can I used PaSD-qc to evaluate my single-cell whole-exome sequencing?*

    In theory, yes, the approach should work on whole-exome sequencing. However, the method to perform the power spectral density estimation should probably be adapted to use information about the targeted capture regions. We are eager to do this. If you have high quality single-cell WES data and are interested in PaSDqc analysis, please contact us. We would be happy to work with you to ensure the results are high-quality. 

10. *I want to use PaSD-qc to analyze the quality of samples aligned to a different genome build. Is this possible*
    
    Yes, but you'll need to a bed file of uniquely mappable positions for the genome. Contact us for help with this.

## Support
__Maxwell Sherman__: maxwell\_sherman {at} hms.harvard.edu

__Peter Park__: peter\_park {at} hms.harvard.edu
