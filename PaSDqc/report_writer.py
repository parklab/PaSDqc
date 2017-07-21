# report_writer.py - methods for writing html report
#
# v 0.1.0
# rev 2017-07-14 (MS: added amplicon size distribution to reports
# Notes:

def writer(df, div_dend, div_psd, div_acf, div_chrom, div_pdf, fout):
    df_html = df.to_html().replace('<table border="1" class="dataframe">','<table class="table table-striped">')

    html_string = report_html(df_html, div_dend, div_psd, div_acf, div_chrom, div_pdf)

    with open(fout, 'w') as f:
        f.write(html_string)

def report_html(df_html, div_dend, div_psd, div_acf, div_chrom, div_amp):
    html_string = '''
<html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
        <style>body{ margin:0 100; background:whitesmoke; }</style>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
        <h1>MDAqc sample report</h1>

        <!-- *** Section 1 *** --->
        <h2>Sample Summary Table</h2>
        ''' + df_html + '''
        <h2>Sample Clustering</h2>
        ''' + div_dend + '''
        <h2>Chromosome Outlier Plots</h2>
        ''' + div_chrom + '''
        <h2>Amplicon size distributions</h2>
        ''' + div_amp + '''
        <h2>Sample Periodograms</h2>
        ''' + div_psd + '''
        <h2>Sample Autocorrelation</h2>
        ''' + div_acf + '''
    </body>
</html>'''

    return html_string
