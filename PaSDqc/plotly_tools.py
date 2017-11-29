# plotly_tools.py - methods for easy plotly plot creation
#
# v 1.1.1
# rev 2017-11-29 (MS: minor bug fix)
# Notes:

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hc
import scipy.signal
import plotly.offline as py
import plotly.graph_objs as go
import plotly.figure_factory as ff

from . import extra_tools

def dendrogram(nd, sample_list, p_cat_spec):

    # load categorical spectra
    f_cat_spec = extra_tools.get_data_file(p_cat_spec)
    cat_spec = pd.read_table(f_cat_spec, index_col=0)
    nd_spec = cat_spec.as_matrix().T
    sample_spec = ['GOLD-{}'.format(c).upper() for c in cat_spec.columns.tolist()]

    nd = np.vstack((nd, nd_spec))
    sample_list.extend(sample_spec)

    def _PSD_dist(nd):
        return hc.distance.pdist(nd, extra_tools.PSD_sym_KL)

    dend = ff.create_dendrogram(nd, labels=sample_list, distfun=_PSD_dist, linkagefun=hc.centroid, orientation='left')
    dend['layout'].update({'width':1500, 'height':800, 
                             'font': dict(size=18), 
                             'margin': go.Margin(l = 450),
                            })
    dend['layout']['xaxis'].update({'title': 'KL Divergence'})

    # div_dend = py.plot(dend, output_type='div')

    return dend

def PSD_plot(freq, nd, sample_list, bulk='bulk_1x.smooth3.spec'):
    # PSD plot
    # f_bulk = extra_tools.get_data_file(bulk)
    # psd_bulk = pd.Series.from_csv(f_bulk, index_col=0, header=None, sep="\t").as_matrix()
    psd_bulk = extra_tools.load_bulk_psd(bulk)
    # psd_bulk = np.loadtxt(bulk)
    # psd_norm = [10*np.log10(row / psd_bulk) for row in nd]
    # psd_smooth = [scipy.signal.savgol_filter(psd, 101, 3) for psd in psd_norm]

    period = 1 / freq
    data = [go.Scatter(x=period, y=10*np.log10(row/psd_bulk), mode='lines', text=s, name=s) for row, s in zip(nd, sample_list)] 
    # data = [go.Scatter(x=freq, y=psd, mode='lines', text=s, name=s) for psd, s in zip(psd_smooth, sample_list)] 
    layout = go.Layout(
        hovermode='closest',
        width=1500, 
        height=800,
        # showlegend=False,
        font=dict(size=18),
        xaxis=dict(
            title='Inverse Genomic Range (1/BP)',
            type='log',
            autorange=True,
            exponentformat='power'
        ),
        yaxis=dict(
            title='Power Spectral Density (dB)',
            # type='log',
            autorange=True,
            exponentformat='power'
        )
    )
    fig = go.Figure(data=data, layout=layout)
    # div_psd = py.plot(fig, output_type='div')

    return fig

# def amplicon_density_plot(amp_list, pdf_list, sample_list):
def amplicon_density_plot(psd_list, sample_list):
    data = [go.Scatter(x=psd.sample_dist['freq'], y=psd.sample_dist['dist'], mode='lines', text=s, name=s) for psd, s in zip(psd_list, sample_list)]
    # data = [go.Scatter(x=amp.freq['dist'], y=pdf, mode='lines', text=s, name=s) for amp, pdf, s in zip(amp_list, pdf_list, sample_list)] 
    # data = [go.Scatter(x=distances, y=pdf, mode='lines', text=s, name=s) for pdf, s in zip(pdf_list, sample_list)] 
    layout = go.Layout(
        hovermode='closest',
        width=1500, 
        height=800,
        # showlegend=False,
        font=dict(size=18),
        xaxis=dict(
            title='Amplicon size',
            type='log',
            autorange=True,
            exponentformat='power'
        ),
        yaxis=dict(
            title='Density',
            # type='log',
            autorange=True,
            exponentformat='power'
        )
    )
    fig = go.Figure(data=data, layout=layout)
    # div_psd = py.plot(fig, output_type='div')

    return fig

def ACF_plot(lags, nd, sample_list):
    # PSD plot
    data = [go.Scatter(x=lags[1:], y=row[1:], mode='lines', text=s, name=s) for row, s in zip(nd, sample_list)] 
    layout = go.Layout(
        hovermode='closest',
        width=1500, 
        height=800,
        # showlegend=False,
        font=dict(size=18),
        xaxis=dict(
            title='Lag (BP)',
            type='log',
            autorange=True,
            exponentformat='power'
        ),
        yaxis=dict(
            title='Autocorrelation',
            # type='log',
            autorange=True,
            exponentformat='power'
        )
    )
    fig = go.Figure(data=data, layout=layout)
    # div_acf = py.plot(fig, output_type='div')

    return fig

# def chrom_KL_plot(j_list, sample_list):
def chrom_KL_plot(psd_list, sample_list):
    # xticks and point labels
    chroms = psd_list[0].chrom_props.index
    text = pd.Series(chroms.tolist(), index=chroms)
    ticks = list(range(1, len(chroms)+1))
    x_pos = pd.Series(ticks, index=chroms)

    # j = j_list[0]
    # chroms = j.index.tolist()
    # text = [c for c in chroms]
    # # chroms[-2:] = ['23', '24']
    # chroms = pd.Series([int(c) for c in chroms], index=j.index)
    # text = pd.Series(text, index=j.index)

    # ticks = sorted(chroms)
    # tick_labs = [t for t in ticks]
    # # tick_labs[-2:] = ['X', 'Y']
    # tick_labs = [str(c) for c in tick_labs]

    data = []
    s_list = []

    # for j, sample in zip(j_list, sample_list):
    for psd, sample in zip(psd_list, sample_list):
        d_tmp, s_tmp = _chrom_trace(psd, sample, x_pos, text)

        data.extend(d_tmp)
        s_list.extend(s_tmp)

    d = dict(args=['visible', [True for s in s_list]],
             label='ALL',
             method='restyle'
        )

    buttons = [d]
    for s in sample_list:
        buttons.append(_chrom_mk_button(s, s_list))

    update_menus = list([dict(x=0.05, 
                              y=1.15,
                              xanchor='left',
                              yanchor='top',
                              active=0,
                              buttons=list(buttons),
                       )])
    layout = dict(width=1500, 
                  height=800,
                  font=dict(size=18),
                  xaxis = dict(tickmode='array',
                               ticktext=text.tolist(),
                               tickvals=ticks
                               title='Chromosome'
                          ),
                  yaxis=dict(title='KL-divergence'),
                  updatemenus = update_menus
             )

    fig = go.Figure(data=data, layout=layout)

    return fig

# def _chrom_trace(j, sample, chroms, text):
def _chrom_trace(psd, sample, x_pos, text):
    # mu = j.median()
    # mad = np.median(np.abs(j-mu))
    # sd_1 = mu + mad
    # sd_2 = mu + 2*mad
    df = psd.chrom_props

    # Passing
    trace0 = go.Scatter(x=x_pos[df.classif == 'Pass'], 
                        y=df.KL[df.classif == 'Pass'], 
                        mode='markers', 
                        marker=dict(color='grey', size=10), 
                        name='pass', 
                        text=sample
                        # text=text[j <= sd_1]
             )
    # Gains
    trace1 = go.Scatter(x=x_pos[df.classif == 'Possible gain'], 
                        y=df.KL[df.classif == 'Possible gain'],
                        mode='markers', 
                        marker=dict(color='green', size=10), 
                        name='Gain?',
                        text=sample
             )
    # Losses
    trace2 = go.Scatter(x=x_pos[df.classif == 'Possible loss'],
                        y=df.KL[df.classif == 'Possible loss'], 
                        mode='markers', 
                        marker=dict(color='blue', size=10), 
                        name='Loss?', 
                        text=sample
                        # text=text[j > sd_2]
                       )
    # Fails
    trace3 = go.Scatter(x=x_pos[df.classif == 'Fail'],
                        y=df.KL[df.classif == 'Fail'], 
                        mode='markers', 
                        marker=dict(color='red', size=10), 
                        name='Fail', 
                        text=sample
                        # text=text[j > sd_2]
                       )

    # data = [trace0, trace1, trace2]

    data = [trace0, trace1, trace2, trace3]
    s_list = [sample for trace in data]

    return data, s_list

def _chrom_mk_button(sample, s_list):
    vis = [True if s == sample else False for s in s_list]

    d = dict(args=['visible', vis],
             label=sample,
             method='restyle'
        )

    return d
