# plotly_tools.py - methods for easy plotly plot creation
#
# v 0.0.8
# rev 2017-03-21 (MS: ACF and chrom plot methods)
# Notes:

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hc
import plotly.offline as py
import plotly.graph_objs as go
import plotly.figure_factory as ff

from . import extra_tools

def dendrogram(nd, sample_list):
 
    def _PSD_dist(nd):
        return hc.distance.pdist(nd, extra_tools.PSD_sym_KL)

    dend = ff.create_dendrogram(nd, labels=sample_list, distfun=_PSD_dist, orientation='left')
    dend['layout'].update({'width':1500, 'height':800, 
                             'font': dict(size=18), 
                             'margin': go.Margin(l = 450),
                            })
    dend['layout']['xaxis'].update({'title': 'KL Divergence'})

    # div_dend = py.plot(dend, output_type='div')

    return dend

def PSD_plot(freq, nd, sample_list):
    # PSD plot
    data = [go.Scatter(x=freq, y=10*np.log(row), mode='lines', text=s, name=s) for row, s in zip(nd, sample_list)] 
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

def chrom_KL_plot(j_list, sample_list):
    # xticks and point labels
    j = j_list[0]
    chroms = j.index.tolist()
    text = [c for c in chroms]
    chroms[-2:] = ['23', '24']
    chroms = pd.Series([int(c) for c in chroms], index=j.index)
    text = pd.Series(text, index=j.index)

    ticks = sorted(chroms)
    tick_labs = [t for t in ticks]
    tick_labs[-2:] = ['X', 'Y']
    tick_labs = [str(c) for c in tick_labs]

    data = []
    s_list = []

    for j, sample in zip(j_list, sample_list):
        d_tmp, s_tmp = _chrom_trace(j, sample, chroms, text)

        data.extend(d_tmp)
        s_list.extend(s_tmp)

    buttons = []
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
                               ticktext=tick_labs,
                               tickvals=ticks
                          ),
                  updatemenus = update_menus
             )

    fig = go.Figure(data=data, layout=layout)

    return fig

def _chrom_trace(j, sample, chroms, text):
    mu = j.median()
    mad = np.median(np.abs(j-mu))
    sd_1 = mu + mad
    sd_2 = mu + 2*mad

    trace0 = go.Scatter(x=chroms[j <= sd_1], 
                        y=j[j <= sd_1], 
                        mode='markers', 
                        marker=dict(color='blue', size=10), 
                        name='pass', 
                        text=text[j <= sd_1]
             )
    trace1 = go.Scatter(x=chroms[(j > sd_1) & (j <= sd_2)], 
                        y=j[(j > sd_1) & (j <= sd_2)],
                        mode='markers', 
                        marker=dict(color='orange', size=10), 
                        name='warn'
             )
    trace2 = go.Scatter(x=chroms[j > sd_2], 
                        y=j[j > sd_2], mode='markers', 
                        marker=dict(color='red', size=10), 
                        name='fail', 
                        text=text[j > sd_2]
                       )

    data = [trace0, trace1, trace2]
    s_list = [sample for trace in data]

    return data, s_list

def _chrom_mk_button(sample, s_list):
    vis = [True if s == sample else False for s in s_list]

    d = dict(args=['visible', vis],
             label=sample,
             method='restyle'
        )

    return d
