# plotly_tools.py - methods for easy plotly plot creation
#
# v 0.0.7
# rev 2017-03-19 (MS: report function creates html report)
# Notes:

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

    div_dend = py.plot(dend, output_type='div')

    return div_dend

def PSD_plot(freq, nd, sample_list):
    # PSD plot
    data = [go.Scatter(x=freq, y=row, mode='lines', text=s, name=s) for row, s in zip(nd, sample_list)] 
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
            type='log',
            autorange=True,
            exponentformat='power'
        )
    )
    fig = go.Figure(data=data, layout=layout)
    div_psd = py.plot(fig, output_type='div')

    return div_psd
