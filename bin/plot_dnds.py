#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Plot the codon usage of various fasta files

This reads an aggregated results.csv file and outputs a single
interactive bar chart (html format, using plotly JS library) that
shows the codon usage patterns in each of the sampled simulated
lineages.

Each row in the csv file represents one sampled lineage.  The dN/dS
ratios is extracted from the row along with the path to the FASTA
sample.  The fasta file is parsed and the distribution of codons is
calculated.  Each sample contributes one set of vertical bars to the
chart.  Samples with dN/dS ratio >900 (indicating failure) are grouped
apart from those with more "reasonable" dN/dS ratios.

You can see an example of the output from this script at 
https://cswarth.github.io/santa-dnds/codon_usage.html

Usage:
    fill in typical command-line usage

'''
from __future__ import print_function

import os
import sys
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

import plotly
from plotly import tools
import plotly.graph_objs as go


# read in results.csv
# plot scater plot xasix is model, y axis is dnds ratio
# Color coded by fitness strentgh


df = pd.read_csv("build/results.csv")
df = df.dropna(axis=0, how='any')
colors = {'0.01':'blue', '0.01025':'red', '0.0105':'green'}

def rand_jitter(arr):
    stdev = .01*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

models = df.model.unique().tolist()
df['model_index'] = rand_jitter(df.model.apply(lambda f: models.index(f)))

# Group by model and add a fitness label.  Within each group there are
# three fitness values, but the set of values if different in each
# group.  Smooth over those differences by assigning labels to the
# fitness values.  After this transform all rows will have one of
# three levels of fitness, ['low', 'medium', 'high'], defined
# according to the values within that group.
#
# it seems like the new categorical values support would be perfect for this....

def mk_fitness_index(df):
    fitness = df.fitness.unique().tolist()
    fitness_label = ['low', 'medium', 'high']
    df['fitness_index'] = df.fitness.apply(lambda f: fitness.index(f))
    df['fitness_label'] = df.fitness_index.apply(lambda i: fitness_label[i])
    return df

df = df.groupby('model').apply(mk_fitness_index)


def hovertext(row):
    return """
mut:{mut:.2E}
gen:{generation}
fit:{fitness}
rep:{replicate}
""".format(**row)

def mktrace(df):
    return go.Scatter(
        x = df.model.tolist(),
        y = df.omega.tolist(),
        mode='markers',
        text=df.apply(hovertext, axis=1).tolist(),
        marker=dict(
            size='14',
            color = df.fitness_index.tolist(),
            colorscale='Viridis'
            ),
            showlegend=False
        )

# create traces, one for each mutation rate in the table.  We will plot these traces next to each other in the same figure.

fig = tools.make_subplots(rows=1, cols=2, subplot_titles=map(lambda m: "mutation rate {:.2E}".format(m), df.mut.unique()))
traces = df.groupby('mut').apply(mktrace).tolist()
for i,t in enumerate(traces):
    fig.append_trace(t, 1, i+1)

# append three more pseudo-traces to fill out the legend
for i,f in enumerate(['low', 'medium', 'high']):
    trace = go.Scatter(
            name=f,
            x = ['noselection'],
            y = [1],
            mode='markers',
            marker=dict(
                size='14',
                color = i+1,
                colorscale='Viridis'
            ),
            showlegend=True
        )
    fig.append_trace(trace, 1, 1)


fig['layout'].update(title='dN/dS ratio by fitness model',
                         autosize=False,
                         width=800,
                         height=800,
                         margin=go.Margin( l=50, r=50, b=150, t=50, pad=4),
                         # paper_bgcolor='#7f7f7f',
                         # plot_bgcolor='#c7c7c7',
                         
                         yaxis=dict(
                             type='log',
                             autorange=True
                             ),
                         hovermode='closest')


url = plotly.offline.plot(fig, show_link=False, auto_open=False)

print(url)
