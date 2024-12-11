# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 08:04:15 2022

@author: hhaque
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FixedLocator


master_fpath = 'main_directory\\'
tasks = ['TSDT', 'GoNoGo']

fig = plt.figure(figsize = (3.7, 1.7))

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

for idx, task in enumerate(tasks):
    fname = 'DFA_mean_ci95_bootstrap_{}_2024_cutoff.csv'.format(task)
    fpath = '{}\\{}'.format(master_fpath, fname)
    
    df = pd.read_csv(fpath, sep = ',')
    df.columns = ['freq', 'group1_mean', 'group2_mean', 'group1_CI1', 'group2_CI1', 'group1_CI2', 'group2_CI2']
    
    colors = sns.color_palette('muted')
    
    ax = fig.add_subplot(1, 2, idx+1)
    
    ax.plot(df.freq, df.group1_mean, color = colors[1], lw = 1)
    plt.fill_between(df.freq, df.group1_CI1, df.group1_CI2, color = colors[1], alpha = 0.25, linewidth = 0.01)
    
    ax.plot(df.freq, df.group2_mean, color = colors[0], lw = 1)
    plt.fill_between(df.freq, df.group2_CI1, df.group2_CI2, color = colors[0], alpha = 0.25, linewidth = 0.01)
    
    ax.tick_params(which='both', top=False, right=False)
    ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
    ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
    ax.set_xscale('log')
    ax.set_xticks([3, 5, 10, 20, 30, 50, 100])
    ax.set_xlim([3, 120])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_minor_locator(FixedLocator([4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 60, 70, 80, 90]))
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(FixedLocator([0.60, 0.64, 0.68]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    plt.title(task, size = 8)
    if task == 'GoNoGo':
        plt.title("Go/NoGo", size = 8)
    plt.ylabel('DFA Neuronal', size = 7)
    plt.xlabel('Frequency (Hz)', size = 7)
    
plt.tight_layout()