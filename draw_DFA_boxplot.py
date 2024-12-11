# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:33:47 2024

@author: hhaque
"""

source_directory    = 'K:\\palva\\Common repos\\OL2015\\source\\Python27\\VWM\\18 - Functions'

import pandas as pd
import seaborn as sns
import sys
sys.path.append(source_directory)
import PtitPrince
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FixedLocator
from scipy import stats
import numpy as np


def cohen_d(group1, group2):
    mean1, mean2 = np.mean(group1), np.mean(group2)
    std1, std2 = np.std(group1, ddof=1), np.std(group2, ddof=1)
    n1, n2 = len(group1), len(group2)
    pooled_std = np.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))
    d = (mean1 - mean2) / pooled_std
    
    return d


sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True


fpath_dfa = 'main_directory\\'

tasks = ['TSDT', 'GoNoGo']
frequency_pairs = [['8.14', '9.02'], ['59.30', '65.30']]
colors = [sns.color_palette('muted')[1], sns.color_palette('muted')[0]]

fig = plt.figure(figsize = (3.7, 1.3))

for idx, task in enumerate(tasks):
    
    low_f = frequency_pairs[idx][0]
    high_f = frequency_pairs[idx][1]
    
    df = pd.read_csv('{}\\{}_DFA_{}_to_{}.csv'.format(fpath_dfa, task, low_f, high_f), sep = ';')
    df.columns = ['Subject', 'DFA']
    
    if task == 'TSDT':
        df = df[df['Subject'] != 'S0120']
    if task == 'GoNoGo':
        df = df[df['Subject'] != 'P0036']
    
    # Create the 'Diagnosis' column
    df['Diagnosis'] = df['Subject'].apply(lambda x: 'NC' if x.startswith('S') else 'ADHD')
    df['variable'] = 'DFA'
    
    ax = fig.add_subplot(1, 2, idx+1)
    PtitPrince.RainCloud(x = 'variable', y = 'DFA', hue = 'Diagnosis', data = df, alpha = .65, dodge = True, 
                         palette = colors, width_viol = .55, box_showfliers = False,
                         width_box = 0.3, box_linewidth = 1, box_whiskerprops = {'linewidth':1}, point_size = 2)
    
    ax.tick_params(which='both', top=False, right=False)
    ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.get_legend().remove()   
    plt.ylabel('DFA', size = 7)
    plt.subplots_adjust(wspace=0.45)
    
    '''Statistical testing'''
    nc_group = df[df['Diagnosis'] == 'NC']['DFA']
    adhd_group = df[df['Diagnosis'] == 'ADHD']['DFA']
    
    t_stat, p_value = stats.ttest_ind(nc_group, adhd_group, equal_var = False)
    
    d = cohen_d(nc_group, adhd_group)
    print('The cohens d for {} is {}'.format(task, d))
    
    if p_value < 0.05:
        x1, x2 = -0.05, 0.05
        col = 'k'
        y = plt.ylim()[1]
        h = y*0.01
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
        plt.text((x1+x2)*.5, y+h, '*', ha='center', va='bottom', color=col, size = 12)