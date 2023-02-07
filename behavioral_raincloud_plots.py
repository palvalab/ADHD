# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 15:37:31 2022

@author: hhaque
"""


import pandas as pd
import seaborn as sns
import PtitPrince
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FixedLocator
from scipy import stats
from decimal import Decimal


def mannwhitney_between_groups(subset_df, task = ''):
    if task != '':
        cat1 = subset_df[(subset_df['group']=='ADHD') & (subset_df['variable']==task)]
        cat2 = subset_df[(subset_df['group']=='Ctrl') & (subset_df['variable']==task)]
    else:
        cat1 = subset_df[subset_df['group']=='ADHD']
        cat2 = subset_df[subset_df['group']=='Ctrl']
    
    adhd_list = list(cat1['value'])
    adhd_list = [x for x in adhd_list if str(x) != 'nan']
    control_list = list(cat2['value'])
    control_list = [x for x in control_list if str(x) != 'nan']
    
    stat, prob = stats.mannwhitneyu(adhd_list, control_list, alternative = 'two-sided')
    
    return stat, prob
    

pd.options.mode.chained_assignment = None  # default='warn'

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

fpath = 'main_directory\\'
all_df = pd.read_csv('{}\\all_subjects_v2.csv'.format(fpath), sep = ';')


new_cols = ['Subject', 'tsdt_hr', 'tsdt_rt', 'tsdt_dfa', 'gng_hr', 'gng_rt', 'gng_dfa', 'gng_fa', 'gng_rtv', 'asrs', 'bis', 'badds',
            'tsdt_task_dfa_8Hz', 'tsdt_task_dfa_16Hz', 'tsdt_rest_dfa_8Hz', 'tsdt_rest_dfa_16Hz', 'gng_task_dfa_8Hz', 'gng_task_dfa_16Hz',
            'gng_rest_dfa_8Hz', 'gng_rest_dfa_16Hz']
all_df.columns = new_cols

patient_df = all_df.iloc[:25,:]
control_df = all_df.iloc[25:,:]

patient_df = pd.melt(patient_df, id_vars='Subject', value_vars=new_cols[1:])
control_df = pd.melt(control_df, id_vars='Subject', value_vars=new_cols[1:])

patient_df['group'] = 'ADHD'
control_df['group'] = 'Ctrl'

all_df = pd.concat([patient_df, control_df], ignore_index=True, axis=0)

fig = plt.figure(figsize=(6, 3))

colors = [sns.color_palette('muted')[1], sns.color_palette('muted')[0]]


#create subset df
all_measures = [['bis'], ['tsdt_hr', 'gng_hr'], ['gng_fa'], ['gng_fa'], ['asrs'], ['tsdt_rt', 'gng_rt'], ['gng_rtv'], ['tsdt_dfa', 'gng_dfa']]
ylabels = ['BIS', 'Hit Rate (%)', 'False Alarm', 'False Alarm', 'ASRS', 'Reaction Time (s)', 'Reaction Time Variability', 'Behavioral DFA']
y_ticks = [[], [40, 60, 80, 100], [], [], [], [0.2, 0.4, 0.6], [0.04, 0.08, 0.12], []]


for idx, measures in enumerate(all_measures):
    subset_df = all_df.loc[all_df['variable'].isin(measures)]
    if len(measures) == 2:
        subset_df.replace(measures[0], 'TSDT', inplace=True)
        subset_df.replace(measures[1], 'Go/NoGo', inplace=True)
    if len(measures) == 1:
        if measures[0][:3] == 'gng':
            subset_df.replace(measures[0], 'Go/NoGo', inplace=True)
        else:
            subset_df.replace(measures[0], '', inplace=True)

    ax = fig.add_subplot(2, 4, idx+1)
    PtitPrince.RainCloud(x = 'variable', y = 'value', hue = 'group', data = subset_df, alpha = .65, dodge = True, 
                         palette = colors, width_viol = .55, box_showfliers = False,
                         width_box = 0.3, box_linewidth = 1, box_whiskerprops = {'linewidth':1}, point_size = 2)
    
    ax.tick_params(which='both', top=False, right=False)
    ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if y_ticks[idx] != []:
        ax.yaxis.set_major_locator(FixedLocator(y_ticks[idx]))
    
    ax.get_legend().remove()
    #plt.legend(frameon = False, prop={'size': 6})
    
    plt.ylabel(ylabels[idx], size = 7)
    plt.xlabel('')
    
    '''statistical tests'''
    if len(measures) == 2:
        stat, prob = mannwhitney_between_groups(subset_df, 'TSDT')
        y = plt.ylim()[1]
        if prob < 0.05:
            x1, x2 = -0.05, 0.05
            col = 'k'
            h = y*0.01
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y+h, '*', ha='center', va='bottom', color=col, size = 12)
        
        stat, prob = mannwhitney_between_groups(subset_df, 'Go/NoGo')
        if prob < 0.05:
            x1, x2 = 0.95, 1.05
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y+h, '*', ha='center', va='bottom', color=col, size = 12)
        
    if len(measures) == 1:
        stat, prob = mannwhitney_between_groups(subset_df)
        if prob < 0.05:
            x1, x2 = -0.05, 0.05
            col = 'k'
            y = plt.ylim()[1]
            h = y*0.01
            plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y+h, '*', ha='center', va='bottom', color=col, size = 12)
    
plt.tight_layout()
