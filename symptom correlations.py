# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 10:42:36 2017

@author: hhaque
"""


import quadfun

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FixedLocator
import matplotlib.lines as mlines
import scipy.stats


def compute_sem(row):
    return scipy.stats.sem(row, nan_policy='omit')


def get_CI_limits(df):
    confidence = 0.95
    
    jackknife_col_names = list(df.columns[2:])
    df_jackknife = df[jackknife_col_names].copy()
    
    jk_n = len(df_jackknife.columns)
    jk_sem = df_jackknife.apply(compute_sem, axis=1)
    h = jk_sem * scipy.stats.t.ppf((1 + confidence) / 2., jk_n-1)
    
    df_jackknife['CI_max'] = df.Property + h
    df_jackknife['CI_min'] = df.Property - h
    
    return df_jackknife


def plot_p_plots(stat_type, master_path, fig, tasks, measures, conditions, y_limits, threshold, count, colors):
    for task in tasks:
        for measure in measures:
            for condition in conditions:
                pos = '{}//{}_python//jackknife//{}_{}_{}_jackknife_pos.csv'.format(master_path, task, measure, stat_type, condition)
                neg = '{}//{}_python//jackknife//{}_{}_{}_jackknife_neg.csv'.format(master_path, task, measure, stat_type, condition)
        
                if condition == 'task':
                    ax = fig.add_subplot(4, 4, count)
                    ax.tick_params(which='both', top=False, right=False)
                    ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
                    ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
                    ax.set_xscale('log')
                    ax.set_xticks([3, 5, 10, 20, 30, 50, 100])
                    ax.set_xlim([3, 120])
                    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
                    ax.xaxis.set_minor_locator(FixedLocator([4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 60, 70, 80, 90]))
                    ax.xaxis.set_minor_formatter(plt.NullFormatter())
                    ax.set_ylim(y_limits)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    
                    plt.title(task, size = 8, color = 'white')
                    
                    y_label = 'P r(DFA, {})'.format(measure)
                    plt.ylabel(y_label, size = 7)
                    plt.xlabel('Frequency (Hz)', size = 7)
                  
                for excel_file in ['pos', 'neg']: 
                    df = pd.read_csv(eval(excel_file), sep = ';')
                    if excel_file == 'neg':
                        df.Property = -df.Property
                    
                    df_jackknife = get_CI_limits(df)
                    if condition == 'rest':
                        ax.plot(df.Frequency, df.Property, lw = 1, color = colors[2])
                        plt.fill_between(df.Frequency, df_jackknife.CI_min, df_jackknife.CI_max, color = colors[2], alpha = 0.25, linewidth = 0.01)
                    if condition == 'task':
                        ax.plot(df.Frequency, df.Property, color = 'black', lw = 1)
                        plt.fill_between(df.Frequency, df_jackknife.CI_min, df_jackknife.CI_max, color = 'black', alpha = 0.25, linewidth = 0.01)
                
                if condition == 'task':
                    count += 1
                    
            plt.fill_between(df.Frequency, -threshold, threshold, color = 'gray', alpha = 0.25, linewidth = 0.01)
            plt.plot(df.Frequency, [0]*len(df.Frequency), color = 'black', linewidth = 0.6)
            
            if count == 2 or count == 10:
                solid_line = mlines.Line2D([], [], lw = 1, color = colors[3], label='Task')
                dashed_line = mlines.Line2D([], [], lw = 1, color = colors[2], label='Rest')
                #plt.legend(handles=[dashed_line, solid_line], frameon = False, prop={'size': 6})


master_path = 'main_directory\\'

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

fig = plt.figure(figsize = (7, 6.5))
colors = [sns.color_palette('muted')[1], sns.color_palette('muted')[0], sns.color_palette('muted')[2], sns.color_palette('muted')[3]]

tasks = ['TSDT', 'GoNoGo']
measures = ['ASRS', 'BIS']
conditions = ['task']

threshold = 0.025

fpath = 'directory\\'
all_df = pd.read_csv('{}\\all_subjects.csv'.format(fpath), sep = ';')

tsdt_exclude, gonogo_exclude = quadfun.exclude_subjects(all_df)

all_df = all_df[['Subject', 'ASRS', 'BIS']]

fpath_dfa = 'directory\\'


'''plotting the linear p plots'''
count = 1
y_limits = [0.0, 0.95]
plot_p_plots('linear', master_path, fig, tasks, measures, conditions, y_limits, threshold, count, colors)


'''plotting the linear scatters'''
count = 5     
stat_test = 'Spearman'
text_coords = [0.85, 0.95]

for task in tasks:
    df_dfa = pd.read_csv('{}\\{}_DFA_46.90_to_52.10.csv'.format(fpath_dfa, task), sep = ';')
    if task == 'TSDT':
        merged_df = pd.merge(all_df, df_dfa[['Subject', task]], on='Subject', how='left')
    if task == 'GoNoGo':
        merged_df = pd.merge(merged_df, df_dfa[['Subject', task]], on='Subject', how='left')
    
for task in tasks:
    if task == 'TSDT':
        cutoff_df = merged_df[~merged_df['Subject'].isin(tsdt_exclude)]
    if task == 'GoNoGo':
        cutoff_df = merged_df[~merged_df['Subject'].isin(gonogo_exclude)]
    
    for measure in measures:
        quadfun.test_within_groups(cutoff_df, [measure, task], stat_test, fig, count, 4, 4, colors, text_coords)
        count += 1   


'''plotting the quadratic p plots'''
y_limits = [-0.16, 0.01]
plot_p_plots('quad', master_path, fig, tasks, measures, conditions, y_limits, threshold, count, colors)


'''plotting the quadratic at the bottom'''
count = 13
stat_test = 'Quadratic'
text_coords = [0.75, 0.95]

for task in tasks:
    df_dfa = pd.read_csv('{}\\{}_DFA_21.60_to_23.73.csv'.format(fpath_dfa, task), sep = ';')
    if task == 'TSDT':
        merged_df = pd.merge(all_df, df_dfa[['Subject', task]], on='Subject', how='left')
    if task == 'GoNoGo':
        merged_df = pd.merge(merged_df, df_dfa[['Subject', task]], on='Subject', how='left')
    
for task in tasks:
    if task == 'TSDT':
        cutoff_df = merged_df[~merged_df['Subject'].isin(tsdt_exclude)]
    if task == 'GoNoGo':
        cutoff_df = merged_df[~merged_df['Subject'].isin(gonogo_exclude)]
    
    for measure in measures:
        quadfun.test_within_groups(cutoff_df, [measure, task], stat_test, fig, count, 4, 4, colors, text_coords)
        count += 1
        
    
plt.tight_layout()