# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 10:42:36 2017

@author: hhaque
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FixedLocator
import matplotlib.lines as mlines

master_path = 'main_directory\\'

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

fig = plt.figure(figsize = (3.7, 1.7))
colors = sns.color_palette('muted')


# y_label = 'P Ctrl > ADHD'
# test = 'Welch'
# tasks = ['TSDT', 'GoNoGo']
# conditions = ['rest', 'task']
# y_limits = [-0.11, 0.04]
# svg_fname = 'Fig_2_Welch'

y_label = 'P Task > Rest'
test = 'Wilcoxon'
tasks = ['TSDT', 'GoNoGo']
conditions = ['control', 'adhd']
y_limits = [-1, 1]
svg_fname = 'Fig_3_Wilcoxon'

# y_label = 'P r(Neural DFA, behav DFA)'
# test = 'Spearman'
# tasks = ['TSDT', 'GoNoGo']
# conditions = ['control', 'adhd']
# y_limits = [-0.5, 0.5]
# svg_fname = 'Fig_4_Spearman'


threshold = 0.025

count = 1
for task in tasks:
    for condition in conditions:
        pos = '{}//{}//{}_{}_pos.xlsx'.format(master_path, task, test, condition)
        neg = '{}//{}//{}_{}_neg.xlsx'.format(master_path, task, test, condition)

        if condition == 'control' or condition == 'rest':
            ax = fig.add_subplot(1, 2, count)
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
            if y_label == 'P Ctrl > ADHD':
                ax.yaxis.set_major_locator(FixedLocator([-0.1, -0.05, 0.00]))
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            
    
            plt.title(task, size = 8)
            if task == 'GoNoGo':
                plt.title("Go/NoGo", size = 8)
            
            #plt.suptitle(task, size = 14)
            #fig.subplots_adjust(hspace= .5)
            
            if count in [1]:
                plt.ylabel(y_label, size = 7)
            plt.xlabel('Frequency (Hz)', size = 7)
            
          
        for excel_file in ['pos', 'neg']: 
            df = pd.read_excel(eval(excel_file))
            df.columns = ['Frequency','Property']
            if excel_file == 'neg':
                df.Property = -df.Property
            
            if condition == 'control' or condition == 'rest':
                ax.plot(df.Frequency, df.Property, lw = 1, color = colors[0])
            if condition == 'adhd' or condition == 'task':
                ax.plot(df.Frequency, df.Property, lw = 1, color = colors[1])
        
        if condition == 'adhd' or condition == 'task':
            count += 1
            
    plt.fill_between(df.Frequency, -threshold, threshold, color = 'gray', alpha = 0.25, linewidth = 0.01)
    
    if y_label == 'P Task > Rest':
        if task == 'TSDT':
            plt.plot([6.62, 19.7], [1.0, 1.0], color = 'black', lw = 1)
        if task == 'GoNoGo':
            plt.plot([31.8, 120], [1.0, 1.0], color = 'black', lw = 1)
            
    if count == 3:
        dashed_line = mlines.Line2D([], [], color=colors[0], lw = 1, label='Ctrl')
        solid_line = mlines.Line2D([], [], color=colors[1], lw = 1, label='ADHD')
        plt.legend(handles=[dashed_line, solid_line], frameon = False, prop={'size': 6})
        
        if y_label == 'P Ctrl > ADHD':
            dashed_line = mlines.Line2D([], [], color=colors[0], lw = 1, label='Rest')
            solid_line = mlines.Line2D([], [], color=colors[1], lw = 1, label='Task')
            plt.legend(handles=[dashed_line, solid_line], frameon = False, prop={'size': 6})
    
plt.tight_layout()
