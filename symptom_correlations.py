# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 10:42:36 2017

@author: hhaque
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FixedLocator
import matplotlib.lines as mlines
import math
import scipy
import scipy.stats
import scipy.special as sc


def quadratic_correlation(x_list, y_list):
    
    fit1 = np.polyfit(x_list, y_list ,1)
    val1 = np.polyval(fit1, x_list)
    
    #quadratic
    fit2 = np.polyfit(x_list, y_list ,2)
    val2 = np.polyval(fit2, x_list)
    
    lin_res = y_list - val1
    
    #partial quadratic
    fit2_part = np.polyfit(x_list, lin_res, 2)
    val2_part = np.polyval(fit2_part, x_list)

    #r-squared
    ybar = np.mean(lin_res)
    ssreg = np.sum((val2_part - ybar)**2)
    sstot = np.sum((lin_res - ybar)**2)
    r2 = ssreg / sstot

    n = len(x_list)
    nt = 3
    dof = len(x_list) - nt

    #get parameter std
    stdind = math.sqrt(np.cov(x_list)*1)
    indepvar_s = np.array(x_list) * (1 / stdind)

    M = np.ones((len(x_list), nt))
    scalefact = np.ones((1, nt))
    modelterms = [2, 1, 0]
    for i in range(nt):
        M[:,i] = M[:,i] * indepvar_s ** modelterms[i]
        scalefact[:,i] = scalefact[:,i] / (stdind ** modelterms[i])
        
    q, r = np.linalg.qr(M)
    coeff = scipy.linalg.solve(r, np.matmul(np.transpose(q), y_list))
    yhat = np.matmul(M, coeff)

    s = np.linalg.norm(y_list - yhat)
    r_inv = scipy.linalg.solve(r, np.identity(nt))
    var = s**2*(np.sum(r_inv**2, 1)) / (n-nt)

    parameter_var = var * (scalefact ** 2)
    parameter_std = np.sqrt(parameter_var)

    t = fit2_part / parameter_std
    p = sc.betainc(dof/2, 0.5, dof/(t**2 + dof))
    
    return r2, p[0][0], lin_res, val2_part, val2


def get_n_of_patient_group(x_list):
    patient_list = x_list[:25]
    patient_list = np.array([x for x in patient_list if str(x) != 'nan'])
    patient_n = len(patient_list)
    
    return patient_n


def test_within_groups(df, measures, stat_test, count, nb_rows, nb_cols, colors, subtitle=''):
    x_list = list(df[measures[0]])
    y_list = list(df[measures[1]])
    
    for i in range(len(x_list)):
        if str(x_list[i]) == 'nan' or str(y_list[i]) == 'nan':
            x_list[i] = np.nan
            y_list[i] = np.nan
    
    patient_n = get_n_of_patient_group(x_list)
    
    x_list = np.array([x for x in x_list if str(x) != 'nan'])
    y_list = np.array([x for x in y_list if str(x) != 'nan'])
    
    if len(x_list) != len(y_list):
        print('sample sizes do not match, sth went wrong')
    else:
        if stat_test == 'Pearson':
            stat, p = scipy.stats.pearsonr(x_list, y_list)
        if stat_test == 'Spearman':
            stat, p = scipy.stats.spearmanr(x_list, y_list)
        if stat_test == 'Quadratic':
            stat, p, lin_res, val2_part, val2 = quadratic_correlation(x_list, y_list)
            
    if p < 0.05:
        color = 'red'
    else:
        color = 'black'
    
    ax = fig.add_subplot(nb_rows, nb_cols, count)
    
    if stat_test == 'Quadratic':
        
        # ax.scatter(x_list[:patient_n], lin_res[:patient_n], color = colors[0])
        # ax.scatter(x_list[patient_n:], lin_res[patient_n:], color = colors[1])
        
        ax.scatter(x_list[:patient_n], y_list[:patient_n], color = colors[0], s = 8)
        ax.scatter(x_list[patient_n:], y_list[patient_n:], color = colors[1], s = 8)
        
        ax.tick_params(which='both', top=False, right=False)
        ax.tick_params(which='major', length = 4, direction='out', labelsize = 6, width=0.6)
        ax.set_ylim([0.54, 0.75])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        plt.xlabel(measures[0].upper(), size = 7)
        plt.ylabel('DFA', size = 7)
        
        textstr = 'R\u00b2={}\np={}'.format(round(stat, 3), round(p, 3))
        ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=6, verticalalignment='top', color = color)
        
        #plt.suptitle(super_title)
        
        # x_list_sorted, val2_part_sorted = zip(*sorted(zip(x_list, val2_part)))
        # plt.plot(x_list_sorted, val2_part_sorted, color = 'black')
        
        x_list_sorted, val2_sorted = zip(*sorted(zip(x_list, val2)))
        plt.plot(x_list_sorted, val2_sorted, color = 'black', lw = 1)
        
    else:
        plt.scatter(x_list, y_list, color = 'black')
        plt.xlabel(measures[0])
        plt.ylabel(measures[1])
        plt.title('rho={}, p={}'.format(round(stat, 3), round(p, 3)), color = color)
        #plt.suptitle(super_title)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        a, b = best_fit(x_list, y_list)
        yfit = [a + b * xi for xi in x_list]
        plt.plot(x_list, yfit, color = 'gray')
        
    
# solve for a and b
def best_fit(X, Y):
    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)
    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    b = numer / denum
    a = ybar - b * xbar
    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))
    return a, b


master_path = 'main_directory\\'

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

fig = plt.figure(figsize = (7, 5))
colors = [sns.color_palette('muted')[1], sns.color_palette('muted')[0]]

groups = ['adhd', 'control']
tasks = ['TSDT', 'GoNoGo']
measures = ['ASRS', 'BIS']
conditions = ['rest', 'task']

y_limits = [-0.8, 0.09]


threshold = 0.025

count = 1
for idg, group in enumerate(groups):
    for task in tasks:
        for measure in measures:
            for condition in conditions:
                pos = '{}//{}//{}_{}_{}_pos.xlsx'.format(master_path, task, measure, group, condition)
                neg = '{}//{}//{}_{}_{}_neg.xlsx'.format(master_path, task, measure, group, condition)
        
                if condition == 'rest':
                    ax = fig.add_subplot(3, 4, count)
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
                    if group == 'control':
                        ax.set_ylim([-0.09, 0.9])
                    #ax.yaxis.set_major_locator(FixedLocator([-0.06, -0.04, -0.02, 0.00, 0.02, 0.04, 0.06]))
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    
            
                    plt.title(task, size = 8, color = 'white')
                    #if task == 'GoNoGo':
                        #plt.title("Go/NoGo", size = 8)
                    
                    #plt.suptitle(task, size = 14)
                    #fig.subplots_adjust(hspace= .5)
                    
                    y_label = 'P r(DFA, {})'.format(measure)
                    plt.ylabel(y_label, size = 7)
                    plt.xlabel('Frequency (Hz)', size = 7)
                    
                  
                for excel_file in ['pos', 'neg']: 
                    df = pd.read_excel(eval(excel_file))
                    df.columns = ['Frequency','Property']
                    if excel_file == 'neg':
                        df.Property = -df.Property
                    
                    if condition == 'rest':
                        ax.plot(df.Frequency, df.Property, lw = 1, linestyle = '--', color = colors[idg])
                    if condition == 'task':
                        ax.plot(df.Frequency, df.Property, lw = 1, color = colors[idg])
                
                if condition == 'task':
                    count += 1
                    
            plt.fill_between(df.Frequency, -threshold, threshold, color = 'gray', alpha = 0.25, linewidth = 0.01)
            
            if count == 2 or count == 6:
                solid_line = mlines.Line2D([], [], color=colors[idg], lw = 1, label='Task')
                dashed_line = mlines.Line2D([], [], color=colors[idg], lw = 1, linestyle = '--', label='Rest')
                plt.legend(handles=[dashed_line, solid_line], frameon = False, prop={'size': 6})
                

'''plotting the quadratic at the bottom'''
        

fpath = 'main_directory\\'
all_df = pd.read_csv('{}\\all_subjects_v2.csv'.format(fpath), sep = ';')

all_df = all_df.drop(columns=['Subject'])
new_cols = ['tsdt_hr', 'tsdt_rt', 'tsdt_dfa', 'gng_hr', 'gng_rt', 'gng_dfa', 'gng_fa', 'gng_rtv', 'asrs', 'bis', 'badds',
            'tsdt_task_dfa_8Hz', 'tsdt_task_dfa_16Hz', 'tsdt_rest_dfa_8Hz', 'tsdt_rest_dfa_16Hz', 'gng_task_dfa_8Hz', 'gng_task_dfa_16Hz',
            'gng_rest_dfa_8Hz', 'gng_rest_dfa_16Hz']
all_df.columns = new_cols

# group = 'All subjects'
# super_title = '{} - TSDT neuronal vs Clinical'.format(group)
# rows = ['asrs', 'bis']
# cols = ['tsdt_task_dfa_8Hz', 'tsdt_task_dfa_16Hz']
tasks = ['tsdt', 'gng']
stat_test = 'Quadratic'


for task in tasks:
    for measure in measures:
        test_within_groups(all_df, [measure.lower(), '{}_task_dfa_16Hz'.format(task)], stat_test, count, 3, 4, colors)
        count += 1        
            
    
plt.tight_layout()

