# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:39:39 2024

@author: hhaque
"""

import numpy as np
import matplotlib.pyplot as plt
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
    
    r2 = r2 * np.sign(fit2_part[0]) #indicate pos or neg relationship

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


def get_n_of_patient_group(x_list, df):
    change_index = df[df['Subject'].str.startswith('S')].index[0]
    
    patient_list = x_list[:change_index]
    patient_list = np.array([x for x in patient_list if str(x) != 'nan'])
    patient_n = len(patient_list)
    
    return patient_n


def test_within_groups(df, measures, stat_test, fig, count, nb_rows, nb_cols, colors, text_coords):
    x_list = list(df[measures[0]])
    y_list = list(df[measures[1]])
    
    for i in range(len(x_list)):
        if str(x_list[i]) == 'nan' or str(y_list[i]) == 'nan':
            x_list[i] = np.nan
            y_list[i] = np.nan
    
    df.reset_index(drop=True, inplace=True)
    patient_n = get_n_of_patient_group(x_list, df)
    
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
        
    ax.scatter(x_list[:patient_n], y_list[:patient_n], color = colors[0], s = 8, label = 'ADHD')
    ax.scatter(x_list[patient_n:], y_list[patient_n:], color = colors[1], s = 8, label = 'NC')
    
    ax.tick_params(which='both', top=False, right=False)
    ax.tick_params(which='major', length = 4, direction='out', labelsize = 6, width=0.6)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    plt.xlabel(measures[0].upper(), size = 7)
    plt.ylabel('DFA', size = 7)
    
    textstr = 'R\u00b2={}\np={}'.format(round(stat, 3), round(p, 3))
    ax.text(text_coords[0], text_coords[1], textstr, transform=ax.transAxes, fontsize=6, verticalalignment='top', color = color)
    #plt.legend(prop = {'size': 6}, ncol = 2)
    
    if stat_test == 'Quadratic':
        x_list_sorted, val2_sorted = zip(*sorted(zip(x_list, val2)))
        plt.plot(x_list_sorted, val2_sorted, color = 'black', lw = 1)
        
    else:
        a, b = best_fit(x_list, y_list)
        yfit = [a + b * xi for xi in x_list]
        plt.plot(x_list, yfit, color = 'black', lw = 1)
        
    
# solve for a and b
def best_fit(X, Y):
    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)
    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    b = numer / denum
    a = ybar - b * xbar
    return a, b


def exclude_subjects(df):
    tsdt_exclude = df.loc[df['HR (tsdt)'] < 25.00, 'Subject'].to_list()
    gonogo_exclude = df.loc[df['HR (gng)'] < 65.00, 'Subject'].to_list()
    
    return tsdt_exclude, gonogo_exclude