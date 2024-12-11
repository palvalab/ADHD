# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:39:39 2024

@author: hhaque
"""

import quadfun
from nptdms import TdmsFile
import numpy as np
from scipy import stats
import pandas as pd
from numpy import genfromtxt
import math
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm


def read_tdms(fpath, TW):
    graphdata = TdmsFile.read(fpath)
    
    group = graphdata[TW]
    channel = group['Vertex Data (CSGreim)__real']
    dfa_values = channel[:]
    
    return dfa_values



def discard_false_discoveries(stats, alpha):
    FDD = round(10 ** -alpha, 2)
    nb_sig_remove = int(len(stats) * FDD)
    
    stats = stats.tolist()
    signs = [1 if x > 0 else -1 for x in stats]
    abs_stats = [abs(x) for x in stats]
    threshold_list = [1 if x > alpha else 0 for x in abs_stats]
    
    sig_data = [a * b * c for a, b, c in zip(signs, abs_stats, threshold_list)]
    
    only_sig_data = []
    for idx, abs_value in enumerate(abs_stats):
        if threshold_list[idx] == 1:
            sig = [abs_stats[idx], idx]
            only_sig_data.append(sig)
            
    only_sig_data_sorted = sorted(only_sig_data, key=lambda x: x[0])
    
    if len(only_sig_data_sorted) > nb_sig_remove:
        for i in range(nb_sig_remove):
            idx = only_sig_data_sorted[i][1]
            sig_data[idx] = 0.0
    else:
        sig_data = [0] * 200
        
    sig_mask = [1 if x != 0 else 0 for x in sig_data]
    sig_mask = np.array(sig_mask)
    
    return sig_mask

    

def morphing_vertex_data(dfa, subject, morphing_op_path):
    morphing_op_path = '{}\\{}.csv'.format(morphing_op_path, subject)
    morphing_op = genfromtxt(morphing_op_path, delimiter = ';')
    
    dfa_morphed = np.dot(dfa, morphing_op)
    
    return dfa_morphed



def collapse_conditions(graphdata_tdms_fpath, freq, conditions, graphfolder):
    dfas = []
    for condition in conditions:
        graphdata_tdms = TdmsFile.read(graphdata_tdms_fpath)
        graphmap_df = graphdata_tdms['Graph Map'].as_dataframe()
        df_row = graphmap_df.loc[(graphmap_df['Condition (STR)'] == condition) & (graphmap_df['Condition Set (STR)'].isin(graphfolder)) & (graphmap_df['f (Hz) (DBL)'] == freq)]
        
        if df_row.shape[0] != 1:
            print("More than one graph data fulfil the condition!")
        
        fpath_condition = df_row['Graph Data Path (STR)'].values[0]
        dfa_condition = read_tdms(fpath_condition, '0')
        dfas.append(dfa_condition)
        
    dfa = dfas[0] - sum(dfas[1:])
    
    return dfa
    


def statsbox(all_vertices, stats_test, pv_values = None):
    if stats_test == 'wilcoxon':
        t, prob = stats.wilcoxon(all_vertices)
        avg_value = np.mean(all_vertices)
    
    if stats_test == 'welch':
        group_a = []
        group_b = []
        
        for idx in range(len(all_vertices)):
            if pv_values[idx] == 1.0:
                group_a.append(all_vertices[idx])
                
            if pv_values[idx] == 0.0:
                group_b.append(all_vertices[idx])
    
        t, prob = stats.ttest_ind(group_a, group_b, equal_var = False)
        avg_value = np.mean(group_a) - np.mean(group_b)
        
    if stats_test == 'spearman':        
        t, prob = stats.spearmanr(all_vertices, pv_values)
        avg_value = t
        
    if stats_test == 'quadratic':
        t, prob, lin_res, val2_part, val2 = quadfun.quadratic_correlation(pv_values, all_vertices) #seems the order of the variables matter
        avg_value = t
        
    if stats_test == 'ancova':
        df = pd.DataFrame(pv_values, columns=['Group', 'Age'])
        df['DFA'] = all_vertices
        
        avg_value = np.mean(df['DFA'][df['Group'] == 1.0]) - np.mean(df['DFA'][df['Group'] == 0.0])
        
        model = sm.formula.ols('DFA ~ C(Group) + Age', data=df).fit()
        ancova_table = anova_lm(model, typ = 2)
        prob = ancova_table.loc['C(Group)', 'PR(>F)']
        
    return avg_value, prob



def gdm_python(subject_list, frequencies, master_folder, morphing_op_path, graphfolder, conditions, stats_test, pv_values = None):
    alpha = 1.3
    
    all_P_positives = [0] * len(frequencies)
    all_P_negatives = [0] * len(frequencies)
    all_avg_values = []
    
    for idx, freq in enumerate(frequencies):
        all_dfa = []
        for subject in subject_list:
            graphdata_tdms_fpath = '{}\\{}\\GraphData.tdms'.format(master_folder, subject)
            
            dfa = collapse_conditions(graphdata_tdms_fpath, freq, conditions, graphfolder)
            dfa = morphing_vertex_data(dfa, subject, morphing_op_path)
            all_dfa.append(dfa)
        
        morphed_parcel_nbs = 200
        avg_values = np.zeros(morphed_parcel_nbs,)
        stats_both = np.zeros(morphed_parcel_nbs,)
        
        for vertex in range(np.shape(stats_both)[0]):
            all_vertices = []
            for subject_dfa in all_dfa:
                all_vertices.append(subject_dfa[vertex])
                
            if pv_values is not None:
                avg_value, prob = statsbox(all_vertices, stats_test, pv_values)
            else:
                avg_value, prob = statsbox(all_vertices, stats_test)
                
            stats_both[vertex] = np.sign(avg_value) * -math.log10(prob)
            avg_values[vertex] = avg_value
                
        sig_mask = discard_false_discoveries(stats_both, alpha)
        avg_values = avg_values * sig_mask
        
        weee_positive = [max(0, num) for num in avg_values]
        weee_negative = [min(0, num) for num in avg_values]
        
        P_positive = np.shape(np.nonzero(weee_positive)[0])[0] / morphed_parcel_nbs
        P_negative = np.shape(np.nonzero(weee_negative)[0])[0] / morphed_parcel_nbs
        
        all_P_positives[idx] = P_positive
        all_P_negatives[idx] = P_negative
        
        all_avg_values.append(avg_values)
        print(freq)
        
    all_avg_values = np.vstack(all_avg_values)
    
    return all_avg_values, all_P_positives, all_P_negatives



def p_to_excel(p_list, output_folder, frequencies, analysis_string, sign):
    fname = '{}\\{}_{}.xlsx'.format(output_folder, analysis_string, sign)
    df = pd.DataFrame({'Frequency': frequencies, 'Property': p_list})
    df.to_excel(fname, index = False)
    
    
    
def avg_value_to_csv(all_avg_values, output_folder, analysis_string):
    fname = '{}\\{}_avg.csv'.format(output_folder, analysis_string)
    np.savetxt(fname, all_avg_values, delimiter = ';')


    
def jackknife_resampling(subject_list, frequencies, master_folder, morphing_op_path, graphfolder, conditions, stats_test, output_folder, analysis_string, pv_values = None):
    all_avg_values, all_P_positives, all_P_negatives = gdm_python(subject_list, frequencies, master_folder, morphing_op_path, graphfolder, conditions, stats_test, pv_values)
    
    df_pos = pd.DataFrame({'Frequency': frequencies, 'Property': all_P_positives})
    df_neg = pd.DataFrame({'Frequency': frequencies, 'Property': all_P_negatives})
    
    for i in range(len(subject_list)):
        jackknife_sample = list(subject_list)
        jackknife_sample.pop(i)
        print('---{}---'.format(subject_list[i]))
        
        if pv_values is not None:
            pv_values_jackknife = list(pv_values)
            pv_values_jackknife.pop(i)
            print('---{}---'.format(pv_values[i]))
        else:
            pv_values_jackknife = None
        
        all_avg_values, jackknife_pos, jackknife_neg = gdm_python(jackknife_sample, frequencies, master_folder, morphing_op_path, graphfolder, conditions, stats_test, pv_values_jackknife)
        
        column_header = 'Jackknife_{}'.format(i)
        df_pos[column_header] = jackknife_pos
        df_neg[column_header] = jackknife_neg

    fname = '{}\\{}_jackknife_pos.csv'.format(output_folder, analysis_string)
    df_pos.to_csv(fname, sep = ';', index = False)
    fname = '{}\\{}_jackknife_neg.csv'.format(output_folder, analysis_string)
    df_neg.to_csv(fname, sep = ';', index = False)
    