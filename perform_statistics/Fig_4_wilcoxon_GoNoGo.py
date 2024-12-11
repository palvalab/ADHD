# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:08:00 2024

@author: hhaque
"""

source_directory    = 'K:\\palva\\Common repos\\OL2015\\source\\Python27\\VWM\\18 - Functions'

import sys
sys.path.append(source_directory)
import dfa_core


frequencies = [3.31, 3.73, 4.15, 4.75, 5.35, 5.93, 6.62, 7.39, 8.14, 9.02, 9.83,
               10.92, 11.89, 13.11, 14.78, 16.3, 17.8, 19.7, 21.6, 23.73, 26.55,
               28.7, 31.8, 34.5, 37.9, 42.5, 46.9, 52.1, 59.3, 65.3, 71, 78, 85.92,
               95.6, 106, 118.5, 131, 143]

output_folder = 'main_directory\\GoNoGo_python\\jackknife'

master_folder = 'main_directory\\'
morphing_op_path = 'main_directory\\Morphing_operators\\GoNoGo'
graphfolder = ['spont', 'Rest_synth']
conditions = ['spont', 'Rest_synth']
stats_test = 'wilcoxon'



'''GoNoGo Wilcoxon Control'''
subject_list = ['list of subjects']

analysis_string = 'Wilcoxon_control'

dfa_core.jackknife_resampling(subject_list, frequencies, master_folder, morphing_op_path, graphfolder, conditions, stats_test, output_folder, analysis_string, pv_values = None)



'''GoNoGo Wilcoxon ADHD'''
subject_list = ['list of subjects']

analysis_string = 'Wilcoxon_adhd'

dfa_core.jackknife_resampling(subject_list, frequencies, master_folder, morphing_op_path, graphfolder, conditions, stats_test, output_folder, analysis_string, pv_values = None)