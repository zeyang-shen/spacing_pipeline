#!/usr/bin/env python

# load data structure packages
import numpy as np
import pandas as pd
import os
import fnmatch
import argparse
from Bio import motifs, SeqIO, Seq

# load package for doing motif-related analysis
import sys
import json
import ast

# load package for doing spacing analysis
from itertools import chain
from scipy import stats

# load package for plotting
import matplotlib.pyplot as plt
import seaborn as sns

def load_motifs(motif_dir, pseudocounts=0.01, key='full'):
    motif_dict = {}
    nuc = ['A', 'C', 'G', 'T']
    for mf in os.listdir(motif_dir):
        with open(motif_dir + mf) as f:
            m = motifs.read(f, 'jaspar')
            counts = np.array([m.counts[n] for n in nuc])
            avg_counts = counts.sum(axis=0).mean()
            m.pseudocounts = avg_counts*pseudocounts
            m.background = None
            if key == 'full':
                motif_dict[m.name+'$'+m.matrix_id] = m
            elif key == 'id':
                motif_dict[m.matrix_id] = m  
    return motif_dict


def getDist_cutoff_list_no_overlap(row):
    peak1id=row[-2]
    peak2id=row[-1]
    for strands in ['++', '--', '+-', '-+']:
        s1=strands[0]
        s2=strands[1]
        col_to_call1="adjusted motif position " + s1 
        col_to_call2="adjusted motif position " + s2 

        tf1_start_pos_list=ast.literal_eval(tf1_peak_df.loc[peak1id][col_to_call1])  # list of positions start 
        tf2_start_pos_list=ast.literal_eval(tf2_peak_df.loc[peak2id][col_to_call2])
        
        if len(tf1_start_pos_list)>0 and len(tf2_start_pos_list)>0:
            for tf1_start_pos in tf1_start_pos_list:
                for tf2_start_pos in tf2_start_pos_list:
                    tf1_end_pos = tf1_start_pos + len(motif_dict[motif_id1])
                    tf2_end_pos = tf2_start_pos + len(motif_dict[motif_id2])
        
                    #compute edge-to-edge motif distance
                    if tf1_start_pos >= tf2_end_pos: #non-overlapping
                        dist = tf1_start_pos - tf2_end_pos
                    elif tf2_start_pos >= tf1_end_pos: #non-overlapping
                        dist = -(tf2_start_pos - tf1_end_pos)
                    elif tf1_start_pos < tf2_start_pos: #overlapping
                        continue
                    elif tf1_start_pos > tf2_start_pos: #overlapping
                        continue
                    elif tf1_start_pos == tf2_start_pos: #overlapping
                        continue
                    spacing_dict[strands].append(dist)

                    
def test_spacing(spacing_array, low=-100, high=100):
    '''
    Statistical testing of spacing relationship
    
    Input:
        spacing_array: an array of spacing
        low: the lower bound of the spacing to consider
        high: the upper bound of the spacing to consider
    
    Outputs:
        pKS: p-value of KS-test on spacing distribution
    '''
    test_list = np.array(spacing_array)
    test_list = test_list[test_list >= low]
    test_list = test_list[test_list <= high]
    low = min(low, min(test_list))
    high = max(high, min(test_list))
    
    bins = np.arange(low,high+1)
    count, div = np.histogram(test_list, bins=bins)
    
    #KS-test
    p_ks_list = []
    for k in range(100):
        bg_spacing_list = np.random.randint(low=low, high=high, size=len(test_list))
        p_ks_list.append(stats.ks_2samp(test_list, bg_spacing_list)[1])
    pKS = np.mean(p_ks_list)
    
#     #base-wise test
#     pBase = []
#     for i in np.arange(low,high+1):
#         p_binom = stats.binom_test(corrected_count[i], n=len(test_list), p=1/(high-low))
#         p_poiss = stats.poisson.pmf(corrected_count[i], len(test_list)*1/(high-low))
#         if p_binom < 0.05/(high-low) or p_poiss < 0.05/(high-low):
#             pBase.append((i, p_binom, p_poiss))
    
    return pKS


def delta(arr, ave=1):
    """ calculate the slop at given point,
    ave: positions to average to compute slope, 
    e.g., if ave = 2, to compute slope of position 8, position 6,7,9,10 will be averaged
    """
    slope=[]
    for i in range(len(arr)):
        if i < ave: # firstt few ones
            current_count = arr[i]
            post = np.average(arr[i+1:i+ave+1])
            slope.append(abs(current_count-post))
        elif len(arr) <= (i+ave): # only use pre
            current_count = arr[i]
            pre=np.average(arr[i-ave:i])
            slope.append(abs(current_count-pre))
        else:
            current_count = arr[i]
            post = np.average(arr[i+1:i+ave+1])
            pre=np.average(arr[i-ave:i])
            ave_slope = (abs(post-current_count) + abs(pre-current_count))/2.0
            slope.append(ave_slope)   
    return slope


def smooth_avg(raw_data, size=2):
    new_data = []
    if size % 2 == 0:
        shift = 0
    else:
        shift = 1
    for i in np.arange(len(raw_data)):
        if i < size//2:
            new_data.append(np.mean(raw_data[:i+size//2+1]))
        elif i+size//2 >= len(raw_data):
            new_data.append(np.mean(raw_data[i-size//2:]))
        else:
            start = max(i-size//2,0)
            end = min(i+size//2+shift,len(raw_data))
            new_data.append(np.mean(raw_data[start:end]))
    return np.array(new_data)


def pval(item):
    #need distribution array in all_slope list
    #compute pvalue given the distribution
    p=1-((stats.percentileofscore(all_slope, item))/100.0)
    return(p)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute spacing between motifs; Example: python characterize_spacing.py ../ENCODE_processed_files/ GATA1 TAL1 --motif_path ../motifs/')
    parser.add_argument("path",
                        help="path to the folder of files from script identify_motif.py",
                        type=str)
    parser.add_argument("tf1",
                        help="TF1",
                        type=str)
    parser.add_argument("tf2",
                        help="TF2",
                        type=str)
    parser.add_argument("--motif_path", 
                        help="path to the folder that stores motif files",
                        type=str)
    args = parser.parse_args()
    
    tf1 = args.tf1
    tf2 = args.tf2
    path = args.path
    motif_path = args.motif_path
    
    original_command = " ".join(["python characterize_spacing.py", path, tf1, tf2, "--motif_path "+motif_path])
    
    #load motif PWMs
    motif_dict = load_motifs(motif_dir=motif_path)
    allmotif = np.array(list(motif_dict.keys()))
    motif_id1 = sorted([i for i in allmotif if i.split('$')[0].upper() == str(tf1).upper()])
    motif_id2 = sorted([i for i in allmotif if i.split('$')[0].upper() == str(tf2).upper()])
    try:
        motif_id1 = motif_id1[0]
        motif_id2 = motif_id2[0]
        print('TF1:', motif_id1, str(len(motif_dict[motif_id1]))+'bp')
        print('TF2:', motif_id2, str(len(motif_dict[motif_id2]))+'bp')
    except:
        sys.exit('ERROR: motif files not found!')
    
    #merge peaks using HOMER
    os.system("mkdir -p "+path+"/merged_files")
    os.system("mkdir -p "+path+"/spacing_files")
    
    try:
        idr1_path = path+fnmatch.filter(os.listdir(path), '*'+tf1+'*_cutoff.tsv')[0]
        idr2_path = path+fnmatch.filter(os.listdir(path), '*'+tf2+'*_cutoff.tsv')[0]
    except:
        sys.exit('ERROR: motif-annotated IDR file not found!')
    
    out_file=path+"/merged_files/merged_"+str(tf1)+"_"+str(tf2)+'_filtered.tsv'
    cmd="mergePeaks -d given " + idr1_path + " " + idr2_path +" > " +out_file
    print('Running HOMER mergePeaks...')
    os.system(cmd)

    #look for co-binding peaks
    merge_df = pd.read_csv(out_file, sep='\t', index_col=0)
    cobind_bools = np.all([np.all(merge_df.iloc[:,-2:].notnull(), axis=1), #peak bound by both TFs
                           merge_df['Total subpeaks'] == 2 #only one peak from each TF to reduce complication
                          ], axis=0)
    tf1_peak_df = pd.read_csv(idr1_path, sep='\t', index_col=0)
    tf2_peak_df = pd.read_csv(idr2_path, sep='\t', index_col=0)
    cobind_peaks = merge_df.loc[cobind_bools]
    
    #compute and save spacings
    spacing_dict = {k:[] for k in ['++', '--', '+-', '-+']}
    cobind_peaks.apply(getDist_cutoff_list_no_overlap,axis=1)
    json_path_np=path+'/spacing_files/'+str(tf1)+'_'+str(tf2)+'_spacing.json'
    with open(json_path_np, "w") as outfile:
        json.dump(spacing_dict, outfile)
    print('Spacings successfully saved to', json_path_np)
    
    #load null distribution
    with open(os.path.dirname(__file__)+'/simul_random_null.json', 'r') as JSON:
        json_dict = json.load(JSON)
    keys = list(json_dict.keys())
    all_slope=[]
    for key in keys:
        s1 = json_dict[key]
        s1 = np.array(s1)
        bins = np.arange(0, 101, 1)
        null_count, _ = np.histogram(s1, bins=bins)
        null_count = null_count/sum(null_count)
        all_slope.append(delta(null_count))
    all_slope = list(chain.from_iterable(all_slope))
    
    #test spacing
    sns.set(style='white', font_scale=1.5)
    plt.figure(figsize=(10,6))
    stats_file = path+'/spacing_files/'+'_'.join([tf1, tf2])+'_spacing_stats.tsv'
    with open(stats_file, 'w') as wfile:
        wfile.write('\t'.join(['orientation', 'sample size', 'constrained spacing', 'KS-test'])+'\n')
        for key in spacing_dict.keys():
            spacing = np.array(spacing_dict[key])
            bins = np.arange(-100, 101)
            count, division = np.histogram(spacing, bins=bins)
            norm_count = count/sum(count)
            pvalue_res=list(map(pval, delta(norm_count)))
            mini_peaks = []
            pv_cutoff = 0.05/200/4
            for i, pv in enumerate(pvalue_res):
                if i == 0 or i+1 == len(pvalue_res):
                    continue
                if pv < pv_cutoff and pvalue_res[i-1] < pv_cutoff and pvalue_res[i+1] < pv_cutoff:
                    mini_peaks.append((division[i], pv))
            mini_peaks = np.array(mini_peaks)
            for p in mini_peaks:
                plt.axvline(p[0], linestyle='--', c='grey')
            plt.plot(division[:-1], smooth_avg(norm_count, size=5), label=key)

            mini_peaks_save = '|'.join([str(int(m[0]))+'('+str(m[1])+')' for m in mini_peaks])
            if mini_peaks_save == '':
                mini_peaks_save = 'nan'
            wfile.write('\t'.join([key, str(len(spacing)), mini_peaks_save, 
                                   str(test_spacing(spacing, low=-100, high=100))])+'\n')
    plt.legend()
    plt.title(' - '.join([tf1, tf2]))
    plt.xlabel('Spacing (bp)')
    plt.savefig(path+'/spacing_files/'+'_'.join([tf1, tf2])+'_spacingDistribution.png', format='png')
    plt.close()
    print('Spacing stats successfully saved to', stats_file)
