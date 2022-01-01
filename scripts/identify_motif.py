#!/usr/bin/env python
import numpy as np
import pandas as pd
import sys
import os
sys.path.append(os.path.dirname(__file__)+'/../')
from Bio import motifs, SeqIO, Seq #biopython package
import argparse

def read_fasta(fasta_file, skip_duplicate=True, fmt='fasta'):
    '''
    Read sequences stored in a FASTA file
    '''
    alphabet = Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA() # need to use this alphabet for motif score calculation
    id_seq_dict = {} # {sequenceID: fastq sequence}
    duplicate_keys = []
    for seq_record in SeqIO.parse(fasta_file, fmt):  
        seq_record.seq.alphabet = alphabet
        if seq_record.id in id_seq_dict.keys():
            duplicate_keys.append(seq_record.id)
        else:
            id_seq_dict[seq_record.id] = seq_record.seq
    # delete duplicate keys
    if skip_duplicate:
        for dk in duplicate_keys:
            del id_seq_dict[dk]
        if len(duplicate_keys) > 0:
            print('Ignore duplicate keys in %s: %s' % (fasta_file, duplicate_keys))
    return id_seq_dict

def load_motifs(motif_dir, pseudocounts=0.01, key='full'):
    '''
    read in motifs; motifs have to be in jaspar format as below:
    
        >MA0002.2       RUNX1
        A  [   287    234    123     57      0     87      0     17     10    131    500 ]
        C  [   496    485   1072      0     75    127      0     42    400    463    158 ]
        G  [   696    467    149      7   1872     70   1987   1848    251     81    289 ]
        T  [   521    814    656   1936     53   1716     13     93   1339   1325   1053 ]
    
    Parameters:
        motif_dir: folder that contains motif files; one file for individual motif
        pseudocounts: w.r.t. position weight matrix, the probability adding to every nucleotide
        key: specify the way to name the motifs in the output dictionary
             options: 'full' (default), 'id'
    '''
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


def find_motif_cutoff(bio_motif, seq_dict, score_cutoff=False, position_cutoff=None, keep_best=False):
    '''
    compute motif score and position for both strand
    input: motif object, sequence dictory
    output: list for score+, pos+, score-, position-
    '''
    
    if score_cutoff:
        #compute background distribution
        background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25} #specify background frequency
        distribution = bio_motif.pssm.distribution(background=background, precision=10**4)
        cutoff = distribution.threshold_fpr(0.001)
    else:
        cutoff = 0
    print('score cutoff:', cutoff)
    
    fwd_pssm = bio_motif.pssm
    rev_pssm = fwd_pssm.reverse_complement()
    scores_p = []
    pos_p = []
    scores_m = []
    pos_m = []
    alphabet = Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA()
    sorted_ids = sorted(seq_dict.keys())
    for sid in sorted_ids:
        seq = seq_dict[sid]
        seq = Seq.Seq(str(seq), alphabet=alphabet)
        if len(seq) < len(bio_motif):
            sys.exit('ERROR: sequence lengths are too short to calculate motif score!')
        fwd_scores = fwd_pssm.calculate(seq) # scores for forward orientation
        rev_scores = rev_pssm.calculate(seq) # scores for reverse orientation
        if type(fwd_scores) == np.float32:
            fwd_scores = np.array([fwd_scores])
        if type(rev_scores) == np.float32:
            rev_scores = np.array([rev_scores])

        #filter by motif scores
        pls_max_scores= fwd_scores[fwd_scores>cutoff]
        pls_max_pos = np.where(fwd_scores>cutoff)[0]

        rev_max_scores= rev_scores[rev_scores>cutoff]
        rev_max_pos = np.where(rev_scores>cutoff)[0]
        
        #keep maximum score(s) only
        if keep_best:
            if len(pls_max_scores) > 0:
                pls_max = np.max(pls_max_scores)
                pls_max_id = np.where(pls_max_scores == pls_max)[0]
                pls_max_scores = pls_max_scores[pls_max_id]
                pls_max_pos = pls_max_pos[pls_max_id]
            
            if len(rev_max_scores) > 0:
                rev_max = np.max(rev_max_scores)
                rev_max_id = np.where(rev_max_scores == rev_max)[0]
                rev_max_scores = rev_max_scores[rev_max_id]
                rev_max_pos = rev_max_pos[rev_max_id]
            
        #filter by motif position
        if position_cutoff is not None:
            if len(pls_max_scores) > 0:
                pls_loc_bools = np.abs(pls_max_pos - (len(fwd_scores)+len(bio_motif)-1)/2) <= position_cutoff
                pls_max_scores = pls_max_scores[pls_loc_bools]
                pls_max_pos = pls_max_pos[pls_loc_bools]
            
            if len(rev_max_scores) > 0:
                rev_loc_bools = np.abs(rev_max_pos - (len(rev_scores)+len(bio_motif)-1)/2) <= position_cutoff
                rev_max_scores = rev_max_scores[rev_loc_bools]
                rev_max_pos = rev_max_pos[rev_loc_bools]
        
        scores_p.append(list(pls_max_scores))
        pos_p.append(list(pls_max_pos))
        scores_m.append(list(rev_max_scores))
        pos_m.append(list(rev_max_pos))

    return scores_p, pos_p, scores_m, pos_m

def filter_to_df(peak_file, seq_dict, bio_motif, size):
    scoreP,posP,scoreM,posM = find_motif_cutoff(bio_motif, seq_dict, 
                                                score_cutoff=score_cutoff, 
                                                position_cutoff=position_cutoff, 
                                                keep_best=keep_best)
    
    try:
        with open(peak_file) as f:
            cnt = 0
            for line in f:
                if line.startswith('#'):
                    cnt += 1
        peak_df = pd.read_csv(peak_file, sep='\t', index_col=0, skiprows=cnt-1)
    except FileNotFoundError:
        with open(peak_file.replace('.tsv', '.txt')) as f:
            cnt = 0
            for line in f:
                if line.startswith('#'):
                    cnt += 1
        peak_df = pd.read_csv(peak_file.replace('.tsv', '.txt'), sep='\t', index_col=0, skiprows=cnt-1)
    peak_df = peak_df.loc[sorted(seq_dict.keys())] #make sure the peak ID order is the same as that from "find_motif"

    peak_df["motif score +"] =scoreP
    peak_df["motif position +"] =posP
    peak_df["motif score -"] =scoreM
    peak_df["motif position -"] =posM

    mids = (peak_df.iloc[:,2]+peak_df.iloc[:,1])//2 #find peak summit/peak center
    peak_df["adjusted motif position +"] = peak_df.apply(lambda x: [y+x["start"] for y in x["motif position +"]],axis=1)
    peak_df["adjusted motif position -"] = peak_df.apply(lambda x: [y+x["start"] for y in x["motif position -"]],axis=1)
    
    #peak_df=peak_df.drop(columns=['motif position +', 'motif position -'])
    peak_df=peak_df[(peak_df["motif score +"].str.len() != 0) | (peak_df["motif score -"].str.len() != 0) ]
    peak_df.iloc[:,1] = mids - size//2
    peak_df.iloc[:,2] = mids + size//2
    output_file = '.'.join(peak_file.split('.')[:-1])+'_cutoff.tsv'
    peak_df.index.name = peak_df.index.name+'('+original_command+')'
    peak_df.to_csv(output_file, sep='\t')
    print('Total valid peaks:', len(peak_df), 'out of', len(seq_dict))
    
    return output_file, peak_df
    
if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description='Scan for motifs; Example: python identify_motif.py ../ENCODE_processed_files/CTCF_idr.fa CTCF --motif_path ../motifs/ --cutoff -d 50')
    parser.add_argument("file", 
                        help="input file of sequences in FASTA format",
                        type=str)
    parser.add_argument("tf", 
                        help="transcription factor symbol",
                        type=str)
    parser.add_argument("--motif_path", 
                        help="path to motif files",
                        type=str, default="./motifs/")
    parser.add_argument("-d", "--distance", 
                        help="distance cutoff to peak center",
                        type=int)
    parser.add_argument("--best",
                        help="Flag to only keep best motifs",
                        action='store_true')
    parser.add_argument("--cutoff", 
                        help="Flag to automatically generate motif score cutoff",
                        action='store_true')
    parser.add_argument("-s", "--size", 
                        help="Re-scaled size of output peaks",
                        type=int, default=400)
    parser.add_argument("--repeat", 
                        help="Repeats annotation file",
                        type=str)
    args = parser.parse_args()
    
    fasta_file = args.file
    tf = args.tf
    motif_path = args.motif_path
    position_cutoff = args.distance
    keep_best = args.best
    score_cutoff = args.cutoff
    scale_size = args.size
    repeat_file = args.repeat #HOMER is required if specified a value instead of "none"
    
    original_command = " ".join(["python identify_motif.py", fasta_file, tf, "--motif_path "+motif_path, 
                                 "--distance "+str(position_cutoff), "--best"*keep_best, "--cutoff"*score_cutoff, 
                                 "--size "+str(scale_size), "--repeat "+repeat_file])
    
    #load motif PWMs
    motif_dict = load_motifs(motif_dir=motif_path)
    allmotif = np.array(list(motif_dict.keys()))
    motif_id_list = sorted([i for i in allmotif if i.split('$')[0].upper() == str(tf).upper()])
    try:
        motif_id = motif_id_list[0]
        print('Scanning for motif:', motif_id)
    except:
        sys.exit('ERROR: no motif file found for '+tf)
    
    #read in sequences
    seq_dict = read_fasta(fasta_file, fmt='fasta')
    if len(seq_dict) < 1:
        seq_dict = read_fasta(fasta_file, fmt='tab')
        if len(seq_dict) < 1:
            sys.exit('No sequence identified!')
    
    #read peak file
    peak_file = fasta_file.replace(".fa",".tsv")
    
    #scan for motifs
    save_file, peak_df = filter_to_df(peak_file, seq_dict, motif_dict[motif_id], scale_size)
    
    #separate peaks in repetitive and nonrepetitive regions
    if repeat_file != None:
        print('Checking repetitive regions...')
        mask_merge_file = save_file.replace('.tsv', '_repeatMerged.tsv')
        cmd="mergePeaks -d given "+save_file+" "+repeat_file+" > " +mask_merge_file #HOMER is required
        os.system(cmd)

        merge_df = pd.read_csv(mask_merge_file, sep='\t', index_col=0, low_memory=False)
        mask_merge_peak = merge_df.iloc[:,-2][np.all([merge_df.iloc[:,-2].notnull(), merge_df.iloc[:,-1].isnull()], axis=0)]
        mask_peak_idx = np.unique([j for i in mask_merge_peak for j in i.split(',')])
        peak_df.loc[mask_peak_idx].to_csv(save_file.replace('.tsv', '_masked.tsv'), sep='\t')
        print('Total peaks outside specified regions:', len(mask_peak_idx))

        mask_merge_peak = merge_df.iloc[:,-2][np.all([merge_df.iloc[:,-2].notnull(), merge_df.iloc[:,-1].notnull()], axis=0)]
        mask_peak_idx = np.unique([j for i in mask_merge_peak for j in i.split(',')])
        peak_df.loc[mask_peak_idx].to_csv(save_file.replace('.tsv', '_inmask.tsv'), sep='\t')
        print('Total peaks in specified regions:', len(mask_peak_idx))
