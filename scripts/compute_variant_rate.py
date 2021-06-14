#!/usr/bin/env python
import numpy as np
import pandas as pd
import os

def extractAF(info):
    if (pd.isna(info) ==True ):
        return (float('nan'))
    else:
        sp=info.split(';')
        af=float([af for af in sp  if "AF=" in af][0].split("=")[1])
        return (af)

def extractAC(info):
    if (pd.isna(info) ==True ):
        
        return (float('nan'))
    else:
        sp=info.split(';')
        ac=float([ac for ac in sp  if "AC=" in ac][0].split("=")[1])
        return (ac)

def get_variant_rate(bed_path, overlap_path, spa=None):
    bed_df = pd.read_csv(bed_path, sep='\t', header=None)
    overlap_df = pd.read_csv(overlap_path, sep='\t', header=None)
    indel_df = pd.merge(bed_df, overlap_df, how='outer', on=[3,3]) #use union of keys from both frames
    
    indel_df['AF']=indel_df.loc[:,11].apply(extractAF).fillna(0)
    indel_df['AC']=indel_df.loc[:,11].apply(extractAC).fillna(0)
    indel_df['delta']=pd.Series(indel_df.loc[:,7].str.len() - indel_df.loc[:,8].str.len()).fillna(0)
    indel_df['delta']=np.absolute(indel_df['delta'])
    indel_df=indel_df[indel_df['delta']<=50]
    
    aggr_df = indel_df.groupby(['0_x','1_x','2_x',3])["AF","AC","delta"].sum().reset_index()
    aggr_df['dist'] = abs(aggr_df['1_x']-aggr_df['2_x'])
    aggr_df['ACfreq'] = aggr_df['AC']/aggr_df['dist']
    if spa is not None:
        print('Restrict spacing at', spa-1, '~', spa+1)
        aggr_df=aggr_df[(aggr_df['dist']>=spa-1) & (aggr_df['dist']<=spa+1)]
    
    return (aggr_df)
    

def normalize_rate(tf,types):
    cols_list=['ACfreq_x', 'ACfreq_y']
    mid=mid_indel_freq(tf)
    sur=sur_indel_freq(tf)
    sur=pd.merge(aggr_df4, aggr_df5, how='outer', on=[3,3])
    sur["AC"]=sur["AC_x"] + sur["AC_y"]
    sur["delta"]=sur["delta_x"] + sur["delta_y"]
    sur["AF"]=sur["AF_x"] + sur["AF_y"]
    sur['ACfreq'] = sur['AC'] /200
    merge_45=merge_45.loc[:,['0_x_x','1_x_x','2_x_x',3,"AF_rare","AC","delta",'deltaXAC',"AF",'ACfreq','deltaXACfreq']]
    
    
    mid=mid.fillna(0)
    sur=sur.fillna(0)
    mid_sur_df=pd.merge(mid,sur,how='inner',on=[3,3])
    select_df=mid_sur_df[cols_list]
    
    return (select_df)

