#!/python
#usage: python3 depth_tTest.py your.control.coverage.file your.HE.coverage.file out.statistics out.txt
#your.control.coverage.file and your.HE.coverage.file are two input files (the results of bedtools coverage)
#out.statistics and out.txt are two output files
import pandas as pd
import numpy as np
import itertools
import sys
from scipy import stats
from collections import Counter
import statsmodels.stats.multitest as smsm
def same_dup_del(ratio,q):###
    if ratio >= 1.5:
        return "dup"
    elif ratio <= 0.5:
        return "del"
    else:
        return "same"
read_len = 150
window_num = 1
control_in = sys.argv[1]#"your.control.coverage.file"
he_in = sys.argv[2]#"your.HE.coverage.file"
control = pd.read_table(control_in,names = ["chr","start","end","control_read_num","total_reads","window_size","control_cov"])#Control.cov
control = control[control["control_cov"] > 0.7]
he = pd.read_table(he_in,names = ["chr","start","end","HE_read_num","total_reads","window_size","HE_cov"])#HE.cov
m = pd.merge(control, he, on = ["chr","start","end","window_size"])
m['control_depth'] = m['control_read_num'] * read_len / (m["window_size"] + 1)
m['HE_depth'] = m['HE_read_num'] * read_len / (m["window_size"] + 1)
depth_factor = (sum(m['control_read_num'])/float(sum(m['window_size'])))/(sum(m['HE_read_num'])/float(sum(m['window_size'])))
m['control_depth_norm'] = m['control_depth']
m['HE_depth_norm'] = m['HE_depth'] * depth_factor
m = m.sort_values(by = ["chr","start"])
j = 0
dic_merge =  {"chr":[],"start":[],"end":[],"type":[]}
for name,group in m.groupby(m["chr"]):
    j += 1
    dic = {"chr":[],"start":[],"end":[],"control_depth_ave":[],"HE_depth_ave":[],"p-value":[]}
    tot_num = len(group)
    #for i in range(0,tot_num - window_num + 1):
        #sub = group.iloc[i:i+window_num,:]
    for i in range(0,tot_num,window_num):
        sub = group.iloc[i:i+window_num,:]
        t,p=stats.ttest_ind(list(sub['control_depth_norm']),list(sub['HE_depth_norm']),equal_var=False)
        chr,start,end,control_depth,HE_depth = list(sub["chr"])[0],min(sub["start"]),max(sub["end"]),np.mean(sub["control_depth_norm"]),np.mean(sub["HE_depth_norm"])
        #print(chr,start,end,control_depth,HE_depth,p,sep = "\t")
        dic["chr"].append(chr)
        dic["start"].append(start)
        dic["end"].append(end)
        dic["control_depth_ave"].append(control_depth)
        dic["HE_depth_ave"].append(HE_depth)
        dic["p-value"].append(p)
    block = pd.DataFrame(dic)
    block = block[['chr','start','end','control_depth_ave','HE_depth_ave','p-value']]
    block = block.sort_values(by = ["chr","start"])
    block['q-value'] = smsm.multipletests(list(block['p-value']),alpha=0.05)[1]
    block['ratio'] = (block['HE_depth_ave'] + 0.01)/(block['control_depth_ave'] + 0.01)
    block["type"] = block.apply(lambda row: same_dup_del(row['ratio'],row["q-value"]),axis = 1)
    #block.to_csv("haha",index = False,sep = "\t")
    if j == 1:
        all = block
    else:
        all = pd.concat([all,block])
    ###merge same types###
    block1 = block.copy(deep=True)
    for k,v in itertools.groupby(list(block["type"])):
        length = len(list(v))
        sub_block = block1.iloc[0:length,:]
        dic_merge["chr"].append(list(sub_block['chr'])[0])
        dic_merge["start"].append(min(sub_block['start']))
        dic_merge["end"].append(max(sub_block['end']))
        dic_merge["type"].append(k)
        block1.drop(block1.head(length).index,inplace = True)
merge = pd.DataFrame(dic_merge)
merge = merge[['chr','start','end','type']]
all.to_csv(sys.argv[3],index = False,sep = "\t")#out.statistics
merge.to_csv(sys.argv[4],index = False,sep = "\t")#out.txt
