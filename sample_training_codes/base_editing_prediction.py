import os 
import re
import sys
import time
import pickle
import itertools
import collections
import pickle
import numpy as np
from scipy.stats import pearsonr

def learn(samples,start=-30,end=10,min_identity=0.8,min_quality=0): 
    mutation_rate_dict = {} 
    for p in range(end-start+1):
        for n, in ["A","T","G","C"]:
             mutation_rate_dict[(p+start,n)] = {"A":[], "T":[], "G":[], "C":[], "N":[]}
    
    mutation_pattern_dict = {} 
    for p,q in itertools.product(range(end-start+1),range(end-start+1)):
        for n,m in itertools.product(["A","T","G","C"],["A","T","G","C"]):
            mutation_pattern_dict[(p+start,q+start),(n,m)] = {"N": {"A":[], "T":[], "G":[], "C":[], "N":[]}, "A": {"A":[], "T":[], "G":[], "C":[], "N":[]}, "T":{"A":[], "T":[], "G":[], "C":[], "N":[]}, "G":{"A":[], "T":[], "G":[], "C":[], "N":[]}, "C":{"A":[], "T":[], "G":[], "C":[], "N":[]}}

    for sample in samples:
        sub_spec, occupancy = sample.get_sub_spectrum(start=start,end=end, min_identity=min_identity, min_quality=min_quality)  
        mut_dict = sample.__comatrixes_dict__[(start,end,"relative",min_identity,min_quality,"all")]["info"]
        astart   = sample.zero_position + start
        aend     = sample.zero_position + end
        for p,m in enumerate(sample.reference[astart:aend+1]):
            for a in ["A","T","G","C","N"]: 
                if p != a:
                    mutation_rate_dict[(p+start,m)][a].append(sub_spec[p] * occupancy[p][a])
                else: 
                    mutation_rate_dict[(p+start,m)][a].append(1-sub_spec[p])

        for p,m in enumerate(sample.reference[astart:aend+1]):
                for q,n in enumerate(sample.reference[astart:aend+1]):
                    for a,b in itertools.product(["A","T","G","C","N"],["A","T","G","C","N"]):
                        edit_freq = mut_dict[(p+start,q+start),(m,n)][a][b] 
                        freq = sum(list(mut_dict[(p+start,q+start),(m,n)][a].values()))
                        if freq != 0:
                            c_prob = edit_freq / freq
                        else:
                            c_prob = None
                        mutation_pattern_dict[(p+start,q+start),(m,n)][a][b].append(c_prob) 
    
    for key1 in mutation_rate_dict:
        for key2 in mutation_rate_dict[key1]:
            if len(mutation_rate_dict[key1][key2]) > 0:
                mutation_rate_dict[key1][key2] = np.mean([value for value in mutation_rate_dict[key1][key2] if value != None]) 
            else:
                mutation_rate_dict[key1][key2] = None 
    for key1 in mutation_pattern_dict:
        for key2 in  mutation_pattern_dict[key1]:
            for key3 in mutation_pattern_dict[key1][key2]:
                if len(mutation_pattern_dict[key1][key2][key3]) > 0:
                    mutation_pattern_dict[key1][key2][key3] = np.mean([value for value in mutation_pattern_dict[key1][key2][key3] if value != None]) 
                else:
                    mutation_pattern_dict[key1][key2][key3] = None

    return mutation_rate_dict, mutation_pattern_dict

def simulation(ref,query,mutation_rate_dict,mutation_pattern_dict,start=-30,end=10):
    pos    = 0 
    value  = 1.0
    values = [] 
    rps    = []
    for pos, (r,q) in enumerate(zip(ref,query)):
        if r != q:
            rps.append(pos)
    
    all_combinations = itertools.product(rps,list(range(end-start+1)))
    for pos in rps:
        p = query[pos] 
        v = mutation_rate_dict[(pos+start,ref[pos])][p]
        if v == None:
            v = 0.0
        value = value * v                     
        values.append(v)
   
    for combi in all_combinations:
        r_combi = (combi[0]+start,combi[1]+start)
        ref_set = (ref[combi[0]],ref[combi[1]]) 
        if combi[0] in rps and combi[1] in rps:
            p = query[combi[0]] 
            q = query[combi[1]]
            v = mutation_pattern_dict[(r_combi,ref_set)][p][q] 
            if v == None:
                v = mutation_rate_dict[(r_combi[1],ref_set[1])][q]
                if v == None:
                    v = 0
                    #v = 1-sum([val for val in mutation_pattern_dict[(r_combi,ref_set)][p].values() if val != None])
            value = value * v
           
        elif combi[0] in rps and combi[1] not in rps:
            p = query[combi[0]]
            v = mutation_pattern_dict[(r_combi,ref_set)][p][ref_set[1]] 
            if v == None:
                v = 1
                #v = 1-sum([val for val in mutation_pattern_dict[(r_combi,ref_set)][p].values() if val != None])
            value = value * v
        
    value = value**(1.0/len(rps))
    if type(value) == complex:
        value 

    if value > min(values):
        value = 1
        for v in values:
            value = value * v
    if np.isnan(value):
        value = 0
    return value

def validation(sample, mutation_rate_dict, mutation_pattern_dict,start=-30,end=10):
    astart   = sample.zero_position + start
    aend     = sample.zero_position + end
    ref      = sample.reference[astart:aend+1] 
    pattern_dict     = sample.get_mut_pattern(start=start, end=end)
    observed_values  = [] 
    predicted_values = [] 
    patterns = [] 
    for pattern in list(pattern_dict.keys())[1:100]:
        if ref != pattern[0] and pattern[1] == "" and "-" not in pattern[0]:
            ovalue = pattern_dict[pattern][0]
            svalue = simulation(ref,pattern[0],mutation_rate_dict,mutation_pattern_dict,start=start,end=end)
            observed_values.append(ovalue)
            predicted_values.append(svalue)
            patterns.append(pattern)
    print(pearsonr(observed_values, predicted_values)) 
    return patterns, observed_values, predicted_values

if __name__ == "__main__":
    start = -25
    end   = -1
    sample = pickle.load(open(sys.argv[1],"rb"))
    m_rate_dict, m_pattern_dict = learn([sample])
    validation(sample,m_rate_dict,m_pattern_dict,start=start,end=end)
