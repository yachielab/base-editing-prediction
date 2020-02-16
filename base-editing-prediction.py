import os 
import re
import sys
import time
import pickle
import itertools
import collections
import pickle
import numpy as np
import argparse
import visualization as vis

def read_csv(fname):
    '''
    Read input model
    '''
    info = []
    with open(fname) as f:
        for line in f:
            if line[0] == "#":
                info.append(line) 
            else: 
                break
        if len(info) > 0:
            name = info[0][1:].rstrip() 
        else:
            name = ".".join(fname.split(".")[:-1])
        m_rate_dict = collections.defaultdict(dict)
        m_pattern_dict = collections.defaultdict(lambda : collections.defaultdict(dict))
        for line in f:
            line = line.rstrip().split(",") 
            if line[0] == "TP": #Save transition probability
                line[2] = line[2].split(":")
                pos  = int(line[2][0]) 
                nuc1 = line[2][1][0] 
                nuc2 = line[2][1][-1] 
                if line[-1] != "N.A.":
                    m_rate_dict[(pos,nuc1)][nuc2] = float(line[-1]) 
                else:
                    m_rate_dict[(pos,nuc1)][nuc2] = None 

            elif line[0] == "CTP": #Save conditional transition probability
                #Conditional transition
                line[1] = line[1].split(":") 
                pos1    = int(line[1][0])
                nuc1    = line[1][1][0] 
                nuc2    = line[1][1][-1] 
                
                #Target transition 
                line[2] = line[2].split(":") 
                pos2    = int(line[2][0]) 
                nuc3    = line[2][1][0] 
                nuc4    = line[2][1][-1] 
                if line[-1] != "N.A.":
                    m_pattern_dict[(pos1,pos2),(nuc1,nuc3)][nuc2][nuc4] = float(line[-1]) 
                else:
                    m_pattern_dict[(pos1,pos2),(nuc1,nuc3)][nuc2][nuc4] = None
            else: 
                pass 
    return m_rate_dict, m_pattern_dict, name

def simulation(ref,query,mutation_rate_dict,mutation_pattern_dict,start=-30,end=10,enzyme=None):
    '''
    Simulation of an editing pattern frequency.
    '''
    pos    = 0 
    value  = 1.0
    values = [] 
    rps    = []
    
    #Search editing position and save into 'rps'.
    for pos, (r,q) in enumerate(zip(ref,query)):
        if r != q:
            rps.append(pos)  
    all_combinations = itertools.product(rps,list(range(end-start+1)))
    
    #Product of base transition probablity at each position.
    for pos in rps:
        p = query[pos] 
        v = mutation_rate_dict[(pos+start,ref[pos])][p]
        if v == None:
            v = 0.0
        value = value * v                     
        values.append(v)
   
    #Product of conditional base transition probablity at each position given a base transition at edited position.
    for combi in all_combinations:
        r_combi = (combi[0]+start,combi[1]+start)
        ref_set = (ref[combi[0]],ref[combi[1]]) 
        if combi[0] in rps: #if conditional base 
            p = query[combi[0]] 
            q = query[combi[1]]
            v = mutation_pattern_dict[(r_combi,ref_set)][p][q] 
            if v == None:
                v = 1
            value = value * v
       
    #Calculation of geometric mean. 
    value = value**(1.0/len(rps)) 
    if np.isnan(value) or value > min(values):
        value = np.prod(values) 
    return value

def simualtion_all(ref,mutation_rate_dict,mutation_pattern_dict,start=-30,end=10,enzyme=None):
    '''
    Simulation of all possible editing pattern frequencies.
    '''
    chars = [] 
    for i, char in enumerate(ref):
        if i + start <= -6:
            if char == "A":
                chars.append(("A","G")) 
            elif char == "C": 
                chars.append(("C","T"))
            else:
                chars.append((char,))
        else:
            chars.append((char,))
    
    #Calculation frequencies of all possible edting pattern 
    pattern_value_dict = {}
    for pattern in itertools.product(*chars): 
        pattern = "".join(pattern)
        if pattern != ref:
            value   = simulation(ref,pattern,mutation_rate_dict,mutation_pattern_dict,start,end) 
            pattern_value_dict[pattern] = value
    
    #Save frequency of each editing pattern
    with open("editing_pattern_freq.txt","w") as o:
        o.write("#Reference sequence:{},Start:{},End:{},Base editor:{}\n".format(ref,start,end,enzyme))
        o.write("Editing_pattern,Editing_frequency\n".format(ref,start,end))
        mono_editing_specs = [0] * len(ref)
        items = list(pattern_value_dict.items()) 
        items.sort(key=lambda x:x[1]) 
        items.reverse()
        for pattern in list(zip(*items))[0]:
            o.write("{},{}\n".format(pattern,pattern_value_dict[pattern])) 
            for i, (r,q) in enumerate(list(zip(ref,pattern))):
                if r != q:
                    mono_editing_specs[i] += pattern_value_dict[pattern]         
    
    #Save values of a editing spectrum
    with open("editing_spectrum.txt","w") as o:
        print("Relative position from the PAM","Reference nucleotide","A-to-G mutation rate","C-to-T mutation rate", sep=",", file=o)
        for i, value in enumerate(mono_editing_specs):
            if ref[i] == "A":
                print(start + i, ref[i], mono_editing_specs[i], "N.A.", sep = ",", file=o)
            elif ref[i] == "C":
                print(start + i, ref[i], "N.A.", mono_editing_specs[i], sep = ",", file=o)
            else: 
                print(start + i, ref[i], "N.A.", "N.A.", sep = ",", file=o)

    vis.main(ref,mono_editing_specs,start,end,enzyme) 

if __name__ == "__main__":
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--input",       type=str,  default="ACACACACACACTCTGATCATACGA")
    p.add_argument("-p","--pattern",      type=str,  default="None")
    p.add_argument("-b","--base_editor", type=str,  default="samples/TargetACEmax.csv")
    p.add_argument("-s","--start",       type=int,  default=-25)
    p.add_argument("-e","--end",         type=int,  default=-1)
    args = p.parse_args() 

    m_rate_dict, m_pattern_dict, name = read_csv(args.base_editor) 
    if args.pattern == "None":
        simualtion_all(args.input,m_rate_dict,m_pattern_dict,start=args.start,end=args.end, enzyme=name)
    else:
        value = simulation(args.input,args.pattern,m_rate_dict,m_pattern_dict,start=args.start,end=args.end, enzyme=name)
        print("Base editor        : {}".format(name))
        print("Reference sequence : {}".format(args.input)) 
        print("Editing pattern    : {}".format(args.pattern))
        print("Editing frequency  : {}".format(value))
