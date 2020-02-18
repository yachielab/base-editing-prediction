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
from tqdm import tqdm

def read_csv(fname):
    '''
    Load input model
    '''
    info = []
    with open(fname) as f:
        for line in f:
            if line[0] == "#":
                info.append(line) 
            else: 
                break
        if len(info) > 0:
            name = info[0].rstrip().split(":")[1].replace(" ","")  
        else:
            name = ".".join(fname.split(".")[:-1])
        m_rate_dict = collections.defaultdict(dict)
        m_pattern_dict = collections.defaultdict(lambda : collections.defaultdict(dict))
        for line in f:
            line = line.rstrip().split(",") 
            #Read conditional transition probability
            if line[0] == "TP": 
                line[2] = line[2].split(":")
                pos  = int(line[2][0]) 
                nuc1 = line[2][1][0] 
                nuc2 = line[2][1][-1] 
                if line[-1] != "N.A.":
                    m_rate_dict[(pos,nuc1)][nuc2] = float(line[-1]) 
                else:
                    m_rate_dict[(pos,nuc1)][nuc2] = None 
            
            #Read conditional transition probability
            elif line[0] == "CTP": 
                #Conditional base transition
                line[1] = line[1].split(":") 
                pos1    = int(line[1][0])
                nuc1    = line[1][1][0] 
                nuc2    = line[1][1][-1] 
                
                #Target base transition 
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

def simulation(ref,query,transition_prob_dict,conditional_transition_prob_dict,start=-30,end=10,model=None):
    '''
    Predicting a frequency of a given editing outcome pattern
    '''
    pos    = 0 
    value  = 1.0
    values = [] 
    edited_positions = []
    
    #Search editing position and save into 'edited_positions'.
    for pos, (r,q) in enumerate(zip(ref,query)):
        if r != q:
            edited_positions.append(pos)  
    all_combinations = itertools.product(edited_positions,list(range(end-start+1)))
    
    #Product of base transition probablity at each position.
    for pos in edited_positions:
        p = query[pos] 
        v = transition_prob_dict[(pos+start,ref[pos])][p]
        if v == None:
            v = 0.0
        value = value * v                     
        values.append(v)
    
    #Product of conditional base transition probablity at each position given a base transition at edited position.
    for combi in all_combinations:
        r_combi = (combi[0]+start,combi[1]+start)
        ref_set = (ref[combi[0]],ref[combi[1]]) 
        p = query[combi[0]] 
        q = query[combi[1]]
        v = conditional_transition_prob_dict[(r_combi,ref_set)][p][q] 
        if v == None:
            v = 1
        value = value * v
    
    #Calculation of geometric mean. 
    value = value**(1.0/len(edited_positions)) 
    if np.isnan(value) or value > min(values):
        value = np.prod(values)
    return value

def simualtion_all(target, transition_prob_dict, conditional_transition_prob_dict, start=-30, end=10, model=None, output=None):
    '''
    Predicting frequencies of the all possible base editing patterns for a given target sequence .
    '''
    chars  = [] 
    values = [] 
    for i, char in enumerate(target):
        values.append(sum(list(transition_prob_dict[i + start,char].values())))
    
    #Generation of all possible base edting pattern 
    for i, char in enumerate(target):
        if values[i] >= 2.0e-3:
            chars.append(("A","G","C","T")) 
            print(i) 
        else:
            chars.append((char,))
            
    #Predicting frequencies of all possible base edting pattern 
    outcome_value_dict = {}
    for outcome in itertools.product(*chars): 
        outcome = "".join(outcome)
        if outcome != target:
            value   = simulation(target,outcome,transition_prob_dict,conditional_transition_prob_dict,start,end) 
            outcome_value_dict[outcome] = value

    #Save frequency of each editing outcome into output file
    with open("{}_allpatterns.csv".format(output),"w") as o:
        o.write("#Model name      : {}\n".format(model))
        o.write("#Target sequence : {}\n".format(target))
        o.write("#Start position  : {}\n".format(start))
        o.write("#End poition     : {}\n".format(end))
        o.write("Editing outcome,Editing frequency\n")
    
        #Calculation of editing spectrum 
        editing_spec = collections.defaultdict(lambda : {"A":0, "T":0, "G":0, "C":0}) 
        items = list(outcome_value_dict.items()) 
        items.sort(key=lambda x:x[1]) 
        items.reverse()
        for outcome in list(zip(*items))[0]:
            o.write("{},{}\n".format(outcome,outcome_value_dict[outcome])) 
            for i, (t,p) in enumerate(list(zip(target,outcome))):
                if t != p:
                    editing_spec[i][p] += outcome_value_dict[outcome]         
     
    #Save values of the editing spectrum into output file
    with open("{}_spectrum.csv".format(output),"w") as o:
        print("#Model name      : {}".format(model), file=o)
        print("#Target sequence : {}".format(args.input), file=o) 
        print("#Start position  : {}".format(start), file=o)
        print("#End poition     : {}".format(end), file=o)
        print("Position to the PAM","Target nucleotide", "Frequency of A","Frquency of T", "Frquency of G", "Frequency of C", sep=",", file=o)
        for i in range(len(target)):
            if   target[i] == "A":
                print(start + i, target[i], 0.0, editing_spec[i]["T"], editing_spec[i]["G"], editing_spec[i]["C"], sep = ",", file=o)
            elif target[i] == "T":
                print(start + i, target[i], editing_spec[i]["A"], 0.0, editing_spec[i]["G"], editing_spec[i]["C"], sep = ",", file=o)
            elif target[i] == "G":
                print(start + i, target[i], editing_spec[i]["A"], editing_spec[i]["T"], 0.0, editing_spec[i]["C"], sep = ",", file=o)
            else:
                print(start + i, target[i], editing_spec[i]["A"], editing_spec[i]["T"], editing_spec[i]["G"], 0.0, sep = ",", file=o)
    vis.main(target, editing_spec, start, end, model) 

if __name__ == "__main__":
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument('-h', "--help", action='store_true', default=False) 
    p.add_argument("-i","--input",          type=str,  default="ACACACACACACTCTGATCATACGA")
    p.add_argument("-o","--outcome",        type=str,  default="None")
    p.add_argument("-m","--training_model", type=str,  default="sample_models/TargetACEmax.csv")
    p.add_argument("-s","--start",          type=int,  default=None)
    p.add_argument("-e","--end",            type=int,  default=None)
    p.add_argument("-f","--filepath",       type=str,  default="./output/None")
    args = p.parse_args() 
    if args.help:
        print()
        with open("help.txt") as f:
            for line in f:
                print(line.rstrip())
        print()
        exit() 

    #Load training model.
    tp_dict, ctp_dict, model = read_csv(args.training_model)
    
    #Make output directory.
    if args.outcome == "None":
        output_path = args.filepath
        if output_path[0:2] != "./" and output_path[0] != "/":
            output_path = "./" +  output_path
        output_dir  = "/".join(output_path.split("/")[:2]) 
        output_file = output_path.split("/")[-1] if output_path.split("/")[-1] != "None" else model
        if output_dir.split("/")[-1] not in os.listdir("/".join(output_dir.split("/")[:-1])):
            os.mkdir(output_dir.split("/")[-1]) 
        os.chdir(output_dir) 

    #Set the position of the input target sequence. 
    target = args.input
    if args.start == None and args.end == None:
        raise ValueError("The argument '-s' or '-e' should be specified.") 
    elif args.start == None:
        end   = args.end
        start = end - len(target) + 1 
    elif args.end == None:
        start = args.start
        end   = start + len(target) - 1 
    else: 
        if (args.end - args.start) + 1 != len(target):
            raise ValueError("The difference between '-s' and '-e' should be equal to legth minus 1 of target sequence") 
        else:
            start = args.start
            end   = args.end
         
    #Execution
    if args.outcome == "None":
        simualtion_all(target, tp_dict, ctp_dict, start=start, end=end, model=model, output=output_file)
    else:
        value = simulation(target, args.outcome, tp_dict, ctp_dict, start=start, end=end, model=model)
        print("Model name        : {}".format(model))
        print("Target sequence   : {}".format(args.input)) 
        print("Editing outcome   : {}".format(args.outcome))
        print("Editing frequency : {}".format(value))
