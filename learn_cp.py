import os 
import re
import sys
import time
import itertools
import argparse
import configparser
import Levenshtein
import collections
import subprocess 
import multiprocessing as mp 
import pickle
import math
import numpy as np
import json

def read_fastq(fastq_name):
    seq_dict = {}
    with open(fastq_name.replace("'","").replace("\\","")) as f:
        n = 0 
        for line in f:
            if line[0] == "@" and n % 4 == 0:
                key = line[1:].rstrip() 
                key = key.split(" ")[0] 
                key = key.replace(":","_")
                seq_dict[key] = [line[1:].rstrip()] 
            else:
                seq_dict[key].append(line.rstrip()) 
            n += 1
    return seq_dict 

def read_fasta(fasta_name):
    seq_dict = {} 
    with open(fasta_name) as f:
        for line in f:
            if line[0] == ">":
                key = line.rstrip()[1:]
                seq_dict[key] = "" 
            else:
                seq_dict[key] += line.rstrip() 
    return seq_dict 

class SAMPLE(object):
    def __init__(self):
        self.reference     = ""
        self.target_dict   = {} 
        self.zero_position = 0 
        self.color_index   = 0 
        self.__subfreq_dict__ = {}
        self.__delfreq_dict__ = {} 
        self.__insfreq_dict__ = {} 
        self.__subspec_dict__ = {} 
        self.__delspec_dict__ = {} 
        self.__insspec_dict__ = {} 
        self.__mutpattern_dict__ = {}
        self.__indelpattern_dict__ = {}
        self.__comatrixes_dict__ = {} 
        self.__readcount_dict__ = {} 
    
    def __getitem__(self, key):
        return self.__dict__[key]
    
    def __setitem__(self, key, value):
        self.__dict__[key] = value 

    def set_reference(self, seq):
        #Refrence sequence direction should be match with Read1 sequence.
        self.reference = seq
            
    def set_zero_position(self, pos, seq=""):
        if self.reference == "":
            raise ValueError("Please set proper reference sequence before setting zero position") 
        if seq == "" or seq == "None":
            self.zero_position = pos
        else:
            seq_position = self.reference.find(seq)  
            if seq_position == -1:
                seq_position = self.reference.translate(str.maketrans("ATGC","TACG"))[::-1].find(seq)  
                if seq_position != -1:
                    self.reference = self.reference.translate(str.maketrans("ATGC","TACG"))[::-1]
                else:
                    raise ValueError("The target sequence is not exist in the reference sequence") 

            self.zero_position = seq_position + pos
    
    def add_element(self, start=0, end=0, seq="", name="", color="None", coordinate="relative"):
        if self.reference == "":
            raise ValueError("Please set proper reference sequence before adding target") 
        
        if name in seq:
            Warning("The target name is already exist in self.target_dict") 

        if (start >= 0 and end > 0) and seq == "":
            if coordinate == "relative":
                seq = self.reference[start + self.zero_position:end + self.zero_position + 1] 

            elif coordinate == "absolute":
                seq = self.reference[start:end] 
        
        elif seq != "":
            if self.reference.count(seq) > 1:
                Warning("The sequence were found two times in the referece sequence") 
            elif self.reference.count(seq) == 0:
                print(self.reference.translate(str.maketrans("ATGC","TACG"))[::-1])
                if self.reference.translate(str.maketrans("ATGC","TACG"))[::-1].count(seq) == 1:
                    self.reference = self.reference.translate(str.maketrans("ATGC","TACG"))[::-1]
                else:
                    raise ValueError("The target sequence is not exist in the reference sequence") 

            if coordinate == "relative":
                start = self.reference.find(seq)  - self.zero_position
                end   = self.reference.find(seq) + len(seq) - self.zero_position

            elif coordinate == "absolute":
                start = self.reference.find(seq)
                end   = self.reference.find(seq) + len(seq) 

        target = {"name":name, "seq":seq, "start":start, "end":end, "color":color, "coordinate":coordinate}
        self.target_dict[name] = target
    
    def merge(self, read1="", read2="", max_overlap=300, allow_outies=True, max_mismatch_density=0.3):        
        print(os.listdir()) 
        if read1 == "":
            R1_dict = read_fastq(self.read1) 
        else:
            R1_dict = read_fastq(read1) 
        if read2 == "":
            R2_dict = read_fastq(self.read2)
        else:
            R2_dict = read_fastq(read2)
            
        if allow_outies == True:
            allow_outies = "--allow-outies"
        else:
            pass

        with open("R1_tmp.fastq","w") as R1_tmp, open("R2_tmp.fastq","w") as R2_tmp:
            for key in list(set(R1_dict.keys()) & set(R2_dict.keys())):
                for i in range(4):
                    R1_tmp.write(R1_dict[key][i] + "\n") 
                    R2_tmp.write(R2_dict[key][i] + "\n") 
        
        os.system("flash --max-overlap {} {} --max-mismatch-density {} {} {} > merge_result.txt".format(max_overlap, allow_outies, max_mismatch_density, "R1_tmp.fastq","R2_tmp.fastq"))
        os.system("mv out.extendedFrags.fastq merged.fastq")
        os.system("cat merged.fastq | awk \'NR % 4 == 1 {{split($0, a, \" \"); gsub(\":\",\"_\",a[1]); print \">\" substr(a[1],2)}} NR % 4 == 2 {{print $0}}\' > merged.fasta")
        os.system("mv out.notCombined_1.fastq notMerged_1.fastq")
        os.system("mv out.notCombined_2.fastq notMerged_2.fastq")
        os.system("rm R1_tmp.fastq R2_tmp.fastq") 
   
    def set_readpath(self,read="merged",path=None):
        if path == None:
            path = os.getcwd()
        self.__read_path__ = path + "/{}_glaln.csv".format(read)
  
    def parse_alignment_result(self, read, start=0, end=0, coordinate="relative"):
        if start == 0 and end == 0:
            start = 0 
            end   = len(self.reference)
            coordinate = "absolute"
        fastq_dict    = read_fastq("{}".format(read))    
        fasta_dict    = collections.OrderedDict() 
        identity_list = [] 
        pair_seq      = [] 
        pair_id       = [] 
        if coordinate == "relative":
            start = start + self.zero_position
            end   = end + self.zero_position
        else: 
            pass
        path = os.getcwd()
        with open(path + "/{}_glaln.csv".format(read.split("/")[-1].replace(".fastq","")),"w") as o:
            count = 0
            for line in open("{}.needle".format(read.split("/")[-1].replace(".fastq",""))):
                if line[0] == ">":
                    if count > 0:
                        pair_id.append(seqid[1:])
                        pair_seq.append(fasta_dict[seqid])

                    if count % 2 == 0 and count > 0:
                        identity = (len(pair_seq[0]) - Levenshtein.hamming(pair_seq[0],pair_seq[1])) * 1.0 / len(pair_seq[0].replace("-",""))
                        identity_list.append(identity * 100)
                        
                        reference = pair_seq[0] 
                        new_reference = reference.rstrip("-") 
                        new_reference = new_reference.lstrip("-") 
                        l = reference.find(new_reference) 
                        r = l + len(new_reference) 
                        reference = new_reference
                        
                        index = 0 
                        new_qualities = [] 
                        qualities = [ord(asc) - 33 for asc in fastq_dict[pair_id[1]][3]]
                        for q in pair_seq[1]:
                            if q == "-":
                                new_qualities.append(100) 
                            else:
                                new_qualities.append(qualities[index]) 
                                index += 1 
                        qualities = new_qualities

                        query     = pair_seq[1][l:r]
                        qualities = qualities[l:r]
                        new_query = ""
                        index     = 0     
                        insert_dict = collections.defaultdict(str) 
                        min_quality = 100
                        for r, q, quality in  zip(reference, query, qualities):
                            if r == "-":
                                insert_dict[index] += q
                            else:
                                new_query += q
                                index += 1
                            if start <= index <= end and quality < min_quality:
                                min_quality = quality
                        
                        insert_cigar = ""       
                        for key in insert_dict:
                            insert_cigar += str(key) 
                            insert_cigar += "".join(insert_dict[key]) 
                        
                        o.write(",".join(list(map(str,[pair_id[1],insert_cigar,new_query,identity,min_quality]))) + "\n") 
                        pair_seq = [] 
                        pair_id  = [] 
                    
                    count += 1 
                    seqid = line.rstrip().replace(":","_").split(" ")[0] 
                    fasta_dict[seqid] = ""

                elif line[0] == "#" or line[0] == "\n":
                        pass     
                else:
                    fasta_dict[seqid] += line.rstrip() 
        self.__read_path__ = path + "/{}_glaln.csv".format(read.split("/")[-1].replace(".fastq",""))
        print(path)  
   
    def alignment(self, read="merged.fastq", num_threads=1, strand="forward"): 
        if "reference" not in os.listdir("."):
            os.mkdir("reference")
        else:
            pass
        
        with open("reference/reference.fasta","w") as o:
            o.write(">reference\n")
            o.write(self.reference + "\n") 

        if read.replace(".fastq",".fasta") not in os.listdir("./"):
            os.system("cat {} | awk \'NR % 4 == 1 {{split($0, a, \" \"); gsub(\":\",\"_\",a[1]); print \">\" substr(a[1],2)}} NR % 4 == 2 {{print $0}}\' > {}".format(read,read.replace("fastq","fasta")))

        if num_threads == 1:
            if strand == "reverse":
                os.system("needleall ./reference/reference.fasta {} -sreverse2 True -gapopen 20.0 -gapextend 2.0 -endopen 0.0 -endextend 0.0 -aformat markx3 -outfile {}".format(read.replace(".fastq",".fasta"),read.replace(".fastq",".needle"))) 
            else:
                os.system("needleall ./reference/reference.fasta {} -gapopen 20.0 -gapextend 2.0 -endopen 0.0 -endextend 0.0 -aformat markx3 -outfile {}".format(read.replace(".fastq",".fasta"),read.replace(".fastq",".needle"))) 
        
        elif num_threads > 1:
            commands  = []
            num_lines = sum([1 for line in open("{}".format(read.replace("fastq","fasta")))])
            os.system("split -l {} {}.fasta tmp_".format(2 * (int((num_lines/2)/num_threads) + 1), target))
            for fasta in [f for f in os.listdir(".") if "tmp_" in f]:
                if read == "R2":
                    commands.append("needleall ./reference/reference.fasta {} -sreverse2 True -gapopen 20.0 -gapextend 2.0 -endopen 0.0 -endextend 0.0 -aformat markx3 -outfile {}.needle".format(fasta,fasta)) 
                else:
                    commands.append("needleall ./reference/reference.fasta {} -gapopen 20.0 -gapextend 2.0 -endopen 0.0 -endextend 0.0 -aformat markx3 -outfile {}.needle".format(fasta,fasta))
            
            os.system("&".join(commands))
            os.system("cat tmp_*.needle > {}.needle".format(read))
            os.system("rm tmp*")
        
    def calc_sub_spectrum(self, min_identity=0.8, min_quality=0):
        sub_list       = [0 for x in range(len(self.reference))]
        occupancy_list = [{"A":0,"T":0,"G":0,"C":0,"N":0} for x in range(len(self.reference))] 
        with open(self.__read_path__) as f:
           count = 0     
           for line in f:
                line     = line.rstrip().split(",") 
                query    = line[2] 
                identity = float(line[3]) 
                quality  = int(line[4])
                if identity >= min_identity and quality >= min_quality:
                    for index, (r,q) in enumerate(zip(self.reference, query)):
                        if r != q and q != "-":
                            sub_list[index]          += 1 
                            occupancy_list[index][q] += 1 
                    count += 1 
        
        for key in range(len(sub_list)):
            sub_list[key] /= (count * 1.0)
        
        for key in range(len(occupancy_list)) :
            sum_value = sum(occupancy_list[key].values())+1 
            for nucl in occupancy_list[key]:
                occupancy_list[key][nucl] /= (sum_value * 1.0)
        
        self.__subspec_dict__[(min_identity,min_quality)] = (sub_list, occupancy_list) 
        return sub_list, occupancy_list,  count
                  
    def calc_coediting_matrix(self, start, end, coordinate="relative", min_identity=0.8, min_quality=0, mutation_type="all"):
        #mode="IoU", "CP", "DP"
        if coordinate == "relative":
            astart = start + self.zero_position
            aend   = end + self.zero_position
        else:
            astart = start 
            aend   = end 
                
        if mutation_type == "all": 
            requirement = lambda x: True
        else: 
            requirement = lambda x: x in mutation_type

        snedit_list, occupancy_list = self.get_sub_spectrum(astart, aend, coordinate="absolute", min_identity=min_identity, min_quality=min_quality)
        diedit_matrix  = [[0] * (aend+1 - astart) for i in range(astart,aend+1)]
        union_matrix   = [[0] * (aend+1 - astart) for i in range(astart,aend+1)]
        
        mut_dict = {} 
        if mutation_type == "all":
            for i1, p1 in enumerate(snedit_list):
                for i2, p2 in enumerate(snedit_list): 
                    union_matrix[i1][i2] = p1 + p2
            
            for p,m in enumerate(self.reference[astart:aend+1]):
                for q,n in enumerate(self.reference[astart:aend+1]):
                    mut_dict[(p+start,q+start),(m,n)] = {"N": {"A":0, "T":0, "G":0, "C":0, "N":0}, "A": {"A":0, "T":0, "G":0, "C":0, "N":0}, "T":{"A":0, "T":0, "G":0, "C":0, "N":0}, "G":{"A":0, "T":0, "G":0, "C":0, "N":0}, "C":{"A":0, "T":0, "G":0, "C":0, "N":0}}
             

            ref = self.reference[astart:aend+1]
            with open(self.__read_path__) as f:
               count = 0     
               for line in f:
                    line      = line.rstrip().split(",") 
                    query     = line[2] 
                    identity  = float(line[3]) 
                    quality   = int(line[4])
                    if identity >= min_identity and quality >= min_quality:
                        for i1, (r1,q1) in enumerate(zip(ref, query[astart:aend+1])):
                            for i2, (r2,q2) in enumerate(zip(ref, query[astart:aend+1])):
                                if q1 != "-" and q2 != "-": 
                                    mut_dict[(i1+start,i2+start),(r1,r2)][q1][q2] += 1
                                    if (r1 != q1) and (r2 != q2):
                                        diedit_matrix[i1][i2] += 1
                        count += 1 
            
            for key1 in mut_dict:
                for key2 in mut_dict[key1]:
                    for key3 in mut_dict[key1][key2]:
                        mut_dict[key1][key2][key3] /= count
                        
            diedit_matrix = np.array(diedit_matrix) / (count * 1.0) 
            IoU_matrix    = np.array(diedit_matrix) / (np.array(union_matrix) - np.array(diedit_matrix))
            coedit_matrix = np.array(diedit_matrix) / np.array(snedit_list) 
        
        else:
            with open(self.__read_path__) as f:
               count = 0     
               single_array = [[0] * len(self.reference)] 
               ref = self.reference[astart:aend+1]
               for line in f:
                    line      = line.rstrip().split(",") 
                    query     = line[2] 
                    identity  = float(line[3]) 
                    quality   = int(line[4])
                    if identity >= min_identity and quality >= min_quality:
                        for index1, (r1,q1) in enumerate(zip(ref, query[astart:aend+1])):
                            if ((r1 != q1) and requirement((r1,q1))):
                                single_array[index1] += 1
                            for index2, (r2,q2) in enumerate(zip(ref, query[astart:aend+1])):
                                if q1 != "-" and q2 != "-": 
                                    if (r1 != q1 and requirement((r1,q1))) or (r2 != q2 and requirement((r2,q2))):
                                        if (r1 != q1 and requirement((r1,q1))) and (r2 != q2 and requirement((r2,q2))):
                                            diedit_matrix[i1][i2] += 1
                                        union_matrix[i1][i2] += 1 
                                        union_matrix[i2][i1] += 1
                        count += 1

            IoU_matrix    = np.array(diedit_matrix) / (np.array(union_matrix) - np.array(diedit_matrix))  
            coedit_matrix = np.array(diedit_matrix) / np.array(single_array)
            diedit_matrix = np.array(diedit_matrix) / (count * 1.0) 
        
        self.__comatrixes_dict__[(start,end,coordinate,min_identity,min_quality,mutation_type)] = {"IU":IoU_matrix, "CP":coedit_matrix, "DP":diedit_matrix, "info":mut_dict}  
        

    def get_reference(self, start, end, coordinate="relative"): 
        if coordinate == "relative":
            astart = start + self.zero_position
            aend   = end + self.zero_position
        else:
            astart = start 
            aend   = end 
        return self.reference[astart:aend+1] 
    
    def get_readcount(self, min_identity=0.8, min_quality=0, coordinate="relative"):
        if (min_identity,min_quality)  not in self.__readcount_dict__.keys():
            with open(self.__read_path__) as f:
               count = 0     
               for line in f:
                    line     = line.rstrip().split(",") 
                    query    = line[2] 
                    identity = float(line[3]) 
                    quality  = int(line[4])
                    if identity >= min_identity and quality >= min_quality:
                        count += 1  
            self.__readcount_dict__[(min_identity,min_quality)] = count
        return self.__readcount_dict__[(min_identity,min_quality)]
    
     
    def get_sub_spectrum(self, start, end, coordinate="relative", min_identity=0.8, min_quality=0, to_csv = False, output_path=None):
        if (min_identity,min_quality) not in self.__subspec_dict__:
            self.calc_sub_spectrum(min_identity=0.8, min_quality=0) 
        else:
            pass
        
        if coordinate == "relative":
            astart = start + self.zero_position
            aend   = end + self.zero_position
        
        else:
            astart = start 
            aend   = end
        
        if to_csv == True:
            sub_spec, occupancy = self.__subspec_dict__[(min_identity,min_quality)][0][astart:aend+1], self.__subspec_dict__[(min_identity,min_quality)][1][astart:aend+1]
            ref    = self.reference[self.zero_position - start: self.zero_position + end + 1]
            header = ["Reference","position","Substitution rate","A","T","G","C","N"]
            rows   = []
            pos    = start
            for r,v in zip(ref,sub_spec):
                rows.append([r,pos,v])
                pos += 1
            for i in range(len(rows)):
                rows[i].append(occupancy[i]["A"])
                rows[i].append(occupancy[i]["T"])
                rows[i].append(occupancy[i]["G"])
                rows[i].append(occupancy[i]["C"])
                rows[i].append(occupancy[i]["N"])
            
            if output == None:
                output = self.__output_path__ + "/substitution_spectrum.csv"
            else:
                output = output_path + "/substitution_spectrum.csv"

            with open(output,"w") as o:
                    o.write(",".join(header) + "\n") 
                    for row in rows:
                        o.write(",".join(list(map(str,row)))  + "\n")

        return self.__subspec_dict__[(min_identity,min_quality)][0][astart:aend+1], self.__subspec_dict__[(min_identity,min_quality)][1][astart:aend+1]
    def calc_mut_pattern(self, start, end, coordinate="relative", min_identity=0.8, min_quality=0): 
        if coordinate == "relative":
            astart = start + self.zero_position
            aend   = end + self.zero_position
        else:
            astart = start 
            aend   = end 
        
        mut_pattern_dict = collections.defaultdict(int)
        with open(self.__read_path__) as f:
           count = 0     
           for line in f:
                line     = line.rstrip().split(",") 
                query    = line[2] 
                identity = float(line[3]) 
                quality  = int(line[4])
                if identity >= min_identity and quality >= min_quality:
                    mut_pattern_dict[(query[astart:aend+1],line[1])] += 1 
                    count += 1 

        new_mut_pattern_dict = collections.OrderedDict() 
        keys = list(mut_pattern_dict.keys()) 
        keys.sort(key=lambda x:mut_pattern_dict[x] * -1.0)
        for key in keys:
            new_mut_pattern_dict[key] = (mut_pattern_dict[key] / (count * 1.0), mut_pattern_dict[key]) 
        self.__mutpattern_dict__[(start,end,coordinate,min_identity,min_quality)] = new_mut_pattern_dict 
        return new_mut_pattern_dict
    

        
    def get_coediting_matrix(self, start, end, coordinate="relative", min_identity=0.8, min_quality=0, mutation_type="all"):
        if (start,end,coordinate,min_identity,min_quality,mutation_type) not in self.__comatrixes_dict__:
            self.calc_coediting_matrix(start,end,coordinate,min_identity,min_quality,mutation_type) 
        else:
            pass
        dp = self.__comatrixes_dict__[(start,end,coordinate,min_identity,min_quality,mutation_type)]["DP"]
        cp = self.__comatrixes_dict__[(start,end,coordinate,min_identity,min_quality,mutation_type)]["CP"]
        iu = self.__comatrixes_dict__[(start,end,coordinate,min_identity,min_quality,mutation_type)]["IU"]
        return {"DP":dp,"CP":cp,"IU":iu}
    
    def get_mut_pattern(self, start, end, coordinate="relative", min_identity=0.8, min_quality=0):
        if (start,end,coordinate,min_identity,min_quality) not in self.__mutpattern_dict__:
            self.calc_mut_pattern(start, end, coordinate, min_identity, min_quality)
        return self.__mutpattern_dict__[(start,end,coordinate,min_identity,min_quality)]

def main(config):
    samples = [] 
    for section in config.sections():
        sample = SAMPLE()
        sample.name = section
        try:
            os.mkdir(config.get(section, "output_dir"))        
        except:
            pass 
        os.chdir(config.get(section, "output_dir")) 
        sample.__output_path__ = os.getcwd() 
        
        #Set reference sequence
        ref    = config.get(section, "reference")   
        sample.set_reference(ref) 
    
        #Set zeroposition
        p = config.get(section, "zero_position")
        p = [sp.split(":") for sp in p.split(",")]
        p = dict(p) 
        sample.set_zero_position(int(p["position"]),seq=p["seq"]) 
        
        #Set elements
        count = 1
        for option in config.options(section):
            if option == "e" + str(count):
                e = config.get(section, option)
                e = [se.split(":") for se in e.split(",")]
                e = dict(e) 
                if e["seq"] == "None":
                    sample.add_element(start=int(e["start"]),end=int(e["end"]),name=e["name"],coordinate=e["coordinate"])   
                else:
                    sample.add_element(seq=e["seq"],name=e["name"],coordinate=e["coordinate"])   
                count += 1
    
        #Alignemnt 
        count = 1
        reads = []
        for option in config.options(section):
            if option == "read" + str(count):
                read = config.get(section, option)
                read = [sr.split(":") for sr in read.split(",")]
                print(read) 
                read = dict(read) 
                reads.append(read) 
                count += 1 

        start, end, coordinate = config.get(section,"range").split(",")
        start = int(start) 
        end   = int(end) 
        if len(reads) == 2:
            sample.merge(read1=reads[0]["path"],read2=reads[1]["path"]) 
            sample.alignment("merged.fastq") 
            sample.parse_alignment_result(read="merged.fastq", start=start, end=end, coordinate=coordinate)
        else:
            sample.alignment(reads[0]["path"],strand=reads[0]["strand"]) 
            sample.parse_alignment_result(read=reads[0]["path"], start=start, end=end, coordinate=coordinate)

        matrix_dict = sample.get_coediting_matrix(start=start, end=end) 
        count = sample.get_readcount() 
        with open("sample.pickle","wb") as o:
            pickle.dump(sample,o) 

        os.chdir("../")
        samples.append(sample) 
    return samples
    
if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-c", "--config", type=str, default=False)
    args   = p.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    main(config) 
