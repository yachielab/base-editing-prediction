import os 
import sys
import pickle

"""
Model name
Target sequence
Editing outcome
Editing frequency
Target sequence
Start position
End position
Editing_outcome
Editing_frequency
Data type
Conditional base transition
Target base transition
Probability
"""

if __name__ == "__main__":
    pickles = [fname for fname in os.listdir("./") if "mutation_pattern_dict.pickle" in fname]
    for pic in pickles:
        with open(pic.replace("_mutation_pattern_dict.pickle",".csv"),"w") as o:
            m_rate_dict    = pickle.load(open(pic.replace("mutation_pattern_dict.pickle","mutation_rate_dict.pickle"),"rb"))
            m_pattern_dict = pickle.load(open(pic,"rb"))
            c_rate_dict    = pickle.load(open("/Users/hideto/Downloads/pickles/EGFP_mutation_rate_dict.pickle","rb"))
            c_pattern_dict = pickle.load(open("/Users/hideto/Downloads/pickles/EGFP_mutation_pattern_dict.pickle","rb"))

            o.write("#Model Name:{}\n".format(pic.replace("_mutation_pattern_dict.pickle","")))
            o.write("#TP:Target base transition probability\n")
            o.write("#CTP:Conditional base transition probability\n")
            o.write("Data type,Conditional transition,Target transition,Probability\n")
            for key in m_rate_dict: 
                for key2 in m_rate_dict[key]:
                    if key2 != "N":
                        value = m_rate_dict[key][key2]
                        if value == None: 
                            value = "N.A."
                        else:
                            value = m_rate_dict[key][key2] - c_rate_dict[key][key2]  
                            value = value if value > 0 else 0  
                        print("TP", "", "{}:{}>{}".format(key[0],key[1],key2), value, sep=",", file=o)
            
            for key in m_pattern_dict: 
                for key2 in m_pattern_dict[key]:
                    for key3 in m_pattern_dict[key][key2]:
                        value = m_pattern_dict[key][key2][key3]
                        if key2 != "N" and key3 != "N":
                            if value == None:
                                value = "N.A."
                            else:
                                value = m_pattern_dict[key][key2][key3]
                                value = value if value > 0 else 0
                            print("CTP", "{}:{}>{}".format(key[0][0],key[1][0],key2), "{}:{}>{}".format(key[0][1],key[1][1],key3), value, sep=",", file=o)
                        else:
                            pass 
