import os 
import sys
import pickle

if __name__ == "__main__":
    pickles = [fname for fname in os.listdir("./") if "mutation_pattern_dict.pickle" in fname]
    for pic in pickles:
        with open(pic.replace("_mutation_pattern_dict.pickle",".csv"),"w") as o:
            m_rate_dict   = pickle.load(open(pic.replace("mutation_pattern_dict.pickle","mutation_rate_dict.pickle"),"rb"))
            m_pattern_dict = pickle.load(open(pic,"rb"))
            o.write("#{}\n".format(pic.replace("mutation_pattern_dict.pickle","")))
            o.write("#TP:Transition probability\n")
            o.write("#CTP:Conditional transition probability\n")
            o.write("Probability type,Conditional base transition,Target base transition,Base transition probability\n")
            for key in m_rate_dict: 
                for key2 in m_rate_dict[key]:
                    if key2 != "N":
                        value = m_rate_dict[key][key2] 
                        if value == None: 
                            value = "N.A."
                        print("TP", "", "{}:{}>{}".format(key[0],key[1],key2), value, sep=",", file=o)
            
            for key in m_pattern_dict: 
                for key2 in m_pattern_dict[key]:
                    for key3 in m_pattern_dict[key][key2]:
                        if key2 != "N" and key3 != "N":
                            value = m_pattern_dict[key][key2][key3]
                            if value == None:
                                value = "N.A."
                            print("CTP", "{}:{}>{}".format(key[0][0],key[1][0],key2), "{}:{}>{}".format(key[0][1],key[1][1],key3), value, sep=",", file=o)
                        else:
                            pass 
