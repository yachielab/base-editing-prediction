import configparser
import argparse
import learn_cp as lcp
from base_editing_prediction import * 

if __name__ == "__main__":
    start = -25
    end   = -1
    p = argparse.ArgumentParser()
    p.add_argument("-c", "--config", type=str, default=False)
    args   = p.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    samples = lcp.main(config)
    m_rate_dict, m_pattern_dict = learn(samples)
    validation(samples[0],m_rate_dict,m_pattern_dict,start=start,end=end)
    
