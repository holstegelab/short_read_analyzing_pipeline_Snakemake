import yaml
import pickle

def read_yaml_config(filename):
    with open(filename,'r', encoding='utf-8') as file:
        yaml_config  = yaml.load(file, Loader=yaml.FullLoader)
    return yaml_config

def save(data,dfilename):
    with open(dfilename,'wb') as df:
        pickle.dump(data,df)

def load(dfilename):
    with open(dfilename,'rb') as df:
        res = pickle.load(df)
    return res


