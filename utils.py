import yaml

def read_yaml_config(filename):
    with open(filename,'r', encoding='utf-8') as file:
        yaml_config  = yaml.load(file, Loader=yaml.FullLoader)
    return yaml_config

