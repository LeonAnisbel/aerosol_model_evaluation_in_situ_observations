import os
from utils_functions import global_vars

os.system('conda env_aer_eval create -f env.yml')
os.system('conda activate env_aer_eval')


try:
    os.makedirs('plots')
except OSError:
    pass

try:
    os.makedirs('outputs')
except OSError:
    pass

