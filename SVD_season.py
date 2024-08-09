##### Import packges
import numpy as np
import pandas as pd
import glob,os
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from sklearn.metrics import mean_squared_error
import math

if __name__ == '__main__':
    # data paths
    data_dir = "pd_files/"

    try:
        os.mkdir('plots')
    except OSError:

        pass

    mix_omf_conc = []
    # for i,v in enumerate(var):

    conc_init = pd.read_pickle(f'{data_dir}prot_conc_SVD18.pkl')
    fig, ax = plt.subplots(figsize=(6, 5))
    print(conc_init)
    conc_pr = conc_init[['Start Date/Time','conc_mod_pro', 'conc_obs_prot_sub', 'conc_obs_ss', 'conc_mod_ss']]   
    #conc_gr = conc.time.dt.floor('1M') # DDTHH:MM --> DD+1T00:00; alternatively, use floor(): DDTHH:MM --> DDT00:00
    print(conc_pr)
    conc_pr = conc_pr.rename(columns={ 'conc_mod_pro':'Model_pr', 'conc_obs_prot_sub': 'Observation_pr','conc_obs_ss': 'Obs_ss', 'conc_mod_ss':'Mod_ss'})

    conc_gr = conc_pr.groupby(['Start Date/Time'], as_index=False).mean()
    print(list(conc_gr.columns), 'columnsss')
    plt.scatter(conc_gr['Start Date/Time'], conc_gr['Model_pr'], color='r', label='Model')
    plt.scatter(conc_gr['Start Date/Time'], conc_gr['Observation_pr'], color='b', label='Observation')
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_ylabel('DCAA | CAA')
    plt.legend()
    ax.set_yscale('log')

    plt.tight_layout()
    plt.savefig('plots/prot_conc_SVD18.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(6, 5))

    plt.scatter(conc_gr['Start Date/Time'], conc_gr['Mod_ss'], color='r', label='Model')
    plt.scatter(conc_gr['Start Date/Time'], conc_gr['Obs_ss'], color='b', label='Observation')
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_ylabel('SS')
    ax.set_yscale('log')

    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/ss_conc_SVD18.png')
