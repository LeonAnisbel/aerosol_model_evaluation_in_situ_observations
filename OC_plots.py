##### Import packges
import numpy as np
import pandas as pd
import glob, os
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from sklearn.metrics import mean_squared_error
import math

import global_vars

def get_stat(obs, mod):
    # stadistics
    MSE = mean_squared_error(obs, mod)
    RMSE = math.sqrt(MSE)
    diff = (mod - obs)
    mean_bias = np.nanmean(diff)
    NMB = np.nanmean(diff/obs)

    R = np.corrcoef(obs, mod)
    return RMSE, mean_bias,NMB, R[0][1]


def plot_fit(conc_mod_obs, obs_col_na, mod_col_na, title, xli, yli, loc):
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=conc_mod_obs, x=obs_col_na, y=mod_col_na, hue="Measurements")
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0, label='1:1 line')
    # print(lims)

    # Fitting poly data
    obs = np.array(conc_mod_obs[obs_col_na].values)
    mod = np.array(conc_mod_obs[mod_col_na].values)

    RMSE, mean_bias, NMB, R = get_stat(obs,mod)
    formatted_pvalues = (f'RMSE= {np.round(RMSE, 2)} \n bias= {np.round(mean_bias, 2)} '
                         f'\n R= {np.round(R, 2)}')
    print(formatted_pvalues)
    ax.text(loc, 0.002, formatted_pvalues, fontsize='14')

    # Customizing axes
    # ax.set_title(title, fontsize='16')

    ax.tick_params(axis='both', labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.xaxis.get_label().set_fontsize(14)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(xli)
    ax.set_ylim(yli)
    plt.legend(loc='upper left')
    ax.grid(linestyle='--', linewidth=0.4)
    fig.tight_layout()

    plt.savefig(f'plots/conc_{title}_{global_vars.exp_name}.png',dpi = 300)
    plt.close()


def plot_oc(conc_ss, title):
    old_cols = 'Aerosol Concentration (µg m$^{-3}$)'

    mod_vals, obs_vals, meas_vals = [], [], []
    list_stat = ['CVAO', 'WAP']
    for m in list_stat:
        conc_loc = conc_ss[conc_ss['Measurements'] == m]
        mod_oc = np.array(conc_loc[conc_loc['']=='Model'][old_cols].values)
        obs_oc = np.array(conc_loc[conc_loc['']=='Observation'][old_cols].values)

        RMSE, mean_bias,NMB, R = get_stat(obs_oc, mod_oc)

        print('MOdel ', mod_oc.mean(), mod_oc.std())
        print('Observations ', obs_oc.mean(), obs_oc.std())
        #
        print('Statistics: ', 'RMSE:', RMSE, 'bias: ', mean_bias, 'NMB: ', NMB, 'R:', R, '\n')

        meas_vals.append(conc_loc['Measurements'].values)
        obs_vals.append(obs_oc)
        mod_vals.append(mod_oc)

    print(len(conc_ss[conc_ss['']=='Observation'][old_cols].values),
          len(conc_ss[conc_ss['']=='Model'][old_cols].values),
          len(conc_ss['Measurements'].values[:len(conc_ss[conc_ss['']=='Model'][old_cols].values)]))

    conc_pd = pd.DataFrame(data={'Measurements': conc_ss['Measurements'].values[:len(conc_ss[conc_ss['']=='Model'][old_cols].values)],
                                    f'Model {title} concentration (µg m$^{-3}$)': conc_ss[conc_ss['']=='Model'][old_cols].values,
                                    'Observations OM concentration (µg m$^{-3}$)': conc_ss[conc_ss['']=='Observation'][old_cols].values*2})
    obs_col_na = 'Observations OM concentration (µg m$^{-3}$)'
    mod_col_na = f'Model {title} concentration (µg m$^{-3}$)'
    print(conc_pd)
    plot_fit(conc_pd,
             obs_col_na, mod_col_na,
             title,
             [1e-5, 1e1],
             [1e-5, 1e1],
             0.000025)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # data paths
    data_dir = "pd_files/"

    try:
        os.mkdir('plots')
    except OSError:

        pass

    mix_omf_conc = []
    all_conc = []
    conc = pd.read_pickle(f'{data_dir}oc_conc_{global_vars.exp_name}.pkl')

    print(conc)

    plot_oc(conc, 'OC')

