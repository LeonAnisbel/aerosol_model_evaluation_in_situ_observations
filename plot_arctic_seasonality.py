import pandas as pd
import global_vars
import matplotlib.pyplot as plt

import read_data_functions


def plot_MC15():
    data = pd.read_pickle(f'pd_files/MH15_conc_{global_vars.exp_name}.pkl')
    print(data.head())
    fig, ax = plt.subplots(2,1)
    ax.flatten()
    obs_MH_15, _, _, _ = read_data_functions.read_data()
    obs_MH_15['seasalt_model'] = data['conc_mod_ss'].values
    obs_MH_15['MOA_model'] = data['conc_mod_tot'].values
    obs_MH_15.rename(columns={'seasalt':'seasalt_obs', 'MOA':'MOA_obs'}, inplace=True)

    data_ss = obs_MH_15[['seasalt_model', 'seasalt_obs']].copy(deep=True)
    data_moa = obs_MH_15[['MOA_model', 'MOA_obs']].copy(deep=True)


    data_ss.plot.bar(ax = ax[0])
    data_moa.plot.bar(ax = ax[1])

    fig.tight_layout()
    plt.savefig(f'plots/MH15_conc_{global_vars.exp_name}.png')

