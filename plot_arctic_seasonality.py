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


def plot_all_arctic_stations():
    data = pd.read_pickle(f'pd_files/all_arctic_stations_conc_{global_vars.exp_name}.pkl')
    print(data.head())
    data_mo_all_stations, dict_metadata = read_data_functions.read_PMOA_all_stations()
    stat_list = list(dict_metadata.keys())
    data_mo_all_stations, dict_metadata = read_data_functions.read_PMOA_all_stations()
    data.drop('conc_obs_tot', axis=1, inplace=True)
    data['conc_obs_tot'] = data_mo_all_stations['PBOA_ug_m3'].values
    data.set_index('months', inplace=True)


    fig, ax = plt.subplots(len(stat_list),1, figsize=(12, 20))
    ax.flatten()
    for idx, sta in enumerate(stat_list):
        #print(sta, data[['ID']==sta])
        data[data['ID']==sta].plot(ax = ax[idx])
        ax[idx].set_title(sta, loc='center')
    fig.tight_layout()
    plt.savefig(f'plots/all_arctic_stations_conc_bar_{global_vars.exp_name}.png')

    fig, ax = plt.subplots(1,1)
    data_mod = data[['ID', 'conc_mod_tot']]
    data_obs = data[['ID', 'conc_obs_tot']]

    m = data_mod.plot(style='-', ax = ax) #'.-'
    o = data_obs.plot(style='--', ax = ax) #'.-'
    fig.legend()
    plt.savefig(f'plots/all_arctic_stations_conc_{global_vars.exp_name}.png')



if __name__ == '__main__':
    plot_all_arctic_stations()