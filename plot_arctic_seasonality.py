import pandas as pd
import global_vars
import matplotlib.pyplot as plt

import read_data_functions


def plot_MC15(yr):
    data = pd.read_pickle(f'pd_files/MH{yr}_conc_{global_vars.exp_name}.pkl')
    print(data.head())
    fig, ax = plt.subplots(2,1)
    ax.flatten()
    obs_MH_15, _, _, _ = read_data_functions.read_data(yr)
    obs_MH_15['seasalt_model'] = data['conc_mod_ss'].values
    obs_MH_15['MOA_model'] = data['conc_mod_tot'].values
    obs_MH_15.rename(columns={'seasalt':'seasalt_obs', 'PMOA':'MOA_obs'}, inplace=True)

    data_ss = obs_MH_15[['seasalt_model', 'seasalt_obs']].copy(deep=True)
    data_moa = obs_MH_15[['MOA_model', 'MOA_obs']].copy(deep=True)


    data_ss.plot(ax = ax[0])
    data_moa.plot(ax = ax[1])

    fig.tight_layout()
    plt.savefig(f'plots/MH{yr}_conc_{global_vars.exp_name}.png')


def plot_all_arctic_stations():
    data = pd.read_pickle(f'pd_files/all_arctic_stations_conc_{global_vars.exp_name}.pkl')
    print(data.head())
    data_mo_all_stations, dict_metadata, data_sel_std = read_data_functions.read_PMOA_all_stations()
    stat_list = list(dict_metadata.keys())

    #data.drop('conc_obs_tot', axis=1, inplace=True)
    #data['conc_obs_tot'] = data_mo_all_stations['PBOA_ug_m3'].values

    data_sel_yr = ((data.groupby(['months', 'ID'])[[ 'conc_obs_tot',  'conc_mod_tot']])
                   .median())
    data_std = ((data.groupby(['months', 'ID'])[[ 'conc_obs_tot',  'conc_mod_tot']])
                   .std())
    data_std['obs_std_fill_min'] = [val-std for val,std in zip(data_sel_yr['conc_obs_tot'].values, data_std['conc_obs_tot'].values)]
    data_std['obs_std_fill_max'] = [val+std for val,std in zip(data_sel_yr['conc_obs_tot'].values, data_std['conc_obs_tot'].values)]
    data_sel_yr = data_sel_yr.reset_index()
    data_std=data_std.reset_index()
    data_sel_yr.set_index('months', inplace=True)
    data_std.set_index('months', inplace=True)

    fig, ax = plt.subplots(len(stat_list),1, figsize=(12, 20))
    ax.flatten()
    for idx, sta in enumerate(stat_list):
        # print(sta, data[['ID']==sta])
        data_sel_yr[data_sel_yr['ID'] == sta].plot(ax=ax[idx])
        data_std_sta = data_std[data_std['ID'] == sta]
        ax[idx].fill_between(data_sel_yr[data_sel_yr['ID'] == sta].index,
                             data_std_sta['obs_std_fill_min'].values,
                             data_std_sta['obs_std_fill_max'].values,
                             facecolor='lightblue',
                             alpha=0.3)
        ax[idx].set_title(sta, loc='center')
    fig.tight_layout()
    plt.savefig(f'plots/all_arctic_stations_conc_bar_yearly_{global_vars.exp_name}.png')
    plt.close()

    date_list = [str(i) + '-' + str(j) for i, j in zip(data['years'].values, data['months'].values)]
    print(date_list)
    data['date'] = date_list
    data.drop('years', axis=1, inplace=True)
    data.drop('months', axis=1, inplace=True)
    data.set_index('date', inplace=True)

    fig, ax = plt.subplots(len(stat_list), 1, figsize=(12, 20))
    ax.flatten()
    for idx, sta in enumerate(stat_list):
        # print(sta, data[['ID']==sta])
        data[data['ID'] == sta].plot(ax=ax[idx])
        ax[idx].set_title(sta, loc='center')
    fig.tight_layout()
    plt.savefig(f'plots/all_arctic_stations_conc_bar_{global_vars.exp_name}.png')
    plt.close()

    fig, ax = plt.subplots(1, 1)
    obs_leg = []
    mod_leg = []
    for idx, sta in enumerate(stat_list):
        data_mod = data[['ID', 'conc_mod_tot']]
        data_obs = data[['ID', 'conc_obs_tot']]
        m = data_mod[data['ID'] == sta].plot(style='-', ax=ax)  # '.-'
        o = data_obs[data['ID'] == sta].plot(style='--', ax=ax)  # '.-'
        obs_leg.append(o)
        mod_leg.append(m)

    fig.legend(handles=obs_leg, labels=stat_list)
    plt.savefig(f'plots/all_arctic_stations_conc_{global_vars.exp_name}.png')



if __name__ == '__main__':
    plot_all_arctic_stations()
    year = '2018'
    print('start seasonality 2018')
    plot_MC15(year)



