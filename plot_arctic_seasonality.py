import numpy as np
import pandas as pd
import global_vars
import matplotlib.pyplot as plt

import read_data_functions


def plot_MC_monthly_seasonality(yr):
    data = pd.read_pickle(f'pd_files/MH{yr}_conc_{global_vars.exp_name}.pkl')
    fig, ax = plt.subplots(3,1, figsize=(6,15))
    ax.flatten()

    _, _, months, _, _ = read_data_functions.read_data(yr)
    obs_MH_15, _, _, _, obs_MH_15_std = read_data_functions.read_data(yr, monthly=True)

    data['months'] = months
    data_month_mean = []
    data_month_std = []
    for i in np.arange(1,13):
        data_mean = data.loc[(data['months']==i), ('conc_mod_ss', 'conc_mod_tot')].mean()
        data_std = data.loc[(data['months']==i), ('conc_mod_ss', 'conc_mod_tot')].std()

        data_month_mean.append(data_mean.values)
        data_month_std.append(data_std.values)

    ss_vals, moa_vals = [i[0] for i in data_month_mean],  [i[1] for i in data_month_mean]
    ss_std, moa_std = [i[0] for i in data_month_std],  [i[1] for i in data_month_std]

    obs_MH_15['Model (SS)'] = ss_vals
    obs_MH_15['Model (PMOA)'] = moa_vals

    obs_MH_15.rename(columns={'seasalt':'Observations (SS)',
                              'PMOA':'Observations (PMOA)',
                              'OMF': 'Observations (OMF)'}, inplace=True)
    data_ss = obs_MH_15[['Model (SS)', 'Observations (SS)']].copy(deep=True)
    data_moa = obs_MH_15[['Model (PMOA)', 'Observations (PMOA)']].copy(deep=True)

    # Read OMF and add it to plot
    data_omf = pd.read_pickle(f'pd_files/tot_omf_MH.pkl')
    times = pd.to_datetime(data_omf['Start Date/Time'], dayfirst=True)
    data_omf_monthly = data_omf.groupby([times.dt.year, times.dt.month])[['OMF_mod_tot']].mean()
    obs_MH_15['Model (OMF)'] = data_omf_monthly['OMF_mod_tot'].values
    data_omf = obs_MH_15[['Model (OMF)', 'Observations (OMF)']].copy(deep=True)

    def fill_between(axs, data, std, c, plot_std=True):
        if plot_std:
            min, max = ([i-j for i,j in zip(data.values, std)],
                        [i+j for i,j in zip(data.values, std)])
            axs.fill_between(np.arange(0,12),
                                min, max,
                               facecolor=c,
                                 alpha=0.1)

    data_ss.plot(ax = ax[0], color = ['b', 'r'])
    fill_between(ax[0],data_ss['Observations (SS)'], obs_MH_15_std['seasalt'].values, 'r')
    fill_between(ax[0], data_ss['Model (SS)'], ss_std, 'b')
    ax[0].set_title('Sea salt')
    ax[0].set_xlabel(' ')
    ax[0].set_ylabel('Concentration (ug m$^{-3}$)')


    data_moa.plot(ax = ax[1], color = ['b', 'r'])
    fill_between(ax[1], data_moa['Observations (PMOA)'], obs_MH_15_std['PMOA'].values, 'r')
    fill_between(ax[1], data_moa['Model (PMOA)'], moa_std, 'b')
    ax[1].set_xlabel(' ')
    ax[1].set_title('PMOA')
    ax[1].set_ylabel('Concentration (ug m$^{-3}$)')

    data_omf.plot(ax = ax[2], color = ['b', 'r'])
    fill_between(ax[2], data_omf['Observations (OMF)'], obs_MH_15_std['OMF'].values, 'r')
    fill_between(ax[2], data_omf['Model (OMF)'], [], 'b', plot_std=False)
    ax[2].set_xlabel(' ')
    ax[2].set_title('OMF')
    ax[2].set_ylabel('OMF')

#    ax[0].semilogy()

    fig.tight_layout()
    plt.savefig(f'plots/MH{yr}_monthly_conc_{global_vars.exp_name}.png')


def plot_MC_daily_seasonality(yr):
    data = pd.read_pickle(f'pd_files/MH{yr}_conc_{global_vars.exp_name}.pkl')
    fig, ax = plt.subplots(2,1, figsize=(10,5))
    ax.flatten()
    obs_MH_15, _, _, _, _ = read_data_functions.read_data(yr)
    obs_MH_15['Model (SS)'] = data['conc_mod_ss'].values
    obs_MH_15['Model (PMOA)'] = data['conc_mod_tot'].values
    obs_MH_15.rename(columns={'seasalt':'Observations (SS)',
                              'PMOA':'Observations (PMOA)'}, inplace=True)

    data_ss = obs_MH_15[['Model (SS)', 'Observations (SS)']].copy(deep=True)
    data_moa = obs_MH_15[['Model (PMOA)', 'Observations (PMOA)']].copy(deep=True)

    data_ss.plot(ax = ax[0], color = ['b', 'r'])
    ax[0].set_title('Sea salt')
    ax[0].set_xlabel(' ')
    ax[0].set_ylabel('Concentration (ug m$^{-3}$)')

    data_moa.plot(ax = ax[1], color = ['b', 'r'])
    ax[1].set_xlabel(' ')
    ax[1].set_title('PMOA')
    ax[1].set_ylabel('Concentration (ug m$^{-3}$)')

    fig.tight_layout()
    plt.savefig(f'plots/MH{yr}_conc_daily_{global_vars.exp_name}.png')



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
    print('start seasonality at MH')
    plot_MC_monthly_seasonality('2018')
    plot_MC_daily_seasonality('2018')

    plot_MC_daily_seasonality('2015')




