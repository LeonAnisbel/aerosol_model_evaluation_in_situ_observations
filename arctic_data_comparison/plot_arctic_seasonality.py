import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils_functions import read_data_functions, global_vars

ff = 12
omf_var = 'OMF_mod_tot'

def concat_omf_mod_obs(data, obs_MH, col_name, yr):
    """
    This function will retrieve observation data and return it as monthly averaged values
    :param data: dataframe with observation and model data
    :param obs_MH: dataframe with all observation data (without temporal averaging)
    :param col_name: string with the name of the column with the data from observations
    :param yr: string indicating the year or year range
    :return: dataframe of monthly averaged values
    """
    if yr == '0209':
        obs_MH[col_name] = data[omf_var].values
    else:
        times = pd.to_datetime(data['Start Date/Time'], dayfirst=True)
        data_omf_monthly = data.groupby([times.dt.year, times.dt.month])[[omf_var]].mean()
        obs_MH[col_name] = data_omf_monthly[omf_var].values

    return obs_MH

def fill_btw(axs, data, std, c, plot_std=True):
    """
    This function is used to plot the standard deviation as shaded area
    :param axs: matplotlib axes
    :var data: monthly data
    :var std: monthly standard deviation
    :param c: color
    :param plot_std: boolean to indicate whether to plot standard deviation
    :returns None
    """
    if plot_std:
        min, max = ([i - j for i, j in zip(data.values, std)],
                    [i + j for i, j in zip(data.values, std)])
        axs.fill_between(np.arange(1, 13),
                         min, max,
                         facecolor=c,
                         alpha=0.1)


def format_func(value, tick_number):
    """
    This function will create the months labels considering that January is month 1
    :returns value
    """
    N = int(value )
    return N


def expand_btwn_start_end_date(data_0209, vars):
    """
    This function will expand monthly data from 0209 to start and end dates
    :param data_0209: dataframe with observation data (no temporal average)
    :param vars: variables to extract from data_0209
    :return: dataframes of monthly data values and std values
    """
    date_end = 'End Date/Time'
    date = 'Start Date/Time'

    data_0209.loc[:, (date_end)] = data_0209[date_end].apply(pd.to_datetime)
    data_0209['dates'] = data_0209.apply(
        lambda row: pd.date_range(row[date], row[date_end], freq='D'),
        axis=1)
    daily = data_0209.explode('dates').rename(columns={'dates': 'date'})
    data_0209_m = (daily
                    .groupby(daily['date'].dt.to_period('M'))[vars]
                    .mean())
    start = []
    for p in data_0209_m.index:
        start.append(p.to_timestamp(how='start'))  # Timestamp('2002-03-01 00:00:00')
        end = p.to_timestamp(how='end')  # Timestamp('2002-03-31 00:00:00')
    data_0209_m['start'] = start
    data_0209_m['start'] = data_0209_m['start'].apply(pd.to_datetime)
    times = data_0209_m['start']

    #compute monthly means and multiannual std
    data_0209_std = (data_0209_m.groupby([times.dt.month])[vars]
                     .std())
    data_0209_mm = (data_0209_m.groupby([times.dt.month])[vars]
                    .mean())
    return data_0209_mm, data_0209_std


def plot_MC_monthly_seasonality(yr, var_names):
    """
    Reads data of observation versus interpolated model values at the station locations
    Creates monthly seasonality plots of OMF at Mace Head and a 3-panel plot of PMOA and SS at Mace Head and PMOA at
    Amsterdam Island
    :param yr: string indicating the year or year range
    :param var_names: names of variables to extract from observational data
    :return: None
    """
    data = pd.read_pickle(f'../outputs/MH{yr}_conc_{global_vars.exp_name}.pkl')
    if yr == '0209':
        stat_exp = 'MH_Rinaldi'
        pmoa_name = 'WIOM'
    else:
        stat_exp = 'MH'
        pmoa_name = 'PMOA'

    _, _, months, _, _ = read_data_functions.read_data(yr, var_names)
    obs_MH_15, _, _, _, obs_MH_15_std = read_data_functions.read_data(yr,
                                                                      var_names,
                                                                      monthly=True)

    data_month_mean = []
    data_month_std = []

    # Read OMF and add it to plot
    data_omf_recom = pd.read_pickle(f'../outputs/tot_omf_{stat_exp}.pkl')
    data_omf_recom_cesm = pd.read_pickle(f'../outputs/fesom_burrows_tot_omf_{stat_exp}.pkl')
    data_omf_cesm = pd.read_pickle(f'../outputs/burrows_tot_omf_{stat_exp}.pkl')

    if yr == '0209':
        pd_data = []
        pd_data_std = []
        for data_0209 in [data_omf_recom, data_omf_recom_cesm, data_omf_cesm]:
            data_0209_omf_hr, data_0209_omf_std = expand_btwn_start_end_date(data_0209,
                                                                             [omf_var])

            pd_data.append(data_0209_omf_hr)
            pd_data_std.append(data_0209_omf_std)

        data_omf_recom = pd_data[0]
        data_omf_recom_cesm = pd_data[1]
        data_omf_cesm = pd_data[2]

        data_0209_conc_hr, data_0209_conc_std = expand_btwn_start_end_date(data,
                                                                           ['conc_mod_ss', 'conc_mod_tot'])
        ss_vals = data_0209_conc_hr['conc_mod_ss']
        ss_std = data_0209_conc_std['conc_mod_ss']
        moa_vals = data_0209_conc_hr['conc_mod_tot']
        moa_std = data_0209_conc_std['conc_mod_tot']
    else:
        for i in np.arange(1, 13):
            data_std = data.loc[(data['months'] == i),
            ('conc_mod_ss', 'conc_mod_tot')].std()
            data_mean = data.loc[(data['months'] == i),
            ('conc_mod_ss', 'conc_mod_tot')].mean()

            data_month_mean.append(data_mean.values)
            data_month_std.append(data_std.values)

        ss_vals, moa_vals = [i[0] for i in data_month_mean], [i[1] for i in data_month_mean]
        ss_std = [i[0] if not np.isnan(i[0]) else 0
                  for i in data_month_std]
        moa_std = [i[1] if not np.isnan(i[1]) else 0
                   for i in data_month_std]

    obs_MH_15['Model (SS)'] = ss_vals
    obs_MH_15['Model (PMOA)'] = moa_vals


    obs_MH_15 = concat_omf_mod_obs(data_omf_recom,
                                   obs_MH_15,
                                   'REcoM',
                                   yr)
    obs_MH_15 = concat_omf_mod_obs(data_omf_recom_cesm,
                                   obs_MH_15,
                                   'REcoM+CESM',
                                   yr)
    obs_MH_15 = concat_omf_mod_obs(data_omf_cesm,
                                   obs_MH_15,
                                   'CESM',
                                   yr)

    obs_MH_15.rename(columns={var_names[1]:'Observations (SS)',
                              var_names[2]:f'Observations ({pmoa_name})',
                              'OMF': 'Observations'},
                     inplace=True)

    data_ss = obs_MH_15[['Model (SS)', 'Observations (SS)']].copy(deep=True)
    data_moa = obs_MH_15[['Model (PMOA)', f'Observations ({pmoa_name})']].copy(deep=True)
    data_omf = obs_MH_15[['REcoM', 'CESM', 'REcoM+CESM', 'Observations']].copy(deep=True)

    # Plot Monthly OMF for Mace Head considering the three cases above. See also global_vars.py for more information
    fig, ax = plt.subplots(1,
                           1,
                           figsize=(5,5))
    cc = ['darkorange','dodgerblue', 'limegreen', 'black']
    data_omf.plot(ax = ax,
                  color =cc)
    fill_btw(ax,
             data_omf['Observations'],
             obs_MH_15_std['OMF'].values,
             'gray')

    plt.legend(loc='upper right')
    ax.set_ylabel('OMF',
                fontsize=ff)
    ax.set_xlabel(' ')
    ax.set_xlabel('Months',
                 fontsize=ff)
    ax.set_xticklabels([])
    ax.set_ylim([0, 0.8])
    xticks = np.arange(1,13)
    ax.set_xticks(xticks)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    ax.xaxis.set_tick_params(labelsize=ff)
    ax.tick_params(axis='both',
                   labelsize=str(ff))
    ax.grid(linestyle='--',
            linewidth=0.4,
            alpha=0.2)
    fig.tight_layout()
    plt.savefig(f'../plots/{stat_exp}{yr}_monthly_OMF_{global_vars.exp_name}.png',
                dpi=300)


    # Plot Monthly PMOA and SS aerosol concentration for Mace Head (FESOM-REcom--ECHAM-HAM results) and a third panel
    # with PMOA at Amsterdam Island
    fig, axs = plt.subplots(2,
                            2,
                            figsize=(10,8))
    axs.flatten()
    cc = ['darkblue', 'black']
    data_ss.plot(ax = axs[0][0], color = cc)
    fill_btw(axs[0][0],
             data_ss['Observations (SS)'],
             obs_MH_15_std[var_names[1]].values,
             'gray')
    fill_btw(axs[0][0],
             data_ss['Model (SS)'],
             ss_std,
             'darkblue')
    axs[0][0].set_ylabel('Concentration (µg m$^{-3}$)',
                     fontsize=ff)

    cc = ['darkred', 'black']
    data_moa.plot(ax = axs[0][1], color = cc)
    fill_btw(axs[0][1],
             data_moa[f'Observations ({pmoa_name})'],
             obs_MH_15_std[var_names[2]].values,
             'gray' )
    fill_btw(axs[0][1],
             data_moa['Model (PMOA)'],
             moa_std,
             'r')
    axs[0][1].set_ylabel('Concentration (µg m$^{-3}$)',
                     fontsize=ff)

    indices = [r'$\bf{(a)}$', r'$\bf{(b)}$']
    for a, i in zip(axs[0], indices):
        a.set_xlabel(' ')
        a.set_title(i, loc='left')
        a.set_xlabel('Months',
                     fontsize=ff)
        a.set_xticklabels([])
        xticks = np.arange(1, 13)
        a.set_xticks(xticks)

        a.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        a.xaxis.set_tick_params(labelsize=ff)

        a.tick_params(axis='both',
                      labelsize=str(ff))
        a.set_title('Mace Head',
                    loc='center')
        a.grid(linestyle='--',
               linewidth=0.4,
               alpha=0.2)
    plot_AI_monthly_seasonality(axs[1][1],
                                '0307',
                                panel_fig=True)
    axs[1][1].grid(linestyle='--',
                   linewidth=0.4,
                   alpha=0.2)
    fig.delaxes(axs[1][0])

    fig.tight_layout()
    plt.savefig(f'../plots/MH{yr}_monthly_conc_{global_vars.exp_name}.png',
                dpi=300)


def plot_MC_daily_seasonality(yr, var_names):
    """
    Reads data of observation versus interpolated model values at the station locations
    Creates plot of daily observations and model interpolated daily values at MAce Head
    :param yr: string indicating the year or year range
    :param var_names: names of variables to extract from observational data
    :return: None
    """
    data = pd.read_pickle(f'../outputs/MH{yr}_conc_{global_vars.exp_name}.pkl')
    fig, ax = plt.subplots(2,
                           1,
                           figsize=(10,5))
    ax.flatten()
    obs_MH_15, _, _, _, _ = read_data_functions.read_data(yr, var_names)
    obs_MH_15['Model (SS)'] = data['conc_mod_ss'].values
    obs_MH_15['Model (PMOA)'] = data['conc_mod_tot'].values
    obs_MH_15.rename(columns={var_names[1]:'Observations (SS)',
                              var_names[2]:'Observations (PMOA)'},
                     inplace=True)

    data_ss = obs_MH_15[['Model (SS)', 'Observations (SS)']].copy(deep=True)
    data_moa = obs_MH_15[['Model (PMOA)', 'Observations (PMOA)']].copy(deep=True)

    cc = ['#4B0082', 'darkred']
    data_ss.plot(ax = ax[0],
                 color = cc)
    ax[0].set_title('Sea salt')
    ax[0].set_xlabel(' ')
    ax[0].set_ylabel('Concentration (ug m$^{-3}$)')

    data_moa.plot(ax = ax[1],
                  color = cc)
    ax[1].set_xlabel(' ')
    ax[1].set_title('PMOA')
    ax[1].set_ylabel('Concentration (ug m$^{-3}$)')

    fig.tight_layout()
    plt.savefig(f'../plots/MH{yr}_conc_daily_{global_vars.exp_name}.png')


def plot_all_arctic_stations():
    """
    Reads data of observation versus interpolated model values for all Arctic stations
    Create plot of all regions with observation and interpolated model values. Find more info about this dataset
    in main_Arctic.py
    :return: None
    """
    data = pd.read_pickle(f'../outputs/all_arctic_stations_conc_{global_vars.exp_name}.pkl')
    data_mo_all_stations, dict_metadata, data_sel_std = read_data_functions.read_PMOA_all_stations()
    stat_list = list(dict_metadata.keys())

    data_sel_yr = ((data.groupby(['months', 'ID'])[[ 'conc_obs_tot',  'conc_mod_tot']])
                   .mean())
    data_std = ((data.groupby(['months', 'ID'])[[ 'conc_obs_tot',  'conc_mod_tot']])
                   .std())
    data_std['obs_std_fill_min'] = [val-std for val,std in zip(data_sel_yr['conc_obs_tot'].values,
                                                               data_std['conc_obs_tot'].values)]
    data_std['obs_std_fill_max'] = [val+std for val,std in zip(data_sel_yr['conc_obs_tot'].values,
                                                               data_std['conc_obs_tot'].values)]
    data_sel_yr = data_sel_yr.reset_index()
    data_std=data_std.reset_index()
    data_sel_yr.set_index('months',
                          inplace=True)
    data_std.set_index('months',
                       inplace=True)

    fig, ax = plt.subplots(len(stat_list),
                           1,
                           figsize=(12, 20))
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
    plt.savefig(f'../plots/all_arctic_stations_conc_bar_yearly_{global_vars.exp_name}.png')
    plt.close()

    date_list = [str(i) + '-' + str(j) for i, j in zip(data['years'].values, data['months'].values)]
    print(date_list)
    data['date'] = date_list
    data.drop('years',
              axis=1,
              inplace=True)
    data.drop('months',
              axis=1,
              inplace=True)
    data.set_index('date',
                   inplace=True)

    fig, ax = plt.subplots(len(stat_list),
                           1,
                           figsize=(12, 20))
    ax.flatten()
    for idx, sta in enumerate(stat_list):
        # print(sta, data[['ID']==sta])
        data[data['ID'] == sta].plot(ax=ax[idx])
        ax[idx].set_title(sta,
                          loc='center')
    fig.tight_layout()
    plt.savefig(f'../plots/all_arctic_stations_conc_bar_{global_vars.exp_name}.png')
    plt.close()

    fig, ax = plt.subplots(1, 1)
    obs_leg = []
    mod_leg = []
    for idx, sta in enumerate(stat_list):
        data_mod = data[['ID', 'conc_mod_tot']]
        data_obs = data[['ID', 'conc_obs_tot']]
        m = data_mod[data['ID'] == sta].plot(style='-',
                                             ax=ax)  # '.-'
        o = data_obs[data['ID'] == sta].plot(style='--',
                                             ax=ax)  # '.-'
        obs_leg.append(o)
        mod_leg.append(m)

    fig.legend(handles=obs_leg,
               labels=stat_list)
    plt.savefig(f'../plots/all_arctic_stations_conc_{global_vars.exp_name}.png')


def plot_AI_monthly_seasonality(ax, yr, panel_fig=False):
    """
    Reads data of observation versus interpolated model values at Amsterdam Island
    Creates monthly seasonality plots of PMOA at Amsterdam Island
    :param ax: matplotlib axes object
    :param yr: string indicating the year or year range
    :param panel_fig: whether to create a figure or plot data on matplotlib axis ax
    :return: None
    """
    data_moa = pd.read_pickle(f'../outputs/AI{yr}_conc_{global_vars.exp_name}.pkl')

    if panel_fig:
        ax = ax
        ax.set_title(r'$\bf{(c)}$', loc='left')
    else:
        fig, ax = plt.subplots(1,1, figsize=(5,5))

    ax.set_title('Amsterdam Island')
    cc = ['darkred', 'black']
    ax.plot(data_moa['months'],
            data_moa['conc_mod_tot'],
            color = cc[0],
            label = 'Model (PMOA)')
    ax.plot(data_moa['months'],
            data_moa['conc_obs_tot'],
            color = cc[1],
            label = 'Observations (WIOM)')

    fill_btw(ax, data_moa[f'conc_obs_tot'],
             data_moa['conc_obs_std'].values,
             'gray' )
    ax.set_ylabel('Concentration (µg m$^{-3}$)',
                     fontsize=ff)

    ax.legend()
    ax.set_xlabel(' ')
    ax.set_xlabel('Months',
                 fontsize=ff)
    ax.set_xticklabels([])
    xticks = np.arange(1, 13)
    ax.set_xticks(xticks)

    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    ax.xaxis.set_tick_params(labelsize=ff)

    ax.tick_params(axis='both', labelsize=str(ff))

    if not panel_fig:
        fig.tight_layout()
        plt.savefig(f'../plots/AI{yr}_monthly_conc_{global_vars.exp_name}.png')
        plt.close()

        # plot also SS
        fig, ax = plt.subplots(1,1, figsize=(5,5))
        ax.plot(data_moa['months'],
                data_moa['conc_mod_ss'],
                color = 'green',
                label = 'Model (SS)')
        plt.legend()

        ax.set_ylabel('Concentration (µg m$^{-3}$)',
                         fontsize=ff)

        ax.set_xlabel(' ')
        ax.set_xlabel('Months',
                     fontsize=ff)
        ax.set_xticklabels([])
        xticks = np.arange(1, 13)
        ax.set_xticks(xticks)

        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        ax.xaxis.set_tick_params(labelsize=ff)

        ax.tick_params(axis='both', labelsize=str(ff))

        fig.tight_layout()
        plt.savefig(f'../plots/AI{yr}_monthly_conc_SS_{global_vars.exp_name}.png')
        plt.close()



if __name__ == '__main__':
    year = '0307'  # from May 2003 to November 2007
    print('start seasonality at AI', year)
    plot_AI_monthly_seasonality(None, year)

    year = '0209' # 2002-2009
    print('start seasonality at MH', year) # Thesis plot (Chapter 5)
    plot_MC_monthly_seasonality(year,
                                global_vars.yr_exp_var_names[year])

    year = '2018'
    print('start seasonality at MH', year)
    plot_MC_monthly_seasonality(year,
                                global_vars.yr_exp_var_names[year])
    plot_MC_daily_seasonality(year,
                              global_vars.yr_exp_var_names[year])

    # year = '2015'
    # print('start seasonality at MH', year)
    # plot_MC_daily_seasonality(year,
    #                             global_vars.yr_exp_var_names[year])

    print('start seasonality all Arctic stations', year)
    plot_all_arctic_stations()





