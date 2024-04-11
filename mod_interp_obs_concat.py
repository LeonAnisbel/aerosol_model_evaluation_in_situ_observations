import pandas as pd
import numpy as np
import read_data_functions
import os
import xarray as xr
import utils_func
import global_vars


def datetime_to_integer(dt_time):
    return pd.to_datetime(dt_time).astype(int) * 1e9  # 10000*dt_time.year + 100*dt_time.month + dt_time.day


def assign_loc_ship(path, exp, ds_btw, ds_sub, ds_lim, ID):
    time_loc = datetime_to_integer(ds_btw['Date/Time'])  # convert datetime to int
    mo_loc = ds_btw['Date/Time'].dt.month.values  # list with months

    time_var_start = datetime_to_integer(ds_lim['Start Date/Time'])  # convert datetime to int
    da_start = ds_lim['Start Date/Time'].dt.day.values  # list with months
    mo_start = ds_lim['Start Date/Time'].dt.month.values  # list with months
    yr_start = ds_lim['Start Date/Time'].dt.year.values  # list with months

    time_var_end = datetime_to_integer(ds_lim['End Date/Time'])  # convert datetime to int
    mo_end = ds_lim['End Date/Time'].dt.month.values  # list with months

    start_4_mod, end_4_mod = [], []
    conc_model_pol, conc_model_pro, conc_model_lip = [], [], []
    conc_model_ss, conc_model_ss_tot, conc_model_tot = [], [], []

    conc_obs_pol_sub, conc_obs_tot_sub = [], []
    conc_obs_pol_sup, conc_obs_tot_sup = [], []

    conc_obs_pol, conc_obs_pro, conc_obs_lip = [], [], []
    conc_obs_ss, conc_obs_ss_tot, conc_obs_tot = [], [], []
    name_lip, name_pro = [], []
    id_camp = []

    for i, (start, end) in enumerate(zip(time_var_start, time_var_end)):

        # times, lats and lons between start and end
        conc_mod_pol, conc_mod_pro, conc_mod_lip, conc_mod_ss, conc_mod_ss_tot = [], [], [], [], []
        interp_btw = []
        interp_btw_pol = []
        interp_btw_pro = []
        interp_btw_lip = []
        interp_btw_ss = []
        interp_btw_ss_tot = []

        # saving lats, longs and months between start and end into a list
        lat_btw, lon_btw, mo_btw = [], [], []
        for l, loc in enumerate(time_loc):
            if loc >= start and loc <= end:  # and \
                #             ds_btw['Latitude'].values[l] >= ds_lim['Start Latitude'].values[i] and \
                #             ds_btw['Latitude'].values[l] <= ds_lim['End Latitude'].values[i]:
                lat_btw.append(ds_btw['Latitude'].values[l])
                lon_btw.append(ds_btw['Longitude'].values[l])

        mod_data, mod_ro_da = utils_func.read_model(path, exp, da_start[i], mo_start[i], yr_start[i], 'tracer')
        print('interpolation data in process for ship trajectories')

        if len(lat_btw) > 1:
            print('ship trajectories between init and end date')
            lon_btw = [lon180 % 360 for lon180 in lon_btw]
            mi_ma_lat = [min(lat_btw), max(lat_btw)]
            mi_ma_lon = [min(lon_btw), max(lon_btw)]

            interp_vars = utils_func.start_interp(mod_data, mod_ro_da, global_vars.variables_names,
                                                  lon_btw, lat_btw, mi_ma_lon,
                                                  mi_ma_lat)
            f_interp_pol = interp_vars['POL']
            f_interp_pro = interp_vars['PRO']
            f_interp_lip = interp_vars['LIP']
            f_interp_ss = interp_vars['SS']
            f_interp_ss_tot = interp_vars['SS_tot']

            interp_btw_pol.append(f_interp_pol)
            interp_btw_pro.append(f_interp_pro)
            interp_btw_lip.append(f_interp_lip)
            interp_btw_ss.append(f_interp_ss)
            interp_btw_ss_tot.append(f_interp_ss_tot)

            conc_mod_pol.append(np.mean(interp_btw_pol))
            conc_mod_pro.append(np.mean(interp_btw_pro))
            conc_mod_lip.append(np.mean(interp_btw_lip))
            conc_mod_ss.append(np.mean(interp_btw_ss))
            conc_mod_ss_tot.append(np.mean(interp_btw_ss_tot))

        else:
            print('ship location at init and end date')

            start_lo = ds_lim['Start Longitude'].values[i] % 360
            start_la = ds_lim['Start Latitude'].values[i]
            mi_ma_lon, mi_ma_lat = [start_lo, start_lo], [start_la, start_la]

            interp_vars = utils_func.start_interp(mod_data, mod_ro_da, global_vars.variables_names,
                                                  start_lo, start_la, mi_ma_lon,
                                                  mi_ma_lat)
            interp_lim_start = interp_vars['POL'][0][0]
            interp_lim_start_pr = interp_vars['PRO'][0][0]
            interp_lim_start_li = interp_vars['LIP'][0][0]
            interp_lim_start_ss = interp_vars['SS'][0][0]
            interp_lim_start_ss_tot = interp_vars['SS_tot'][0][0]

            end_lo = ds_lim['End Longitude'].values[i] % 360
            end_la = ds_lim['End Latitude'].values[i]
            mi_ma_lon, mi_ma_lat = [end_lo, end_lo], [end_la, end_la]

            interp_vars = utils_func.start_interp(mod_data, mod_ro_da, global_vars.variables_names,
                                                  end_lo, end_la, mi_ma_lon,
                                                  mi_ma_lat)
            interp_lim_end = interp_vars['POL'][0][0]
            interp_lim_end_pr = interp_vars['PRO'][0][0]
            interp_lim_end_li = interp_vars['LIP'][0][0]
            interp_lim_end_ss = interp_vars['SS'][0][0]
            interp_lim_end_ss_tot = interp_vars['SS_tot'][0][0]

            conc_mod_pol.append((interp_lim_start + interp_lim_end) / 2)
            conc_mod_pro.append((interp_lim_start_pr + interp_lim_end_pr) / 2)
            conc_mod_lip.append((interp_lim_start_li + interp_lim_end_li) / 2)
            conc_mod_ss.append((interp_lim_start_ss + interp_lim_end_ss) / 2)
            conc_mod_ss_tot.append((interp_lim_start_ss_tot + interp_lim_end_ss_tot) / 2)

        # if ds_lim['conc_pol'].values[i] > 0.:
        start_4_mod.append(ds_lim['Start Date/Time'].values[i])
        end_4_mod.append(ds_lim['End Date/Time'].values[i])

        conc_obs_pol_sub.append(ds_sub['CCHO_µg_per_m3'].values[i])
        conc_obs_tot_sub.append(ds_sub['OC_µg_per_m3'].values[i])
        conc_obs_pro.append(ds_sub['FAA_µg_per_m3'].values[i])
        conc_obs_lip.append(ds_sub['FFA_PG_µg_per_m3'].values[i])
        conc_obs_ss.append(ds_sub['SS_µg_per_m3'].values[i])
        if i < len(ds_lim['SS_µg_per_m3'].values):
            conc_obs_ss_tot.append(ds_lim['SS_µg_per_m3'].values[i])
        else:
            conc_obs_ss_tot.append(np.nan)

        interp_lim_start_pol = sum(conc_mod_pol) / len(conc_mod_pol)
        interp_lim_start_pro = sum(conc_mod_pro) / len(conc_mod_pro)
        interp_lim_start_lip = sum(conc_mod_lip) / len(conc_mod_lip)
        interp_lim_start_ss = sum(conc_mod_ss) / len(conc_mod_ss)
        interp_lim_start_ss_tot = sum(conc_mod_ss_tot) / len(conc_mod_ss_tot)

        conc_model_pol.append(interp_lim_start_pol)
        conc_model_pro.append(interp_lim_start_pro)
        conc_model_lip.append(interp_lim_start_lip)
        conc_model_ss.append(interp_lim_start_ss)
        conc_model_ss_tot.append(interp_lim_start_ss_tot)

        conc_model_tot.append(interp_lim_start_pol +
                              interp_lim_start_pro +
                              interp_lim_start_lip)

        name_pro.append('DCAA   DFAA')
        name_lip.append('Lipids  FFA_PG')

        id_camp.append(ID)

    # create new dataframe to store the data after filtering all lat, lons and Date/Time
    pd_da = pd.DataFrame({'ID': id_camp, 'Start Date/Time': start_4_mod, 'End Date/Time': end_4_mod,
                          'conc_mod_poly': conc_model_pol, 'conc_mod_pro': conc_model_pro,
                          'conc_mod_lip': conc_model_lip,
                          'conc_obs_poly_sub': conc_obs_pol_sub, 'conc_obs_prot_sub': conc_obs_pro,
                          'conc_obs_lipi_sub': conc_obs_lip,
                          'conc_obs_tot_sub': conc_obs_tot_sub, 'conc_mod_tot': conc_model_tot,
                          'conc_obs_ss': conc_obs_ss, 'conc_mod_ss': conc_model_ss,
                          'conc_obs_ss_tot': conc_obs_ss_tot, 'conc_mod_ss_tot': conc_model_ss_tot,
                          'Name_pro': name_pro, 'Name_lip': name_lip})

    pd_da['Start Date/Time'] = pd_da['Start Date/Time'].apply(pd.to_datetime)
    pd_da['End Date/Time'] = pd_da['End Date/Time'].apply(pd.to_datetime)

    return pd_da


# interpolation function
def interp_conc_stations(path, exp, obs, obs_tot, ID):
    start_4_mod, end_4_mod = [], []
    id_camp = []
    conc_model_pol, conc_model_pro, conc_model_lip, conc_model_tot, conc_model_ss, conc_model_ss_tot = [], [], [], [], [], []
    conc_obs_pol, conc_obs_pro, conc_obs_lip, conc_obs_tot, conc_obs_ss, conc_obs_ss_tot = [], [], [], [], [], []
    name_lip, name_pro = [], []
    months = obs['Start Date/Time'].dt.month.values
    years = obs['Start Date/Time'].dt.year.values
    days = obs['Start Date/Time'].dt.day.values

    for m, mo in enumerate(months):
        lat, lon = obs['Start Latitude'].values[0], obs['Start Longitude'].values[0] % 360
        mod_data, mod_ro_da = utils_func.read_model(path, exp, days[m], mo, years[m], 'tracer')

        mi_ma_lon, mi_ma_lat = [lon, lon], [lat, lat]
        # interpolation function
        interp_vars = utils_func.start_interp(mod_data, mod_ro_da, global_vars.variables_names,
                                              lon, lat, mi_ma_lon,
                                              mi_ma_lat)
        interp_lim_start_pol = interp_vars['POL'][0][0]
        interp_lim_start_pro = interp_vars['PRO'][0][0]
        interp_lim_start_lip = interp_vars['LIP'][0][0]
        interp_lim_start_ss = interp_vars['SS'][0][0]
        interp_lim_start_ss_tot = interp_vars['SS_tot'][0][0]

        start_4_mod.append(obs['Start Date/Time'].values[m])
        end_4_mod.append(obs['End Date/Time'].values[m])

        conc_model_pol.append(interp_lim_start_pol)
        conc_model_pro.append(interp_lim_start_pro)
        conc_model_lip.append(interp_lim_start_lip)
        conc_model_tot.append(interp_lim_start_pol +
                              interp_lim_start_pro +
                              interp_lim_start_lip)
        conc_model_ss.append(interp_lim_start_ss)
        conc_model_ss_tot.append(interp_lim_start_ss_tot)

        conc_obs_pol.append(obs['CCHO_µg_per_m3'].values[m])
        conc_obs_tot.append(obs['OC_µg_per_m3'].values[m])
        conc_obs_pro.append(obs['FAA_µg_per_m3'].values[m])
        conc_obs_lip.append(obs['FFA_PG_µg_per_m3'].values[m])
        conc_obs_ss.append(obs['SS_µg_per_m3'].values[m])
        if m < len(obs_tot['SS_µg_per_m3'].values):
            conc_obs_ss_tot.append(obs_tot['SS_µg_per_m3'].values[m])
        else:
            conc_obs_ss_tot.append(np.nan)

        name_pro.append('DCAA   DFAA')
        name_lip.append('Lipids  FFA_PG')

        id_camp.append(ID)

    # create new dataframe to store the model data after interpolation together with obs.
    pd_da = pd.DataFrame({'ID': id_camp, 'Start Date/Time': start_4_mod, 'End Date/Time': end_4_mod,
                          'conc_mod_poly': conc_model_pol, 'conc_mod_pro': conc_model_pro,
                          'conc_mod_lip': conc_model_lip,
                          'conc_obs_poly_sub': conc_obs_pol, 'conc_obs_prot_sub': conc_obs_pro,
                          'conc_obs_lipi_sub': conc_obs_lip,
                          'conc_obs_tot_sub': conc_obs_tot, 'conc_mod_tot': conc_model_tot,
                          'conc_obs_ss': conc_obs_ss, 'conc_mod_ss': conc_model_ss,
                          'conc_obs_ss_tot': conc_obs_ss_tot, 'conc_mod_ss_tot': conc_model_ss_tot,
                          'Name_pro': name_pro, 'Name_lip': name_lip})

    pd_da['Start Date/Time'] = pd_da['Start Date/Time'].apply(pd.to_datetime)
    pd_da['End Date/Time'] = pd_da['End Date/Time'].apply(pd.to_datetime)

    return pd_da
