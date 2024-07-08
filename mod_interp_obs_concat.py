import pandas as pd
import numpy as np
import read_data_functions
import os
import xarray as xr
import utils_func
import global_vars
from datetime import datetime


def datetime_to_integer(dt_time):
    return pd.to_datetime(dt_time).astype(int) * 1e9  # 10000*dt_time.year + 100*dt_time.month + dt_time.day


def get_interp_fixed_loc(lo, la, mod_da, mod_ro_da):
    mi_ma_lon, mi_ma_lat = [lo, lo], [la, la]

    interp_vars = utils_func.start_interp(mod_da,
                                          mod_ro_da,
                                          global_vars.variables_names,
                                          lo,
                                          la,
                                          mi_ma_lon,
                                          mi_ma_lat)
    return (interp_vars['POL'][0][0],
            interp_vars['PRO'][0][0],
            interp_vars['LIP'][0][0],
            interp_vars['SS'][0][0],
            interp_vars['SS_tot'][0][0])


def concat_lists(multi_lists):
    if len(multi_lists) > 1:
        new_concat_list = []
        for i in range(len(multi_lists)):
            new_concat_list.append(multi_lists[i][0])
    else:
        new_concat_list = multi_lists[0]
    return np.concatenate(new_concat_list)


def assign_loc_ship(path, exp, ds_btw, ds_sub, ds_lim, ID):
    global conc_mod_pol, conc_mod_pro, conc_mod_lip, conc_mod_ss, conc_mod_ss_tot
    time_loc = datetime_to_integer(ds_btw['Date/Time'])  # convert datetime to int
    da_loc = ds_btw['Date/Time'].dt.day.values  # list with months
    mo_loc = ds_btw['Date/Time'].dt.month.values  # list with months
    yr_loc = ds_btw['Date/Time'].dt.year.values  # list with months

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
        lat_btw, lon_btw, date_btw = [], [], []
        init_date_btw = [ds_btw['Date/Time'].dt.day.values[0], mo_loc[0], ds_btw['Date/Time'].dt.year.values[0]]

        lat_btw, lon_btw, mo_btw = [], [], []
        dates_btw = []
        dates_btw_btw = []
        lat_btw_btw = []
        lon_btw_btw = []
        n = 0
        for l, loc in enumerate(time_loc):
            if start <= loc <= end:
                if n == 0:
                    ref_date_btw = [da_loc[l], mo_loc[l], yr_loc[l]]
                    n += 1
                if ref_date_btw == [da_loc[l], mo_loc[l], yr_loc[l]]:
                    dates_btw_btw.append([da_loc[l], mo_loc[l], yr_loc[l]])
                    lat_btw_btw.append(ds_btw['Latitude'].values[l])
                    lon_btw_btw.append(ds_btw['Longitude'].values[l])
                else:
                    dates_btw.append(dates_btw_btw)
                    lat_btw.append(lat_btw_btw)
                    lon_btw.append(lon_btw_btw)

                    ref_date_btw = [da_loc[l], mo_loc[l], yr_loc[l]]
                    # print(ref_date_btw,'ref date')

                    dates_btw_btw = []
                    lat_btw_btw = []
                    lon_btw_btw = []

        print('interpolation data in process for ship trajectories')

        if len(dates_btw) > 1:
            print(len(dates_btw[0]))
            for id_btw, da_btw_list in enumerate(dates_btw):
                print((da_btw_list[0], 'date for ship trajectories'))
                mod_data, mod_ro_da = utils_func.read_model(path,
                                                            exp,
                                                            da_btw_list[0][0],
                                                            da_btw_list[0][1],
                                                            da_btw_list[0][2],
                                                            'tracer')

                print('ship trajectories between init and end date')
                lon_btw_360 = [lon180 % 360 for lon180 in lon_btw[id_btw]]
                mi_ma_lat = [min(lat_btw[id_btw]), max(lat_btw[id_btw])]
                mi_ma_lon = [min(lon_btw_360), max(lon_btw_360)]

                interp_vars = utils_func.start_interp(mod_data,
                                                      mod_ro_da,
                                                      global_vars.variables_names,
                                                      lon_btw_360,
                                                      lat_btw[id_btw],
                                                      mi_ma_lon,
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

            interp_btw_pol = concat_lists(interp_btw_pol)
            interp_btw_pro = concat_lists(interp_btw_pro)
            interp_btw_lip = concat_lists(interp_btw_lip)
            interp_btw_ss = concat_lists(interp_btw_ss)
            interp_btw_ss_tot = concat_lists(interp_btw_ss_tot)

            conc_mod_pol.append(np.nanmean(interp_btw_pol))
            conc_mod_pro.append(np.nanmean(interp_btw_pro))
            conc_mod_lip.append(np.nanmean(interp_btw_lip))
            conc_mod_ss.append(np.nanmean(interp_btw_ss))
            conc_mod_ss_tot.append(np.nanmean(interp_btw_ss_tot))

        else:
            print('ship location at init and end date')
            mod_data, mod_ro_da = utils_func.read_model(path,
                                                        exp,
                                                        da_start[i],
                                                        mo_start[i],
                                                        yr_start[i], 'tracer')
            start_lo = ds_lim['Start Longitude'].values[i] % 360
            start_la = ds_lim['Start Latitude'].values[i]
            interp_lim_start = get_interp_fixed_loc(start_lo,
                                                    start_la,
                                                    mod_data,
                                                    mod_ro_da)

            end_lo = ds_lim['End Longitude'].values[i] % 360
            end_la = ds_lim['End Latitude'].values[i]
            interp_lim_end = get_interp_fixed_loc(end_lo,
                                                  end_la,
                                                  mod_data,
                                                  mod_ro_da)

            print('start end interp ', interp_lim_start[0], interp_lim_start[1])
            conc_mod_pol.append(np.nanmean([interp_lim_start[0], interp_lim_end[0]]))
            conc_mod_pro.append(np.nanmean([interp_lim_start[1], interp_lim_end[1]]))
            conc_mod_lip.append(np.nanmean([interp_lim_start[2], interp_lim_end[2]]))
            conc_mod_ss.append(np.nanmean([interp_lim_start[3], interp_lim_end[3]]))
            conc_mod_ss_tot.append(np.nanmean([interp_lim_start[4], interp_lim_end[4]]))

        start_4_mod.append(ds_lim['Start Date/Time'].values[i])
        end_4_mod.append(ds_lim['End Date/Time'].values[i])

        conc_obs_pol_sub.append(ds_sub['CCHO_µg_per_m3'].values[i])
        conc_obs_tot_sub.append(ds_sub['OM_µg_per_m3'].values[i])
        conc_obs_pro.append(ds_sub['CAA_µg_per_m3'].values[i])
        conc_obs_lip.append(ds_sub['PG_µg_per_m3'].values[i])
        conc_obs_ss.append(ds_sub['SS_µg_per_m3'].values[i])
        if i < len(ds_lim['SS_µg_per_m3'].values):
            conc_obs_ss_tot.append(ds_lim['SS_µg_per_m3'].values[i])
        else:
            conc_obs_ss_tot.append(np.nan)

        interp_lim_start_pol = np.nanmean(conc_mod_pol)
        interp_lim_start_pro = np.nanmean(conc_mod_pro)
        interp_lim_start_lip = np.nanmean(conc_mod_lip)
        interp_lim_start_ss = np.nanmean(conc_mod_ss)
        interp_lim_start_ss_tot = np.nanmean(conc_mod_ss_tot)

        conc_model_pol.append(interp_lim_start_pol)
        conc_model_pro.append(interp_lim_start_pro)
        conc_model_lip.append(interp_lim_start_lip)
        conc_model_ss.append(interp_lim_start_ss)
        conc_model_ss_tot.append(interp_lim_start_ss_tot)

        conc_model_tot.append(interp_lim_start_pol +
                              interp_lim_start_pro +
                              interp_lim_start_lip)

        name_pro.append('DCAA   CAA')
        name_lip.append('Lipids  PG')

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


def interp_conc_stations(path, exp, obs, obs_tot, ID):
    start_4_mod, end_4_mod = [], []
    id_camp = []
    conc_model_pol, conc_model_pro, conc_model_lip, conc_model_tot, conc_model_ss, conc_model_ss_tot = [], [], [], [], [], []
    conc_obs_pol, conc_obs_pro, conc_obs_lip, conc_obs_tot, conc_obs_ss, conc_obs_ss_tot = [], [], [], [], [], []
    svd_data = [[], [], []]
    name_lip, name_pro = [], []
    dates_start = obs['Start Date/Time']
    dates_end = obs['End Date/Time']

    months = dates_start.dt.month.values
    years = dates_start.dt.year.values
    days = dates_start.dt.day.values
    lat, lon = obs['Start Latitude'].values[0], obs['Start Longitude'].values[0] % 360

    for m, da_st in enumerate(dates_start.values):
        if ID[:3] == 'SVD' or ID[:2] == 'RS':
            mod_da_btw, mod_ro_btw = [], []
            btw_start_end = np.arange(da_st,
                                      dates_end.values[m],
                                      dtype="datetime64[D]").astype(datetime)
            btw_start_end = [dt.strftime("%Y-%m-%d") for dt in pd.to_datetime(btw_start_end)]
            btw_start_end_object = [datetime.strptime(data_btw, '%Y-%m-%d').date() for data_btw in btw_start_end]
            months_btw_start_end = [i.month for i in btw_start_end_object]
            days_btw_start_end = [i.day for i in btw_start_end_object]
            years_btw_start_end = [i.year for i in btw_start_end_object]
            print(days_btw_start_end, months_btw_start_end)

            for m_id, mo_btw in enumerate(months_btw_start_end):
                mod_da, mod_ro = utils_func.read_model(path, exp, days_btw_start_end[m_id], mo_btw,
                                                       years_btw_start_end[m_id], 'tracer')
                mod_da_btw.append(mod_da)
                mod_ro_btw.append(mod_ro)
            mod_da_btw_ds = xr.concat(mod_da_btw, dim='time')
            mod_ro_btw_ds = xr.concat(mod_ro_btw, dim='time')

            interp_lim_start = get_interp_fixed_loc(lon,
                                                    lat,
                                                    mod_da_btw_ds,
                                                    mod_ro_btw_ds)
        else:
            mod_data, mod_ro_da = utils_func.read_model(path, exp, days[m], months[m], years[m], 'tracer')

            interp_lim_start = get_interp_fixed_loc(lon,
                                                    lat,
                                                    mod_data,
                                                    mod_ro_da)

        # get_dataframe()

        start_4_mod.append(obs['Start Date/Time'].values[m])
        end_4_mod.append(obs['End Date/Time'].values[m])

        conc_model_pol.append(interp_lim_start[0])
        conc_model_pro.append(interp_lim_start[1])
        conc_model_lip.append(interp_lim_start[2])
        conc_model_tot.append(np.sum(interp_lim_start[0] +
                                     interp_lim_start[1] +
                                     interp_lim_start[2]))
        conc_model_ss.append(interp_lim_start[3])
        conc_model_ss_tot.append(interp_lim_start[4])

        conc_obs_pol.append(obs['CCHO_µg_per_m3'].values[m])
        conc_obs_tot.append(obs['OM_µg_per_m3'].values[m])
        conc_obs_pro.append(obs['CAA_µg_per_m3'].values[m])
        conc_obs_lip.append(obs['PG_µg_per_m3'].values[m])
        conc_obs_ss.append(obs['SS_µg_per_m3'].values[m])
        if m < len(obs_tot['SS_µg_per_m3'].values):
            conc_obs_ss_tot.append(obs_tot['SS_µg_per_m3'].values[m])
        else:
            conc_obs_ss_tot.append(np.nan)

        name_pro.append('DCAA   CAA')
        name_lip.append('Lipids  PG')

        id_camp.append(ID)

    print(len(conc_obs_pol), len(conc_model_pol), '\n', conc_obs_pol, conc_model_pol)

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

    if ID == 'SVD18':
        pd_da_18 = pd.DataFrame({'ID': id_camp, 'Start Date/Time': start_4_mod, 'End Date/Time': end_4_mod,
                                 'conc_mod_pro': conc_model_pro, 'conc_obs_prot_sub': conc_obs_pro,
                                 'conc_obs_ss': conc_obs_ss, 'conc_mod_ss': conc_model_ss, })

        pd_da_18.to_pickle(f'pd_files/prot_conc_{ID}.pkl')

    return pd_da
