import pandas as pd
import numpy as np
import xarray as xr
from utils_functions import utils_func, global_vars
from datetime import datetime


def datetime_to_integer(dt_time):
    """
    Converts a datetime object to an integer.
    """
    return pd.to_datetime(dt_time).astype(int) * 1e9  # 10000*dt_time.year + 100*dt_time.month + dt_time.day


def get_interp_fixed_loc(lo, la, mod_da, mod_ro_da):
    """
    Calls interpolation function and returns interpolated model data for point locations.
    :param lo: station location longitude
    :param la: station location latitude
    :param mod_da: model dataset
    :param mod_ro_da: model dataset with air density
    :return: Interpolated values for SS, OC and SS all sizes
    """
    mi_ma_lon, mi_ma_lat = [lo, lo], [la, la]

    interp_vars = utils_func.start_interp(mod_da,
                                          mod_ro_da,
                                          global_vars.variables_names,
                                          lo,
                                          la,
                                          mi_ma_lon,
                                          mi_ma_lat)
    return (interp_vars['SS'][0][0],
            interp_vars['OC'][0][0],
            interp_vars['SS_tot'][0][0])


def concat_lists(multi_lists):
    """
    Concatenate multiple lists into one numpy array.
    :param multi_lists: list of interpolated values
    :return: array of concatenated values
    """
    if len(multi_lists) > 1:
        new_concat_list = []
        for i in range(len(multi_lists)):
            new_concat_list.append(multi_lists[i][0])
    else:
        new_concat_list = multi_lists[0]
    return np.concatenate(new_concat_list)


def assign_loc_ship(path, exp, ds_btw, ds_sub, ds_lim, ID):
    """
    Iterates through the available observational data to perform model interpolation for ship data (changing location).
    Interpolates model results to each location and then average according to limits
    Use for ship-based campaigns (PASCAL and PI_ICE)
    :var path: path to model data
    :var exp: model experiment name
    :var ds_btw: dataframe with ship trajectory
    :var ds_sub: dataframe with submicron aerosol concentration
    :var ds_lim: dataframe with all sizes aerosol concentration
    :param ID: Campaign ID
    :return: Dataframe with observation, interpolated model values, date, species and campaign name for polysaccharides
    """
    global conc_mod_ss, conc_mod_ss_tot
    time_loc = datetime_to_integer(ds_btw['Date/Time'])  # convert datetime to int
    da_loc = ds_btw['Date/Time'].dt.day.values  # list with months
    mo_loc = ds_btw['Date/Time'].dt.month.values  # list with months
    yr_loc = ds_btw['Date/Time'].dt.year.values  # list with months

    time_var_start = datetime_to_integer(ds_lim['Start Date/Time'])  # convert datetime to int
    da_start = ds_lim['Start Date/Time'].dt.day.values  # list with months
    mo_start = ds_lim['Start Date/Time'].dt.month.values  # list with months
    yr_start = ds_lim['Start Date/Time'].dt.year.values  # list with months
    time_var_end = datetime_to_integer(ds_lim['End Date/Time'])  # convert datetime to int

    start_4_mod, end_4_mod = [], []
    conc_model_oc = []
    conc_model_ss, conc_model_ss_tot = [], []


    conc_obs_oc = []
    conc_obs_ss, conc_obs_ss_tot = [], []
    id_camp = []

    for i, (start, end) in enumerate(zip(time_var_start, time_var_end)):

        # times, lats and lons between start and end
        conc_mod_ss, conc_mod_ss_tot, conc_mod_oc = [], [], []
        interp_btw_oc, interp_btw_ss, interp_btw_ss_tot = [], [], []

        # saving lats, longs and months between start and end into a list
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

                f_interp_ss = interp_vars['SS']
                f_interp_oc = interp_vars['OC']
                f_interp_ss_tot = interp_vars['SS_tot']


                interp_btw_ss.append(f_interp_ss)
                interp_btw_oc.append(f_interp_oc)
                interp_btw_ss_tot.append(f_interp_ss_tot)


            interp_btw_ss = concat_lists(interp_btw_ss)
            interp_btw_oc = concat_lists(interp_btw_oc)
            interp_btw_ss_tot = concat_lists(interp_btw_ss_tot)


            conc_mod_ss.append(np.nanmean(interp_btw_ss))
            conc_mod_oc.append(np.nanmean(interp_btw_oc))
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
            conc_mod_ss.append(np.nanmean([interp_lim_start[0], interp_lim_end[0]]))
            conc_mod_oc.append(np.nanmean([interp_lim_start[1], interp_lim_end[1]]))
            conc_mod_ss_tot.append(np.nanmean([interp_lim_start[2], interp_lim_end[2]]))

        start_4_mod.append(ds_lim['Start Date/Time'].values[i])
        end_4_mod.append(ds_lim['End Date/Time'].values[i])

        conc_obs_oc.append(ds_sub['OC_µg_per_m3'].values[i])
        conc_obs_ss.append(ds_sub['SS_µg_per_m3'].values[i])
        if i < len(ds_lim['SS_µg_per_m3'].values):
            conc_obs_ss_tot.append(ds_lim['SS_µg_per_m3'].values[i])
        else:
            conc_obs_ss_tot.append(np.nan)

        interp_lim_start_ss = np.nanmean(conc_mod_ss)
        interp_lim_start_oc = np.nanmean(conc_mod_oc)
        interp_lim_start_ss_tot = np.nanmean(conc_mod_ss_tot)

        conc_model_ss.append(interp_lim_start_ss)
        conc_model_oc.append(interp_lim_start_oc)
        conc_model_ss_tot.append(interp_lim_start_ss_tot)

        id_camp.append(ID)

    # create new dataframe to store the data after filtering all lat, lons and Date/Time
    pd_da = pd.DataFrame({'ID': id_camp, 'Start Date/Time': start_4_mod, 'End Date/Time': end_4_mod,
                          'conc_obs_ss': conc_obs_ss, 'conc_mod_ss': conc_model_ss,
                          'conc_obs_oc': conc_obs_oc, 'conc_mod_oc': conc_model_oc,
                          'conc_obs_ss_tot': conc_obs_ss_tot, 'conc_mod_ss_tot': conc_model_ss_tot})

    pd_da['Start Date/Time'] = pd_da['Start Date/Time'].apply(pd.to_datetime)
    pd_da['End Date/Time'] = pd_da['End Date/Time'].apply(pd.to_datetime)

    return pd_da


def interp_conc_stations(path, exp, obs, obs_tot, ID):
    """
    Iterates through the available observational data to perform model interpolation for fixed location stations.
    Also saves dataframe in a pickle file
    :param path: Path to the data directory where model data is stored
    :param exp: model experiment ID
    :param obs: observation data submicron size (aerosol concentration)
    :param obs_tot: observation data all sizes (aerosol concentration)
    :param ID: Id of the observation campaign
    :return: dataframe with observational data and interpolated model results for lipids and DCAA
    """
    start_4_mod, end_4_mod = [], []
    id_camp = []
    conc_model_ss, conc_model_oc, conc_model_ss_tot = [], [], []
    conc_obs_ss, conc_obs_oc, conc_obs_ss_tot = [], [], []
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


        start_4_mod.append(obs['Start Date/Time'].values[m])
        end_4_mod.append(obs['End Date/Time'].values[m])


        conc_model_ss.append(interp_lim_start[0])
        conc_model_oc.append(interp_lim_start[1])
        conc_model_ss_tot.append(interp_lim_start[2])

        conc_obs_ss.append(obs['SS_µg_per_m3'].values[m])
        conc_obs_oc.append(obs['OC_µg_per_m3'].values[m])

        if m < len(obs_tot['SS_µg_per_m3'].values):
            conc_obs_ss_tot.append(obs_tot['SS_µg_per_m3'].values[m])
        else:
            conc_obs_ss_tot.append(np.nan)


        id_camp.append(ID)

    pd_da = pd.DataFrame({'ID': id_camp, 'Start Date/Time': start_4_mod, 'End Date/Time': end_4_mod,
                          'conc_obs_ss': conc_obs_ss, 'conc_mod_ss': conc_model_ss,
                          'conc_obs_oc': conc_obs_oc, 'conc_mod_oc': conc_model_oc,
                          'conc_obs_ss_tot': conc_obs_ss_tot, 'conc_mod_ss_tot': conc_model_ss_tot,})

    pd_da['Start Date/Time'] = pd_da['Start Date/Time'].apply(pd.to_datetime)
    pd_da['End Date/Time'] = pd_da['End Date/Time'].apply(pd.to_datetime)

    if ID == 'SVD18':
        pd_da_18 = pd.DataFrame({'ID': id_camp, 'Start Date/Time': start_4_mod, 'End Date/Time': end_4_mod,
                                 'conc_obs_ss': conc_obs_ss, 'conc_mod_ss': conc_model_ss, })

        pd_da_18.to_pickle(f'../outputs/prot_conc_{ID}.pkl')

    return pd_da
