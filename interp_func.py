import pandas as pd
import numpy as np
from scipy.interpolate import RectBivariateSpline
import read_data_functions
import os
import xarray as xr
import xesmf as xe

## interpolation functions

def regrid(ds_in,dr):
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], np.arange(np.nanmin(ds_in.lat.values),
                                       np.nanmax(ds_in.lat.values), 1.5),
                    {"units": "degrees_north"}),
            "lon": (["lon"], np.arange(0, 360, 1.5), {"units": "degrees_east"}),
        }
    )
    regridder = xe.Regridder(ds_in, ds_out, "conservative")
    dr_out = regridder(dr, keep_attrs=True)

    return dr_out

def get_array(da):
    da = da.to_dataset(name='VAR')
    da = da.to_array()
    da = np.nan_to_num(da, copy=True, nan=0.0)
    return da[0]

def interp(path, exp,dy,mo, yr,varb,modes,ext):


    unit_factor = 1e9
    
    data_dir = f"{path}"
    if mo < 10:
        mo_str = f'0{mo}'
    else:
        mo_str = f'{mo}'

    files = f'{data_dir}{exp}_{yr}{mo_str}.01_{ext}.nc'
    file_ro = f'{data_dir}{exp}_{yr}{mo_str}.01_vphysc.nc'

#    file = f'{data_dir}{exp}_201701.01_{ext}.nc'
#    file_uv = f'{data_dir}{exp}_201701.01_echam.nc'


    f_interp = []
    if os.path.exists(files):
        print('reading ',files , 'for interpolation with data month = ', mo)
        data = read_data_functions.read_model_spec_data(files)
        data_ro = read_data_functions.read_model_spec_data(file_ro)
        da_m_ro = data_ro['rhoam1'].isel(time = dy -1 ).isel(lev=46).compute()

        # only accumulation soluble mode
        for var in varb:
            da_m = 0
            da_m = data[f'{var}_AS'].isel(time = dy -1 ).isel(lev=46).compute()

#ICE test
            #da_m_ice = data_uv['seaice'].isel(time = dy -1 ).compute()
            #da_m_regrid_ice = get_array(regrid(data_uv,da_m_ice))
# ICE test

            da_m_new = regrid(data, da_m)
            da_m_regrid = get_array(regrid(data,da_m))
            da_m_regrid_ro = get_array(regrid(data_ro,da_m_ro))
            

            # calculating wind velocity vector for aerosol concentration
            #da_m_regrid_uv = np.sqrt(da_m_regrid_u ** 2 + da_m_regrid_v ** 2)
            # print(da_m_new,len(da_m_regrid),len(da_m_regrid[0]))
            # interpolation
            f_interp.append(
                RectBivariateSpline(da_m_new.lat.values, da_m_new.lon.values, da_m_regrid*da_m_regrid_ro*unit_factor))

        # creating interpolation grid for SS accumulation and coarse mode together
        da_m = 0
        for mo in modes:
            suma = da_m + data[f'SS_{mo}'].isel(time=dy - 1).isel(lev=46).compute()  # .to_dataset(name = 'VAR')
            da_m = suma
        da_m_new = regrid(data, da_m)
        da_m_regrid = get_array(regrid(data, da_m))
        da_m_regrid_ro = get_array(regrid(data_ro, da_m_ro))
        f_interp.append(
            RectBivariateSpline(da_m_new.lat.values, da_m_new.lon.values, da_m_regrid * da_m_regrid_ro * unit_factor))

    # ICE test
            #f_interp.append(RectBivariateSpline(da_m_new.lat.values, da_m_new.lon.values, da_m_regrid_ice))
# ICE test
    return f_interp


# save all lat, lon and date/Time in locations files between start and end in conc files
def assign_loc_ship(path, exp,ds_btw,ds_sub,ds_lim,ID,modes):
    def datetime_to_integer(dt_time):
        return pd.to_datetime(dt_time).astype(int) / 10 ** 9  # 10000*dt_time.year + 100*dt_time.month + dt_time.day

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
    conc_model_ss, conc_model_ss_tot , conc_model_tot = [], [], []

    conc_obs_pol_sub, conc_obs_tot_sub = [], []
    conc_obs_pol_sup, conc_obs_tot_sup = [], []

    conc_obs_pol, conc_obs_pro, conc_obs_lip =  [], [], []
    conc_obs_ss,conc_obs_ss_tot, conc_obs_tot = [], [], []
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

        # defining interpolation function
        f_int = interp(path, exp,da_start[i],mo_start[i], yr_start[i], ['POL', 'PRO', 'LIP','SS'],modes,'tracer')
        f_interp_pol = f_int[0]
        f_interp_pro = f_int[1]
        f_interp_lip = f_int[2]
        f_interp_ss = f_int[3]
        f_interp_ss_tot = f_int[4]



        if len(lat_btw) > 1:
            for la, lo in zip(lat_btw, lon_btw):
                interp_btw_pol.append(f_interp_pol(la,lo)[0][0])
                interp_btw_pro.append(f_interp_pro(la,lo)[0][0])
                interp_btw_lip.append(f_interp_lip(la,lo)[0][0])
                interp_btw_ss.append(f_interp_ss(la,lo)[0][0])
                interp_btw_ss_tot.append(f_interp_ss_tot(la,lo)[0][0])

            conc_mod_pol.append(sum(interp_btw_pol) / len(interp_btw_pol))
            conc_mod_pro.append(sum(interp_btw_pro) / len(interp_btw_pro))
            conc_mod_lip.append(sum(interp_btw_lip) / len(interp_btw_lip))
            conc_mod_ss.append(sum(interp_btw_ss) / len(interp_btw_ss))
            conc_mod_ss_tot.append(sum(interp_btw_ss_tot) / len(interp_btw_ss_tot))

        else:
            interp_lim_start = f_interp_pol(ds_lim['Start Latitude'].values[i],
                                            ds_lim['Start Longitude'].values[i])[0][0]
            interp_lim_end = f_interp_pol(ds_lim['End Latitude'].values[i],
                                          ds_lim['End Longitude'].values[i])[0][0]
            conc_mod_pol.append((interp_lim_start + interp_lim_end) / 2)

            interp_lim_start_pr = f_interp_pro(ds_lim['Start Latitude'].values[i],
                                               ds_lim['Start Longitude'].values[i])[0][0]
            interp_lim_end_pr = f_interp_pro(ds_lim['End Latitude'].values[i],
                                             ds_lim['End Longitude'].values[i])[0][0]
            conc_mod_pro.append((interp_lim_start_pr + interp_lim_end_pr) / 2)

            interp_lim_start_li = f_interp_lip(ds_lim['Start Latitude'].values[i],
                                               ds_lim['Start Longitude'].values[i])[0][0]
            interp_lim_end_li = f_interp_lip(ds_lim['End Latitude'].values[i],
                                             ds_lim['End Longitude'].values[i])[0][0]
            conc_mod_lip.append((interp_lim_start_li + interp_lim_end_li) / 2)

            interp_lim_start_ss = f_interp_ss(ds_lim['Start Latitude'].values[i],
                                               ds_lim['Start Longitude'].values[i])[0][0]
            interp_lim_end_ss = f_interp_ss(ds_lim['End Latitude'].values[i],
                                             ds_lim['End Longitude'].values[i])[0][0]
            conc_mod_ss.append((interp_lim_start_ss + interp_lim_end_ss) / 2)

            interp_lim_start_ss_tot = f_interp_ss_tot(ds_lim['Start Latitude'].values[i],
                                               ds_lim['Start Longitude'].values[i])[0][0]
            interp_lim_end_ss_tot = f_interp_ss_tot(ds_lim['End Latitude'].values[i],
                                             ds_lim['End Longitude'].values[i])[0][0]
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

        name_pro.append('DHAA   DFAA')
        name_lip.append('Lipids  FFA_PG')


        #
        # conc_obs_pol_sup.append(ds_sup['conc_pol'].values[i])
        # conc_obs_tot_sup.append(ds_sup['conc_tot'].values[i])

        id_camp.append(ID)

    # print('PROOOOO')
    # print(id_camp)
    # print(conc_obs_pol_sub)
    # print(conc_model_pol)
    #
    # print('PROOOOO TOTAL')
    # print(id_camp)
    # print(conc_obs_tot_sub)
    # print(conc_model_tot)

    # create new dataframe to store the data after filtering all lat, lons and Date/Time

    pd_da = pd.DataFrame({'ID': id_camp, 'Start Date/Time': start_4_mod, 'End Date/Time': end_4_mod,
                          'conc_mod_poly': conc_model_pol, 'conc_mod_pro': conc_model_pro, 'conc_mod_lip': conc_model_lip,
                          'conc_obs_poly_sub': conc_obs_pol_sub, 'conc_obs_prot_sub': conc_obs_pro,
                          'conc_obs_lipi_sub': conc_obs_lip,
                          'conc_obs_tot_sub': conc_obs_tot_sub,'conc_mod_tot': conc_model_tot,
                          'conc_obs_ss': conc_obs_ss,'conc_mod_ss': conc_model_ss,
                          'conc_obs_ss_tot': conc_obs_ss_tot, 'conc_mod_ss_tot': conc_model_ss_tot,
                          'Name_pro':name_pro,'Name_lip':name_lip})

    pd_da['Start Date/Time'] = pd_da['Start Date/Time'].apply(pd.to_datetime)
    pd_da['End Date/Time'] = pd_da['End Date/Time'].apply(pd.to_datetime)

    return pd_da


# interpolation function
def interp_conc_stations(path, exp,obs,obs_tot, ID,modes):
    from scipy.interpolate import interp2d

    start_4_mod, end_4_mod = [], []
    id_camp = []
    conc_model_pol, conc_model_pro, conc_model_lip, conc_model_tot,conc_model_ss,conc_model_ss_tot = [], [], [], [],[], []
    conc_obs_pol, conc_obs_pro, conc_obs_lip, conc_obs_tot,conc_obs_ss,conc_obs_ss_tot = [], [], [], [], [],[]
    name_lip, name_pro = [], []
    months = obs['Start Date/Time'].dt.month.values
    years = obs['Start Date/Time'].dt.year.values
    days = obs['Start Date/Time'].dt.day.values

    for m, mo in enumerate(months):
        lat, lon = obs['Start Latitude'].values[0], obs['Start Longitude'].values[0]% 360

        # interpolation function
        f_int = interp(path,exp,days[m],mo, years[m], ['POL', 'PRO', 'LIP','SS'],modes,'tracer')
        f_interp_pol = f_int[0]
        f_interp_pro = f_int[1]
        f_interp_lip = f_int[2]
        f_interp_ss = f_int[3]
        f_interp_ss_tot = f_int[4]


        # interpolate obs coordinates
        interp_lim_start_pol = f_interp_pol(lat,lon)[0][0]
        interp_lim_start_pro = f_interp_pro(lat,lon)[0][0]
        interp_lim_start_lip = f_interp_lip(lat,lon)[0][0]
        interp_lim_start_ss = f_interp_ss(lat,lon)[0][0]
        interp_lim_start_ss_tot = f_interp_ss_tot(lat,lon)[0][0]


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

        name_pro.append('DHAA   DFAA')
        name_lip.append('Lipids  FFA_PG')

        id_camp.append(ID)


    # create new dataframe to store the model data after interpolation together with obs.
    pd_da = pd.DataFrame({'ID': id_camp, 'Start Date/Time': start_4_mod, 'End Date/Time': end_4_mod,
                          'conc_mod_poly': conc_model_pol, 'conc_mod_pro': conc_model_pro, 'conc_mod_lip': conc_model_lip,
                          'conc_obs_poly_sub': conc_obs_pol, 'conc_obs_prot_sub': conc_obs_pro,
                          'conc_obs_lipi_sub': conc_obs_lip,
                          'conc_obs_tot_sub': conc_obs_tot,'conc_mod_tot': conc_model_tot,
                          'conc_obs_ss': conc_obs_ss, 'conc_mod_ss': conc_model_ss,
                          'conc_obs_ss_tot': conc_obs_ss_tot, 'conc_mod_ss_tot': conc_model_ss_tot,
                          'Name_pro':name_pro,'Name_lip':name_lip})

    pd_da['Start Date/Time'] = pd_da['Start Date/Time'].apply(pd.to_datetime)
    pd_da['End Date/Time'] = pd_da['End Date/Time'].apply(pd.to_datetime)

    return pd_da

