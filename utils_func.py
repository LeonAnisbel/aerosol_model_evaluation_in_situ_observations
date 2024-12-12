import numpy as np
import xarray as xr
import os
import read_data_functions
from scipy.interpolate import griddata


def read_model(path, exp, day, mo, yr, ext):
    data_dir = f"{path}"
    if mo < 10:
        mo_str = f'0{mo}'
    else:
        mo_str = f'{mo}'
    
    print(yr, mo_str, day)
    files = f'{data_dir}{exp}_{yr}{mo_str}.01_{ext}.nc'
    file_ro = f'{data_dir}{exp}_{yr}{mo_str}.01_vphysc.nc'

    if os.path.exists(files):
        print('reading ', files, 'for interpolation with data month = ', mo)
        data = read_data_functions.read_model_spec_data(files)
        data_ro = read_data_functions.read_model_spec_data(file_ro)
        ti_sel = [day*2-1, day*2]
        print(ti_sel)
        da_ro, da_ds = [], []
        for ti in ti_sel:
            da_ro.append(data_ro['rhoam1'].isel(time=ti).isel(lev=46))
            da_ds.append(data.isel(time=ti).isel(lev=46))
        da_m_ro = xr.concat(da_ro, dim='time')
        da_m_ds = xr.concat(da_ds, dim='time')

    return da_m_ds, da_m_ro


def def_box(ds, lat_obs, lon_obs):
    bx_size = 5
    ds_bx = ds.where((ds.lat < lat_obs[1] + bx_size) & (ds.lat > lat_obs[0] - bx_size) &
                     (ds.lon < lon_obs[1] + bx_size) & (ds.lon > lon_obs[0] - bx_size), drop=True)
    return ds_bx


def get_mod_box(dr_aer, lat_obs, lon_obs):
    dr_aer_bx = def_box(dr_aer, lat_obs, lon_obs)

    # variables for the interpolation
    lat_mod = dr_aer_bx.lat.values
    lon_mod = dr_aer_bx.lon.values
    model_lo_la = []
    model_data = []
    for la in range(len(dr_aer_bx.values)):
        for lo in range(len(dr_aer_bx[0].values)):
            model_lo_la.append([lon_mod[lo], lat_mod[la]])
            model_data.append(float(dr_aer_bx[la][lo].values))

    # clean nans
    model_lo_la_notnan = []
    model_data_notnan = []
    for idx in range(len(model_data)):
        if model_data[idx] > 0:
            model_lo_la_notnan.append(model_lo_la[idx])
            model_data_notnan.append(model_data[idx])

    return model_lo_la_notnan, model_data_notnan


def start_interp(mod_data, mod_ro_da, var_names, lon_btw, lat_btw, mi_ma_lon, mi_ma_lat):
    interp_var_list = {}
    for va_na in var_names:
        interp_var_list[va_na] = interp_func(mod_data, mod_ro_da, va_na,
                                             lon_btw, lat_btw, mi_ma_lon,
                                             mi_ma_lat)

    interp_var_list['SS_tot'] = interp_func(mod_data, mod_ro_da, 'SS',
                                       lon_btw, lat_btw, mi_ma_lon,
                                       mi_ma_lat, both_modes=True)
    return interp_var_list


def interp_func(mod_ds, mod_dr_ro, var, obs_lon, obs_lat, obs_lon_mi_ma, obs_lat_mi_ma, both_modes=False):
    if both_modes:
        dr_var = (mod_ds[f'{var}_AS'] + mod_ds[f'{var}_CS'])
    else:
        dr_var = mod_ds[f'{var}_AS']

    unit_factor = 1e9
    dr_aer = dr_var * mod_dr_ro * unit_factor
    dr_aer = dr_aer.mean(dim='time', skipna=True)

    points, values = get_mod_box(dr_aer, obs_lat_mi_ma, obs_lon_mi_ma)
    grid_lon, grid_lat = np.meshgrid(obs_lon, obs_lat)
    f = griddata(points, values, (grid_lon, grid_lat), method='cubic')
    return f
