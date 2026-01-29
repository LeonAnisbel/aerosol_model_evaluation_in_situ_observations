import numpy as np
import xarray as xr
from utils_functions import read_data_functions
from scipy.interpolate import griddata


def read_model(path, exp, day, mo, yr, ext, monthly=False, multiyear=False):
    """
    Reads the model data for the specific month and year of the observational data
    :param path: path to model data
    :param exp: experiment name
    :param day: day of the observation
    :param mo: month of the observation
    :param yr: year of the observation
    :param ext: model file extension id
    :param monthly: whether to use monthly data
    :param multiyear: whether to use multiyear data
    :return: datasets of aerosol concentration and air density
    """
    data_dir = f"{path}"
    if mo < 10:
        mo_str = f'0{mo}'
    else:
        mo_str = f'{mo}'

    files = f'{data_dir}{exp}_{yr}{mo_str}.01_{ext}.nc'
    file_ro = f'{data_dir}{exp}_{yr}{mo_str}.01_vphysc.nc'

    if multiyear:
        files = [f'{data_dir}{exp}_{y}{mo_str}.01_{ext}.nc' for y in yr]
        file_ro = [f'{data_dir}{exp}_{y}{mo_str}.01_vphysc.nc' for y in yr]

    # if os.path.exists(files):
    print('reading ', files, 'for interpolation with data month = ', mo)
    data = read_data_functions.read_model_spec_data(files)
    data_ro = read_data_functions.read_model_spec_data(file_ro)

    da_ro, da_ds = [], []
    if monthly or multiyear:
        da_m_ro = data_ro['rhoam1'].isel(lev=46)
        da_m_ds = data.isel(lev=46)
    else:
        ti_sel = [day * 2 - 1, day * 2 - 2]  # based on a 12h output
        print(data.time.values[ti_sel], ti_sel)
        for ti in ti_sel:
            da_ro.append(data_ro['rhoam1'].isel(time=ti).isel(lev=46))
            da_ds.append(data.isel(time=ti).isel(lev=46))
        da_m_ro = xr.concat(da_ro, dim='time')
        da_m_ds = xr.concat(da_ds, dim='time')

    return da_m_ds, da_m_ro


def def_box(ds, lat_obs, lon_obs):
    """
    Defines a box around the station location (lat_obs, lon_obs) to consider for the interpolation
    :var ds: xarray dataset
    :var lat_obs: latitude of station location
    :var lon_obs: longitude of station location
    :return: a dataset after selecting the box considered for the interpolation
    """
    bx_size = 5
    ds_bx = ds.where((ds.lat < lat_obs[1] + bx_size) & (ds.lat > lat_obs[0] - bx_size) &
                     (ds.lon < lon_obs[1] + bx_size) & (ds.lon > lon_obs[0] - bx_size), drop=True)
    return ds_bx


def get_mod_box(dr_aer, lat_obs, lon_obs):
    """
    Calls the function to create a box around the station location (lat_obs, lon_obs) and extract the data within this
    box from the model dataArray
    :var dr_aer: xarray dataset
    :var lat_obs: latitude of station location
    :var lon_obs: longitude of station location
    :return: latitude and longitude values defined as a box around the station location leaving out nans
    """

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


def start_interp(mod_data, mod_ro_da, var_names, lon_btw, lat_btw, mi_ma_lon, mi_ma_lat, all_modes=False):
    """
    This function will iterate through all variable (aerosol tracers) to perform the interpolation
    :var mod_data: xarray dataset with aerosol mass
    :var mod_ro_da: xarray dataset with air density
    :param var_names: list of variable names (aerosol species names in ECHAM-HAM)
    :var lat_btw: latitude of station location
    :var lon_btw: longitude of station location
    :var mi_ma_lon: longitude range (relevant for ship-based campaigns)
    :var mi_ma_lat: latitude range (relevant for ship-based campaigns)
    :param all_modes: boolean to determine if only Accumulation modes should be considered or also coarse mode
    :return: dictionary containing interpolated value with variable names as keys
    """
    interp_var_list = {}
    for va_na in var_names:
        interp_var_list[va_na] = interp_func(mod_data, mod_ro_da, va_na,
                                             lon_btw, lat_btw, mi_ma_lon,
                                             mi_ma_lat, all_modes=all_modes)

    interp_var_list['SS_tot'] = interp_func(mod_data, mod_ro_da, 'SS',
                                       lon_btw, lat_btw, mi_ma_lon,
                                       mi_ma_lat, all_modes=True)
    return interp_var_list


def interp_func(mod_ds, mod_dr_ro, var, obs_lon, obs_lat, obs_lon_mi_ma, obs_lat_mi_ma, all_modes=False):
    """
    This function will create and adapt the required grids for the interpolation (griddata)
    :var mod_ds: xarray dataset with aerosol mass
    :var mod_dr_ro: xarray dataArray with air density
    :param var: variable name (aerosol species name in ECHAM-HAM)
    :var obs_lat: latitude of station location
    :var obs_lon: longitude of station location
    :var obs_lon_mi_ma: longitude range (relevant for ship-based campaigns)
    :var obs_lat_mi_ma: latitude range (relevant for ship-based campaigns)
    :param all_modes: boolean to determine if only Accumulation modes should be considered or also coarse mode
    :return: interpolated value
    """
    if all_modes:
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
