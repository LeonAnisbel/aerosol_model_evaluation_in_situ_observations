import xarray as xr
import glob, os
from matplotlib import cm
import mod_interp_obs_concat
import plot_functions
import numpy as np
import codecs
import pandas as pd
import itertools

################################
def read_obs_data_loc():
    #### Reading station locations
    main_dir = '/work/bb1005/b381361/echam_postproc/'
    plot_dir = '../OMF_plots/'
    loc_dir = 'Aerosol_sample_coordinates/'
    files = ['PASCAL_lat_lon_aer.csv', 'PI_ICE_lat_lon_aer.csv']
    pol_data = 'AER_OMF_pol_lip_pro_all_sizes.csv'

    doc = codecs.open(main_dir + loc_dir + files[0], 'r', 'UTF-8')  #open for reading with "universal" type set
    PASCAL_loc = pd.read_csv(doc, sep=',')

    doc = codecs.open(main_dir + loc_dir + files[1], 'r', 'UTF-8')  #open for reading with "universal" type set
    PI_ICE_loc = pd.read_csv(doc, sep=',')

    # converting 'Date/Time' column to datetime data type
    PASCAL_loc['Date/Time'] = PASCAL_loc['Date/Time'].apply(pd.to_datetime)
    PI_ICE_loc['Date/Time'] = PI_ICE_loc['Date/Time'].apply(pd.to_datetime)

    doc = codecs.open(main_dir + pol_data, 'r', 'UTF-8')  # open for reading with "universal" type set
    data = pd.read_csv(doc, sep=',')

    # converting 'Date/Time' column to datetime data type
    data['Start Date/Time'] = data['Start Date/Time'].apply(pd.to_datetime)
    data['End Date/Time'] = data['End Date/Time'].apply(pd.to_datetime)

    # Save data for PI_ICE
    data_PI_ICE = data[:24]
    data_PI_ICE_subm = data_PI_ICE[:8]  # data for submicron aerosol size
    data_PI_ICE_super = data_PI_ICE[8:16]  # data for supermicron aerosol size
    data_PI_ICE_tot = data_PI_ICE[16:]  # data for all aerosol sizes

    # Save data for PASCAL
    data_PASCAL = data[24:33]
    data_PASCAL_subm = data_PASCAL[:3]  # data for submicron aerosol size
    data_PASCAL_super = data_PASCAL[3:6]  # data for supermicron aerosol size
    data_PASCAL_tot = data_PASCAL[6:]  # data for all aerosol sizes

    data_CVAO = data[33:69]  # data for PM1 aerosol sizes
    data_CVAO_subm = data_CVAO[:28]
    # data_CVAO_super = data_CVAO[28:36]
    data_CVAO_tot = data_CVAO[28:]

    data_SVAL = data[69:207]
    data_SVAL_14_15_subm = data_SVAL[:23]
    data_SVAL_18_19_subm = data_SVAL[23:45]
    data_SVAL_14_15_tot = data_SVAL[45:69]
    data_SVAL_18_19_tot = data_SVAL[69:]

    data_RS = data[207:]

    dates = []
    for i, d in enumerate(data['Start Date/Time'].dt.month):
        dates.append([d, data['Start Date/Time'].dt.year[i]])
        dates.append([data['End Date/Time'].dt.month[i], data['End Date/Time'].dt.year[i]])
    dates.sort()
    dates = list(k for k, _ in itertools.groupby(dates))

    PI_ICE = [PI_ICE_loc, data_PI_ICE, data_PI_ICE_subm, data_PI_ICE_super, data_PI_ICE_tot]
    PASCAL = [PASCAL_loc, data_PASCAL, data_PASCAL_subm, data_PASCAL_super, data_PASCAL_tot]
    CVAO = [data_CVAO, data_CVAO_subm, data_CVAO_tot]
    SVAL_14_15 = [data_SVAL, data_SVAL_14_15_subm, data_SVAL_14_15_tot]
    SVAL_18_19 = [data_SVAL, data_SVAL_18_19_subm, data_SVAL_18_19_tot]

    return dates, PASCAL, PI_ICE, CVAO, SVAL_14_15, SVAL_18_19, data_RS


def read_model_spec_data(file):
    return xr.open_mfdataset(file, concat_dim='time', combine='nested')


#reading data
def read_nc_ds(files, path):
    data_model = []
    data_month = [[] for i in range(12)]
    da_month_mean = [[] for i in range(12)]

    pp = len(path)

    for fi in files:
        if int(fi[pp + 11:pp + 15]) == 2015 or int(
                fi[pp + 11:pp + 15]) == 2010:  # or int(fi[-16:-12]) == 2015:# or int(fi[-16:-12]) == 2019:
            da = xr.open_mfdataset(fi)
            data = da.mean(dim='time')
            data_model.append(data)
            for i in range(12):
                if int(fi[pp + 15:pp + 17]) - 1 == i:
                    data_month[i].append(data)

    for l, da in enumerate(data_month):
        if len(da) > 0:
            dada = xr.concat(da, dim='time')
            da_month_mean[l] = dada.mean(dim='time')

    return xr.concat(data_model, dim='time'), data_month, da_month_mean
