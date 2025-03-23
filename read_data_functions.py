import xarray as xr
import global_vars
import codecs
import pandas as pd
import itertools


################################
def read_obs_data_loc(main_dir, loc_dir):
    #### Reading station locations
    files = ['PASCAL_lat_lon_aer.csv', 'PI_ICE_lat_lon_aer.csv']
    pol_data = 'AER_OMF_pol_lip_pro_all_sizes.csv'

    doc = codecs.open(main_dir + loc_dir + files[0], 'r', 'UTF-8')  #open for reading with "universal" type set
    PASCAL_loc = pd.read_csv(doc, sep=',')

    doc = codecs.open(main_dir + loc_dir + files[1], 'r', 'UTF-8')  #open for reading with "universal" type set
    PI_ICE_loc = pd.read_csv(doc, sep=',')

    # converting 'Date/Time' column to datetime data type
    print(PASCAL_loc['Date/Time'].values)
    PASCAL_loc['Date/Time'] = [pd.Timestamp(row).to_pydatetime() for row in PASCAL_loc['Date/Time'].values]
    PI_ICE_loc['Date/Time'] = [pd.Timestamp(row).to_pydatetime() for row in PI_ICE_loc['Date/Time'].values]

    print(PI_ICE_loc['Date/Time'])

    doc = codecs.open(main_dir + pol_data, 'r', 'UTF-8')  # open for reading with "universal" type set
    data = pd.read_csv(doc, sep=',')

    # converting 'Date/Time' column to datetime data type
    data['Start Date/Time'] = data['Start Date/Time'].apply(pd.to_datetime)
    data['End Date/Time'] = data['End Date/Time'].apply(pd.to_datetime)
    print(data['Event'])
    # Save data for PI_ICE
    data_PI_ICE_subm_1 = data[data['Station_sizes'] == 'PI-ICE_0.05_1.2']  # data for submicron aerosol size
    data_PI_ICE_subm_2 = data[data['Station_sizes'] == 'PI-ICE_0.14_1.2']  # data for submicron aerosol size
    data_PI_ICE_super = data[data['Station_sizes'] == 'PI-ICE_1.2_10']  # data for supermicron aerosol size
    data_PI_ICE_tot = data[data['Station_sizes'] == 'PI-ICE_0.05_10']  # data for all aerosol sizes

    # Save data for PASCAL
    # data_PASCAL = data[24:33]
    data_PASCAL_subm_1 = data[data['Station_sizes'] == 'PASCAL_0.05_1.2']  # data for submicron aerosol size
    data_PASCAL_subm_2 = data[data['Station_sizes'] == 'PASCAL_0.14_1.2']  # data for submicron aerosol size
    data_PASCAL_super = data[data['Station_sizes'] == 'PASCAL_1.2_10']  # data for supermicron aerosol size
    data_PASCAL_tot = data[data['Station_sizes'] == 'PASCAL_0.05_10']  # data for all aerosol sizes

    data_CVAO = data[33:69]  # data for PM1 aerosol sizes
    data_CVAO_subm_1 = data[data['Station_sizes'] == 'CVAO_0.05-1.2 µm']
    data_CVAO_subm_2 = data[data['Station_sizes'] == 'CVAO_PM1']

    # data_CVAO_super = data_CVAO[28:36]
    data_CVAO_tot = data[data['Station_sizes'] == 'CVAO_0.05-10 µm']

    data_SVAL = data[69:207]
    data_SVAL_14_subm = data[data['Station_sizes'] == 'Sval_0-1.5µm_14']
    data_SVAL_15_subm_1 = data[data['Station_sizes'] == 'Sval_0-1.5µm_15']
    data_SVAL_15_subm_2 = data[data['Station_sizes'] == 'Sval_0.95-1.5µm_15']

    data_SVAL_18_subm = data[data['Station_sizes'] == 'Sval_0-1.5µm_18']

    data_SVAL_14_tot = data[data['Station_sizes'] == 'Sval_0-10µm_14']
    data_SVAL_15_tot = data[data['Station_sizes'] == 'Sval_0-10µm_15']

    data_SVAL_18_tot = data[data['Station_sizes'] == 'Sval_0-10µm_18']

    data_RS = data[data['Station_sizes'] == 'Ross_sea_0-1.5µm']


    dates = []
    for i, d in enumerate(data['Start Date/Time'].dt.month):
        dates.append([d, data['Start Date/Time'].dt.year[i]])
        dates.append([data['End Date/Time'].dt.month[i], data['End Date/Time'].dt.year[i]])
    dates.sort()
    dates = list(k for k, _ in itertools.groupby(dates))

    PI_ICE = [PI_ICE_loc, data_PI_ICE_subm_2, data_PI_ICE_super, data_PI_ICE_tot]
    PASCAL = [PASCAL_loc, data_PASCAL_subm_2, data_PASCAL_super, data_PASCAL_tot]
    CVAO = [data_CVAO, data_CVAO_subm_2, data_CVAO_tot]
    SVAL_14 = [data_SVAL, data_SVAL_14_subm, data_SVAL_14_tot]
    SVAL_15 = [data_SVAL, data_SVAL_15_subm_2, data_SVAL_15_tot]

    # SVAL_18_19 = [data_SVAL, data_SVAL_18_19_subm, data_SVAL_18_19_tot]
    SVAL_18 = [data_SVAL, data_SVAL_18_subm, data_SVAL_18_tot]

    return dates, PASCAL, PI_ICE, CVAO, SVAL_14, SVAL_15, SVAL_18, data_RS



def read_data(yr, monthly=False):
    da_dir = global_vars.main_data_dir+'MH_PMOAseasalt_'+yr+'.csv'
    #da_dir = '/home/manuel/Downloads/'+'MH_PMOAseasalt_'+yr+'.csv'

    data_all = codecs.open(da_dir,
                           'r')
    data = pd.read_csv(data_all, sep=',')
    data_15 = data[['date', 'seasalt', 'PMOA']].copy(deep=True)

    data_15['OMF'] = (data_15['PMOA'].values /
                      (data_15['PMOA'].values +
                       data_15['seasalt'].values))    # data_15['OMF'] = data_15['MOA'].values/(data_15['MOA'].values+data_15['seasalt'].values) #*(1/0.3061)

    data_15.loc[:, ('date')] = data_15['date'].apply(pd.to_datetime, dayfirst=True)

    times = pd.to_datetime(data_15['date'], dayfirst=True)

    if monthly:
        days = []
        data_15_hr = (data_15.groupby([times.dt.year, times.dt.month])[['seasalt', 'PMOA', 'OMF']]
                      .mean())
        data_15_std = (data_15.groupby([times.dt.year, times.dt.month])[['seasalt', 'PMOA', 'OMF']]
                      .std())
    else:
        data_15_hr = (data_15.groupby([times.dt.year, times.dt.month, times.dt.day])[['seasalt', 'PMOA', 'OMF']]
                  .median())
        data_15_std = (data_15.groupby([times.dt.year, times.dt.month])[['seasalt', 'PMOA', 'OMF']]
                      .std())
        days = [i[2] for i in data_15_hr.index]

    months = [i[1] for i in data_15_hr.index]
    years = [i[0] for i in data_15_hr.index]
    print(data_15_hr)

    return data_15_hr, days,months,years, data_15_std


def read_PMOA_all_stations():
    da_dir = global_vars.main_data_dir+'ArcticOA_dataset_PBOA.csv'
    data = codecs.open(da_dir, 'r')
    data = pd.read_csv(data, sep=',')
    data_sel = data.copy(deep=True)
    data_sel['PBOA_ug_m3'] = data_sel['PBOA_ng_m-3']*1e-3 # convert to ug/m3
    data_sel['Median date'] = data_sel['Median date'].apply(pd.to_datetime, dayfirst=False)
    times = pd.to_datetime(data_sel['Median date'], dayfirst=True)

    data_sel_yr = (
        (data_sel.groupby(['Station', 'Latitude', 'Longitude', times.dt.year, times.dt.month])[['PBOA_ug_m3']])
        .median())
    data_sel_std = ((data_sel.groupby(['Station','Latitude','Longitude', times.dt.year, times.dt.month])[['PBOA_ug_m3']])
                   .std())

    # create dictionary with stations metadata
    stations = list(set([i[0] for i in data_sel_yr.index]))
    stations.sort()
    dict_yr = {}
    for sta in stations:
        dict_yr[sta] = {}
        dict_yr[sta]['location'] = list(set([(i[1], i[2]) for i in data_sel_yr.index if i[0] == sta]))
        dict_yr[sta]['years'] = list([i[3] for i in data_sel_yr.index if i[0]==sta])
        dict_yr[sta]['months'] = list([i[4] for i in data_sel_yr.index if i[0]==sta])
    return data_sel_yr, dict_yr, data_sel_std


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
