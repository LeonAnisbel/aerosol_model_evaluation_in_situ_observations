import xarray as xr
import global_vars
import codecs
import pandas as pd
import itertools

def read_obs_data_loc(main_dir, loc_dir):
    """
    Read observational data from various campaigns
    :returns dates: list of dates
             PASCAL: PASCAL campaign data list (len=4) of dataframes
                    [ship location, sub-micron size, super-micron size, total size]
             PI_ICE: PI_ICE campaign data list (len=4) of dataframes
                    [ship location, sub-micron size, super-micron size, total size]
             CVAO: CVAO data list (len=2)  of dataframes [sub-micron size, total size]
             SVAL_15: Svalbard data list (len=2) of dataframes [sub-micron size, total size]
    """

    #### Reading station locations
    files = ['PASCAL_lat_lon_aer.csv', 'PI_ICE_lat_lon_aer.csv']
    pol_data = 'AER_OMF_pol_lip_pro_all_sizes.csv'

    doc = codecs.open(main_dir + loc_dir + files[0],
                      'r',
                      'UTF-8')  #open for reading with "universal" type set
    PASCAL_loc = pd.read_csv(doc, sep=',')

    doc = codecs.open(main_dir + loc_dir + files[1],
                      'r',
                      'UTF-8')  #open for reading with "universal" type set
    PI_ICE_loc = pd.read_csv(doc, sep=',')

    # converting 'Date/Time' column to datetime data type
    PASCAL_loc['Date/Time'] = [pd.Timestamp(row).to_pydatetime() for row in PASCAL_loc['Date/Time'].values]
    PI_ICE_loc['Date/Time'] = [pd.Timestamp(row).to_pydatetime() for row in PI_ICE_loc['Date/Time'].values]

    doc = codecs.open(main_dir + pol_data, 'r', 'UTF-8')  # open for reading with "universal" type set
    data = pd.read_csv(doc, sep=',')

    # converting 'Date/Time' column to datetime data type
    data['Start Date/Time'] = data['Start Date/Time'].apply(pd.to_datetime)
    data['End Date/Time'] = data['End Date/Time'].apply(pd.to_datetime)
    # Save data for PI_ICE
    data_PI_ICE_subm_1 = data[data['Station_sizes'] == 'PI-ICE_0.05_1.2']  # data for submicron aerosol size
    data_PI_ICE_subm_2 = data[data['Station_sizes'] == 'PI-ICE_0.14_1.2']  # data for submicron aerosol size
    data_PI_ICE_super = data[data['Station_sizes'] == 'PI-ICE_1.2_10']  # data for supermicron aerosol size
    data_PI_ICE_tot = data[data['Station_sizes'] == 'PI-ICE_0.05_10']  # data for all aerosol sizes

    # Save data for PASCAL
    data_PASCAL_subm_1 = data[data['Station_sizes'] == 'PASCAL_0.05_1.2']  # data for submicron aerosol size
    data_PASCAL_subm_2 = data[data['Station_sizes'] == 'PASCAL_0.14_1.2']  # data for submicron aerosol size
    data_PASCAL_super = data[data['Station_sizes'] == 'PASCAL_1.2_10']  # data for supermicron aerosol size
    data_PASCAL_tot = data[data['Station_sizes'] == 'PASCAL_0.05_10']  # data for all aerosol sizes

    data_CVAO_subm_1 = data[data['Station_sizes'] == 'CVAO_0.05-1.2 µm']
    data_CVAO_subm_2 = data[data['Station_sizes'] == 'CVAO_PM1']

    # data_CVAO_super = data_CVAO[28:36]
    data_CVAO_tot = data[data['Station_sizes'] == 'CVAO_0.05-10 µm']

    data_SVAL_15_subm_1 = data[data['Station_sizes'] == 'Sval_0-1.5µm_15']
    data_SVAL_15_subm_2 = data[data['Station_sizes'] == 'Sval_<0.49_0.95µm_15']

    data_SVAL_15_tot = data[data['Station_sizes'] == 'Sval_0-10µm_15']


    dates = []
    for i, d in enumerate(data['Start Date/Time'].dt.month):
        dates.append([d, data['Start Date/Time'].dt.year[i]])
        dates.append([data['End Date/Time'].dt.month[i], data['End Date/Time'].dt.year[i]])
    dates.sort()
    dates = list(k for k, _ in itertools.groupby(dates))

    PI_ICE = [PI_ICE_loc, data_PI_ICE_subm_2, data_PI_ICE_super, data_PI_ICE_tot]
    PASCAL = [PASCAL_loc, data_PASCAL_subm_2, data_PASCAL_super, data_PASCAL_tot]
    CVAO = [data_CVAO_subm_2, data_CVAO_tot]
    SVAL_15 = [data_SVAL_15_subm_2, data_SVAL_15_tot]

    return dates, PASCAL, PI_ICE, CVAO, SVAL_15



def read_data(yr, var_names, monthly=False):
    """
    Read observational data from Mace Head for a certain year.
    :param yr: string with year id
    :param var_names: list of variable to read according to the year selected
    :param monthly: Boolean, whether to read compute monthly or daily data
    :return: dataframe with data, lists of days, months and years, and dataframe with standard deviation
    """
    data_15_hr = []
    da_dir = global_vars.main_data_dir+'MH_PMOAseasalt_'+yr+'.csv'
    print('Reading data from '+da_dir)
    #da_dir = '/home/manuel/Downloads/'+'MH_PMOAseasalt_'+yr+'.csv'

    date = var_names[0]
    ss = var_names[1]
    moa = var_names[2]

    data_all = codecs.open(da_dir,
                           'r')
    data = pd.read_csv(data_all, sep=',')

    if yr == '0209':
        dayfirst = False
        date_end = var_names[3]
    else:
        dayfirst=True

    data_15 = data[var_names].copy(deep=True)

    data_15['OMF'] = (data_15[moa].values /
                      (data_15[moa].values +
                       data_15[ss].values))    #compute OMF from observations
    data_15.loc[:, (date)] = data_15[date].apply(pd.to_datetime,
                                                 dayfirst=dayfirst)
    times = pd.to_datetime(data_15[date],
                           dayfirst=dayfirst) # extract times


    if yr == '0209': # years 2002-2009
        if monthly: # compute multiannual monthly mean
            data_15.loc[:, (date_end)] = data_15[date_end].apply(pd.to_datetime,
                                                                 dayfirst=dayfirst)
            data_15['dates'] = data_15.apply(
                                lambda row: pd.date_range(row[date],
                                                          row[date_end],
                                                          freq='D'),
                                axis=1)

            daily = data_15.explode('dates').rename(columns={'dates': 'date'})
            data_15_hr = (daily
                       .groupby(daily['date'].dt.to_period('M'))[[ss, moa, 'OMF']]
                       .mean())
            start = []

            for p in data_15_hr.index:
                start.append(p.to_timestamp(how='start'))  # Timestamp('2002-03-01 00:00:00')
            data_15_hr['start'] = start
            data_15_hr['start'] = data_15_hr['start'].apply(pd.to_datetime)
            times = data_15_hr['start']
            data_15_std = (data_15_hr.groupby([times.dt.month])[[ss, moa, 'OMF']]
                          .std())
            data_15_hr = (data_15_hr.groupby([times.dt.month])[[ss, moa, 'OMF']]
                          .mean())

        else:
            data_15_hr = data_15
            data_15_std = None
            data_15_hr[date] = data_15_hr[date].apply(pd.to_datetime)
            data_15_hr[date_end] = data_15_hr[date_end].apply(pd.to_datetime)

    if yr != '0209':
        if monthly:
            days = []
            data_15_hr = (data_15.groupby([times.dt.year, times.dt.month])[[ss, moa, 'OMF']]
                          .mean())
            data_15_std = (data_15.groupby([times.dt.year, times.dt.month])[[ss, moa, 'OMF']]
                          .std())
        else:
            data_15_hr = (data_15.groupby([times.dt.year, times.dt.month, times.dt.day])[[ss, moa, 'OMF']]
                      .mean())
            data_15_std = (data_15.groupby([times.dt.year, times.dt.month])[[ss, moa, 'OMF']]
                          .std())
            days = [i[2] for i in data_15_hr.index]


    # save months, years, days as lists
    if yr == '0209':
        if monthly:
            months = data_15_hr.index.to_list
            years = None
            days = None
        else:
            months = data_15_hr[date].dt.month.values
            days = data_15_hr[date].dt.day.values
            years = data_15_hr[date].dt.year.values
    else:
        months = [i[1] for i in data_15_hr.index.to_numpy()]
        years = [i[0] for i in data_15_hr.index.to_numpy()]

    return data_15_hr, days, months, years, data_15_std


def read_PMOA_all_stations():
    """
    Read organic aerosols from multiple Arctic stations
    :return: dataframe with monthly averaged values, dictionary with stations metadata and dataframe with monthly std
    """
    da_dir = global_vars.main_data_dir+'ArcticOA_dataset_PBOA.csv'
    data = codecs.open(da_dir, 'r')
    data = pd.read_csv(data, sep=',')
    data_sel = data.copy(deep=True)
    data_sel['PBOA_ug_m3'] = data_sel['PBOA_ng_m-3']*1e-3 # convert to ug/m3
    data_sel['Median date'] = data_sel['Median date'].apply(pd.to_datetime, dayfirst=False)
    times = pd.to_datetime(data_sel['Median date'], dayfirst=True)

    data_sel_yr = (
        (data_sel.groupby(['Station', 'Latitude', 'Longitude', times.dt.year, times.dt.month])[['PBOA_ug_m3']])
        .mean())
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
    """
    This function reads (with dask) the data from a certain file type group
    :param file: files to read
    :return: dataset
    """
    return xr.open_mfdataset(file,
                             concat_dim='time',
                             combine='nested')


def read_nc_ds(files, path):
    """
    Reads model netcdf files from a certain path
    :param files: files or file type to read
    :param path: path to access files
    :return: Dataset with model data, list of months and list of datasets as monthly averages
    """
    data_model = []
    data_month = [[] for i in range(12)]
    da_month_mean = [[] for i in range(12)]

    pp = len(path)

    for fi in files: # Select and read specific time steps from model output
        if int(fi[pp + 11:pp + 15]) == 2015 or int(
                fi[pp + 11:pp + 15]) == 2010:
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
