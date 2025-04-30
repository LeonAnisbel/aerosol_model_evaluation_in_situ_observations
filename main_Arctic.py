import numpy as np

import global_vars
import mod_interp_obs_concat
import read_data_functions
import pandas as pd
import sys

if sys.argv[1] == 'all':
    data_mo_all_stations, dict_metadata, _ = read_data_functions.read_PMOA_all_stations()
    data_sta_list = []
    for sta in list(dict_metadata.keys()):
        dict_metadata[sta]['model_PMOA'] = mod_interp_obs_concat.interp_all_arctic_stations('ac3_arctic',
                                                                    data_mo_all_stations,
                                                                    sta,
                                                                    float(dict_metadata[sta]['location'][0][0]),
                                                                    float(dict_metadata[sta]['location'][0][1])%360,
                                                                    [],
                                                                    dict_metadata[sta]['months'],
                                                                    dict_metadata[sta]['years'])
        data_sta_list.append(dict_metadata[sta]['model_PMOA'])

    pd_sta_all = pd.concat(data_sta_list)
    pd_sta_all.to_pickle(f'pd_files/all_arctic_stations_conc_{global_vars.exp_name}.pkl')


if sys.argv[1] == 'MH':
    lat, lon = 53.3333, -9.9 % 360
    year = '0209'
    data_MH, days,months,years, _ = read_data_functions.read_data(year,
                                                                  global_vars.yr_exp_var_names[year],
                                                                  monthly=False)
    pd_MH0209 = mod_interp_obs_concat.interp_conc_arctic_stations('ac3_arctic',
                                                      data_MH,
                                                      'MH'+year,
                                                      lat,
                                                      lon,
                                                      days,
                                                      months,
                                                      years,
                                                      global_vars.yr_exp_var_names[year][1:])
    pd_MH0209.to_pickle(f'pd_files/MH{year}_conc_{global_vars.exp_name}.pkl')


    year = '2018'
    data_MH, days,months,years, _ = read_data_functions.read_data(year,
                                                                  global_vars.yr_exp_var_names[year],
                                                                  monthly=False)
    pd_MH15 = mod_interp_obs_concat.interp_conc_arctic_stations('ac3_arctic',
                                                      data_MH,
                                                      'MH18',
                                                      lat,
                                                      lon,
                                                      days,
                                                      months,
                                                      years,
                                                      global_vars.yr_exp_var_names[year][1:])
    pd_MH15.to_pickle(f'pd_files/MH{year}_conc_{global_vars.exp_name}.pkl')


if sys.argv[1] == 'AI':
    # from May 2003 to November 2007
    months_yr = 4*[np.arange(2004,2008)] + 7*[np.arange(2003, 2008)]
    months_yr.append(np.arange(2003, 2007))
    lat, lon = -37.8, 77.5667 % 360

    factor_OC_OM = 1.4 # conversion factor suggested by Facchini et al. 2008 for WIOC
    factor_unit =  1000 # convert from ng to ug

    # data from Sciare et al. 2009
    data_AI0307 = np.array([189, 122,  96,  68,  64,  41,  49,  47,  64,  54,  64, 130]) * factor_OC_OM / factor_unit
    data_AI0307_std = np.array([62, 53,  6, 15, 32, 28, 10, 17, 43, 30, 27, 39]) * factor_OC_OM / factor_unit
    obs_AI=pd.DataFrame({'WIOM_ug_m3': data_AI0307,
                       'WIOM_ug_m3_std': data_AI0307_std})

    pd_AI = mod_interp_obs_concat.interp_all_arctic_stations('ac3_arctic',
                                                             obs_AI,
                                                             'AI0307',
                                                              lat,
                                                              lon,
                                                             None,
                                                             np.arange(1, 13),
                                                             months_yr)

    pd_AI.to_pickle(f'pd_files/AI0307_conc_{global_vars.exp_name}.pkl')

