import global_vars
import mod_interp_obs_concat
import read_data_functions
import pandas as pd
import sys

if sys.argv[1] == 'all':
    data_mo_all_stations, dict_metadata = read_data_functions.read_PMOA_all_stations()
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
    data_MH_15, days,months,years = read_data_functions.read_data()
    lat, lon = 53.3333, -9.9 % 360
    pd_MH15 = mod_interp_obs_concat.interp_conc_arctic_stations('ac3_arctic',
                                                      data_MH_15,
                                                      'MH15',
                                                      lat,
                                                      lon,
                                                      days,
                                                      months,
                                                      years)
    pd_MH15.to_pickle(f'pd_files/MH15_conc_{global_vars.exp_name}.pkl')
