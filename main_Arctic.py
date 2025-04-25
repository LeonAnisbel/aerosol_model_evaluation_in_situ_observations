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

    exit()


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


