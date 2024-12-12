import global_vars
import mod_interp_obs_concat
import read_data_functions

data_MH_15, days,months,years = read_data_functions.read_data()
lat, lon = 53.3333, 9.9 % 360
pd_MH15 = mod_interp_obs_concat.interp_conc_arctic_stations('ac3_arctic',
                                                  data_MH_15,
                                                  'MH15',
                                                  lat,
                                                  lon,
                                                  days,
                                                  months,
                                                  years)
pd_MH15.to_pickle(f'pd_files/MH15_conc_{global_vars.exp_name}.pkl')


