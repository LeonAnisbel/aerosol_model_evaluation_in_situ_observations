import mod_interp_obs_concat, global_vars, combine_pd
import os
import pandas as pd
exp = global_vars.exp_name
data_dir = global_vars.data_directory

try:
    os.mkdir('plots')
except OSError:

    pass
try:
    os.mkdir('plots/Sfc_conc_plots')
except OSError:
    pass

plot_dir = './plots/'
plot_dir_sfc = plot_dir + 'Sfc_conc_plots/global_'


def interp_piice_pascal(PI_ICE, PASCAL, pd_obs_mo_cvao):
    pd_obs_mo_pi_ice = mod_interp_obs_concat.assign_loc_ship(data_dir,
                                                             exp,
                                                             PI_ICE[0],
                                                             PI_ICE[1],
                                                             PI_ICE[3],
                                                             'WAP')
    pd_obs_mo_pascal = mod_interp_obs_concat.assign_loc_ship(data_dir,
                                                             exp,
                                                             PASCAL[0],
                                                             PASCAL[1],
                                                             PASCAL[3],
                                                             'NAO')
    print('Finished interpolation WAP and NAO')
    conc_all_po = pd.concat([pd_obs_mo_pascal, pd_obs_mo_cvao, pd_obs_mo_pi_ice])
    combine_pd.pd_combine_group(conc_all_po,
                                'PCHO|CCHO',
                                'conc_mod_poly',
                                'conc_obs_poly_sub',
                                'poly')

    combine_pd.pd_combine_group(conc_all_po,
                                '(PCHO+DCAA+PL)|OM',
                                'conc_mod_tot',
                                'conc_obs_tot_sub',
                                'tot')


def interp_svd_rs(SVAL_14, SVAL_15, SVAL_18, RS_18_20, pd_obs_mo_cvao):
    # pd_obs_mo_sval_18_19 = mod_interp_obs_concat.interp_conc_stations(data_dir, exp, SVAL_18_19[1], SVAL_18_19[2],
    #                                                                   'SVD18')

    pd_obs_mo_sval_15 = mod_interp_obs_concat.interp_conc_stations(data_dir,
                                                                   exp,
                                                                   SVAL_15[1],
                                                                   SVAL_15[2],
                                                                   'SVD15')
    pd_obs_mo_sval_14 = mod_interp_obs_concat.interp_conc_stations(data_dir,
                                                                      exp,
                                                                      SVAL_14[1],
                                                                      SVAL_14[2],
                                                                      'SVD14')
    pd_obs_mo_sval_18 = mod_interp_obs_concat.interp_conc_stations(data_dir,
                                                                   exp,
                                                                   SVAL_18[1],
                                                                   SVAL_18[2],
                                                                   'SVD18')
    pd_obs_mo_rs_18_20 = mod_interp_obs_concat.interp_conc_stations(data_dir,
                                                                    exp,
                                                                    RS_18_20,
                                                                    RS_18_20,
                                                                    'RS')



    print('Finished interpolation SVD')
    conc_all_pr = pd.concat([pd_obs_mo_cvao, pd_obs_mo_sval_14, pd_obs_mo_sval_15, pd_obs_mo_sval_18, pd_obs_mo_rs_18_20])
    combine_pd.pd_combine_group(conc_all_pr,
                                'DCAA|FAA',
                                'conc_mod_pro',
                                'conc_obs_prot_sub',
                                'prot')

def interp_cvao(CVAO):
    pd_obs_mo_cvao = mod_interp_obs_concat.interp_conc_stations(data_dir,
                                                                exp,
                                                                CVAO[1],
                                                                CVAO[2],
                                                                'CVAO')

    combine_pd.pd_combine_group(pd_obs_mo_cvao,
                                'PL|\n(PG)',
                                'conc_mod_lip',
                                'conc_obs_lipi_sub',
                                'lipi')

    print('Finished interpolation CVAO ')
    return pd_obs_mo_cvao