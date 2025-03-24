import global_vars, combine_pd
import os
import pandas as pd

exp = global_vars.exp_name
data_dir = global_vars.data_directory
if global_vars.exp_name != 'echam_base':
    import mod_interp_obs_concat as mod_inter_obs
else:
    import mod_interp_obs_concat_echam_base as mod_inter_obs


def interp_piice_pascal(PI_ICE, PASCAL, pd_obs_mo_cvao):
    pd_obs_mo_pi_ice = mod_inter_obs.assign_loc_ship(data_dir,
                                                     exp,
                                                     PI_ICE[0],
                                                     PI_ICE[1],
                                                     PI_ICE[3],
                                                     'WAP')
    pd_obs_mo_pascal = mod_inter_obs.assign_loc_ship(data_dir,
                                                     exp,
                                                     PASCAL[0],
                                                     PASCAL[1],
                                                     PASCAL[3],
                                                     'NAO')
    print('Finished interpolation WAP and NAO')
    conc_all_po = pd.concat([pd_obs_mo_pascal, pd_obs_mo_cvao, pd_obs_mo_pi_ice])
    if global_vars.exp_name != 'echam_base':
        combine_pd.pd_combine_group(conc_all_po,
                                    'PCHO|CCHO',
                                    'conc_mod_poly',
                                    'conc_obs_poly_sub',
                                    'poly')

        # combine_pd.pd_combine_group(conc_all_po,
        #                             '(PCHO+DCAA+PL)|OM',
        #                             'conc_mod_tot',
        #                             'conc_obs_tot_sub',
        #                             'tot')

    combine_pd.pd_combine_group(conc_all_po,
                                'OC',
                                'conc_mod_oc',
                                'conc_obs_oc',
                                'oc')


def interp_svd_rs(SVAL_15, pd_obs_mo_cvao):
    if global_vars.exp_name != 'echam_base':
        import mod_interp_obs_concat as mod_inter_obs
    pd_obs_mo_sval_15 = mod_inter_obs.interp_conc_stations(data_dir,
                                                           exp,
                                                           SVAL_15[0],
                                                           SVAL_15[1],
                                                           'SVD')


    print('Finished interpolation SVD')
    conc_all_pr = pd.concat([pd_obs_mo_cvao,
                             pd_obs_mo_sval_15,
                             ])
    combine_pd.pd_combine_group(conc_all_pr,
                                'DCAA|CAA',
                                'conc_mod_pro',
                                'conc_obs_prot_sub',
                                'prot')


def interp_cvao(CVAO):
    pd_obs_mo_cvao = mod_inter_obs.interp_conc_stations(data_dir,
                                                        exp,
                                                        CVAO[0],
                                                        CVAO[1],
                                                        'CVAO')

    if global_vars.exp_name != 'echam_base':
        combine_pd.pd_combine_group(pd_obs_mo_cvao,
                                    'PL|\n(PG)',
                                    'conc_mod_lip',
                                    'conc_obs_lipi_sub',
                                    'lipi')
    # combine_pd.pd_combine_group(pd_obs_mo_cvao,
    #                             'OC',
    #                             'conc_mod_oc',
    #                             'conc_obs_oc',
    #                             'oc')

    print('Finished interpolation CVAO ')
    return pd_obs_mo_cvao
