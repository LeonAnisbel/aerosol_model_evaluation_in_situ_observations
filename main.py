##### Import packges
import numpy as np
import pandas as pd
import read_data_functions
import mod_interp_obs_concat
import glob,os
import combine_pd
import map_with_loc
import global_vars

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # data paths
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
    plot_dir_sfc = plot_dir+'Sfc_conc_plots/global_'

    #read_data_functions.read_do_mean_plot(data_dir, plot_dir_sfc,poles = False,make_plots = True, get_aver = True)
    #
    #read_data_functions.read_do_monthly_mean_plot(data_dir, plot_dir_sfc,poles = False,plot_month = True, plot_season = True)

    dates, PASCAL, PI_ICE, CVAO, SVAL_14_15, SVAL_18_19  = (
        read_data_functions.read_obs_data_loc())
    #PI_ICE_loc,data_PI_ICE,data_PI_ICE_subm,data_PI_ICE_super,data_PI_ICE_tot

    pd_obs_mo_pascal = mod_interp_obs_concat.assign_loc_ship(data_dir, exp,PASCAL[0], PASCAL[2], PASCAL[4],'NAO')
    pd_obs_mo_pi_ice = mod_interp_obs_concat.assign_loc_ship(data_dir, exp,PI_ICE[0],PI_ICE[2],PI_ICE[4],'WAP')
    print('Finished interpolation WAP and NAO')

    # #
    pd_obs_mo_cvao = mod_interp_obs_concat.interp_conc_stations(data_dir, exp,CVAO[1],CVAO[2],'CVAO')
    pd_obs_mo_sval_14_15 = mod_interp_obs_concat.interp_conc_stations(data_dir, exp,SVAL_14_15[1],SVAL_14_15[2],'SVD14')
    pd_obs_mo_sval_18_19 = mod_interp_obs_concat.interp_conc_stations(data_dir, exp,SVAL_18_19[1],SVAL_18_19[2],'SVD18')

    print('Finished interpolation CVAO and SVD')
    # #
    # # combine data of poly for all stations and generate box plot with statistics
    conc_pd_pol_tot = pd.concat([pd_obs_mo_pascal,pd_obs_mo_cvao,pd_obs_mo_pi_ice,pd_obs_mo_sval_14_15,pd_obs_mo_sval_18_19])
    combine_pd.pd_combine_pol(exp,conc_pd_pol_tot,plot_dir)
    print('Finished plotting SS and POL')
    #
    #
    # # combine data of lip and prot for all stations and generate box plot with statistics
    conc_pd_pr_li_tot = pd.concat([pd_obs_mo_cvao,pd_obs_mo_sval_14_15,pd_obs_mo_sval_14_15])
    combine_pd.pd_combine_pro_lip(conc_pd_pr_li_tot, plot_dir)
    #
    #
    # # generate box plot for lip and proteins over CVAO
    combine_pd.pd_combine_pr_li(pd_obs_mo_cvao, plot_dir)


    conc_all_po = pd.concat([pd_obs_mo_pascal, pd_obs_mo_cvao, pd_obs_mo_pi_ice])
    conc_all_pr = pd.concat([pd_obs_mo_cvao, pd_obs_mo_sval_14_15, pd_obs_mo_sval_18_19])
    conc_all_li = pd_obs_mo_cvao

    combine_pd.pd_combine_all(conc_all_po, conc_all_pr, conc_all_li, plot_dir)



    # map_with_loc.loc_map_plot(5,2017,PASCAL[0], PI_ICE[0], data_CVAO,plot_dir)

