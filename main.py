##### Import packges
import numpy as np
import pandas as pd
import read_data_functions
import mod_interp_obs_concat
import glob, os
import combine_pd
import map_with_loc
import global_vars, var_executor

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
    plot_dir_sfc = plot_dir + 'Sfc_conc_plots/global_'

    #read_data_functions.read_do_mean_plot(data_dir, plot_dir_sfc,poles = False,make_plots = True, get_aver = True)
    #
    #read_data_functions.read_do_monthly_mean_plot(data_dir, plot_dir_sfc,poles = False,plot_month = True, plot_season = True)

    dates, PASCAL, PI_ICE, CVAO, SVAL_14, SVAL_15, SVAL_18, RS_18_20 = (
        read_data_functions.read_obs_data_loc())


    pd_obs_mo_cvao = var_executor.interp_cvao(CVAO)
    var_executor.interp_piice_pascal(PI_ICE, PASCAL, pd_obs_mo_cvao)
    var_executor.interp_svd_rs(SVAL_14, SVAL_15, SVAL_18, RS_18_20, pd_obs_mo_cvao)

     # # combine data of poly for all stations and generate box plot with statistics
    # conc_pd_pol_tot = pd.concat(
    #      [pd_obs_mo_pascal, pd_obs_mo_cvao, pd_obs_mo_pi_ice, pd_obs_mo_sval_14_15, pd_obs_mo_sval_18, pd_obs_mo_sval_19])
    # combine_pd.pd_combine_pol(exp, conc_pd_pol_tot, plot_dir)
    # #print('Finished plotting SS and POL')
    # #
    # #
    # # # combine data of lip and prot for all stations and generate box plot with statistics
    # conc_pd_pr_li_tot = pd.concat([pd_obs_mo_cvao, pd_obs_mo_sval_14_15, pd_obs_mo_sval_14_15])
    # combine_pd.pd_combine_pro_lip(conc_pd_pr_li_tot, plot_dir)
    # #
    # #
    # # # generate box plot for lip and proteins over CVAO
    # combine_pd.pd_combine_pr_li(pd_obs_mo_cvao, plot_dir)
    #
    # conc_all_po = pd.concat([pd_obs_mo_pascal, pd_obs_mo_cvao, pd_obs_mo_pi_ice])
    # conc_all_pr = pd.concat([pd_obs_mo_cvao, pd_obs_mo_sval_14_15, pd_obs_mo_sval_18, pd_obs_mo_sval_19, pd_obs_mo_rs_18_20])
    # conc_all_li = pd_obs_mo_cvao
    # #
    # combine_pd.pd_combine_all(conc_all_po, conc_all_pr, conc_all_li, plot_dir)

    #map_with_loc.loc_map_plot(5,2017,PASCAL[0], PI_ICE[0], data_CVAO,plot_dir)
