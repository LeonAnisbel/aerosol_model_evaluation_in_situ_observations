import plot_functions
import pandas as pd
import numpy as np


def pd_combine_pol(exp, conc_pd, plot_dir):
    # combining all conc obs and model into a datframe (new_conc_pol) for box plot
    obs_conc_pol = (conc_pd['conc_mod_poly'].to_list() +
                    conc_pd['conc_obs_poly_sub'].to_list())

    obs_conc_all = (conc_pd['conc_mod_tot'].to_list() +
                    conc_pd['conc_obs_tot_sub'].to_list())

    obs_conc_ss = (conc_pd['conc_mod_ss'].to_list() +
                   conc_pd['conc_obs_ss'].to_list())

    obs_conc_ss_tot = (conc_pd['conc_mod_ss_tot'].to_list() +
                       conc_pd['conc_obs_ss_tot'].to_list())

    ID_all = (conc_pd['ID'].to_list() +
              conc_pd['ID'].to_list())

    mod_obs = (len(conc_pd['ID'].values) * ['Model'] +
               len(conc_pd['ID'].values) * ['Measurements'])

    new_conc_pd = pd.DataFrame({'PCHO Concentration (µg m$^{-3}$)': obs_conc_pol,
                                'Measurement': ID_all, '': mod_obs,
                                'Total aerosol Concentration (µg m$^{-3}$)': obs_conc_all,
                                'SS Concentration (µg m$^{-3}$)': obs_conc_ss,
                                'Total SS Concentration (µg m$^{-3}$)': obs_conc_ss_tot})

    new_conc_pd_tot_nan = new_conc_pd.dropna(subset=['Total aerosol Concentration (µg m$^{-3}$)'])
    new_conc_pd_pol_nan = new_conc_pd.dropna(subset=['PCHO Concentration (µg m$^{-3}$)'])
    new_conc_pd_ss_nan = new_conc_pd.dropna(subset=['SS Concentration (µg m$^{-3}$)'])
    new_conc_pd_ss_tot_nan = new_conc_pd.dropna(subset=['Total SS Concentration (µg m$^{-3}$)'])

    new_conc_pd_ss_nan.to_pickle('pd_files/ss_acc_conc.pkl')
    new_conc_pd_ss_tot_nan.to_pickle('pd_files/ss_tot_conc.pkl')


    # plot_functions.box_plot(new_conc_pd_tot_nan, plot_dir,
    #                         "Total aerosol Concentration (µg m$^{-3}$)",
    #                         "",
    #                         f'tot_conc_box_{exp}')
    # plot_functions.box_plot(new_conc_pd_pol_nan, plot_dir,
    #                         "PCHO Concentration (µg m$^{-3}$)",
    #                         "PCHO",
    #                         f'CCHO_conc_box_{exp}')
    # plot_functions.box_plot(new_conc_pd_ss_nan, plot_dir,
    #                         "SS Concentration (µg m$^{-3}$)",
    #                         "Sea Salt",
    #                         f'SS_conc_box_{exp}')
    # plot_functions.box_plot(new_conc_pd_ss_tot_nan, plot_dir,
    #                         "Total SS Concentration (µg m$^{-3}$)",
    #                         "Sea Salt",
    #                         f'SS_conc_box_{exp}')

    # Plotting box plots with statistics
    names = ['NAO', 'CVAO', 'WAP', 'SVD14', 'SVD15', 'SVD18']
    dict_comp = dict((name, {'pol_mod': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_poly_sub'])[
                                 'conc_mod_poly'].to_list(),
                             'pol_obs': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_poly_sub'])[
                                 'conc_obs_poly_sub'].to_list(),
                             'tot_mod': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_tot_sub'])[
                                 'conc_mod_tot'].to_list(),
                             'tot_obs': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_tot_sub'])[
                                 'conc_obs_tot_sub'].to_list(),
                             'ss_mod': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_ss'])[
                                 'conc_mod_ss'].to_list(),
                             'ss_obs': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_ss'])[
                                 'conc_obs_ss'].to_list(),
                             'ssto_mod': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_ss_tot'])[
                                 'conc_mod_ss'].to_list(),
                             'ssto_obs': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_ss_tot'])[
                                 'conc_obs_ss_tot'].to_list()
                             }) for name in names)

    plot_functions.box_plot_stat(new_conc_pd_pol_nan, dict_comp,
                                 'pol',
                                 f'pcho_stat_{exp}',
                                 [1e-7, 1e0],
                                 "PCHO Concentration (µg m$^{-3}$)",
                                 plot_dir,
                                 0.5,
                                 0.05)
    plot_functions.box_plot_stat(new_conc_pd_tot_nan, dict_comp,
                                 'tot',
                                 f'total_stat_{exp}',
                                 [1e-4, 1e1],
                                 "Total aerosol Concentration (µg m$^{-3}$)",
                                 plot_dir,
                                 5,
                                 0.05)
    plot_functions.box_plot_stat(new_conc_pd_ss_nan, dict_comp,
                                 'ss',
                                 f'ss_stat_{exp}',
                                 [1e-3, 1e1],
                                 "SS Concentration (µg m$^{-3}$)",
                                 plot_dir,
                                 5,
                                 0.05)
    plot_functions.box_plot_stat(new_conc_pd_ss_tot_nan, dict_comp,
                                 'ssto',
                                 f'ss_stat_{exp}',
                                 [1e-3, 1e2],
                                 "Total SS Concentration (µg m$^{-3}$)",
                                 plot_dir,
                                 50,
                                 0.05)


def pd_combine_pr_li(pd_obs_mo_cvao, plot_dir):
    # Compiling all conc obs and model into a datframe (new_conc_pro_lip) for box plot
    obs_conc_pro_lip = pd_obs_mo_cvao['conc_mod_pro'].to_list() + pd_obs_mo_cvao['conc_mod_lip'].to_list() + \
                       pd_obs_mo_cvao['conc_obs_prot_sub'].to_list() + pd_obs_mo_cvao['conc_obs_lipi_sub'].to_list()

    ID_all = (pd_obs_mo_cvao['Name_pro'].to_list() + pd_obs_mo_cvao['Name_lip'].to_list() +
              pd_obs_mo_cvao['Name_pro'].to_list() + pd_obs_mo_cvao['Name_lip'].to_list())  # Prot or lip

    mod_obs = len(pd_obs_mo_cvao['Name_pro'].values) * ['Model'] + len(pd_obs_mo_cvao['Name_pro'].values) * ['Model'] + \
              len(pd_obs_mo_cvao['Name_pro'].values) * ['Measurements'] + len(pd_obs_mo_cvao['Name_pro'].values) * [
                  'Measurements']

    new_conc_pd_pl = pd.DataFrame({'Aerosol Concentration': obs_conc_pro_lip, 'Species': ID_all, '': mod_obs})

    plot_functions.box_plot_lip(new_conc_pd_pl, plot_dir)


def pd_combine_pro_lip(conc_pd, plot_dir):
    # combining all conc obs and model into a datframe (new_conc_pol) for box plot
    conc_pd_pro = conc_pd.dropna(subset=['conc_obs_prot_sub'])
    conc_pd_lip = conc_pd.dropna(subset=['conc_obs_lipi_sub'])

    obs_conc_pro = (conc_pd_pro['conc_mod_pro'].to_list() +
                    list(conc_pd_pro['conc_obs_prot_sub'].values))

    obs_conc_lip = (conc_pd_lip['conc_mod_lip'].to_list() +
                    conc_pd_lip['conc_obs_lipi_sub'].to_list())

    ID_pro = (conc_pd_pro['ID'].to_list() +
              conc_pd_pro['ID'].to_list())

    ID_lip = (conc_pd_lip['ID'].to_list() +
              conc_pd_lip['ID'].to_list())

    mod_obs_pro = (len(conc_pd_pro['ID'].values) * ['Model'] +
                   len(conc_pd_pro['ID'].values) * ['Measurements'])

    mod_obs_lip = (len(conc_pd_lip['ID'].values) * ['Model'] +
                   len(conc_pd_lip['ID'].values) * ['Measurements'])

    new_conc_pd_pro = pd.DataFrame({'Proteins Concentration (µg m$^{-3}$)': obs_conc_pro,
                                    'Measurement': ID_pro,
                                    '': mod_obs_pro,
                                    })

    new_conc_pd_lip = pd.DataFrame({'Lipids Concentration (µg m$^{-3}$)': obs_conc_lip,
                                    'Measurement': ID_lip,
                                    '': mod_obs_lip,
                                    })

    new_conc_pd_lip_nan = new_conc_pd_lip.dropna(subset=['Lipids Concentration (µg m$^{-3}$)'])
    new_conc_pd_pro_nan = new_conc_pd_pro.dropna(subset=['Proteins Concentration (µg m$^{-3}$)'])

    # Plotting box plots with statistics
    names = ['CVAO', 'SVD14', 'SVD15', 'SVD18']
    dict_comp_pro = dict((name, {
        'pro_mod': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_prot_sub'])['conc_mod_pro'].to_list(),
        'pro_obs': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_prot_sub'])['conc_obs_prot_sub'].to_list(),
        }) for name in names)

    name = 'CVAO'
    dict_comp_lip = {name: {
        'lip_mod': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_lipi_sub'])['conc_mod_lip'].to_list(),
        'lip_obs': conc_pd[conc_pd['ID'] == name].dropna(subset=['conc_obs_lipi_sub'])['conc_obs_lipi_sub'].to_list()}}

    # print(dict_comp_lip)
    plot_functions.box_plot_stat(new_conc_pd_pro_nan,
                                 dict_comp_pro,
                                 'pro',
                                 'proteins_stat',
                                 [1e-6, 1e-1],
                                 "Proteins Concentration (µg m$^{-3}$)",
                                 plot_dir,
                                 0.05,
                                 -0.15)
    plot_functions.box_plot_stat(new_conc_pd_lip_nan,
                                 dict_comp_lip,
                                 'lip',
                                 'lipids_stat',
                                 [1e-4, 1e0],
                                 "Lipids Concentration (µg m$^{-3}$)",
                                 plot_dir,
                                 0.5,
                                 -0.4)


def create_dataframe(da_pd_nan, mod, obs, mod_ss, obs_ss, mac_na):
    da_pd = da_pd_nan.dropna(subset=[obs])
    print(da_pd)
    # combining all omf obs and model into a datframe (new_omf_pol) for box plot
    obs_pd = (list(da_pd[mod].values) +
              list(da_pd[obs].values))

    obs_ss_pd = (list(da_pd[mod_ss].values) +
              list(da_pd[obs_ss].values))

    id_all = (da_pd['ID'].to_list() +
              da_pd['ID'].to_list())

    mod_obs = (len(da_pd['ID'].values) * ['Model'] +
               len(da_pd['ID'].values) * ['Observation'])

    mac = len(mod_obs) * [mac_na]

    new_omf_pd = pd.DataFrame({'Aerosol Concentration (µg m$^{-3}$)': obs_pd,
                               'SS Concentration (µg m$^{-3}$)': obs_ss_pd,
                               'Measurements': id_all,
                               '': mod_obs,
                               'Macromolecules': mac})

    # new_omf_pd_tot_nan = new_omf_pd.dropna(subset=['Aerosol omf'])
    # new_omf_pd_po_nan = new_omf_po.dropna(subset=['Aerosol Concentration (µg m$^{-3}$)'])

    return new_omf_pd


def rename_func(data_pd, na, new_na):
    pd_new = data_pd
    list_col = data_pd['Measurements'].to_list()
    for i in range(len(list_col)):
        if list_col[i] == na:
            list_col[i] = new_na
    pd_new = pd_new.drop(columns='Measurements')
    pd_new['Measurements'] = list_col

    return pd_new


def pd_combine_all(dicc_po_to, dicc_pr, dicc_li, plot_dir):
    mac_names = ['PCHO|CCHO', 'DCAA|FAA', 'PL|PG', '(PCHO+DCAA+PL)|OM']
    conc_pd_po = create_dataframe(dicc_po_to, 'conc_mod_poly', 'conc_obs_poly_sub', mac_names[0])
    conc_pd_pr = create_dataframe(dicc_pr, 'conc_mod_pro', 'conc_obs_prot_sub', mac_names[1])
    conc_pd_li = create_dataframe(dicc_li, 'conc_mod_lip', 'conc_obs_lipi_sub', mac_names[2])
    conc_pd_to = create_dataframe(dicc_po_to, 'conc_mod_tot', 'conc_obs_tot_sub', mac_names[3])

    conc_pd_li_new = rename_func(conc_pd_li, 'CVAO', 'CVAO ')
    conc_pd_pr_new = rename_func(conc_pd_pr, 'CVAO', 'CVAO  ')
    conc_pd_to1 = rename_func(conc_pd_to, 'CVAO', 'CVAO   ')
    conc_pd_to2 = rename_func(conc_pd_to1, 'NAO', 'NAO ')
    conc_pd_to_new = rename_func(conc_pd_to2, 'WAP', 'WAP ')

    plot_functions.box_plot_vert(pd.concat([conc_pd_po, conc_pd_pr_new, conc_pd_li_new, conc_pd_to_new]),
                                 mac_names, ['pol', 'pro', 'lip'],
                                 'All_groups', [1e-6, 1e1])

    conc_pd_po.to_pickle('pd_files/poly_conc.pkl')
    conc_pd_pr_new.to_pickle('pd_files/prot_conc.pkl')
    conc_pd_li_new.to_pickle('pd_files/lipi_conc.pkl')
    conc_pd_to_new.to_pickle('pd_files/tot_conc.pkl')


def pd_combine_group(dicc_va, mac_names, mod_key_na, obs_key_na, id_na):
    conc_pd = create_dataframe(dicc_va,
                               mod_key_na,
                               obs_key_na,
                               'conc_mod_ss',
                               'conc_obs_ss',
                               mac_names)

    if id_na == 'lipi':
        conc_pd_new = rename_func(conc_pd, 'CVAO', 'CVAO ')
    if id_na == 'prot':
        conc_pd_new = rename_func(conc_pd, 'CVAO', 'CVAO  ')
    if id_na == 'tot':
        conc_pd_to1 = rename_func(conc_pd, 'CVAO', 'CVAO   ')
        conc_pd_to2 = rename_func(conc_pd_to1, 'NAO', 'NAO ')
        conc_pd_new = rename_func(conc_pd_to2, 'WAP', 'WAP ')
    if id_na == 'poly':
        conc_pd_new = conc_pd

    # plot_functions.box_plot_vert(pd.concat([conc_pd_po, conc_pd_pr_new, conc_pd_li_new, conc_pd_to_new]),
    #                              mac_names, ['pol', 'pro', 'lip'],
    #                              'All_groups', [1e-6, 1e1])

    conc_pd_new.to_pickle(f'pd_files/{id_na}_conc.pkl')
