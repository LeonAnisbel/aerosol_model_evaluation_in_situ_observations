import global_vars
import plot_functions
import pandas as pd
import numpy as np

def create_dataframe(da_pd_nan, mod, obs, mod_ss, obs_ss, mod_oc, obs_oc, mac_na):
    da_pd = da_pd_nan.dropna(subset=[obs])
    print(da_pd)
    # combining all omf obs and model into a datframe (new_omf_pol) for box plot
    obs_pd = (list(da_pd[mod].values) +
              list(da_pd[obs].values))

    ss_pd = (list(da_pd[mod_ss].values) +
              list(da_pd[obs_ss].values))
    oc_pd = (list(da_pd[mod_oc].values) +
              list(da_pd[obs_oc].values))

    id_all = (da_pd['ID'].to_list() +
              da_pd['ID'].to_list())

    mod_obs = (len(da_pd['ID'].values) * ['Model'] +
               len(da_pd['ID'].values) * ['Observation'])

    mac = len(mod_obs) * [mac_na]

    new_omf_pd = pd.DataFrame({'Aerosol Concentration (µg m$^{-3}$)': obs_pd,
                               'SS Concentration (µg m$^{-3}$)': ss_pd,
                               'OC Concentration (µg m$^{-3}$)': oc_pd,
                               'Measurements': id_all,
                               '': mod_obs,
                               'Macromolecules': mac})

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

def pd_combine_group(dicc_va, mac_names, mod_key_na, obs_key_na, id_na):
    conc_pd = create_dataframe(dicc_va,
                               mod_key_na,
                               obs_key_na,
                               'conc_mod_ss',
                               'conc_obs_ss',
                               'conc_mod_oc',
                               'conc_obs_oc',
                               mac_names)

    if id_na == 'lipi':
        conc_pd_new = rename_func(conc_pd, 'CVAO', 'CVAO ')
    if id_na == 'prot':
        conc_pd_new = rename_func(conc_pd, 'CVAO', 'CVAO  ')
    if id_na == 'tot':
        conc_pd_to1 = rename_func(conc_pd, 'CVAO', 'CVAO   ')
        conc_pd_to2 = rename_func(conc_pd_to1, 'NAO', 'NAO ')
        conc_pd_new = rename_func(conc_pd_to2, 'WAP', 'WAP ')
    if id_na == 'poly' or id_na == 'oc':
        conc_pd_new = conc_pd


    # plot_functions.box_plot_vert(pd.concat([conc_pd_po, conc_pd_pr_new, conc_pd_li_new, conc_pd_to_new]),
    #                              mac_names, ['pol', 'pro', 'lip'],
    #                              'All_groups', [1e-6, 1e1])

    conc_pd_new.to_pickle(f'pd_files/{id_na}_conc_{global_vars.exp_name}.pkl')
