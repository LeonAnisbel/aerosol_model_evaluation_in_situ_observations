import pandas as pd
import os
import matplotlib.pyplot as plt
import OC_plots, SS_plots, PMOA_plots
from utils_functions import read_data_functions, utils_plots, global_vars
import cartopy.crs as ccrs
import codecs

if __name__ == '__main__':
    data_dir = "../outputs/"
    obs_dir = '../Observations/'

    with_oc = True   #if False only PMOA should be considered, if True it will compute PMOA+OC

    try:
        os.mkdir('../plots')
    except OSError:

        pass

    var = ['poly', 'prot', 'lipi', 'tot']  #
    mac_names = ['PCHO$_{aer}$|CCHO$_{aer}$', 'DCAA$_{aer}$|CAA$_{aer}$', 'PL$_{aer}$|PG$_{aer}$'
                ,'(PCHO$_{aer}$+DCAA$_{aer}$+PL$_{aer}$)|OM$_{aer}$']  #
    station_names = [['NAO', 'CVAO', 'WAP'],  #
                     ['SVD', 'CVAO  '],  #
                     ['CVAO '],
                     ['NAO ', 'CVAO   ', 'WAP ']]  #
    fig_title = 'All_groups'
    mix_omf_conc = []
    all_conc = []

    for i, v in enumerate(var): # iterate though biomolecules
        # Read OMF data
        omf = pd.read_pickle(f'{data_dir}{v}_omf.pkl')

        if v == 'prot':
            stat_names = []
            for stat in omf["Measurements"]:
                stat_new = stat
                if stat == 'SVD15':
                    stat_new = 'SVD'
                stat_names.append(stat_new)
            omf.drop("Measurements",
                     axis=1,
                     inplace=True)
            omf['Measurements'] = stat_names


        omf_rename = utils_plots.rename_func(omf,
                                             '',
                                             'Model',
                                             'Model offline OMF') # reformat column names
        omf_rename = utils_plots.rename_func(omf_rename,
                                             '',
                                             'Observation',
                                             'Observation OMF') # reformat column names
        omf_rename = omf_rename.rename(columns={'Aerosol OMF':
                                                    'Aerosol OMF'}) # reformat column names

        # Read aerosol concentration data
        conc = pd.read_pickle(f'{data_dir}{v}_conc_{global_vars.exp_name}.pkl')

        # To consider PMOA+OC, "with_oc" must be True, otherwise only PMOA will be considered
        if v == 'tot' and with_oc:
            loc_dir = ''

            # get original observation data
            _, PASCAL, PI_ICE, CVAO, _ = read_data_functions.read_obs_data_loc(obs_dir,
                                                                               loc_dir)
            pd_concat = pd.concat([PASCAL[1]['OC_µg_per_m3'],
                                   CVAO[0]['OC_µg_per_m3'],
                                   PI_ICE[1]['OC_µg_per_m3']]).to_list() # concatenate data from all stations
            pd_concat = [i*2 for i in pd_concat] # consider factor of 2 to convert OC to OM
            obs_oc_conc = pd.DataFrame(data = {'OC observation': pd_concat})
            obs_oc_conc.dropna(inplace=True)

            # redefine names
            mac_names = ['PCHO$_{aer}$|CCHO$_{aer}$',
                         'DCAA$_{aer}$|CAA$_{aer}$',
                         'PL$_{aer}$|PG$_{aer}$',
                         '(PCHO$_{aer}$+DCAA$_{aer}$+PL$_{aer}$+OC)|OM$_{aer}$']
            fig_title = 'With_OC_All_groups'

            # read OC interpolated model results
            conc_oc = pd.read_pickle(f'{data_dir}oc_conc_{global_vars.exp_name}.pkl')
            conc_oc_mod = conc_oc[conc_oc[''] == 'Model']['OC Concentration (µg m$^{-3}$)'].to_list()

            # sum OC and PMOA to get OM from model (in model OM = OC+PMOA)
            conc_oc_pmoa_mod = [i + j for i, j in zip(conc_oc_mod,conc[conc[''] == 'Model']['Aerosol Concentration (µg m$^{-3}$)'].values)]
            meas_name_mod = ['Model' for i in range(len(conc_oc_mod))]
            meas_name_obs = ['Observation' for i in range(len(obs_oc_conc['OC observation'].to_list()))]

            name_sta = [sta  for sta in conc_oc[conc_oc[''] == 'Model']['Measurements']] # create list of

            # concatenate model and observations into a new dataframe with OC and OM
            conc_oc_om = pd.DataFrame(data={'Measurements': name_sta+name_sta,
                                      '':meas_name_mod+ meas_name_obs,
                                      'OC Aerosol concentration':conc_oc_mod +  obs_oc_conc['OC observation'].to_list(),
                                      'Aerosol Concentration (µg m$^{-3}$)': conc_oc_pmoa_mod+obs_oc_conc['OC observation'].to_list()}
                                      )

            # # rename columns (avoid confusion)
            # conc_oc['PMOA Aerosol concentration'] = conc['Aerosol Concentration (µg m$^{-3}$)'].values
            # conc_oc['OC Aerosol concentration'] = conc_oc['Aerosol Concentration (µg m$^{-3}$)'].values
            #
            # new_conc_om = (conc_oc[conc_oc['']=='Model']['Aerosol Concentration (µg m$^{-3}$)'].values +
            #                conc[conc_oc['']=='Model']['Aerosol Concentration (µg m$^{-3}$)'].values)
            #
            # conc.drop(columns = ['Aerosol Concentration (µg m$^{-3}$)'])
            # conc['Aerosol Concentration (µg m$^{-3}$)'] = (list(new_conc_om) +
            #                                                conc[conc[''] == 'Observation'][
            #                                                    'Aerosol Concentration (µg m$^{-3}$)'].to_list())
            #
            # conc_oc.drop(columns = ['Aerosol Concentration (µg m$^{-3}$)'])
            # conc_oc['Aerosol Concentration (µg m$^{-3}$)'] = (list(new_conc_om) +
            #                                                conc_oc[conc_oc[''] == 'Observation'][
            #                                                    'Aerosol Concentration (µg m$^{-3}$)'].to_list())

            OC_plots.plot_oc(conc_oc_om,
                             'OC_PMOA_') # Leon-Marcos et al. 2025 plot
            OC_plots.plot_oc(conc_oc_om,
                             'OC_PMOA_thesis_',
                             no_stat=True) # Thesis plot

############################################################################
        conc_rename = utils_plots.rename_func(conc,
                                              '',
                                              'Model',
                                              'ECHAM-HAM aerosol concentration') # reformat column names
        conc_rename = utils_plots.rename_func(conc_rename,
                                              '',
                                              'Observation',
                                              'Observation aerosol concentration') # reformat column names
        all_conc.append(conc_rename)
        conc_rename = conc_rename.rename(columns={'Aerosol Concentration (µg m$^{-3}$)':
                                                      'Aerosol OMF'})  #'Aerosol Concentration (µg m$^{-3}$)

        # concatenate OMF and aerosol concentration for all stations and species
        mix = pd.concat([omf_rename, conc_rename])
        mix_omf_conc.append(mix)

    # Create box plot fig with Sea salt aerosol concentration for all species and stations (Leon-Marcos et al. 2025)
    fig, ax = plt.subplots(figsize=(5, 4))
    SS_plots.plot_ss(ax,
                     pd.concat(all_conc[:-1]))
    fig.tight_layout()
    plt.savefig(f'../plots/SS_conc_SS submicrom_box_{global_vars.exp_name}.png',
                dpi = 300)
    plt.close()

    # Create box plot with aerosol concentration comparison only (for PhD Defense)
    PMOA_plots.define_df_plot_conc(mix_omf_conc, mac_names, fig_title)

    # Create box plot fig with omf and organic aerosol concentration for all species and stations (Leon-Marcos et al. 2025)
    fig, ax = plt.subplots(figsize=(13,8))#15, 8
    PMOA_plots.box_plot_vert(ax,
                             pd.concat([mix_omf_conc[0],
                                        mix_omf_conc[1],
                                        mix_omf_conc[2]]),
                             mac_names,
                             [1e-7, 1e1])
    plt.savefig(f'../plots/mixed_omf_conc_{fig_title}_box_{global_vars.exp_name}.png',
                dpi=300,
                bbox_inches="tight")
    plt.close()

############################################################################
#### Thesis plots
    # Create box plot fig with Sea salt aerosol concentration for all species and stations with map of station locations
    fig1, ax = plt.subplots(figsize=(5,4))
    SS_plots.plot_ss(ax, pd.concat(all_conc[:-1]),
                     with_map=True)

    file_water = obs_dir+"SEAWATER_data.csv"
    doc = codecs.open(file_water,
                      'r',
                      'UTF-8')  # open for reading with "universal" type set 'rU'
    data_water = pd.read_csv(doc, sep=',')

    WAP = data_water[data_water['Event'] == 'WAP']
    NAO = data_water[data_water['Event'] == 'PASCAL']
    CVAO = data_water[data_water['Event'] == 'CVAO']
    var_list = [NAO, WAP, CVAO]
    loc_list = []
    for loc in var_list:
        l_list = [[loc['Latitude'].values[m], loc['Longitude'].values[m]]
                  for m in range(len(loc['Longitude'].values))]
        loc_list.append(l_list)
    loc_list.append([[78.9175, 11.89417]])
    names = ['NAO', 'WAP', 'CVAO', 'SVD']

    proj = ccrs.PlateCarree()
    ax1 = fig1.add_axes([0.64, 0.55, 0.4, 0.4],
                       projection=proj)
    PMOA_plots.function_plot_two_pannel(ax1,
                                        loc_list,
                                        names,
                                        5)

    plt.savefig(f'../plots/SS_conc_SS submicrom_box_map_{global_vars.exp_name}.png',
                dpi = 300)
    plt.close()

