##### Import packges
import numpy as np
import pandas as pd
import glob, os
import seaborn as sns
import matplotlib.pyplot as plt

import OC_plots
import read_data_functions
import SS_plots
import global_vars
import utils_plots
import cartopy.crs as ccrs
import codecs

def function_plot_two_pannel(ax, location, names):
    ax.coastlines(resolution='110m', color='gray')
    ax.set_extent([-80, 30, -90, 90])

    ind = 0
    for loc, n in zip(location, names):
        if len(loc) == 2:
            pass
        else:
            ind = 1
            for l in loc:
                ax.scatter(x=l[1], y=l[0],
                           color='k',
                           label=n,
                           marker='o',
                           alpha=1,
                           transform=ccrs.PlateCarree())
            if n=='NAO': s1, s2 = 3,13
            else: s1, s2 = 10, 3

            plt.text(l[1] + s1, l[0] - s2, n,
                     horizontalalignment='right',
                     transform=ccrs.Geodetic(),
                     bbox={'facecolor': 'lightgray',
                           'boxstyle': 'round',
                           'pad': 0.2,
                           'alpha': 0.8})

    gl = ax.gridlines(draw_labels=True,
                      x_inline=False,
                      y_inline=False)  # adding grid lines with labels
    gl.top_labels = False
    gl.right_labels = False
    return ax


def box_plot_vert(ax, dict_df, mol_name, ID, title, lim):

    states_palette = sns.color_palette("seismic", n_colors=4)
    my_pal = {"Model offline OMF": "skyblue",
              "ECHAM-HAM aerosol concentration": "pink",
              "Observation OMF": "royalblue",
              "Observation aerosol concentration": "palevioletred"}
    l_extra = list(my_pal.keys())
    my_pal_keys_list = [[l_extra[0], l_extra[2]],
                        [l_extra[1], l_extra[3]]]
    # Plot with seaborn
    # df = df.drop(df[df.score < 50].index)
    col_na = "Aerosol OMF"
    # df = df.drop(df[df.score < 50].index)
    dict_df = dict_df.drop(dict_df[dict_df[col_na] < 0].index)
    # print(svd["Aerosol OMF & Concentration (µg m$^{-3}$)"].values)
    flier = utils_plots.get_marker_flier()
    bx = sns.boxplot(data=dict_df,
                     x="Measurements",
                     y="Aerosol OMF",
                     hue="",
                     palette=my_pal,
                     flierprops=flier,
                     width=.7)

    # The box shows the quartiles of the
    # dataset while the whiskers extend to
    # show the rest of the distribution,
    # except for points that are determined
    # to be “outliers” using a method that
    # is a function of the inter-quartile range.

    utils_plots.add_species_text_name(ax, mol_name[0], 0.3, 4, 14)
    utils_plots.add_species_text_name(ax, mol_name[1], 2.9, 4, 14)
    utils_plots.add_species_text_name(ax, mol_name[2], 4.63, 4, 14)
    # utils_plots.add_species_text_name(ax, mol_name[3], 5.8, 4, 14)


    stat_pos= [-0.1, 0.9, 1.9, 2.9, 3.9,  4.9]
    stat_names = ['NAO', 'CVAO', 'WAP',  'CVAO  ','SVD', 'CVAO ']
    for i, stat in enumerate(stat_names):
        print('\n', '\n')
        print(stat, '\n')
        df_stat = dict_df[dict_df["Measurements"] == stat]
        n = len(df_stat[df_stat[""]==my_pal_keys_list[-1][-1]]["Aerosol OMF"])
        formatted_pvalues = (f'n = {np.round(n, 2)} ')
        ax.text(stat_pos[i], 0.000008, formatted_pvalues, fontsize=8)

        for var in my_pal_keys_list:
            print(stat, var, '\n')
            obs = df_stat[df_stat[""]==var[1]]["Aerosol OMF"].values
            mod = df_stat[df_stat[""]==var[0]]["Aerosol OMF"].values

            RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs, mod)
            formatted_pvalues = (f'RMSE= {np.round(RMSE, 4)}, bias= {np.round(mean_bias, 4)} '
                                 ,f'R= {np.round(R, 4)},  NMB= {np.round(NMB, 4)} '
                                 ,f' pval= {np.round(pval, 2)}')
            print(formatted_pvalues, '\n')


    # Customizing axes
    ax.tick_params(axis='both', labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.set_xlabel('', fontsize=14)
    ax.set_yscale('log')
    ax.grid(linestyle='--', linewidth=0.4)
    ax.set_ylim(lim)

    plt.legend(loc="lower left", fontsize='14')  # bbox_to_anchor=(1.04, 1),

    # create secondary axis
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    ax2.set_ylim(lim)
    ax2.set_ylim(ax.get_ylim())
    # ax2.set_yticklabels(ax.get_yticks())
    ax2.tick_params(axis='both', labelsize='14')
    ax2.set_ylabel('Concentration (µg m$^{-3}$)', fontsize=14)

    # color axis
    # ax2.spines['right'].set_color('palevioletred')
    ax2.yaxis.label.set_color('palevioletred')
    ax2.tick_params(axis='y', colors='palevioletred')

    # ax.spines['left'].set_color('royalblue')
    ax.yaxis.label.set_color('royalblue')
    ax.tick_params(axis='y', colors='royalblue')

    # dotted lines to separate groups
    ax.axvline(2.5, color=".3", dashes=(2, 2))
    ax.axvline(4.5, color=".3", dashes=(2, 2))
    # ax.axvline(5.55, color=".3", dashes=(2, 2))


def box_plot_vert_OM(dict_df, mol_name, ID, title, lim):
    font = 10
    fig, ax = plt.subplots(figsize=(5, 4))

    new_meas_vals = []
    for i in dict_df[''].values:
        if i == 'ECHAM-HAM aerosol concentration':
            i_new = 'Model'
        if i == 'Observation aerosol concentration':
            i_new = 'Observations'
        new_meas_vals.append(i_new)
    dict_df.drop(columns=[''])
    dict_df[''] = new_meas_vals

    my_pal = {"Model": "pink",
              "Observations": "palevioletred"}

    col_na = "Aerosol Concentration (µg m$^{-3}$)"
    dict_df = dict_df.drop(dict_df[dict_df[col_na] < 0].index)
    bx = sns.boxplot(data=dict_df,
                     x="Measurements",
                     y=col_na,
                     hue="",
                     palette=my_pal,
                     width=.5)

    utils_plots.add_species_text_name(ax, mol_name, 0.02, 4, font)

    # Customizing axes
    ax.tick_params(axis='both', labelsize=font)
    ax.yaxis.get_label().set_fontsize(font)
    ax.set_xlabel('', fontsize=font)
    ax.set_yscale('log')
    ax.grid(linestyle='--', linewidth=0.4)
    ax.set_ylim(lim)
    ax.set_ylabel('Aerosol concentration (µg m$^{-3}$)', fontsize=font)

    plt.legend(loc="lower left", fontsize=font)  # bbox_to_anchor=(1.04, 1),

    ax.spines['left'].set_color('palevioletred')
    ax.yaxis.label.set_color('palevioletred')
    ax.tick_params(axis='y', colors='palevioletred')

    plt.tight_layout()
    plt.savefig(f'plots/mixed_conc_{title}_OC_box_{global_vars.exp_name}.png', dpi=300, bbox_inches="tight")
    plt.show()

if __name__ == '__main__':
    # data paths
    data_dir = "pd_files/"
    with_oc = True   # False only PMOA should be considered, True it will compute PMOA+OC

    try:
        os.mkdir('plots')
    except OSError:

        pass

    var = ['poly', 'prot', 'lipi', 'tot']  #
    mac_names = ['PCHO$_{aer}$|CCHO$_{aer}$', 'DCAA$_{aer}$|CAA$_{aer}$', 'PL$_{aer}$|PG$_{aer}$'
                ,'(PCHO$_{aer}$+DCAA$_{aer}$+PL$_{aer}$)|OM$_{aer}$']  #
    station_names = [['NAO', 'CVAO', 'WAP'],  #
                     ['SVD14', 'SVD', 'SVD18', 'CVAO  ', 'RS'],  #
                     ['CVAO '],
                     ['NAO ', 'CVAO   ', 'WAP ']]  #
    fig_title = 'All_groups'
    mix_omf_conc = []
    all_conc = []
    for i, v in enumerate(var):
        omf = pd.read_pickle(f'{data_dir}{v}_omf.pkl')

        if v == 'prot':
            stat_names = []
            for stat in omf["Measurements"]:
                stat_new = stat
                if stat == 'SVD15':
                    stat_new = 'SVD'
                stat_names.append(stat_new)
            omf.drop("Measurements", axis=1, inplace=True)
            omf['Measurements'] = stat_names


        omf_rename = utils_plots.rename_func(omf, '', 'Model', 'Model offline OMF')
        omf_rename = utils_plots.rename_func(omf_rename, '', 'Observation', 'Observation OMF')
        omf_rename = omf_rename.rename(columns={'Aerosol OMF':
                                                    'Aerosol OMF'})

        conc = pd.read_pickle(f'{data_dir}{v}_conc_{global_vars.exp_name}.pkl')

# To consider PMOA+OC, with_oc must be True, otherwise only PMOA should be considered
############################################################################
        if v == 'tot' and with_oc:
            main_dir = './'
            loc_dir = ''

            _, PASCAL, PI_ICE, CVAO, _ = read_data_functions.read_obs_data_loc(main_dir, loc_dir)
            pd_concat = pd.concat([PASCAL[1]['OC_µg_per_m3'], CVAO[0]['OC_µg_per_m3'], PI_ICE[1]['OC_µg_per_m3']]).to_list()
            pd_concat = [i*2 for i in pd_concat]
            obs_oc_conc = pd.DataFrame(data = {'OC observation': pd_concat})
            obs_oc_conc.dropna(inplace=True)

            mac_names = ['PCHO$_{aer}$|CCHO$_{aer}$', 'DCAA$_{aer}$|CAA$_{aer}$', 'PL$_{aer}$|PG$_{aer}$',
                         '(PCHO$_{aer}$+DCAA$_{aer}$+PL$_{aer}$+OC)|OM$_{aer}$']
            fig_title = 'With_OC_All_groups'

            conc_oc = pd.read_pickle(f'{data_dir}oc_conc_{global_vars.exp_name}.pkl')
            conc_oc_mod = conc_oc[conc_oc[''] == 'Model']['OC Concentration (µg m$^{-3}$)'].to_list()
            #conc_oc_mod_25 = [i-i*0.25 for i in conc_oc[conc_oc[''] == 'Model']['OC Concentration (µg m$^{-3}$)'].to_list()]
            #conc_oc_mod = conc_oc_mod_25

            conc_oc_pmoa_mod = [i + j for i, j in zip(conc_oc_mod,conc[conc[''] == 'Model']['Aerosol Concentration (µg m$^{-3}$)'].values)]
            meas_name_mod = ['Model' for i in range(len(conc_oc_mod))]
            meas_name_obs = ['Observation' for i in range(len(obs_oc_conc['OC observation'].to_list()))]

            name_sta = [sta  for sta in conc_oc[conc_oc[''] == 'Model']['Measurements']]

            conc_oc_om = pd.DataFrame(data={'Measurements': name_sta+name_sta,
                                      '':meas_name_mod+ meas_name_obs,
                                      'OC Aerosol concentration':conc_oc_mod +  obs_oc_conc['OC observation'].to_list(),
                                      'Aerosol Concentration (µg m$^{-3}$)': conc_oc_pmoa_mod+obs_oc_conc['OC observation'].to_list()}
                                      )

            conc_oc['PMOA Aerosol concentration'] = conc['Aerosol Concentration (µg m$^{-3}$)'].values
            conc_oc['OC Aerosol concentration'] = conc_oc['Aerosol Concentration (µg m$^{-3}$)'].values

            # print(conc_oc[conc_oc[''] == 'Model']['Aerosol Concentration (µg m$^{-3}$)'],
            #       conc_oc[conc_oc[''] == 'Model']['OC Concentration (µg m$^{-3}$)'], )
            # # OC_plots.plot_oc(conc_oc, 'OC')

            new_conc_om = (conc_oc[conc_oc['']=='Model']['Aerosol Concentration (µg m$^{-3}$)'].values +
                           conc[conc_oc['']=='Model']['Aerosol Concentration (µg m$^{-3}$)'].values)

            conc.drop(columns = ['Aerosol Concentration (µg m$^{-3}$)'])
            conc['Aerosol Concentration (µg m$^{-3}$)'] = (list(new_conc_om) +
                                                           conc[conc[''] == 'Observation'][
                                                               'Aerosol Concentration (µg m$^{-3}$)'].to_list())

            conc_oc.drop(columns = ['Aerosol Concentration (µg m$^{-3}$)'])
            conc_oc['Aerosol Concentration (µg m$^{-3}$)'] = (list(new_conc_om) +
                                                           conc_oc[conc_oc[''] == 'Observation'][
                                                               'Aerosol Concentration (µg m$^{-3}$)'].to_list())
            # print(conc_oc[conc_oc[''] == 'Observation']['OC Aerosol concentration'].values)
            # print(conc_oc[conc_oc[''] == 'Observation']['Measurements'].values)
            OC_plots.plot_oc(conc_oc_om, 'OC+PMOA')
############################################################################

        conc_rename = utils_plots.rename_func(conc, '', 'Model', 'ECHAM-HAM aerosol concentration')
        conc_rename = utils_plots.rename_func(conc_rename, '', 'Observation', 'Observation aerosol concentration')
        all_conc.append(conc_rename)
        conc_rename = conc_rename.rename(columns={'Aerosol Concentration (µg m$^{-3}$)':
                                                      'Aerosol OMF'})  #'Aerosol Concentration (µg m$^{-3}$)

        mix = pd.concat([omf_rename, conc_rename])
        mix_omf_conc.append(mix)

        # compute_aver_std(conc, 'Aerosol Concentration (µg m$^{-3}$)', i, v)
        print('\n', '\n')
        # compute_aver_std(omf, 'Aerosol OMF', i, v)



    SS_plots.plot_ss(pd.concat(all_conc[:-1]))

    fig, ax = plt.subplots(figsize=(12,8))#15, 8
    box_plot_vert(ax, pd.concat([mix_omf_conc[0], mix_omf_conc[1],
                             mix_omf_conc[2]]),
                              mac_names, ['pol', 'pro', 'lip'],
                              fig_title, [1e-7, 1e1])
    plt.savefig(f'plots/mixed_omf_conc_{fig_title}_box_{global_vars.exp_name}.png', dpi=300, bbox_inches="tight")


#################
#### Thesis
    fig, ax = plt.subplots(figsize=(15,8))
    box_plot_vert(ax, pd.concat([mix_omf_conc[0], mix_omf_conc[1],
                             mix_omf_conc[2]]),
                              mac_names, ['pol', 'pro', 'lip'],
                              fig_title, [1e-7, 1e1])
    file_water = "/home/manuel/Downloads/Observ_sites_maps/SEAWATER_data.csv"
    doc = codecs.open(file_water, 'r', 'UTF-8')  # open for reading with "universal" type set 'rU'
    data_water = pd.read_csv(doc, sep=',')

    WAP = data_water[data_water['Event'] == 'WAP']
    NAO = data_water[data_water['Event'] == 'PASCAL']
    CVAO = data_water[data_water['Event'] == 'CVAO']
    var_list = [NAO, WAP, CVAO]
    loc_list = []
    for loc in var_list:
        l_list = [[loc['Latitude'].values[m], loc['Longitude'].values[m]] for m in range(len(loc['Longitude'].values))]
        loc_list.append(l_list)
    loc_list.append([[78.9175, 11.89417]])
    names = ['NAO', 'WAP', 'CVAO', 'SVD']  # ,, 'RS' 'Mace \n Head']

    proj = ccrs.PlateCarree()
    ax1 = fig.add_axes([0.88, 0.33, 0.4, 0.55], projection=proj)
    function_plot_two_pannel(ax1, loc_list, names)



    plt.savefig(f'plots/mixed_omf_conc_map_{fig_title}_box_{global_vars.exp_name}.png', dpi=300, bbox_inches="tight")

    # if with_oc:
    #     box_plot_vert_OM(all_conc[-1],
    #                   mac_names[-1],
    #                   station_names[-1],
    #                   fig_title,
    #                   [1e-3, 1e1])

