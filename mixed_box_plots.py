##### Import packges
import numpy as np
import pandas as pd
import glob, os
import seaborn as sns
import matplotlib.pyplot as plt
import SS_plots

def rename_func(data_pd, col, na, new_na):
    pd_new = data_pd
    list_col = data_pd[col].to_list()
    for i in range(len(list_col)):
        if list_col[i] == na:
            list_col[i] = new_na
    pd_new = pd_new.drop(columns=col)
    pd_new[col] = list_col

    return pd_new


def add_text(ax, mol_name, loc1, loc2):
    ax.text(loc1, loc2,
            mol_name,
            fontsize='14',
            weight='bold',
            bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})


def box_plot_vert(dict_df, mol_name, ID, title, lim):
    # Create new plot
    #print(ID, list(dict_df["Measurements"]))
    fig, ax = plt.subplots(figsize=(15, 8))
    # YlGnBu
    states_palette = sns.color_palette("seismic", n_colors=4)
    my_pal = {"Model offline OMF": "skyblue",
              "ECHAM-HAM aerosol concentration": "pink",
              "Observation OMF": "royalblue",
              "Observation aerosol concentration": "palevioletred"}
    # Plot with seaborn
    # df = df.drop(df[df.score < 50].index)
    col_na = "Aerosol OMF"
    # df = df.drop(df[df.score < 50].index)
    dict_df = dict_df.drop(dict_df[dict_df[col_na] < 0].index)
    # print(svd["Aerosol OMF & Concentration (µg m$^{-3}$)"].values)
    bx = sns.boxplot(data=dict_df, x="Measurements",
                     y="Aerosol OMF", hue="",
                     palette=my_pal,
                     width=.7)
    # The box shows the quartiles of the
    # dataset while the whiskers extend to
    # show the rest of the distribution,
    # except for points that are determined
    # to be “outliers” using a method that
    # is a function of the inter-quartile range.

    add_text(ax, mol_name[0], 0.3, 4)
    add_text(ax, mol_name[1], 2.9, 4)
    add_text(ax, mol_name[2], 4.6, 4)
    add_text(ax, mol_name[3], 5.8, 4)

    # plot_text(dict_macrom[0], ax, ID[0], 0.1, 2, 50, 'PCHO')
    # plot_text(dict_macrom[1], ax, ID[1], 5, 7, 50, 'DCAA')
    # # plot_text(dict_macrom[2], ax, ID[2], 10, 9.6, 20, 'PL')

    # Customizing axes
    ax.tick_params(axis='both', labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.set_xlabel('', fontsize=14)
    ax.set_yscale('log')
    ax.grid(linestyle='--', linewidth=0.4)
    ax.set_ylim(lim)

    plt.legend(loc="lower right", fontsize='14')  # bbox_to_anchor=(1.04, 1),

    # create secondary axis
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    ax2.set_ylim(lim)
    ax2.set_ylim(ax.get_ylim())
    # ax2.set_yticklabels(ax.get_yticks())
    ax2.tick_params(axis='both', labelsize='14')
    ax2.set_ylabel('Concentration (µg m$^{-3}$)', fontsize=14)

    # color axis
    ax.spines['left'].set_color('royalblue')
    ax.yaxis.label.set_color('royalblue')
    ax.tick_params(axis='y', colors='royalblue')

    ax2.spines['right'].set_color('palevioletred')
    ax2.yaxis.label.set_color('palevioletred')
    ax2.tick_params(axis='y', colors='palevioletred')

    # dotted lines to separate groups
    ax.axvline(2.5, color=".3", dashes=(2, 2))
    ax.axvline(4.5, color=".3", dashes=(2, 2))
    ax.axvline(5.55, color=".3", dashes=(2, 2))


    plt.savefig(f'plots/mixed_omf_conc_{title}_box.png', dpi=300, bbox_inches="tight")
    plt.show()


# Press the green button in the gutter to run the script.
def compute_aver_std(vals, val_na, i, id_biom):
    print('\n', val_na)
    cases = ['Model', 'Observation']
    for sta in station_names[i]:
        sel_vals = vals[vals['Measurements'] == sta]
        mod_sel_vals = sel_vals[sel_vals[""] == cases[0]][val_na]
        aver = mod_sel_vals.mean()
        std = mod_sel_vals.std()
        print('Model', id_biom, sta, '\n', len(mod_sel_vals), aver, std, '\n', '\n')

        obs_sel_vals = sel_vals[sel_vals[""] == cases[1]][val_na]
        aver = obs_sel_vals.mean()
        std = obs_sel_vals.std()
        print('Observation', id_biom, sta, '\n', len(obs_sel_vals), aver, std, '\n', '\n')

        print("Model mean bias")
        print((np.array(mod_sel_vals.values)[:len(obs_sel_vals.values)] - np.array(obs_sel_vals.values)).mean())


if __name__ == '__main__':
    # data paths
    data_dir = "pd_files/"

    try:
        os.mkdir('plots')
    except OSError:

        pass

    var = ['poly', 'prot', 'lipi', 'tot']  #
    mac_names = ['PCHO$_{aer}$|CCHO$_{aer}$', 'DCAA$_{aer}$|CAA$_{aer}$', 'PL$_{aer}$|PG$_{aer}$',
                 '(PCHO$_{aer}$+DCAA$_{aer}$+PL$_{aer}$)|OM']  #
    station_names = [['NAO', 'CVAO', 'WAP'],  #
                     ['SVD14', 'SVD15', 'SVD18', 'CVAO  ', 'RS'],  #
                     ['CVAO '],
                     ['NAO ', 'CVAO   ', 'WAP ']]  #
    mix_omf_conc = []
    all_conc = []
    for i, v in enumerate(var):
        omf = pd.read_pickle(f'{data_dir}{v}_omf.pkl')

        omf_rename = rename_func(omf, '', 'Model', 'Model offline OMF')
        omf_rename = rename_func(omf_rename, '', 'Observation', 'Observation OMF')
        omf_rename = omf_rename.rename(columns={'Aerosol OMF':
                                                    'Aerosol OMF'})

        conc = pd.read_pickle(f'{data_dir}{v}_conc.pkl')
        conc_rename = rename_func(conc, '', 'Model', 'ECHAM-HAM aerosol concentration')
        conc_rename = rename_func(conc_rename, '', 'Observation', 'Observation aerosol concentration')
        all_conc.append(conc_rename)
        conc_rename = conc_rename.rename(columns={'Aerosol Concentration (µg m$^{-3}$)':
                                                      'Aerosol OMF'})  #'Aerosol Concentration (µg m$^{-3}$)

        mix = pd.concat([omf_rename, conc_rename])
        mix_omf_conc.append(mix)

        compute_aver_std(conc, 'Aerosol Concentration (µg m$^{-3}$)', i, v)
        print('\n', '\n')
        compute_aver_std(omf, 'Aerosol OMF', i, v)



    SS_plots.plot_ss(pd.concat(all_conc[:-1]))

    box_plot_vert(pd.concat([mix_omf_conc[0], mix_omf_conc[1],
                             mix_omf_conc[2], mix_omf_conc[3]]),
                  mac_names, ['pol', 'pro', 'lip'],
                  'All_groups', [1e-6, 1e1])

