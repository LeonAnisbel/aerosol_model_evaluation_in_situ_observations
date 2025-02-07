##### Import packges
import numpy as np
import pandas as pd
import  os
import seaborn as sns
import matplotlib.pyplot as plt
import global_vars
import utils_plots

def plot_fit(ax, conc_mod_obs, obs_col_na, mod_col_na, title, xli, yli, loc):
    font = 10

    pl = sns.scatterplot(data=conc_mod_obs, x=obs_col_na, y=mod_col_na, hue="Measurements", ax = ax)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    # ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0, label='1:1 line')
    # print(lims)

    # Fitting poly data
    obs = np.array(conc_mod_obs[obs_col_na].values)
    mod = np.array(conc_mod_obs[mod_col_na].values)

    RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs,mod)
    formatted_pvalues = (f'RMSE= {np.round(RMSE, 2)} \n bias= {np.round(mean_bias, 2)} '
                         f'\n R= {np.round(R, 2)} \n NMB= {np.round(NMB, 2)} '
                         f'\n pval= {np.round(pval, 2)}')
    ax.text(loc, 1.2, formatted_pvalues, fontsize=font)

    # Customizing axes
    # ax.set_title(title, fontsize='16')

    ax.tick_params(axis='both', labelsize=font)
    ax.set_ylabel(mod_col_na+' (µg m$^{-3}$)')
    ax.yaxis.get_label().set_fontsize(font)
    ax.xaxis.get_label().set_fontsize(font)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(xli)
    ax.set_ylim(yli)

    ax.legend_.remove()
        # pl.legend(loc='upper center',
        #           bbox_to_anchor=(0.5, 1.5),
        #           ncol=4,
        #           fontsize=font)

    ax.grid(linestyle='--', linewidth=0.4)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1-0.04,
                         box.width, box.height ])

    k = 10
    utils_plots.set_log_ax(ax,
               [10 ** -3, 10 ** 1],
               [10 ** -3, 10 ** 1],
              "-",)
    utils_plots.set_log_ax(ax,
               [10 ** -3 , 10 ** 1],
               [10 ** -3*k, 10 ** 1 * k],
           ":")
    utils_plots.set_log_ax(ax,
              [10 ** -3, 10 ** 1] ,
              [10 ** -3/k, 10 ** 1 / k],
                ":")

    return pl


def plot_box(ax, conc_mod_obs, col_na, yli, title, label, pos, right_pannel = False):
    font = 10
    my_pal = {"Model": "pink",
              "Observation": "palevioletred"}

    flier = utils_plots.get_marker_flier()
    pl = sns.boxplot(data=conc_mod_obs,
                     x='Measurements',
                     y=col_na,
                     hue="",
                     palette=my_pal,
                     ax = ax,
                     flierprops=flier,
                     width=.5)
    utils_plots.add_species_text_name(ax, title + '|OM$_{aer}$', pos, 4, 8)

    ax.tick_params(axis='both', labelsize=font)
    ax.set_ylabel( 'Aerosol concentration (µg m$^{-3}$)')
    ax.set_xlabel( ' ')

    ax.yaxis.get_label().set_fontsize(font)
    ax.xaxis.get_label().set_fontsize(font)
    ax.set_yscale('log')
    ax.set_ylim(yli)
    ax.set_title(label, loc='left', fontsize=font )

    if right_pannel:
        ax.set_yticklabels([])
        ax.set_ylabel(' ')

    ax.legend_.remove()

    ax.grid(linestyle='--', linewidth=0.4)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1-0.04,
                         box.width, box.height ])

    return pl


def plot_oc(conc_ss, ti):
    conc_total = conc_ss.copy()
    old_cols_names = ['OC Aerosol concentration',
                'Aerosol Concentration (µg m$^{-3}$)']
    title_list = ['OC', '(OC+PMOA)']
    list_stat = ['NAO','CVAO', 'WAP']

    fig1, axs1 = plt.subplots(1, 2, figsize=(6, 5))
    axes1 = axs1.flatten()

    fig2, axs2 = plt.subplots(1, 2, figsize=(6, 5))
    axes2 = axs2.flatten()

    bool_panel = [False, True]
    pos_box_title = [0.6, 0.3]
    stat_pos = [-0.3, 0.85, 1.85]
    labels = [r'$\bf{(a)}$', r'$\bf{(b)}$']

    for old_cols, title, ax1, ax2, panel, pos, lab  in zip(old_cols_names, title_list, axes1, axes2, bool_panel, pos_box_title, labels):
        pl1 = plot_box(ax1,
                      conc_total,
                      old_cols,
                      [1e-5, 1e1],
                      title,
                       lab,
                      pos,
                      right_pannel=panel)

        mod_vals, obs_vals, meas_vals = [], [], []
        for i, m in enumerate(list_stat):
            conc_loc = conc_ss[conc_ss['Measurements'] == m]
            mod_oc = np.array(conc_loc[conc_loc['']=='Model'][old_cols].values)
            obs_oc = np.array(conc_loc[conc_loc['']=='Observation'][old_cols].values)

            RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs_oc,mod_oc)
            if i == 0:
                formatted_pvalues = (f'RMSE= {np.round(RMSE, 2)} '
                                     f'\n MB= {np.round(mean_bias, 2)} '
                                     f'\n NMB= {np.round(NMB, 2)} '
                                     f'\n    n= {np.round(len(obs_oc), 2)} ')
            else:
                formatted_pvalues = (f'{np.round(RMSE, 2)} '
                                     f'\n {np.round(mean_bias, 2)} '
                                     f'\n {np.round(NMB, 2)} ' 
                                     f'\n    {np.round(len(obs_oc), 2)} ')
            ax1.text(stat_pos[i], 0.0000105,  formatted_pvalues, fontsize=6)

        conc_pd = pd.DataFrame(data={'Measurements': conc_ss['Measurements'].values[:len(conc_ss[conc_ss['']=='Model'][old_cols].values)],
                                        f'Model {title} concentration': conc_ss[conc_ss['']=='Model'][old_cols].values,
                                        'Observations OM concentration (µg m$^{-3}$)': conc_ss[conc_ss['']=='Observation'][old_cols].values*2})
        obs_col_na = 'Observations OM concentration (µg m$^{-3}$)'
        mod_col_na = f'Model {title} concentration'


        pl2 = plot_fit(ax2,
                 conc_pd,
                 obs_col_na, mod_col_na,
                 title,
                 [1e-5, 1e1],
                 [1e-5, 1e1],
                 0.002)

    utils_plots.customize_legend_sv_fig(pl1, fig1, 'OC_PMOA_scatter_panel_plot')
    utils_plots.customize_legend_sv_fig(pl2, fig2, 'OC_PMOA_box_panel_plot')



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # data paths
    data_dir = "pd_files/"

    try:
        os.mkdir('plots')
    except OSError:

        pass

    mix_omf_conc = []
    all_conc = []
    conc = pd.read_pickle(f'{data_dir}oc_conc_{global_vars.exp_name}.pkl')

    print(conc)

    plot_oc(conc, 'OC')

