##### Import packges
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import utils_plots

def plot_text_all_stat(ax, conc_mod_obs, obs_col_na, mod_col_na, font, loc):
    """
    Adds to plot the text with statistics
    :var ax: matplotlib axes
    :var conc_mod_obs: dataframe with observation and model values
    :param obs_col_na: observation column name
    :param mod_col_na: model column name
    :param font: font size
    :param loc: position to add text with species names
    :return : matplotlib object
    """
    obs = np.array(conc_mod_obs[obs_col_na].values)
    mod = np.array(conc_mod_obs[mod_col_na].values)

    RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs, mod)
    formatted_pvalues = (f'Pearson corr. = {np.round(R, 2)} \n'
                         f'RMSE= {np.round(RMSE, 2)} \n'
                         f'MB= {np.round(mean_bias, 2)} \n'
                         f'NMB= {np.round(NMB, 2)} ')
    ax.text(loc[0], loc[1], formatted_pvalues, fontsize=font)


def plot_fit(ax, conc_mod_obs, obs_col_na, mod_col_na, title, xli, yli, loc):
    """
    Creates a two-panel scatter plot of modelled OC (SPMOAoff) and OC+PMOA (SPMOAon) .vs. observations of OM
    :var ax: matplotlib axes
    :var conc_mod_obs: dataframe with observation and model values
    :param obs_col_na: observation column name
    :param mod_col_na: model column name
    :param title: subplot title
    :param xli: lower and upper x-axis limits
    :param yli: lower and upper y-axis limits
    :param loc: position to add text with species names
    :return : matplotlib object
    """
    font = 10

    pl = sns.scatterplot(data=conc_mod_obs, x=obs_col_na, y=mod_col_na, hue="Measurements", ax = ax)
    plot_text_all_stat(ax, conc_mod_obs, obs_col_na, mod_col_na, font, loc)

    ax.tick_params(axis='both', labelsize=font)
    ax.set_ylabel(mod_col_na+' (µg m$^{-3}$)')
    ax.yaxis.get_label().set_fontsize(font)
    ax.xaxis.get_label().set_fontsize(font)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(xli)
    ax.set_ylim(yli)

    ax.legend_.remove()

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


def plot_box(ax, conc_mod_obs, col_na, yli, title, label, pos, font, right_pannel = False):
    """
    Creates a two-panel box plot of modelled OC (SPMOAoff) and OC+PMOA (SPMOAon) .vs. observations of OM
    :var ax: matplotlib axes
    :var conc_mod_obs: dataframe with observation and model values
    :param col_na: column name
    :param yli: y-axis limits
    :param title: subplot title
    :param label: y-axis label
    :param pos: position to add text with species names
    :param font: font size
    :param right_pannel: do not plot y ticks or labels if right pannel
    :return : matplotlib object
    """
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
    utils_plots.add_species_text_name(ax,
                                      title + '|OM$_{aer}$',
                                      pos,
                                      4,
                                      font-2)

    ax.tick_params(axis='both',
                   labelsize=font)
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


def plot_oc(conc, ti, no_stat=False):
    """
    Prepares the data to plot each station for both experiments of modelled OC concentration (SPMOAoff) and
    OC+PMOA (SPMOAon) .vs. observations of OM . It also computes the statistics and calls the function that creates
    the figures
    :var conc: dataframe containing observations of OM and model interpolated results
    :var ti: use in figure title
    :param no_stat: if True, adds statistics to plot
    :return : None
     """

    f = 8
    conc_total = conc.copy()
    old_cols_names = ['OC Aerosol concentration',
                'Aerosol Concentration (µg m$^{-3}$)']
    title_list = ['OC', '(OC+PMOA)']
    list_stat = ['NAO','CVAO', 'WAP']

    if no_stat: w = 4
    else: w = 5
    fig1, axs1 = plt.subplots(1,
                              2,
                              figsize=(6, w))
    axes1 = axs1.flatten()

    bool_panel = [False, True]
    pos_box_title = [0.6, 0.3]
    stat_pos = [-0.3, 0.85, 1.85]
    labels = [r'$\bf{(a)}$', r'$\bf{(b)}$']

    for old_cols, title, ax1, panel, pos, lab  in zip(old_cols_names,
                                                      title_list,
                                                      axes1,
                                                      bool_panel,
                                                      pos_box_title,
                                                      labels):
        pl1 = plot_box(ax1,
                      conc_total,
                      old_cols,
                      [1e-5, 1e1],
                      title,
                       lab,
                      pos,
                       f+2,
                      right_pannel=panel)

        if not no_stat:
            for i, m in enumerate(list_stat):
                conc_loc = conc[conc['Measurements'] == m]
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
                ax1.text(stat_pos[i], 0.0000105,  formatted_pvalues, fontsize=f-2)

        conc_pd = pd.DataFrame(data={'Measurements': conc['Measurements'].values[:len(conc[conc['']=='Model'][old_cols].values)],
                                        f'Model {title} concentration': conc[conc['']=='Model'][old_cols].values,
                                        'Observations OM concentration (µg m$^{-3}$)': conc[conc['']=='Observation'][old_cols].values*2})
        obs_col_na = 'Observations OM concentration (µg m$^{-3}$)'
        mod_col_na = f'Model {title} concentration'

        if no_stat:
            plot_text_all_stat(ax1,
                               conc_pd,
                               obs_col_na,
                               mod_col_na,
                               f,
                               [0.002, 0.00002])

    utils_plots.customize_legend_sv_fig(pl1,
                                        fig1,
                                        ti+'box_panel_plot')


    fig2, axs2 = plt.subplots(1,
                              2,
                              figsize=(6, 5))
    axes2 = axs2.flatten()
    for old_cols, title, ax2, in zip(old_cols_names, title_list, axes2):

        conc_pd = pd.DataFrame(data={'Measurements': conc['Measurements'].values[:len(conc[conc['']=='Model'][old_cols].values)],
                                        f'Model {title} concentration': conc[conc['']=='Model'][old_cols].values,
                                        'Observations OM concentration (µg m$^{-3}$)': conc[conc['']=='Observation'][old_cols].values*2})
        obs_col_na = 'Observations OM concentration (µg m$^{-3}$)'
        mod_col_na = f'Model {title} concentration'
        pl2 = plot_fit(ax2,
                 conc_pd,
                 obs_col_na, mod_col_na,
                 title,
                 [1e-5, 1e1],
                 [1e-5, 1e1],
                 [0.0002, 1.2])
    utils_plots.customize_legend_sv_fig(pl2,
                                        fig2,
                                        ti+'scatter_panel_plot')

