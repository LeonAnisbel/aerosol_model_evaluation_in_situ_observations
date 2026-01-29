import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from utils_functions import utils_plots, global_vars


def plot_fit(conc_mod_obs, obs_col_na, mod_col_na, title, xli, yli, loc):
    """
    creates Old SS scatter plot of model vs. observations
    :var conc_mod_obs: dataframe with SS concentrations from observations and model interpolated results
    :var obs_col_na: column name of observations
    :var mod_col_na: column name of model values
    :param title: fig title
    :param xli: x-axis max. limit
    :param yli: y-axis max. limit
    :param loc:
    :return: None
    """
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=conc_mod_obs,
                    x=obs_col_na,
                    y=mod_col_na,
                    hue="Measurements")
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    ax.plot(lims,
            lims,
            'k--',
            alpha=0.75,
            zorder=0,
            label='1:1 line')

    # Fitting poly data
    obs = np.array(conc_mod_obs[obs_col_na].values)
    mod = np.array(conc_mod_obs[mod_col_na].values)

    RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs, mod)
    formatted_pvalues = (f'RMSE= {np.round(RMSE, 2)} \n bias= {np.round(mean_bias, 2)} '
                         f'\n R= {np.round(R, 2)}', f'\n pval= {np.round(pval, 2)}')
    print('SEA SALT STATISTICS', formatted_pvalues)

    ax.tick_params(axis='both',
                   labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.xaxis.get_label().set_fontsize(14)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(xli)
    ax.set_ylim(yli)
    plt.legend(loc='upper left')
    ax.grid(linestyle='--',
            linewidth=0.4)
    fig.tight_layout()

    plt.savefig(f'../plots/SS_conc_{title}_scatter_{global_vars.exp_name}.png',
                dpi = 300)
    plt.close()


def plot_box_ss(ax, conc_mod_obs, mod_pd, obs_pd, with_map=False):
    """
    Creates box plot figure of modelled Sea salt concentration and observations for all available stations
    :var ax: matplotlib axes object
    :var conc_mod_obs: dataframe with SS concentrations from observations and model interpolated results
    :var mod_pd: list with model interpolated SS concentrations
    :var obs_pd: list with observations of SS concentrations
    :param with_map: boolean to determine whether to plot a map in addition to the SS concentration plot
    :return: None
    """

    font = 10
    my_pal = {"Model": "gray",
              "Observations": "black"}
    conc_mod_obs = conc_mod_obs.reset_index(drop=True)

    pl1 = sns.boxplot(data=conc_mod_obs,
                      x="Measurements_renamed",
                      y='SS Concentration (µg m$^{-3}$)',
                      hue="",
                      palette=my_pal,
                      fill=False,
                      gap=0.1)
    pl2 = sns.stripplot(data=conc_mod_obs,
                        x="Measurements_renamed",
                        y='SS Concentration (µg m$^{-3}$)',
                        hue="",
                        dodge=True,
                        palette=my_pal,)
    ax.legend_.remove()
    ax.yaxis.get_label().set_fontsize(font)
    ax.tick_params(axis='both',
                   labelsize=font)

    handles, labels = pl1.get_legend_handles_labels()
    plt.legend(handles=handles[:2],
                       labels= labels[:2],
                       ncol=1,
                       # bbox_to_anchor=(0.18, 0.1),
                       loc='lower left',
                       fontsize=font)

    ax.set_yscale('log')
    ax.set_xlabel('')
    ax.set_ylim([10**-3, 10**1])

    ax.grid(linestyle='--',
            linewidth=0.4)

    list_stat = ['NAO','CVAO', 'WAP', 'SVD']
    stat_pos = [-0.3, 0.85, 1.85, 2.85]

    if not with_map:
        for i, stat in enumerate(list_stat):
            mod = mod_pd[mod_pd['Measurements_renamed']==stat]
            obs = obs_pd[obs_pd['Measurements_renamed']==stat]

            RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs['Aerosol Concentration (µg m$^{-3}$)'].values,
                                                                 mod['Aerosol Concentration (µg m$^{-3}$)'].values)

            print( 'SS statistics \n',   stat, 'RMSE',  RMSE,'MB',  mean_bias, 'NMB',  NMB)
            if i == 0:
                formatted_values = (f'\n NMB= {np.round(NMB, 2)} '
                                     f'\n    n = {np.round(len(obs), 2)} ')
            else:
                formatted_values = (f'\n {np.round(NMB, 2)} ' 
                                     f'\n    {np.round(len(obs), 2)} ')

            ax.text(stat_pos[i],
                    4,
                    formatted_values,
                    fontsize=6)


def plot_ss(ax, conc_ss, with_map=False):
    """
    Prepares the data to plot each station with an additional identifier to indicate to
    which dataset are the SS samples associated. It also computes the statistics and calls
    the function that creates the SS concentration scatter plot
    :var ax: matplotlib axes object
    :var conc_ss: dataframe with SS concentrations from observations and model interpolated results
    :param with_map: boolean to determine whether to plot a map in addition to the SS concentration plot
    :return : None
    """
    new_meas_vals = []
    for i in conc_ss['Measurements'].values:
        i_new = i
        if i[:4] == 'CVAO':
            i_new = 'CVAO'
        new_meas_vals.append(i_new)
    conc_ss['Measurements_renamed'] = new_meas_vals

    new_meas_vals = []
    for i in conc_ss[''].values:
        if i == 'ECHAM-HAM aerosol concentration':
            i_new = 'Model'
        if i == 'Observation aerosol concentration':
            i_new = 'Observations'
        new_meas_vals.append(i_new)
    conc_ss.drop(columns=[''])
    conc_ss[''] = new_meas_vals

    mod_pd = conc_ss[conc_ss[''] == 'Model']
    obs_pd = conc_ss[conc_ss[''] == 'Observations']


    plot_box_ss(ax,
                conc_ss,
                mod_pd, obs_pd,
                with_map=with_map)
