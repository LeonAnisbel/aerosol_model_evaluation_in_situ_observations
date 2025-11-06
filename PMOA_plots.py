import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import utils_plots
import cartopy.crs as ccrs

def function_plot_two_pannel(ax, location, names, ff):
    """Map of station location """
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
                           'alpha': 0.8},
                     fontsize=ff)
    #
    # gl = ax.gridlines(draw_labels=True,
    #                   x_inline=False,
    #                   y_inline=False)  # adding grid lines with labels
    # gl.top_labels = False
    # gl.right_labels = False
    return ax


def box_plot_vert(ax, dict_df, mol_name, ID, title, lim):
    """Plotting model values and ovservations of OMF and aerosol mass concentration"""
    f = 18
    states_palette = sns.color_palette("seismic", n_colors=4)
    my_pal = {"Model offline OMF": "skyblue",
              "ECHAM-HAM aerosol concentration": "pink",
              "Observation OMF": "royalblue",
              "Observation aerosol concentration": "palevioletred"}
    l_extra = list(my_pal.keys())
    my_pal_keys_list = [[l_extra[0], l_extra[2]],
                        [l_extra[1], l_extra[3]]]
    # Plot with seaborn
    col_na = "Aerosol OMF"
    dict_df = dict_df.drop(dict_df[dict_df[col_na] < 0].index)
    flier = utils_plots.get_marker_flier()
    bx = sns.boxplot(data=dict_df,
                     x="Measurements",
                     y="Aerosol OMF",
                     hue="",
                     palette=my_pal,
                     flierprops=flier,
                     width=.7)

    # The box shows the quartiles of the dataset while the whiskers extend to show the rest of the distribution, except
    # for points that are determined to be “outliers” using a method that is a function of the inter-quartile range.

    utils_plots.add_species_text_name(ax, mol_name[0], 0.3, 4, f-2)
    utils_plots.add_species_text_name(ax, mol_name[1], 2.9, 4, f-2)
    utils_plots.add_species_text_name(ax, mol_name[2], 4.63, 4, f-2)


    stat_pos= [-0.1, 0.9, 1.9, 2.9, 3.85,  4.9]
    stat_names = ['NAO', 'CVAO', 'WAP',  'CVAO  ','SVD', 'CVAO ']
    for i, stat in enumerate(stat_names):
        print('\n', '\n')
        print(stat, '\n')
        df_stat = dict_df[dict_df["Measurements"] == stat]
        n = len(df_stat[df_stat[""]==my_pal_keys_list[-1][-1]]["Aerosol OMF"])
        formatted_pvalues = (f'n = {np.round(n, 2)} ')
        ax.text(stat_pos[i], 0.000008, formatted_pvalues, fontsize=f-2)

        for var in my_pal_keys_list:
            print(stat, var, '\n')
            obs = df_stat[df_stat[""]==var[1]]["Aerosol OMF"].values
            mod = df_stat[df_stat[""]==var[0]]["Aerosol OMF"].values

            RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs, mod)
            formatted_pvalues = (f'RMSE= {np.round(RMSE, 4)}, bias= {np.round(mean_bias, 4)} '
                                 ,f'R= {np.round(R, 4)},  NMB= {np.round(NMB, 4)} '
                                 ,f' pval= {np.round(pval, 2)}')


    # Customizing axes
    ax.tick_params(axis='both', labelsize=f'{f}')
    ax.yaxis.get_label().set_fontsize(f)
    ax.set_xlabel('', fontsize=f)
    ax.set_yscale('log')
    ax.grid(linestyle='--', linewidth=0.4)
    ax.set_ylim(lim)

    plt.legend(loc="lower left", fontsize=f'{f-3}')

    # create secondary axis
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    ax2.set_ylim(lim)
    ax2.set_ylim(ax.get_ylim())
    ax2.tick_params(axis='both', labelsize=f'{f}')
    ax2.set_ylabel('Concentration (µg m$^{-3}$)', fontsize=f)

    # color axis
    ax2.yaxis.label.set_color('palevioletred')
    ax2.tick_params(axis='y', colors='palevioletred')

    ax.yaxis.label.set_color('royalblue')
    ax.tick_params(axis='y', colors='royalblue')

    # dotted lines to separate groups
    ax.axvline(2.5, color=".3", dashes=(2, 2))
    ax.axvline(4.5, color=".3", dashes=(2, 2))
