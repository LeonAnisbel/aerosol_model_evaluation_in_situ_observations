import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy as cart
import matplotlib.colors as mcolors
import numpy as np
from matplotlib import ticker, cm
from matplotlib.ticker import FuncFormatter
import seaborn as sns
from sklearn.metrics import mean_squared_error
import matplotlib.path as mpath


def function_plot_mult(x, y, z_pol, z_pro, z_lip, ICE, label, title):
    # Define the figure and each axis for the 3 rows and 3 columns
    fig, axs = plt.subplots(nrows=1, ncols=3,
                            subplot_kw={'projection': ccrs.NorthPolarStereo()},  # Robinson
                            figsize=(12, 6))

    clevs = np.array([0, 0.001, 0.005, 0.05, 0.1, 0.5, 1, 2, 3, 5, 7])

    for i in range(len(axs)):
        axs[i].set_extent([-180, 180, 50, 90], ccrs.PlateCarree())

    colors = (
        '#3000e2', '#2548ff', '#27a7ff', '#70e4ff', '#cbfdff',
        '#fdf9d2', '#fec783', '#ff5d46', '#fe0229', '#b40018')
    cmap1 = mcolors.ListedColormap(colors)
    norm1 = mcolors.BoundaryNorm(clevs, cmap1.N)

    ic = axs[0].contourf(x, y, ICE[:, :], levels=4, cmap='Greys_r',
                         transform=ccrs.PlateCarree())
    m = axs[0].pcolormesh(x, y, z_lip[:, :], norm=norm1, cmap=cmap1,
                          transform=ccrs.PlateCarree())
    axs[0].set_title('Lipids', fontsize=18)
    # axs[0].text(.81, 0.05, "(a)", transform=axs[0].transAxes, fontsize=18, fontweight='bold')

    ic = axs[1].contourf(x, y, ICE[:, :], levels=4, cmap='Greys_r',
                         transform=ccrs.PlateCarree())
    m = axs[1].pcolormesh(x, y, z_pro[:, :], norm=norm1, cmap=cmap1,
                          transform=ccrs.PlateCarree())
    axs[1].set_title('Proteins', fontsize=18)
    # axs[1].text(.81, 0.05, "(b)", transform=axs[1].transAxes, fontsize=18, fontweight='bold')

    ic = axs[2].contourf(x, y, ICE[:, :], levels=4, cmap='Greys_r',
                         transform=ccrs.PlateCarree())
    m = axs[2].pcolormesh(x, y, z_pol[:, :], norm=norm1, cmap=cmap1,
                          transform=ccrs.PlateCarree())
    axs[2].set_title('Polysaccharides', fontsize=18)
    # axs[2].text(.81, 0.05, "(c)", transform=axs[2].transAxes, fontsize=18, fontweight='bold')

    # ic = axs[2].contourf(x, y, ICE[:, :],levels=4, cmap='Greens_r',
    #    transform=ccrs.PlateCarree())
    # m = axs[2].pcolormesh(x, y, z_lip[:, :]+z_pro[:, :]+z_pol[:, :], norm = norm1, cmap = cmap1,
    #             transform=ccrs.PlateCarree())
    # axs[2].set_title('MOA', fontsize=18)

    for ax in axs.flat:
        ax.gridlines()
        ax.add_feature(
            cart.feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgray'))
        ax.coastlines(resolution='50m')

        # ic_bar = plt.colorbar(ic, extendfrac='auto',
        #                         shrink=0.6,orientation='horizontal')
        # ic_bar.set_label('Sea ice cover (%)',fontsize = '12')
        # labels = np.arange(20,100,15)
        # ic_bar.set_ticklabels(labels)
        #
        # cbar = plt.colorbar(m, norm=norm1, cmap=cmap1, ticks=clevs,
        #                     boundaries=clevs, extend='max', extendfrac='auto',
        #                     shrink=0.9,orientation='horizontal')
        # cbar.set_label(label,fontsize = '16')

    # Add a colorbar axis at the bottom of the graph
    cax = fig.add_axes([0.3, 0.2, 0.5, 0.04])
    # cax = fig.add_axes([axs[2].get_position().x1 + 0.01, axs[2].get_position().y0, 0.02, axs[2].get_position().height])
    # Draw the colorbar
    cbar = fig.colorbar(m, norm=norm1, cmap=cmap1, ticks=clevs,
                        boundaries=clevs, extend='max', extendfrac='auto',
                        shrink=0.9, cax=cax, orientation='horizontal')  # ,cax=cax)
    cbar.set_label(label, fontsize='12')

    # Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.4, 0.07, 0.3, 0.01])
    #Draw the colorbar
    ic_bar = fig.colorbar(ic, extendfrac='auto', shrink=0.6,
                          cax=cbar_ax, orientation='horizontal')
    ic_bar.set_label('Sea ice cover (%)', fontsize='12')
    labels = np.arange(20, 100, 15)
    ic_bar.set_ticklabels(labels)
    #

    #     plt.suptitle(title)
    plt.savefig(title + '.png', dpi=300, bbox_inches="tight")
    plt.show()


def plot_map_mean(C, name, id_var, step, ma, outdir_polts, factor, unit, cmap, poles=True, total=True):
    if poles:
        fig, ax = plt.subplots(1, 1,  # define figure with cartopy
                               subplot_kw={'projection': ccrs.NorthPolarStereo()}, figsize=(5, 6))
        ax.set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
    else:
        fig, ax = plt.subplots(1, 1,  # define figure with cartopy
                               subplot_kw={'projection': ccrs.Robinson()}, figsize=(5, 6))

    fig.suptitle(name, fontsize=20)

    Z = (C[id_var[0]] + C[id_var[1]] + C[id_var[2]]) / factor
    id_var_a = id_var
    Z_da = np.ma.masked_where(Z <= 0, Z)

    levels = np.linspace(0.0, ma, 20)
    im = ax.contourf(Z.lon, Z.lat, Z_da, levels=levels, extend='max',
                     cmap=cmap, transform=ccrs.PlateCarree())  # , levels=np.arange(-10,10.1, 0.1))
    ax.add_feature(
        cart.feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgray'))
    ax.coastlines()
    ax.gridlines()

    if ma >= 6:
        fmt = lambda x, pos: '{:.0f}'.format(x)
    if ma < 6 and ma >= 1:
        fmt = lambda x, pos: '{:.1f}'.format(x)
    if ma < 0.05:
        fmt = lambda x, pos: '{:.4f}'.format(x)

    if ma >= 0.05 and ma <= 0.9: fmt = lambda x, pos: '{:.2f}'.format(x)

    cbar = fig.colorbar(im, format=FuncFormatter(fmt), orientation="horizontal", label=unit)
    cbar.set_label(unit, fontsize='16')
    cbar.ax.tick_params(rotation=45)

    plt.savefig(outdir_polts + name + '_mean.png', dpi=300, bbox_inches="tight")
    plt.show()


def plot_map(C, name, id_var, step, ma, outdir_polts, factor, unit, cmap, poles=True, total=True):
    if poles:
        fig, axes = plt.subplots(4, 1,  # define figure with cartopy
                                 subplot_kw={'projection': ccrs.NorthPolarStereo()}, figsize=(5, 12))
        for i in range(len(axes)):
            axes[i].set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
    else:
        fig, axes = plt.subplots(4, 1,  # define figure with cartopy
                                 subplot_kw={'projection': ccrs.Robinson()}, figsize=(5, 12))

    axflat = axes.flatten()

    plt.subplots_adjust(left=0, right=1, top=0.9, bottom=0)

    fig.suptitle(name, fontsize=20)
    months = [1, 4, 8, 10]
    titles = ['February', 'May', 'September', 'November']

    for idx, var in enumerate(months):
        if total:
            Z = (C[var][id_var[0]] + C[var][id_var[1]] + C[var][id_var[2]]) / factor
            id_var_a = id_var[0][:3]
        else:
            Z = C[var][id_var] / factor
            id_var_a = id_var

        Z_da = np.ma.masked_where(Z <= 0, Z)
        # Z_ma = ma
        # Z_mi = mi

        # np.logspace(np.log10(Z_mi), np.log10(Z_ma), 30),locator=ticker.LogLocator()
        levels = np.linspace(0.0, ma, 20)
        im = axflat[idx].contourf(Z.lon, Z.lat, Z_da, levels=levels, extend='max',
                                  cmap=cmap, transform=ccrs.PlateCarree())  # , levels=np.arange(-10,10.1, 0.1))
        axflat[idx].set_title(titles[idx], fontsize='14')
        axflat[idx].add_feature(
            cart.feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgray'))
        axflat[idx].coastlines()
        axflat[idx].gridlines()

    cbar_ax = fig.add_axes([0.1, -0.05, 0.8, 0.03])

    if ma >= 6:
        fmt = lambda x, pos: '{:.0f}'.format(x)
    if ma < 6 and ma >= 1:
        fmt = lambda x, pos: '{:.1f}'.format(x)
    if ma < 0.05:
        fmt = lambda x, pos: '{:.4f}'.format(x)

    if ma >= 0.05 and ma <= 0.9: fmt = lambda x, pos: '{:.2f}'.format(x)

    cbar = fig.colorbar(im, format=FuncFormatter(fmt), cax=cbar_ax, orientation="horizontal", label=unit)
    # step = ma/10
    # clevs = [x for x in np.arange(0, ma + step, step)]
    cbar.set_label(unit, fontsize='16')
    cbar.ax.tick_params(rotation=45)
    #cbar.locator = ticker.LogLocator()
    # cbar.set_ticks(cbar.locator.tick_values(Z_mi, Z_ma))
    # cbar.minorticks_off()

    plt.savefig(outdir_polts + name + '.png', dpi=300, bbox_inches="tight")
    plt.show()


def box_plot(data_pd, plot_dir, yaxis, title, name):
    # box plot using seaborn
    fig, ax = plt.subplots(1, 1, figsize=(9, 8))  # creating figure
    # The box shows the quartiles of the dataset while the whiskers extend to show the rest of the distribution,
    # except for points that are determined to be “outliers” using a method that is a function of the inter-quartile
    # range.
    sns.boxplot(data=data_pd, x="Measurement", y=yaxis, hue="")

    # Customizing axes
    ax.tick_params(axis='both', labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.set_xlabel('')
    ax.set_yscale('log')
    ax.set_title(title, fontsize='16')
    ax.grid(linestyle='--', linewidth=0.4)
    ax.set_ylim([1e-5, 1e1])

    plt.savefig(f'{plot_dir}{name}.png', dpi=300)
    plt.show()


def box_plot_lip(new_pd_pl, plot_dir):
    # box plot using seaborn
    import matplotlib.ticker as ticker

    fig, ax = plt.subplots(figsize=(9, 8))
    bx = sns.boxplot(data=new_pd_pl, x="Species", y="Aerosol Concentration", hue="",
                     showfliers=False)  # The box shows the quartiles of the
    # dataset while the whiskers extend to
    # show the rest of the distribution,
    # except for points that are determined
    # to be “outliers” using a method that
    # is a function of the inter-quartile range.
    #
    # bx.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    # bx.xaxis.set_major_formatter(ticker.ScalarFormatter())

    # Customizing axes
    ax.tick_params(axis='both', labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.set_xlabel('', fontsize=14)
    ax.set_yscale('log')
    ax.set_ylim([1e-5, 1e1])
    ax.grid(linestyle='--', linewidth=0.4)
    ax.set_title('Aerosol concentration (Proteins and Lipids) over CVAO', fontsize='16')
    plt.savefig(plot_dir + 'Pro_Lip_omf_box.png', dpi=300)


# Statistical indicators
def equat_stat(model, observ):
    # Root Mean Squared Error
    rmse = np.sqrt(mean_squared_error(model, observ))

    # Correlation coefficient (R)
    corrcoef = np.corrcoef(model, observ)[0][1]

    # Normalised Mean Bias (NMB)
    obs_mean = sum(observ) / len(observ)
    mod_mean = sum(model) / len(model)
    bias = mod_mean - obs_mean
    nmb = bias / obs_mean

    # Normalised Mean Standard Deviation (NMSD)
    std_obs = np.std(observ)
    std_mod = np.std(model)
    nmsd = (std_mod - std_obs) / std_obs

    return rmse, corrcoef, nmb, nmsd


def stat(da, id_mod, id_obs):
    rmse_val = []
    corrcoef = []
    nmb_val = []
    nmsd_val = []

    mod_all = []
    obs_all = []
    for na, _ in da.items():
        if len(da[na][id_mod]) > 0:  # this won't be necessary when the data is complete

            model = da[na][id_mod]
            observ = da[na][id_obs]

            for i, m in enumerate(model):  # For global statistics
                mod_all.append(m)  # appending values to a list of modeled data
                obs_all.append(observ[i])  # appending values to a list of observed data

            rmse, corr, nmb, nmsd = equat_stat(model, observ)
            #             da[na]['RMSE'].append(rmse)

            rmse_val.append(rmse)
            corrcoef.append(corr)
            nmb_val.append(nmb)
            nmsd_val.append(nmsd)

    rmse_tot, corr_tot, nmb_tot, nmsd_tot = equat_stat(mod_all, obs_all)

    stat_all = [rmse_tot, corr_tot, nmb_tot, nmsd_tot]

    return rmse_val, corrcoef, nmb_val, nmsd_val, stat_all


def box_plot_stat(dict_df, dict_macrom, ID, title, lim, yaxis, plot_dir, loc_high, loc_left):
    # Create new plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot with seaborn
    bx = sns.boxplot(data=dict_df, x="Measurement", y=yaxis, hue="")

    rmse, corr, nmb, nmsd, stat_all = stat(dict_macrom, f'{ID}_mod', f'{ID}_obs')
    #     print('RMSE', 'R','NMB','nmsd','\n',stat_all)

    if ID == 'pro' or ID == 'lip':
        # Correlation is a measure of how well two vectors track with each other as they change.
        # You can't track mutual change when one vector doesn't change
        # (https://stackoverflow.com/questions/45897003/python-numpy-corrcoef-runtimewarning-invalid-value-encountered-in-true-divide).
        formatted_pvalues = f' NMB= {np.round(stat_all[2], 2)}   NMSD= {np.round(stat_all[3], 2)}'
    else:
        formatted_pvalues = f'RMSE = {np.round(stat_all[0], 3)}    R= {np.round(stat_all[1], 2)}   NMB= {np.round(stat_all[2], 2)}   NMSD= {np.round(stat_all[3], 2)}'

    ax.text(loc_left, loc_high, formatted_pvalues, fontsize='14', weight='bold',
            bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

    # Customizing axes
    ax.tick_params(axis='both', labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.set_xlabel('', fontsize=14)
    ax.set_yscale('log')
    ax.grid(linestyle='--', linewidth=0.4)
    ax.set_ylim(lim)

    plt.legend(loc='lower left')

    # plt.legend(loc="lower right", fontsize='14')  # bbox_to_anchor=(1.04, 1),

    plt.savefig(plot_dir + f'{title}_{ID}_box.png', dpi=300)


def box_plot_vert(dict_df, mol_name, ID, title, lim):
    # Create new plot
    fig, ax = plt.subplots(figsize=(15, 8))
    # YlGnBu
    states_palette = sns.color_palette("BuPu", n_colors=2)
    # Plot with seaborn
    bx = sns.boxplot(data=dict_df, x="Measurements",
                     y="Aerosol Concentration (µg m$^{-3}$)", hue="", palette=states_palette,
                     width=.7)
    # The box shows the quartiles of the
    # dataset while the whiskers extend to
    # show the rest of the distribution,
    # except for points that are determined
    # to be “outliers” using a method that
    # is a function of the inter-quartile range.

    ax.text(0.1, 4, mol_name[0], fontsize='14', weight='bold', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    ax.text(3.2, 4, mol_name[1], fontsize='14', weight='bold', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    ax.text(5.7, 4, mol_name[2], fontsize='14', weight='bold', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})
    ax.text(7, 4, mol_name[3], fontsize='14', weight='bold', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

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

    # dotted lines to separate groups
    ax.axvline(2.5, color=".3", dashes=(2, 2))
    ax.axvline(6.5, color=".3", dashes=(2, 2))
    ax.axvline(7.5, color=".3", dashes=(2, 2))

    plt.legend(loc="lower left", fontsize='14')  # bbox_to_anchor=(1.04, 1),

    plt.savefig(f'plots/all_final_conc_{title}_box.png', dpi=300, bbox_inches="tight")
