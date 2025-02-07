##### Import packges
import numpy as np
import pandas as pd
import glob, os
import seaborn as sns
import matplotlib.pyplot as plt
import global_vars
import utils_plots


def plot_fit(conc_mod_obs, obs_col_na, mod_col_na, title, xli, yli, loc):
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=conc_mod_obs, x=obs_col_na, y=mod_col_na, hue="Measurements")
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0, label='1:1 line')
    # print(lims)

    # Fitting poly data
    obs = np.array(conc_mod_obs[obs_col_na].values)
    mod = np.array(conc_mod_obs[mod_col_na].values)

    RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs,mod)
    formatted_pvalues = (f'RMSE= {np.round(RMSE, 2)} \n bias= {np.round(mean_bias, 2)} '
                         f'\n R= {np.round(R, 2)}', f'\n pval= {np.round(pval, 2)}')
    print('SEA SALT STATISTICS', formatted_pvalues)
    #    ax.text(loc, 0.0015, formatted_pvalues, fontsize='14')

    # Customizing axes
    # ax.set_title(title, fontsize='16')

    ax.tick_params(axis='both', labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.xaxis.get_label().set_fontsize(14)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(xli)
    ax.set_ylim(yli)
    plt.legend(loc='upper left')
    ax.grid(linestyle='--', linewidth=0.4)
    fig.tight_layout()

    plt.savefig(f'plots/SS_conc_{title}_scatter_{global_vars.exp_name}.png',dpi = 300)
    plt.close()


def plot_violin(conc_mod_obs, mod_pd, obs_pd, obs_col_na, mod_col_na, title, xli, yli, loc):
    font = 10
    fig, ax = plt.subplots(figsize=(5, 4))
    my_pal = {"Model": "gray",
              "Observations": "black"}
    conc_mod_obs = conc_mod_obs.reset_index(drop=True)

    pl1 = sns.boxplot(data=conc_mod_obs,
                      x="Measurements_renamed",
                      y='SS Concentration (µg m$^{-3}$)',
                      hue="",
                      #saturation=0.5, TRY TO RUN it locally to display empty boxes with the edges in the pallete colors
                      palette=my_pal,fill=False, gap=0.1)
    pl2 = sns.stripplot(data=conc_mod_obs,
                        x="Measurements_renamed",
                        y='SS Concentration (µg m$^{-3}$)',
                        hue="",
                        dodge=True,
                        palette=my_pal,
                        ax = ax)
    ax.legend_.remove()
    ax.yaxis.get_label().set_fontsize(font)
    ax.tick_params(axis='both', labelsize=font)

    handles, labels = pl1.get_legend_handles_labels()
    fig.legend(handles=handles[:2],
                       labels= labels[:2],
                       ncol=1,
                       bbox_to_anchor=(0.18, 0.1),
                       loc='lower left',
                       fontsize=font)

    ax.set_yscale('log')
    ax.set_xlabel('')
    ax.set_ylim([10**-3, 10**1])

    ax.grid(linestyle='--', linewidth=0.4)

    list_stat = ['NAO','CVAO', 'WAP', 'SVD']
    stat_pos = [-0.3, 0.85, 1.85, 2.85]

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
        ax.text(stat_pos[i], 4, formatted_values, fontsize=6)

    fig.tight_layout()
    plt.savefig(f'plots/SS_conc_{title}_box_{global_vars.exp_name}.png',dpi = 300)
    plt.close()

def transform_pd(conc, old_cols, obs_col, mod_col):
    conc_model = conc.loc[(conc[''] == 'Model')][:-1]
    conc_obs = conc.loc[(conc[''] == 'Measurements')][:-1]
    conc_model_rename = conc_model.rename(columns={old_cols:
                                                       mod_col})
    conc_obs_rename = conc_obs.rename(columns={old_cols:
                                                   obs_col})

    conc_mod_obs = pd.DataFrame()
    conc_mod_obs[mod_col] = conc_model_rename[mod_col].values
    conc_mod_obs[obs_col] = conc_obs_rename[obs_col].values
    conc_mod_obs["Measurements"] = conc_model_rename["Measurements"]
    return conc_mod_obs


def plot_ss(conc_ss):
    mac_names = 'SS submicrom'
    old_cols = 'SS Concentration (µg m$^{-3}$)'
    obs_col_na = 'Observations SS concentration (µg m$^{-3}$)'
    mod_col_na = 'Model SS concentration (µg m$^{-3}$)'

    new_meas_vals = []
    for i in conc_ss['Measurements'].values:
        # len(i) > 3 and i[:3] == 'CVAO        else:
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


    plot_violin(conc_ss,
                mod_pd, obs_pd,
             obs_col_na, mod_col_na,
             mac_names,
             [1e-3, 1e1],
             [1e-3, 1e1],
             0.7)


    conc_ss = conc_ss[(conc_ss['Measurements'] == 'CVAO ') |
                      (conc_ss['Measurements'] == 'CVAO  ') |
                      (conc_ss['Measurements'] == 'CVAO') |
                      (conc_ss['Measurements'] == 'WAP') |
                      (conc_ss['Measurements'] == 'NAO') |
                      #(conc_ss['Measurements'] == 'SVD14') |
                      (conc_ss['Measurements'] == 'SVD')]
                      #(conc_ss['Measurements'] == 'SVD18')


    #print(conc_ss[conc_ss['Measurements']=='CVAO  '])


    new_meas_vals = []
    for i in obs_pd['Measurements'].values:
        # len(i) > 3 and i[:3] == 'CVAO        else:
        i_new = i
        if i == 'CVAO ':
            i_new = 'CVAO_PG$_{aer}$'
        if i == 'CVAO':
            i_new = 'CVAO_CCHO$_{aer}$'
        if i == 'CVAO  ':
            i_new = 'CVAO_AA$_{aer}$'
        if i == 'SVD':
            i_new = 'SVD'
        new_meas_vals.append(i_new)

    print(len(new_meas_vals), len(obs_pd[old_cols].values), len(mod_pd[old_cols].values))
    conc_ss_pd = pd.DataFrame(data={'Measurements': new_meas_vals,
                                    'Observations SS concentration (µg m$^{-3}$)': obs_pd[old_cols].values,
                                    'Model SS concentration (µg m$^{-3}$)': mod_pd[old_cols].values})

    new_conc_ss = conc_ss_pd.dropna(subset=[obs_col_na])
    new_conc_ss = new_conc_ss.dropna(subset=[mod_col_na])

    #print(new_conc_ss[new_conc_ss['Measurements'] == 'SVD14'][mod_col_na].min(),
     #     new_conc_ss[new_conc_ss['Measurements'] == 'SVD14'][mod_col_na].max())
    # new_meas_vals.append('NO_STAT')
    # print(new_conc_ss)
    for m in range(len(new_meas_vals)):
        mod_ss = new_conc_ss[new_conc_ss['Measurements'] == new_meas_vals[m]][mod_col_na]
        print('SS', new_meas_vals[m])

        # print(mod_ss)

        obs_ss = new_conc_ss[new_conc_ss['Measurements'] == new_meas_vals[m]][obs_col_na]
        RMSE, mean_bias, NMB, R, pval = utils_plots.get_stat(obs_ss, mod_ss)

        # # if new_meas_vals[m] != new_meas_vals[m + 1]:
        #print('SS', new_meas_vals[m])
        print('MOdel ', mod_ss.median(), mod_ss.std())
        print('Observations ', obs_ss.median(), obs_ss.std())
        #
        print('Statistics: ', 'RMSE:', RMSE, 'bias: ', mean_bias, 'NMB: ', NMB, 'R:', R, '\n')
    # conc_mod_obs_ss_acc = transform_pd(new_conc_ss, old_cols, obs_col_na, mod_col_na)
    plot_fit(new_conc_ss,
             obs_col_na, mod_col_na,
             mac_names,
             [1e-3, 1e1],
             [1e-3, 1e1],
             0.7)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # data paths
    data_dir = "pd_files/"

    try:
        os.mkdir('plots')
    except OSError:

        pass

    var = ['ss_acc', 'ss_tot']
    mac_names = ['SS submicrom', 'SS supermicron']  #
    old_cols = ['SS Concentration (µg m$^{-3}$)', 'Total SS Concentration (µg m$^{-3}$)']  #
    obs_col_na = 'Observations SS concentration (µg m$^{-3}$)'
    mod_col_na = 'Model SS concentration (µg m$^{-3}$)'

    mix_omf_conc = []
    # for i,v in enumerate(var):

    conc_ss = pd.read_pickle(f'{data_dir}{var[0]}_conc.pkl')
    conc_mod_obs_ss_acc = transform_pd(conc_ss, old_cols[0], obs_col_na, mod_col_na)
    conc_mod_obs_ss_tot = transform_pd(conc_ss, old_cols[1], obs_col_na, mod_col_na)
    conc_mod_obs_ss_tot = conc_mod_obs_ss_tot.dropna(subset=obs_col_na)
    plot_fit(conc_mod_obs_ss_acc,
             obs_col_na, mod_col_na,
             mac_names[0],
             [1e-3, 1e1],
             [1e-3, 1e1],
             0.5)
    plot_fit(conc_mod_obs_ss_tot,
             obs_col_na, mod_col_na,
             mac_names[1],
             [1e-3, 1e2],
             [1e-3, 1e2],
             4.)

    # print(conc_mod_obs)
