##### Import packges
import numpy as np
import pandas as pd
import glob,os
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from sklearn.metrics import mean_squared_error
import math

def rename_func(data_pd, col, na, new_na):
    pd_new = data_pd
    list_col = data_pd[col].to_list()
    for i in range(len(list_col)):
        if list_col[i] == na:
            list_col[i] = new_na
    pd_new = pd_new.drop(columns=col)
    pd_new[col] = list_col

    return pd_new


def plot_fit(conc_mod_obs, obs_col_na,mod_col_na,title,xli,yli,loc):
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=conc_mod_obs, x=obs_col_na, y=mod_col_na, hue="Measurement")
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0, label='1:1 line')
    # print(lims)

    # Fitting poly data
    obs = np.array(conc_mod_obs[obs_col_na].values)
    mod = np.array(conc_mod_obs[mod_col_na].values)

    # stadistics
    MSE = mean_squared_error(obs, mod)
    RMSE = math.sqrt(MSE)
    mean_bias = (mod - obs).mean()
    R = np.corrcoef(obs, mod)
    print(RMSE,mean_bias, R[0][1])
    formatted_pvalues = (f'RMSE= {np.round(RMSE,2)} \n bias= {np.round(mean_bias,2)} '
                         f'\n R= {np.round(R[0][1],2)}')
    ax.text(loc,0.0015,formatted_pvalues,fontsize = '14')


    # Customizing axes
    # ax.set_title(title, fontsize='16')

    ax.tick_params(axis='both', labelsize='14')
    ax.yaxis.get_label().set_fontsize(14)
    ax.xaxis.get_label().set_fontsize(14)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(xli)
    ax.set_ylim(yli)
    plt.legend()
    ax.grid(linestyle='--', linewidth=0.4)
    fig.tight_layout()

    plt.savefig(f'plots/SS_conc_{title}_box.png')
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
    conc_mod_obs["Measurement"] = conc_model_rename["Measurement"]
    return conc_mod_obs

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # data paths
    data_dir = "pd_files/"

    try:
        os.mkdir('plots')
    except OSError:

        pass



    var = ['ss_acc', 'ss_tot']
    mac_names = ['SS submicrom', 'SS supermicron']
    old_cols = ['SS Concentration (µg m$^{-3}$)','Total SS Concentration (µg m$^{-3}$)']
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


