from matplotlib import pyplot as plt
from sklearn.metrics import mean_squared_error
import math
from scipy.stats import linregress
import numpy as np

def rename_func(data_pd, col, na, new_na):
    """Renames columns of the data_pd data frame from the old name "na" to "new_na" """
    pd_new = data_pd
    list_col = data_pd[col].to_list()
    for i in range(len(list_col)):
        if list_col[i] == na:
            list_col[i] = new_na
    pd_new = pd_new.drop(columns=col)
    pd_new[col] = list_col

    return pd_new


def get_stat(obs, mod):
    """Computes statistical indexes"""
    MSE = mean_squared_error(obs, mod)

    RMSE = math.sqrt(MSE)
    diff = np.subtract(mod, obs)
    mean_bias = np.nanmean(diff)
    NMB = np.nansum(diff)/np.nansum(obs)

    # correlation coefficients
    res_lin_reg = linregress(obs, mod)
    pearsons_coeff = res_lin_reg.rvalue
    pval_corr = res_lin_reg.pvalue

    return RMSE, mean_bias,NMB, pearsons_coeff, pval_corr


def set_log_ax(axis, x, y, style):
    axis.loglog(x, y,
              color="black",
              linestyle=style,
              alpha=0.7,
              linewidth=0.5)


def add_species_text_name(ax, mol_name, loc1, loc2, f):
    ax.text(loc1, loc2,
            mol_name,
            fontsize= f,
            weight='bold',
            bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})



def customize_legend_sv_fig(pl, fig, title):
    """ Customize box plots  """
    handles, labels = pl.get_legend_handles_labels()
    fig.legend(handles=handles,
                       labels= labels,
                       ncol=4,
                       bbox_to_anchor=(0.53, 1.01),
                       loc='upper center',
                       fontsize=12)
    fig.tight_layout()

    plt.savefig(f'plots/{title}.png',dpi = 300)#{title}_{global_vars.exp_name}
    plt.close()


def get_marker_flier():
    return dict(marker='D',
              markerfacecolor='gray',
              markersize=2.5,
              linestyle='none')
