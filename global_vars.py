import os

exp_name = 'ac3_arctic'#'echam_base'#ac3_longno'#'ac3_gongst' #'ac3_arctic'
data_directory = f"/scratch/b/b381361/{exp_name}/"#f"/work/bb1005/b381361/my_experiments/{exp_name}/"#f"/scratch/b/b381361/{exp_name}/"
if exp_name == 'echam_base':
    variables_names = ['SS', 'OC']
else:
    variables_names = ['POL', 'PRO', 'LIP', 'SS', 'OC']

# plot_dir = './plots/'
# plot_dir_sfc = plot_dir + 'Sfc_conc_plots/global_'
# try:
#     os.mkdir('plots')
# except OSError:
#
#     pass
# try:
#     os.mkdir('plots/Sfc_conc_plots')
# except OSError:
#     pass
