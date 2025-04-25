import os

main_data_dir = '../../'#'/work/bb1005/b381361/echam_postproc/'
exp_name = 'ac3_arctic'
data_directory = f"/work/bb1005/b381361/my_experiments/{exp_name}/"
yr_exp_var_names = {'0209': ['Start Date/Time', 'SS_µg_per_m3', 'WIOC_µg_per_m3', 'End Date/Time'],
                    '2015':['date', 'seasalt', 'PMOA'],
                    '2018':['date', 'seasalt', 'PMOA']}
if exp_name == 'echam_base':
    variables_names = ['SS', 'OC']
else:
    variables_names = ['POL', 'PRO', 'LIP', 'SS', 'OC']

try:
    os.mkdir('plots')
except OSError:

    pass

