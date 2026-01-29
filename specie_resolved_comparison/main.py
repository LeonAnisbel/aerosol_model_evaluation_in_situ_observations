##### Import packges
import sys
from utils_functions import read_data_functions, global_vars
import var_executor


if __name__ == '__main__':
    # data paths
    exp = global_vars.exp_name
    data_dir = global_vars.data_directory
    main_dir = global_vars.main_data_dir
    loc_dir = 'Aerosol_sample_coordinates/'

    # Read all observational data (station or ship trajectory and aerosol concentration)
    dates, PASCAL, PI_ICE, CVAO, SVAL_15 = (
        read_data_functions.read_obs_data_loc(main_dir, loc_dir))
    pd_obs_mo_cvao = var_executor.interp_cvao(CVAO)

    # The argument given in the run command determines whether to interpolate model results at Svalbard ('SVD') or
    # other stations
    if sys.argv[1] == 'SVD':
        var_executor.interp_svd_rs(SVAL_15, pd_obs_mo_cvao)
    else:
        var_executor.interp_piice_pascal(PI_ICE, PASCAL, pd_obs_mo_cvao)
