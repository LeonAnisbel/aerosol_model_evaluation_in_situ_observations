##### Import packges
import sys
import read_data_functions
import global_vars, var_executor

# Press the green button in the gutter to run the script.

if __name__ == '__main__':
    # data paths
    exp = global_vars.exp_name
    data_dir = global_vars.data_directory

    dates, PASCAL, PI_ICE, CVAO, SVAL_14, SVAL_15, SVAL_18, RS_18_20 = (
        read_data_functions.read_obs_data_loc())
    pd_obs_mo_cvao = var_executor.interp_cvao(CVAO)

    if sys.argv[1] == 'SVD':
        var_executor.interp_svd_rs(SVAL_15, pd_obs_mo_cvao)
    else:
        var_executor.interp_piice_pascal(PI_ICE, PASCAL, pd_obs_mo_cvao)
