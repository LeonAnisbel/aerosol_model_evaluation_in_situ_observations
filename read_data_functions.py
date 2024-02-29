import xarray as xr
import glob, os
from matplotlib import cm
import mod_interp_obs_concat
import plot_functions
import numpy as np
import codecs
import pandas as pd
import itertools



def read_do_mean_plot(path,plot_dir_sfc,poles = False,make_plots = False, get_aver = False):
    da = xr.open_mfdataset(path + '*.01_emi.nc', concat_dim='time', combine='nested')
    da_emi_mean = da.mean(dim='time')
    da = xr.open_mfdataset(path + '*.01_burden.nc', concat_dim='time', combine='nested')
    da_burden_mean = da.mean(dim='time')

    if get_aver:
        # print('fraction')
        num = (da_emi_mean['emi_LIP'] +
               da_emi_mean['emi_POL'].rename('emi_LIP') +
               da_emi_mean['emi_PRO'].rename('emi_LIP'))
        num = num.compute()
        num_regrid = interp_func.regrid(da_emi_mean, num)
        num_regrid_mean = sum(num_regrid) * (111.1*1.5)**2 *  1e-9 / 3.171e-8 / 1e-6
        print('Global MOA emission (Tg/yr)')
        print(np.nanmean(num_regrid_mean.mean))
        print('')

        burden = (da_burden_mean['burden_LIP'] +
                  da_burden_mean['burden_POL'].rename('burden_LIP')  +
                  da_burden_mean['burden_PRO'].rename('burden_LIP'))
        burden = burden.compute()
        burden_regrid = interp_func.regrid(da_burden_mean, burden)
        burden_regrid_mean = sum(burden_regrid) * (111.1*1.5)**2  * 1e-9 / 1e-6
        print('Global BURDEN of MOA (Tg)')
        print(np.nanmean(burden_regrid_mean))
        print('')


        burden_ss = da_burden_mean['burden_SS']
        burden_ss = burden_ss.compute()
        burden_ss_regrid = interp_func.regrid(da_burden_mean, burden_ss)
        burden_ss_regrid_mean = sum(burden_ss_regrid) * (111.1*1.5)**2  * 1e-9 / 1e-6
        print('Global BURDEN of SS (Tg)')
        print(np.nanmean(burden_ss_regrid_mean))
        print('')


        den = da_emi_mean['emi_SS'].rename('emi_LIP')
        den = den.compute()
        den_regrid = interp_func.regrid(da_emi_mean, den)
        den_regrid_mean = sum(den_regrid) * (111.1*1.5)**2 *  1e-9 / 3.171e-8 / 1e-6
        print('Global SS emission (Tg/yr)')
        print(np.nanmean(den_regrid_mean))
        print('')


        div = num_regrid/den_regrid
        # weights = np.cos(np.deg2rad(num.lat))
        # print(div.weighted(weights).mean(dim={"lat", "lon"})* 100)
        print('fraction oc organics in SS')
        print(div.mean()* 100)


    if make_plots:
        print('Plotting mean emission maps')
        unit: str = 'µg m$^{-2}$ s$^{-1}$'
        cmap = cm.rainbow
        title_id = 'Emission mass flux'
        plot_functions.plot_map_mean(da_emi_mean, f'MOA {title_id}',
                                     ['emi_POL','emi_LIP','emi_PRO'], 0,0.8,
                                    plot_dir_sfc,1e-11,'10$^{-2}$'+unit,cmap,
                                     poles = poles,total = True)



        print('Plotting burden maps')
        unit: str = 'µg m$^{-2}$'
        cmap = cm.rainbow
        title_id = 'Atmospheric burden of '
        print('Plotting mean burden maps')
        plot_functions.plot_map_mean(da_burden_mean, f'{title_id} MOA',
                                     ['burden_POL', 'burden_LIP', 'burden_PRO'], 0, 1.5,
                                    plot_dir_sfc, 1e-6, 'mg m$^{-2}$',cmap, poles = poles, total=True)




def read_do_monthly_mean_plot(path,plot_dir_sfc,poles = False,plot_month = False, plot_season = False):

    da_month_emi = []
    da_month_burden = []
    da_month_ice = []


    for l in range(12):
        if l < 9:
            mo_str = f'0{l+1}'
        else:
            mo_str = f'{l+1}'

        print('reading ', path + '*' + mo_str + '.01_emi.nc')

        da = xr.open_mfdataset(path + '*' + mo_str + '.01_emi.nc', concat_dim='time', combine='nested')
        da_month_emi.append(da.mean(dim='time'))

        if plot_month:

            da3 = xr.open_mfdataset(path + '*' + mo_str + '.01_burden.nc', concat_dim='time', combine='nested')
            da_month_burden.append(da3.mean(dim='time'))

        if plot_season:
            da4 = xr.open_mfdataset(path + '*' + mo_str + '.01_echam.nc', concat_dim='time', combine='nested')
            da_month_ice.append(da4.mean(dim='time'))

   ### -------------------------------------------------------------------------------------------------------------------
    if plot_month:

        print('Plotting emi maps')
        unit: str = 'µg m$^{-2}$ s$^{-1}$'
        cmap = cm.rainbow
        title_id = 'Emission mass flux'


        plot_functions.plot_map(da_month_emi, f'Polysaccharides {title_id}', 'emi_POL', 0,10,
                 plot_dir_sfc,1e-14,'10$^{-4}$'+unit,cmap,poles = poles,total = False)
        plot_functions.plot_map(da_month_emi, f'Proteins {title_id}', 'emi_PRO', 0,5,
                 plot_dir_sfc,1e-13,'10$^{-4}$'+unit,cmap,poles = poles,total = False)
        plot_functions.plot_map(da_month_emi, f'Lipids  {title_id}', 'emi_LIP', 0,1.5,
                 plot_dir_sfc,1e-11,'10$^{-2}$'+unit,cmap,poles = poles,total = False)
        plot_functions.plot_map(da_month_emi, f'Sea Salt {title_id}', 'emi_SS', 0,30,
                 plot_dir_sfc,1e-11,'10$^{-2}$'+unit,cmap,poles = poles,total = False)

        plot_functions.plot_map(da_month_emi, f'MOA {title_id}', ['emi_POL','emi_LIP','emi_PRO'], 0,1.5,
                 plot_dir_sfc,1e-11,'10$^{-2}$'+unit,cmap,poles = poles,total = True)

        print('Finished Plotting emi')

        ### -------------------------------------------------------------------------------------------------------------------

        print('Plotting burden maps')
        unit: str = 'µg m$^{-2}$'
        cmap = cm.rainbow
        title_id = 'Atmospheric burden of '

        plot_functions.plot_map(da_month_burden, f'{title_id} Polysaccharides', 'burden_POL', 0, 20,
                 plot_dir_sfc, 1e-9,  unit,cmap, poles = poles, total=False)
        plot_functions.plot_map(da_month_burden, f'{title_id} proteins', 'burden_PRO', 0, 100,
                 plot_dir_sfc, 1e-9, unit,cmap, poles = poles, total=False)
        plot_functions.plot_map(da_month_burden, f'{title_id} lipids', 'burden_LIP', 0, 2,
                 plot_dir_sfc, 1e-6,  'mg m$^{-2}$',cmap, poles = poles, total=False)
        plot_functions.plot_map(da_month_burden, f'{title_id} Sea Salt', 'burden_SS', 0, 50,
                 plot_dir_sfc, 1e-6, 'mg m$^{-2}$',cmap, poles = poles, total=False)

        plot_functions.plot_map(da_month_burden, f'{title_id} MOA', ['burden_POL', 'burden_LIP', 'burden_PRO'], 0, 2,
                 plot_dir_sfc, 1e-6, 'mg m$^{-2}$',cmap, poles = poles, total=True)

        print('Finished Plotting burden')

    if plot_season:
    # Plot emiss flux for summer and winter
        print('Plotting maps for seasons')
        DEF_ice = (da_month_ice[11]['seaice'].compute() +
                   da_month_ice[0]['seaice'].compute() +
                   da_month_ice[1]['seaice']).compute() * 100 / 3

        JJA_ice = (da_month_ice[5]['seaice'].compute() +
                   da_month_ice[6]['seaice'].compute() +
                   da_month_ice[7]['seaice']).compute() * 100 / 3

        DEF_ice = DEF_ice.to_dataset(name = 'var')
        JJA_ice = JJA_ice.to_dataset(name = 'var')

        DEF_ice = DEF_ice.to_array()[0]
        JJA_ice = JJA_ice.to_array()[0]


        # DEF_ice = DEF_ice.where(DEF_ice < 20, 0.)
        # JJA_ice = JJA_ice.where(JJA_ice < 20, 0.)
        #
        for p in range(len(DEF_ice[0])):
            for s in range(len(DEF_ice)):
                if DEF_ice[s][p] < 20:
                    DEF_ice[s][p] = 0.
                if JJA_ice[s][p] < 20:
                    JJA_ice[s][p] = 0.

        T_vars_summer = []
        T_vars_winter = []
        lat = np.array(da_month_emi[11].lat.values)
        lon = np.array(da_month_emi[11].lon.values)

        vars = ['emi_POL','emi_PRO','emi_LIP']
        for i in vars:
            DEF = (da_month_emi[11][i].compute() +
                   da_month_emi[0][i].compute() +
                   da_month_emi[11][i].compute())
            JJA = (da_month_emi[5][i].compute() +
                   da_month_emi[6][i].compute() +
                   da_month_emi[7][i].compute())

            # DEF = DEF.where(DEF_ice != 0.,  np.nan)
            # JJA = JJA.where(JJA_ice != 0.,  np.nan)
            DEF = DEF.to_dataset(name = 'var')
            JJA = JJA.to_dataset(name = 'var')

            DEF = DEF.to_array()[0]
            JJA = JJA.to_array()[0]

            for p in range(len(DEF_ice[0])):
                for s in range(len(DEF_ice)):
                    if DEF_ice[s][p] != 0.:
                        DEF[s][p] = np.nan
                    if JJA_ice[s][p] != 0.:
                        JJA[s][p] = np.nan

            T_vars_winter.append(DEF/3)
            T_vars_summer.append(JJA/3)

        x, y = np.meshgrid(lon, lat)

        print('Plotting maps for seasons')
        title = 'Multiannual Emission Mass mass flux (JJA)'
        factor = 1e12
        units = '10$^{-3}$ µg m$^{-2}$ s$^{-1}$'

        plot_functions.function_plot_mult(x, y, T_vars_summer[0] * factor ,
                                          T_vars_summer[1] * factor ,
                                          T_vars_summer[2] * factor,
                                          JJA_ice,
                                          units,
                                        title)

        title = 'Multiannual Emission Mass mass flux (DJF)'
        plot_functions.function_plot_mult(x, y, T_vars_winter[0] * factor,
                                          T_vars_winter[1] * factor ,
                                          T_vars_winter[2] * factor,
                                          DEF_ice,
                                          units,
                                         title)


        print('Finish plots maps for seasons')


################################
def read_obs_data_loc():
    #### Reading station locations
    main_dir = '/work/bb1005/b381361/echam_postproc/'
    plot_dir = '../OMF_plots/'
    loc_dir = 'Aerosol_sample_coordinates/'
    files = ['PASCAL_lat_lon_aer.csv','PI_ICE_lat_lon_aer.csv']
    pol_data = 'AER_OMF_pol_lip_pro_all_sizes.csv'
    pro_li_data = 'AER_OMF_pro_lip.csv'

    doc = codecs.open(main_dir+loc_dir+files[0],'r','UTF-8') #open for reading with "universal" type set
    PASCAL_loc = pd.read_csv(doc, sep=',')

    doc = codecs.open(main_dir+loc_dir+files[1],'r','UTF-8') #open for reading with "universal" type set
    PI_ICE_loc = pd.read_csv(doc, sep=',')

    # converting 'Date/Time' column to datetime data type
    PASCAL_loc['Date/Time'] =  PASCAL_loc['Date/Time'].apply(pd.to_datetime)
    PI_ICE_loc['Date/Time'] =  PI_ICE_loc['Date/Time'].apply(pd.to_datetime)

    doc = codecs.open(main_dir + pol_data, 'r', 'UTF-8')  # open for reading with "universal" type set
    data = pd.read_csv(doc, sep=',')

    # converting 'Date/Time' column to datetime data type
    data['Start Date/Time'] = data['Start Date/Time'].apply(pd.to_datetime)
    data['End Date/Time'] = data['End Date/Time'].apply(pd.to_datetime)

    # Save data for PI_ICE
    data_PI_ICE = data[:24]
    data_PI_ICE_subm = data_PI_ICE[:8]  # data for submicron aerosol size
    data_PI_ICE_super = data_PI_ICE[8:16]  # data for supermicron aerosol size
    data_PI_ICE_tot = data_PI_ICE[16:]  # data for all aerosol sizes

    # Save data for PASCAL
    data_PASCAL = data[24:33]
    data_PASCAL_subm = data_PASCAL[:3]  # data for submicron aerosol size
    data_PASCAL_super = data_PASCAL[3:6]  # data for supermicron aerosol size
    data_PASCAL_tot = data_PASCAL[6:]  # data for all aerosol sizes

    data_CVAO = data[33:69]  # data for PM1 aerosol sizes
    data_CVAO_subm = data_CVAO[:28]
    # data_CVAO_super = data_CVAO[28:36]
    data_CVAO_tot = data_CVAO[28:]

    data_SVAL = data[69:]
    data_SVAL_14_15_subm = data_SVAL[:23]
    data_SVAL_18_19_subm = data_SVAL[23:45]
    data_SVAL_14_15_tot = data_SVAL[45:69]
    data_SVAL_18_19_tot = data_SVAL[69:]
    # print(data_CVAO_subm,data_SVAL_14_15_subm)


    dates = []
    for i,d in enumerate(data['Start Date/Time'].dt.month):
        dates.append([d,data['Start Date/Time'].dt.year[i]])
        dates.append([data['End Date/Time'].dt.month[i],data['End Date/Time'].dt.year[i]])
    dates.sort()
    dates = list(k for k,_ in itertools.groupby(dates))

    PI_ICE = [PI_ICE_loc,data_PI_ICE,data_PI_ICE_subm,data_PI_ICE_super,data_PI_ICE_tot]
    PASCAL = [PASCAL_loc,data_PASCAL,data_PASCAL_subm,data_PASCAL_super,data_PASCAL_tot]
    CVAO = [data_CVAO,data_CVAO_subm, data_CVAO_tot]
    SVAL_14_15 = [data_SVAL,data_SVAL_14_15_subm,data_SVAL_14_15_tot]
    SVAL_18_19 = [data_SVAL,data_SVAL_18_19_subm,data_SVAL_18_19_tot]

    return dates, PASCAL, PI_ICE, CVAO, SVAL_14_15, SVAL_18_19




def read_model_spec_data(file):
    return xr.open_mfdataset(file, concat_dim='time', combine='nested')




#reading data
def read_nc_ds(files,path):
    data_model = []
    data_month = [[] for i in range(12)]
    da_month_mean = [[] for i in range(12)]

    pp = len(path)

    for fi in files:
        if int(fi[pp+11:pp+15]) == 2015 or int(fi[pp+11:pp+15]) == 2010:# or int(fi[-16:-12]) == 2015:# or int(fi[-16:-12]) == 2019:
            da = xr.open_mfdataset(fi)
            data = da.mean(dim = 'time')
            data_model.append(data)
            for i in range(12):
                if int(fi[pp+15:pp+17]) -1 == i:
                    data_month[i].append(data)

    for l,da in enumerate(data_month):
        if len(da) > 0:
            dada = xr.concat(da,dim = 'time')
            da_month_mean[l] = dada.mean(dim = 'time')

    return xr.concat(data_model,dim = 'time'), data_month, da_month_mean



