# Interpolate aerosol concentration from ECHAM-HAM model grid to the station locations
> This project is used to interpolate ECHAM6.3-HAM2.3 model grid to the station locations where marine aerosol observations were available
> 
> The scripts use the aerosol-climate model output (this data is archived on Levante HPC system). See more information and results in [Leon-Marcos et al. (2025a)](https://doi.org/10.5194/egusphere-2025-2829) and Anisbel Leon Marcos' Doctoral Thesis (Chapter 5) at the University of Leipzig.
> 
> The conda environment for this project is contained in env.yml file. Run "conda activate env.yml" to set up the environment for this project.

* Run main.py to perform the interpolation adding "SVD" or "NAO" as argument words to perform the interpolation of lipids and proteins or polysaccharides, respectively. 
 Also "sbatch run_python_main.sh" can be used to run it in levante with "SVD" as default option. 
 This script will create pickle files with the observations and interpolated model data ([outputs](outputs)/*_conc.pkl) for each biomolecule that can later be used for the plots.


* global_vars.py contains all global variables
  1. Directories to model and observational data
  2. model experiment
  3. model tracers considered in the comparison


* mixed_box_plots.py will create all relevant figures (as box plots):\
    -- Box plot of observed and interpolated values of all organic marine aerosol tracers concentration and OMF (PMOA_plots.py). OMF values from observations and interpolated model results ([outputs](outputs)/*_omf.pkl) are generated within the project [OMF_model_evaluation](https://github.com/LeonAnisbel/OMF_model_evaluation.git).\
    -- Box plot of observed and interpolated values of sea salt concentration (SS_plots.py file)\
    -- Box plot of observed (OM) and interpolated values of OC (SPMOAoff experiment) concentration and OC+PMOA (SPMOAon experiment) (OC_plots.py file)\
    -- Statistics of this comparison is file statistics.log

### Analysis of Arctic stations, Mace Head and Amsterdam Island
* Run main_Arctic.py to perform the interpolation of PMOA at inner Arctic stations, Mace Head and Amsterdam Island. Add "all", "MH" or "AI" as argument in the python execution command line to interpolate model values, respectively.
  Additionally, "sbatch run_python_main_arctic.sh" can be used to run it in levante


* plot_arctic_seasonaltiy.py creates figures of monthly observed and interpolated values at Mace Head (August 2015 and 2018) and other stations within the Arctic
