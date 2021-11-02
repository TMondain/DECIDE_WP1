
# DECIDE Species Distribution Modelling workflow

This repository contains all the code necessary to run the species distribution models (SDMs) for the DECIDE project. These models are currently being run for both day-flying moths and butterflies, with plans to run orthoptera and all moths at a later stage. All of the methods are coded in 'R', but some of the scripts need to be run on the Lotus HPC, because of the computational requirements.

# Directory structure

Files are organised like so:

-   `scripts/`

    -   `1_data_preparation/` - Processing the raw environmental data and species data
        -  `environmental/` - processing the landcover map and cliamte data ready for use in SDMs
        -  `species/` - (subdivided again into `butterfly` and `moth`) processing the species occurence data ready for SDMs
    -   `2_effort_mapping/` - creating the effort layer that is used to downweight the decide score
    -   `3_run_sdms/` - running the SDMs from the processed data
    -   `4_create_decide_score/` - combining the SDM uncertainty and 
    -   `x_exploratory_analyses/` - a place for small standalone scripts used to testing or diagnosing

-   `data/` - Inputs for the SDMs
    -   `raw_data/` - Containing the raw data files as it was provided/downloaded from source
        -  `environmental/` - eg. raw landcover map as downloaded
        -  `species/` - eg. species records as downloaded from iRecord
    -   `derived_data/` - Containing data that has been manipulated
        -  `environmental/` - Prepared environmental data eg. landcover map adapted to 100x100m resolution ready for SDM 
        -  `species/` - Prepared species occurence data ready for SDM

-   `outputs/` - Outputs from the SDMs

The `scripts/` folder contains the bulk of the work. All functions are found in the subdirectory `scripts/functions/`. All of the modelling relies on the functions in the script: `scripts/functions/edited_rob_functions.R`; running the models on Lotus relies on the scripts prefixed with 'lotus'. All of the SDM functions are based on **Rob Boyd's** `soaR` package (https://github.com/robboyd/soaR).

# Workflow and functions

The workflow for the SDMs are broadly divided into two sections, sorting out the data and running the models. Sorting out the data is all done on my own computer (except for maybe the digital elevation map layer - see below) and all the modelling is done on the lotus HPC.

## 1. Data sorting

### Environmental data

First, there are several data layers that need to be downloaded:

-   UKCEH Modal Landcover map, any year, 25m resolution raster
-   HADUK monthly total rainfall and minimum and maximum temperature for the UK and years of interest. (Currently 2010-2019 inclusive)
-   Copernicus digital elevation map, 25m resolution, tiles X and Y for GB

These need to be converted to the desired 100m resolution maps and reprojected onto the same grid. The code to get the environmental data processed into the final format are contained in the scripts:

`scripts/main_scripts/environmental/Convert CEHLCM_25_to_100.R` `scripts/main_scripts/Envrinomental_data_sort_final.Rmd`

The first converts the UKCEH LCM from the 25m resolution to the 100m percentage cover of each of the 21 landcover classifications (be warned, this script isn't very well commented or tidy - needs improving). The second, covers all the rest of the data processing (except elevation):

-   converting HAD-UK raw data into the 19 bioclim variables, and reducing these based on correlations between them
-   getting slope and aspect from the elevation data
-   masking the three layers to include only Great Britain
-   reprojecting all layers to the LCM projection

In the `Envrinomental_data_sort_final.Rmd` script, there is also code to combine the elevation rasters and convert them from 25m to 100m. This is currently commented out because it kept crashing my computer. However, it **has** worked in the past and so might work for you. I ended up having to run it on Datalabs, but this is all documented in the `Envrinomental_data_sort_final.Rmd` script.

### Species data

The only two taxa that have been analysed so far (08/10/2021) are butterflies and day-flying moths. The processing from raw data to analysable datasets is documented in the script: `scripts/main_scripts/Species_data_sort_final.Rmd`. For day-flying moths there is an extra script to explore the list length that individuals were found on, here: `scripts/main_scripts/moth/Explore_list_lengths_for_pseudoabsences_moths.Rmd`. This was to ensure that only day-flying moth species were included in the final analyses.

### Transfer to Lotus

The files generated from these two scripts, the environmental raster stack and the species datasets, need to be transferred to Lotus for the rest of the workflow. Directories will need to be changed before running the scripts.

## 2. The modelling

All of the scripts for doing the modelling, except for step 1, through to creating the DECIDE score are found in the `scripts/functions/` directory, with the prefix 'lotus' because all are run on the Lotus HPC. Step 1, generating the pseudoabsences, is still run on Lotus but I haven't yet created an automated script for this yet like I have done for the other steps. I won't go into too much detail here about how each of the scripts work, but will explain the general workflow. The way the models are currently set up, the scripts need to be run on your PC and then the folders created need to be transferred to your directory on Lotus.

The modelling workflow follows five steps (with scripts):

1.  Creating pseudoabsences (scripts discussed below)
2.  Running the SDMs and producing predictions (`automated_lotus_sdms_predictions.R`)
3.  Error check and rerun script (`error_check.R`)
4.  Combining the models for each species (`Lotus_combine_models.R`)
5.  Creating the DECIDE score (`lotus_decide_score.R`)

### Step 1 - Pseudoabsences

The pseudoabsence function is called `cpa()` and is found in the `scripts/functions/Edited_Rob_Functions.R` script. The scripts to run the pseudoabsences are found here `scripts/lotus/moth/Pseudoabsence_script.R` for moths and here `scripts/lotus/moth/pseudoabsence_scripts/Pseudoabsence_script_butterflies.R` for butterflies. For the moths, the code needs to be copy and pasted into the R terminal. For the butterflies, the script needs to be run on your own computer and then transferred to lotus to run.

#### Outputs

Each of the two scripts produces an rdata file which contains a list of presences and pseudoabsences for each species. So is a list of lists.

### Step 2 - Running the SDMs and producing predictions for GB

The functions for the modeling workflow are all found here `scripts/functions/Edited_Rob_Functions.R`. The models are run using the function `fsdm()` and the predictions are obtained from the outputs using `get_predictions()`. These functions have been coded to work with five different model types (codes in the modelling process and the package they're from in brackets): general linear model (glm, base R), general additive model ('gam', package `mgcv`), random forest ('rf', package `randomForest`), maxent ('me', package `dismo`) and lasso regression ('lrReg', package `glmnet`). Currently, GAMs fail when there are very few data, although this is expected given the number of covariates we are using. Lasso regression also fails because of the scale we are working at. Get an error to do with sparse matrices probably because of the large number of 0s in the matrices for carrying out the models at 100m resolution. Lasso regression works on smaller scales and would probably work at 1km resolution GB wide.

#### Outputs

Each run for each species produces the following output:

-   GB wide raster of mean predictions for each species across all the bootstrapped models
-   GB wide raster of standard deviation for each species across all the bootstrapped models
-   .csv of the AUCs for all model runs and the mean AUC
-   the model objects (containing each bootstrapped model)
-   the data used to run the models (does not have the folds the data were split into)

### Step 3 - Error check and failed SDMs

This script is run in the R terminal on Lotus. It checks models that have been run and identifies those that failed. It also finds and loads the error scripts to see why the models failed. Then you need to select the models that need rerunning (i.e. the models that failed because of node errors rather than too few data) and run the SDMs. The code to do this is all contained in the script, it resubmits jobs straight from the R terminal on Lotus. All the rerun model outputs are stored in the same place as the original models.

The error log retrieval script does **NOT** work for models that have been rerun, only for models that were run the first time. This is because of where the error scripts for the rerun models are stored.

### Step 4 - Combining the models for each species

Combines the model predictions for each species, taking the mean across the different modelling types weighted by the mean AUC across all bootstrapped models (from the step above). Before taking the mean, the code drops models whose AUC falls below a given threshold (which can be altered). If the AUCs of all models fall below the threshold then all models are used in the mean and flagged for later.

#### Outputs

-   GB wide mean predictions across all model types
-   GB wide mean standard deviation across all model types
-   a csv storing the names of all the models that failed

### Step 5 - Create the DECIDE score

Combine the models for all species within each taxon to create a DECIDE score. Has a bunch of options to create the score which can be changed and added to as needed. This also contains the code to downweight the score by recording activity. It automatically loads in the raw data files and creates the weighting layer within the script.

#### Outputs

-   DECIDE score raster layer for each taxa

## Where does everything go?

The DECIDE score layers for each taxa can be downloaded to our PC and then transferred to the appdev folder. The species-level files, if you need to look at them, can be transferred to the object store using s4cmd. Files on the object store can be read in by a notebook on datalabs that has been granted permission.

To transfer files to the object store, need to use the terminal in mobaxterm or other client. The code to run is:

    #############################
    ####     One time setup
    #############################

    # load jaspy
    module load jaspy

    mkdir convert

    cd convert

    # Activate a virtualenvironment so packages can be installed locally
    virtualenv venv

    # install s4cmd
    pip install s4cmd

    #############################


    ### do this every time
    cd /location/where/you/installed/convert
    . venv/bin/activate

    # then change working directory to location of .s3cmd configuration file:
    cd /location/of/s3cmd/configuration/file/

    # then run the code to transfer the files
    # this code works recursively because of -r
    s4cmd --endpoint=datalabs_url dsync -r -c 15 /location/of/files/to/transfer/ location/of/file/destination/on/object/store

# Fin.
