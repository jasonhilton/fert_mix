Forecasting fertility with Bayesian parametric mixtures
================
Jason Hilton
13 March 2019

This repository provides the code needed to produce the results for the paper "Forecasting Fertility using Bayesian Parametric Mixture Models" by Jason Hilton, Erengul Dodd, Jonathan J. Forster, Peter W.F. Smith and Jakub Bijak.

The code is design so that the results can be replicated by running a series of `R` and shell scripts.

Requirements
------------

Running the code requires a recent version of `R` (&gt;3.5 is probably best) and of `rstan`, the `R` interface to (see here)\[<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>\]. To run the assessment across a range of models and countries, a High Performance Computing (HPC) environment with the same setup and which accepts jobs with Portable Batch System (PBS) is needed.

A range of other R packages is needed. These can be installed by running the `install_packages.R` script from same directory as this read-me file (this will make copy files to your machine):

``` bash
Rscript scripts/install_packages.R
```

Get the data
------------

Before downloading the data, it is necessary to obtain a [Human Fertility Database](https://www.humanfertility.org/cgi-bin/main.php) account, and to set environment variables containing your username and password. You can do this by adding the following lines to your `.Renviron.site` file (create one if you do not already have one), which will be in the directory specified by running the command `R.home()` in the R console, in the `etc` sub-folder.

``` bash
HFD_user=user
HFD_pass=pass
```

The data can then be obtained running the following shell commands:

``` bash
Rscript scripts/get_data.R VV
Rscript scripts/get_data.R RR
```

The two arguments to the scripts, "VV" and "RR", refer to the two different data configurations, age-cohort and age-period respectively.

Run the model
-------------

A single model can be run by running another script, passing in a reference to a config file that contains the model setup, together with specifying the country for which you want to run the model. The default settings are provided in the file `config/base_hb.yaml`, but you can write your own configuration by editing this file and re-saving, if desired. To run the model with these settings for England and Wales, run:

``` bash
Rscript scripts/run_model.R GBRTENW base_hb
```

GBRTENW is the Human Fertility Database code for England and Wales. The other possibilities are FRATNP, USA, and SWE for France, USA and Sweden respectively.

Run a batch of models
---------------------

**Requires HPC environment**

A set of model runs are specified and save to csv through the file `script/make_csv_files.R`. The `bash` script `send_csv_jobs.sh` runs this R script and then submits a PBS job array to run every model specified therein, running assessments on each model and saving out the results. All that is required is to run the shell script from within the HPC environment:

``` bash
send_csv_jobs.sh
```

This will save results to a sub-folder of the `mcmc_results` folder named `run_` plus the date and time at which the model was run ( format `yyyymmdd_hhmmss`). After all runs are completed, this folder will contain sub-folders for each country with one `.rds` file for each model specification, each contain a list with the data used to fit the model, the fitted stan model object, and the git hash identifying the version of the code used. There will also be one summary `assessment_` file for each model, indexed by the line in the csv file it relates to.

Recalculate LOOIC values
------------------------

Pareto Smoothed Importance Sampling approximations of Leave-One-Out predictive densities (Vehtari et al. 2016) are used to assess the model's ability to predict a single left out point within the range of the observed data. However, many models have high values of the Pareto-K diagnostic for some observations, indicating that the calculated values of the Leave-one-out Information Criterion are unreliable. After the batch of the stan model sampling runs triggered by the script `send_jobs_csv.sh` is complete, an additional script is called (`scripts/compute_loos.R`) that identifies problematic points for each model. This script also specifies a new set of model runs to refit the model with groups of such points left out, so that the 'true' value of the predictive score can be calculated. These additional model runs are detailed in a csv file in the folder `run_specs/lko`.

These jobs can be executed with the HPC environment through another script.

``` bash
send_jobs_lko_csv.sh
```

At this point, it is recommended to copy the mcmc\_results and processed folders from your HPC environment to your local machine.

Note that the generated results folder has a time-stamp in its file name. `run_yyyymmdd_hhmmss`.

Finally, the corrected values of the LOOIC metric is constructed from the results of the refitting exercise with the following script.

``` bash
Rscript scripts/recompute_loo.R
```

Run comparator models
---------------------

Results can be compared against some competitor models. A minimalist R package written by this author `bfcf` (Bayesian Forecasting of Cohort Fertility) is required to implement the method of Schmertmann et al. This is also provided for reviewers. This package is based entirely on the code kindly made publicly available at <http://schmert.net/cohort-fertility/>. The code for fitting the Myrskylla model is adapted from that generously provided in the supplemental material to the Bohk-Ewald et al. 2018.

The results of these can be obtained (and saved in the `processed` folder) by running the scripts below.

``` bash
Rscript scripts/freeze_rates_eval.R
Rscript scripts/process_myrskylla.R
Rscripts scripts/process_schmertmann.R
```

Create Plot data
----------------

The plot of the posterior distribution of rates requires a bit of simple preprocessing. This is done using the command:

``` bash
Rscript scripts/get_rate_plot_data.R
```

Create the paper
----------------

Finally, the paper can be created by running the script below.

``` bash
Rscript scripts/process_rmd.R
```

System info
-----------

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 7 x64 (build 7601) Service Pack 1
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] splines   stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] LexisPlotR_0.3        ggfan_0.1.3           scales_1.0.0         
    ##  [4] abind_1.4-5           matrixStats_0.54.0    loo_2.1.0            
    ##  [7] rstan_2.19.2          ggplot2_3.2.1         StanHeaders_2.18.1-10
    ## [10] yaml_2.2.0            knitr_1.24            rmarkdown_1.14       
    ## [13] stringr_1.4.0         stringi_1.4.3         readr_1.3.1          
    ## [16] magrittr_1.5          tibble_2.1.3          tidyr_0.8.3          
    ## [19] dplyr_0.8.3           purrr_0.3.2          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.2        pillar_1.4.2      compiler_3.5.1   
    ##  [4] prettyunits_1.0.2 tools_3.5.1       pkgbuild_1.0.4   
    ##  [7] zeallot_0.1.0     digest_0.6.20     evaluate_0.14    
    ## [10] gtable_0.3.0      pkgconfig_2.0.2   rlang_0.4.0      
    ## [13] cli_1.1.0         parallel_3.5.1    xfun_0.8         
    ## [16] gridExtra_2.3     withr_2.1.2       vctrs_0.2.0      
    ## [19] hms_0.5.0         stats4_3.5.1      grid_3.5.1       
    ## [22] tidyselect_0.2.5  glue_1.3.1        inline_0.3.15    
    ## [25] R6_2.4.0          processx_3.4.1    callr_3.3.1      
    ## [28] ps_1.3.0          backports_1.1.4   htmltools_0.3.6  
    ## [31] assertthat_0.2.1  colorspace_1.4-1  lazyeval_0.2.2   
    ## [34] munsell_0.5.0     crayon_1.3.4

Iridis
------

``` r
readRDS("processed/hpc_sesh.Rds")
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux Server release 6.10 (Santiago)
    ## 
    ## Matrix products: default
    ## BLAS: /home/local/software/R/3.5.1/lib64/R/lib/libRblas.so
    ## LAPACK: /home/local/software/R/3.5.1/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ## [1] C
    ## 
    ## attached base packages:
    ## [1] splines   stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] abind_1.4-5        matrixStats_0.54.0 loo_2.0.0         
    ##  [4] rstan_2.18.2       StanHeaders_2.18.1 ggplot2_3.0.0     
    ##  [7] git2r_0.23.0       yaml_2.1.19        stringr_1.3.1     
    ## [10] stringi_1.2.3      readr_1.1.1        magrittr_1.5      
    ## [13] tibble_1.4.2       tidyr_0.8.1        dplyr_0.7.6       
    ## [16] purrr_0.2.5       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.17      pillar_1.1.0      compiler_3.5.1   
    ##  [4] plyr_1.8.4        bindr_0.1.1       prettyunits_1.0.2
    ##  [7] tools_3.5.1       pkgbuild_1.0.2    gtable_0.2.0     
    ## [10] pkgconfig_2.0.1   rlang_0.2.1       cli_1.0.0        
    ## [13] parallel_3.5.1    bindrcpp_0.2.2    gridExtra_2.2.1  
    ## [16] withr_2.1.2       hms_0.3           stats4_3.5.1     
    ## [19] grid_3.5.1        tidyselect_0.2.4  glue_1.2.0       
    ## [22] inline_0.3.14     R6_2.2.2          processx_3.2.1   
    ## [25] callr_3.1.1       scales_0.5.0      ps_1.3.0         
    ## [28] assertthat_0.2.0  colorspace_1.4-0  lazyeval_0.2.1   
    ## [31] munsell_0.5.0     crayon_1.3.4
