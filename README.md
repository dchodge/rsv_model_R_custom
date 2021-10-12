# R package implementation of the cpp code from Hodgson et al. 2020
## Summary

This package provides an R implementation of the cpp code for the RSV intervention programmes which were evaluated in [Hodgson et al. 2020]( https://doi.org/10.1186/s12916-020-01802-8). In addition, this code provides an interface for custom RSV intervention programmes to be evaluated which are not explored in the paper. 

## Installation guide

*Quick-install instructions*

Clone the repository and look through the RMarkdown files in the /vingettes folder.

## Outline of model structure

After loading the relevant packages, a Rcpp module is called which contains the cpp code for the intervention model (see `src/RunInterventions.h` and `src/RunInterventions.cpp`). The Rcpp module must be paramterised with a large set of data (see `RunInteventionsClass.R` for full list). 

Once this is done, the paramters associated with the intervention programme calendarrs must be defined. There are seven possible different vaccination calendars including:
* pal, palivizumab programme
* mAB_VHR, mabs programmes aimed at very-high-risk infants (< 9 months old)
* mAB_HR, mabs programmes aimed at high-risk infants (< 12 months old)
* mAB_LR, mabs programmes aimed at low-risk persons (everyone who's not high- or very-high-risk)
* mat, maternal programme
* LAV_HR, vaccination agaisnt high-risk infants (<12 months old)
* LAV_LR, vaccination agaisnt low-risk persons.

Each calendar requires a list of parameters, these are:
* id, indicator stating if the calendar should be implemented. (If calendar is not defined it is not implemented.)
* age_id, age group which should be targeted (binary vector of values).
* t_start, start time in weeks (week 1 starting Jul 1st)
* t_end, end time in weeks, (week 1 starting Jul 1st) 
* eff, efficacy value of the programme, usually needs a seed value.
* cov, coverage of the programme (between 0 and 1).
There are some parameters specific to a programme type:
* catchup (mAB_HR, mAB_LR), whether a catchup programme should be included during the first four weeks of the programme.
* age_id_catchup (mAB_HR, mAB_LR), age group which should be targeted by catchup (binary vector of values).
* uptake_type (LAV_HR, LAV_LR), if = 2 then uptake is assumed to be constant.  If equal to 1 then a uptake vector must be provided (See remake_hodgson.R for an examples of this).
* eff_inf/eff_mat (mat), efficay of vaccination in prevevnt disease in infants (eff_inf), and in mothers (eff_mat)

In addition, a list called `vac_par_info` is needed, where:
* om_mab, is the average rate of loss of immunity due to long-acting mabs in days (default = 1 / 250)
* direct, boolean indicating if the direct effects of vaccination should be calcualted only (default = false)
* xi_boost, a value indicating the proportional increase of decrease in duration of immunity due to maternal vaccination relative to normal maternal immunity.

Please see both `rsv_inter_custom.Rmd` and `rsv_inter_hodgson.Rmd` for examples of how to implement this practically. 

## Overview of files
### /data folder
* `cpp_model_output/outcomes_all.RDS`, is dataframe which contains the annual incidence for each age group of each healthcare outcome under each intervention programme from the cpp code.
* `cpp_model_output/inci_strat.RDS`, is dataframe which contains the annual incidence of each age group, social group, and risk group under each intervention programme  from the cpp code.
* `inter_model_input/rsv_data_uk.RData`, a list of UK-specific demographic data used to define the ODEs of the transmission model (see `inter_model_input/rsv_data_desc.R` for a full description).
* `inter_model_input/inter_data_uk.RData`, a further list of UK-specific demographic data used to define the ODEs of the transmission model. (See `R/RunInterventionsClass.R` for a full description).
* `inter_model_input/seed_samples.csv`, a vector of integers corresponding to the seed values in the random number generation
* `inter_model_input/hodgson_programmes.RData`, .Rdata file of the inputs required to run the 15 intervention programmes, created from the `R/remake_hodgson.R` rile
* `inter_model_output/all_prog_sim.RData`, output, (incidence, costs, qaly) from the R implementation of the 15 intervention programmes from Hodgson et al. (see `rsv_inter_hodgson.Rmd`)

 ### /R folder
 * `RunInterventionsClass.R`, a function to initialised the RunInterventions class
 * `vac_cal.R`, functions to calculate the daily vaccination calendar, given a set of parameters values assocaited with the vaccination programmes. 
 * `cal_outcomes.R`, Converting the incidence from the RunInterventions model into health and economic outcomes
 * `remake_hodgson.R`, parameter definitions for the 15 programmes outlined in Hodgson et al. 2020

### /src folder
* The Rcpp/cpp code for the RSV intervention model is contained here. 

 ### /vignettes
 * `rsv_inter_hodgson.Rmd`, a walktrhough how to evaluate the 15 intervention programmes outlined in Hodgson et al.
 * `rsv_inter_custom.Rmd`, a walkthrough how to evaluate your own intervention programme and output the health and economic outcomes. 
 * `plt_outcomes.Rmd`, some plotting functions for the 15 intervention programmes outlined in Hodgson et al. 2020
 * `r_cpp_compare.Rmd`, Comparing the outputs of the cpp model and the R model.

## Linked publications

If you wish to use any part of this code, please reference

Hodgson, D., Pebody, R., Panovska-Griffiths, J. et al. Evaluating the next generation of RSV intervention strategies: a mathematical modelling study and cost-effectiveness analysis. BMC Med 18, 348 (2020). https://doi.org/10.1186/s12916-020-01802-8

## Contact details

Please email david.hodgson@lshtm.ac.uk with any queries relating to this code.
