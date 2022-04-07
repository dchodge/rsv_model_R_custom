
make_data_list <- function() {
  source("R/vac_cal.R") # generates the vaccination calendars
  source("R/calc_outcomes.R")  # calculate the outcomes
  # Load posteriors
  load(here("data", "inter_model_output", "posteriors.Rda"))  # posteriors from fitting in Hodgson et al. 2020
  # Load seeds
  seeds <- read.csv(here("data", "inter_model_input", "seed_samples.csv"), header = FALSE)[, 1] + 1

  # Ex.1 Long-acting monoclonal antibodies at HR, LR, and VHR, given seasonally with 90% coverage.
  ind_pal <- c(rep(1, 9), rep(0, 16)) # VHR (<8 months)
  ind_mabs <- c(rep(1, 1), rep(0, 24)) # Birth only (0 months only)
  G_plus <- c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

  require(triangle)

  pal_eff <- rweibull(1250, 12.4311, 0.772251) 
  mab_eff <- rtriangle(1250, 0.496, 0.871, 0.745) # efficacy of long-acting mabs

  # VHR, HR and LR at birth during the winter
  make_vac_program_info_none <- function(seed) {
    list(
        pal = list(id = TRUE, age_id = ind_pal, t_start = 0, t_end = 0, eff = pal_eff[seed], cov = 0.0)
      )
  }

  make_vac_program_info_pal <- function(seed) {
    list(
        pal = list(id = TRUE, age_id = ind_pal, t_start = 15 * 7, t_end = 32 * 7, eff = pal_eff[seed], cov = 0.9)
      )
  }

  make_vac_program_info_mabs_vhr_seasonal <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 15 * 7, t_end = 32 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }

  make_vac_program_info_mabs_seasonal <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }

  make_vac_program_info_mabs_year_round <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 0 * 7, t_end = 52 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 0, t_end = 52 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 0, t_end = 52 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }

  make_vac_program_info_mabs_seasonal_catchup <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_HR =  list(id = TRUE, catchup = TRUE, catchupnip = FALSE, age_id_catchup = G_plus, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_LR =  list(id = TRUE, catchup = TRUE, catchupnip = FALSE, age_id_catchup = G_plus, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }


  make_vac_program_info_mabs_seasonal_catchup_nip <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = TRUE, age_id = ind_mabs, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = TRUE, age_id = ind_mabs, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }

  ## Different assumptions about the duration of protection
  vac_par_info <- list(om_mab = 1 / 150, direct = FALSE, xi_boost = 1)

  data_list <- 
    list(
      none = make_vac_program_info_none, 
      pal = make_vac_program_info_pal, 
      vhr_seasonal = make_vac_program_info_mabs_vhr_seasonal,
      seasonal = make_vac_program_info_mabs_seasonal,
      year_round = make_vac_program_info_mabs_year_round, 
      seasonal_catchup = make_vac_program_info_mabs_seasonal_catchup, 
      seasonal_catchup_nip = make_vac_program_info_mabs_seasonal_catchup_nip,
      vac_par_info = vac_par_info,
      pal_eff = pal_eff,
      mab_eff = mab_eff,
      post = post,
      seeds = seeds,
      S = 100
  )
  save(data_list, file = here::here("outputs", "nirsevimab", "data_list.RData"))
}


run_model_base <- function(data_list) {
  seeds <- data_list$seeds
  post <- data_list$post
  vac_par_info <- data_list$vac_par_info
  S <- data_list$S

  output_default_none <- run_sample_custom(seeds[1:S], data_list$none, vac_par_info, 0, post)
  output_default_base <- run_sample_custom(seeds[1:S], data_list$pal, vac_par_info, 0, post)
  output_season_vhr_base <- run_sample_custom(seeds[1:S], data_list$vhr_seasonal, vac_par_info, 0, post)
  output_season_base <- run_sample_custom(seeds[1:S], data_list$seasonal, vac_par_info, 0, post)
  output_yr_base <- run_sample_custom(seeds[1:S], data_list$year_round, vac_par_info, 0, post)
  output_season_catchup_base <- run_sample_custom(seeds[1:S], data_list$seasonal_catchup, vac_par_info, 0, post)
  output_season_catchup_nip_base <- run_sample_custom(seeds[1:S], data_list$seasonal_catchup_nip, vac_par_info, 0, post)


  save(output_default_none, file = here("outputs", "nirsevimab", "impact", "base", "none.RData"))
  save(output_default_base, file = here("outputs", "nirsevimab", "impact", "base", "status_quo_base.RData"))
  save(output_season_vhr_base, file = here("outputs", "nirsevimab", "impact", "base", "output_season_vhr_base.RData"))
  save(output_season_base, file = here("outputs", "nirsevimab", "impact", "base", "output_season_base.RData"))
  save(output_yr_base, file = here("outputs", "nirsevimab", "impact", "base", "output_yr_base.RData"))
  save(output_season_catchup_base, file = here("outputs", "nirsevimab", "impact", "base",  "output_season_catchup_base.RData"))
  save(output_season_catchup_nip_base, file = here("outputs", "nirsevimab", "impact", "base", "output_season_catchup_nip_base.RData"))
}

run_model_base_direct <- function(data_list) {
  seeds <- data_list$seeds
  post <- data_list$post
  vac_par_info <- data_list$vac_par_info
  S <- data_list$S
  vac_par_info$direct <- TRUE 

  output_default_none_d <- run_sample_custom(seeds[1:S], data_list$none, vac_par_info, 0, post)
  output_default_base_d <- run_sample_custom(seeds[1:S], data_list$pal, vac_par_info, 0, post)
  output_season_vhr_base_d <- run_sample_custom(seeds[1:S], data_list$vhr_seasonal, vac_par_info, 0, post)
  output_season_base_d <- run_sample_custom(seeds[1:S], data_list$seasonal, vac_par_info, 0, post)
  output_yr_base_d <- run_sample_custom(seeds[1:S], data_list$year_round, vac_par_info, 0, post)
  output_season_catchup_base_d <- run_sample_custom(seeds[1:S], data_list$seasonal_catchup, vac_par_info, 0, post)
  output_season_catchup_nip_base_d <- run_sample_custom(seeds[1:S], data_list$seasonal_catchup_nip, vac_par_info, 0, post)

  save(output_default_none_d, file = here("outputs", "nirsevimab", "impact", "base", "none_d.RData"))
  save(output_default_base_d, file = here("outputs", "nirsevimab", "impact", "base", "status_quo_base_d.RData"))
  save(output_season_vhr_base_d, file = here("outputs", "nirsevimab", "impact", "base", "output_season_vhr_base_d.RData"))
  save(output_season_base_d, file = here("outputs", "nirsevimab", "impact", "base", "output_season_base_d.RData"))
  save(output_yr_base_d, file = here("outputs", "nirsevimab", "impact", "base", "output_yr_base_d.RData"))
  save(output_season_catchup_base_d, file = here("outputs", "nirsevimab", "impact", "base",  "output_season_catchup_base_d.RData"))
  save(output_season_catchup_nip_base_d, file = here("outputs", "nirsevimab", "impact", "base", "output_season_catchup_nip_base_d.RData"))
}



run_icer_base <- function(icer_threshold, filename) {
  source("R/cea.R")

  load(file = here("outputs", "nirsevimab", "impact", "base", "none.RData")) # output_default_none
  load(file = here("outputs", "nirsevimab", "impact", "base", "status_quo_base.RData")) # output_default_base
  load(file = here("outputs", "nirsevimab", "impact", "base", "output_season_vhr_base.RData")) # output_season_vhr_base
  load(file = here("outputs", "nirsevimab", "impact", "base", "output_season_base.RData")) # output_season_base
  load(file = here("outputs", "nirsevimab", "impact", "base", "output_yr_base.RData")) # output_yr_base
  load(file = here("outputs", "nirsevimab", "impact", "base", "output_season_catchup_base.RData")) # output_season_catchup_base
  load(file = here("outputs", "nirsevimab", "impact", "base", "output_season_catchup_nip_base.RData")) # output_season_catchup_nip_base

  ## Look at estimating the dominance from these programmes

  list_outputs <- list(output_default_base, output_season_vhr_base, output_season_base,
    output_season_catchup_base, output_season_catchup_nip_base)
  list_names <- list("vhr", "seasonal",
    "seasonal_and_catchup", "seasonal_and_catchup_nip")

  ppd_base <- cal_impact_standard(list_outputs, list_names, icer_threshold)
  save(ppd_base, file = here::here("outputs", "nirsevimab", "icer", paste0(filename, ".RData")))
}


run_icer_mabs_dur <- function(datalist, threshold) {
  source("R/cea.R")

  load(file = here("outputs", "nirsevimab", "impact", "base", "status_quo_base.RData")) # output_default_base

  seeds <- datalist$seeds
  post <- datalist$post
  S <- datalist$S

  vac_par_info_d_250 <- list(om_mab = 1 / 250, direct = FALSE, xi_boost = 1)

  output_season_vhr_d_250 <- run_sample_custom(seeds[1:S], datalist$vhr_seasonal, vac_par_info_d_250, 0, post)
  output_season_d_250 <- run_sample_custom(seeds[1:S], datalist$seasonal, vac_par_info_d_250, 0, post)
  output_yr_d_250 <- run_sample_custom(seeds[1:S], datalist$year_round, vac_par_info_d_250, 0, post)
  output_season_catchup_d_250 <- run_sample_custom(seeds[1:S], datalist$seasonal_catchup, vac_par_info_d_250, 0, post)
  output_season_catchup_nip_d_250 <- run_sample_custom(seeds[1:S], datalist$seasonal_catchup_nip, vac_par_info_d_250, 0, post)

  save(output_season_vhr_d_250, file = here("outputs", "nirsevimab", "impact", "mabs_dur", "output_season_vhr_d_250.RData"))
  save(output_season_d_250, file = here("outputs", "nirsevimab", "impact", "mabs_dur", "output_season_d_250.RData"))
  save(output_yr_d_250, file = here("outputs", "nirsevimab", "impact", "mabs_dur", "output_yr_d_250.RData"))
  save(output_season_catchup_d_250, file = here("outputs", "nirsevimab", "impact", "mabs_dur", "output_season_catchup_d_250.RData"))
  save(output_season_catchup_nip_d_250, file = here("outputs", "nirsevimab", "impact", "mabs_dur", "output_season_catchup_nip_d_250.RData"))

  list_outputs <- list(output_default_base, output_season_vhr_d_250, output_season_d_250,
    output_season_catchup_d_250, output_season_catchup_nip_d_250)
  list_names <-  list("vhr", "seasonal",
    "seasonal_and_catchup", "seasonal_and_catchup_nip")
  ppd_d_250 <- cal_impact_standard(list_outputs, list_names, threshold)
  save(ppd_d_250, file = here::here("outputs", "nirsevimab", "icer", "d_250_icer.RData"))

  vac_par_info_d_360 <- list(om_mab = 1 / 360, direct = FALSE, xi_boost = 1)

  output_season_vhr_d_360 <- run_sample_custom(seeds[1:S],  datalist$vhr_seasonal, vac_par_info_d_360, 0, post)
  output_season_d_360 <- run_sample_custom(seeds[1:S],  datalist$seasonal, vac_par_info_d_360, 0, post)
  output_yr_d_360 <- run_sample_custom(seeds[1:S], datalist$year_round, vac_par_info_d_360, 0, post)
  output_season_catchup_d_360 <- run_sample_custom(seeds[1:S], datalist$seasonal_catchup, vac_par_info_d_360, 0, post)
  output_season_catchup_nip_d_360 <- run_sample_custom(seeds[1:S],  datalist$seasonal_catchup_nip, vac_par_info_d_360, 0, post)

  save(output_season_vhr_d_360, file = here("outputs", "nirsevimab", "impact", "mabs_dur","output_season_vhr_d_360.RData"))
  save(output_season_d_360, file = here("outputs", "nirsevimab", "impact", "mabs_dur","output_season_d_360.RData"))
  save(output_yr_d_360, file = here("outputs", "nirsevimab", "impact", "mabs_dur","output_yr_d_360.RData"))
  save(output_season_catchup_d_360, file = here("outputs", "nirsevimab", "impact", "mabs_dur", "output_season_catchup_d_360.RData"))
  save(output_season_catchup_nip_d_360, file = here("outputs", "nirsevimab", "impact", "mabs_dur","output_season_catchup_nip_d_360.RData"))


  list_outputs <- list(output_default_base, output_season_vhr_d_360, output_season_d_360,
    output_season_catchup_d_360, output_season_catchup_nip_d_360)
  list_names <-  list("vhr", "seasonal",
    "seasonal_and_catchup", "seasonal_and_catchup_nip")
  ppd_d_360 <- cal_impact_standard(list_outputs, list_names, threshold)
  save(ppd_d_360, file = here::here("outputs", "nirsevimab", "icer", "d_360_icer.RData"))

}

run_icer_coverage <- function(datalist, threshold, coverage) {

    seeds <- datalist$seeds
    post <- datalist$post
    S <- datalist$S
    vac_par_info <- datalist$vac_par_info

    pal_eff <- datalist$pal_eff
    mab_eff <- datalist$mab_eff

    ind_pal <- c(rep(1, 9), rep(0, 16)) # VHR (<8 months)
    ind_mabs <- c(rep(1, 1), rep(0, 24)) # Birth only (0 months only)
    G_plus <- c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
    source("R/cea.R")

    load(file = here("outputs", "nirsevimab", "impact", "base", "status_quo_base.RData")) # output_default_base

    make_vac_program_info_mabs_vhr_seasonal <- function(seed) {
        list(
            mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 15 * 7, t_end = 32 * 7, eff = mab_eff[seed], cov = 0.9)
        )
    }

    make_vac_program_info_mabs_seasonal_low_cov <- function(seed) {
        list(
            mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
            mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = coverage),
            mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = coverage)
        )
    }

    make_vac_program_info_mabs_year_round_low_cov <- function(seed) {
        list(
            mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 0 * 7, t_end = 52 * 7, eff = mab_eff[seed], cov = 0.9),
            mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 0, t_end = 52 * 7, eff = mab_eff[seed], cov = coverage),
            mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 0, t_end = 52 * 7, eff = mab_eff[seed], cov = coverage)
        )
    }

    make_vac_program_info_mabs_seasonal_catchup_low_cov <- function(seed) {
        list(
            mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
            mAB_HR =  list(id = TRUE, catchup = TRUE, catchupnip = FALSE, age_id_catchup = G_plus, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = coverage),
            mAB_LR =  list(id = TRUE, catchup = TRUE, catchupnip = FALSE, age_id_catchup = G_plus, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = coverage)
        )
    }


    make_vac_program_info_mabs_seasonal_catchup_nip_low_cov <- function(seed) {
        list(
            mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
            mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = TRUE, age_id = ind_mabs, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = coverage),
            mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = TRUE, age_id = ind_mabs, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = coverage)
        )
    }

    output_season_vhr_low_cov <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_vhr_seasonal, vac_par_info, 0, post)
    output_season_low_cov <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_seasonal_low_cov, vac_par_info, 0, post)
    output_yr_low_cov <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_year_round_low_cov, vac_par_info, 0, post)
    output_season_catchup_low_cov <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_seasonal_catchup_low_cov, vac_par_info, 0, post)
    output_season_catchup_nip_low_cov <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_seasonal_catchup_nip_low_cov, vac_par_info, 0, post)

    save(output_season_vhr_low_cov, file = here("outputs", "nirsevimab", "impact", "coverage", "output_season_vhr_low_cov.RData"))
    save(output_season_low_cov, file = here("outputs", "nirsevimab", "impact","coverage",  "output_season_low_cov.RData"))
    save(output_yr_low_cov, file = here("outputs", "nirsevimab", "impact", "coverage", "output_yr_low_cov.RData"))
    save(output_season_catchup_low_cov, file = here("outputs", "nirsevimab", "impact", "coverage", "output_season_catchup_low_cov.RData"))
    save(output_season_catchup_nip_low_cov, file = here("outputs", "nirsevimab", "impact", "coverage", "output_season_catchup_nip_low_cov.RData"))

    list_outputs <- list(output_default_base, output_season_vhr_low_cov, output_season_low_cov,
        output_season_catchup_low_cov, output_season_catchup_nip_low_cov)
    list_names <-  list("vhr", "seasonal",
    "seasonal_and_catchup", "seasonal_and_catchup_nip")
    ppd_low_cov <- cal_impact_standard(list_outputs, list_names, threshold)
    save(ppd_low_cov, file = here::here("outputs", "nirsevimab", "icer", "low_cov_icer.RData"))
}

run_icer_admin_year <- function(datalist, threshold, ind_mabs) {

  seeds <- datalist$seeds
  post <- datalist$post
  S <- datalist$S
  vac_par_info <- datalist$vac_par_info

  pal_eff <- datalist$pal_eff
  mab_eff <- datalist$mab_eff

  ind_pal <- c(rep(1, 9), rep(0, 16)) # VHR (<8 months)
  G_plus <- c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

  source("R/cea.R")

  load(file = here("outputs", "nirsevimab", "impact", "base", "status_quo_base.RData")) # output_default_base


  make_vac_program_info_mabs_vhr_seasonal <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 15 * 7, t_end = 32 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }

  make_vac_program_info_mabs_seasonal_2mo <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9)
        )
  }

  make_vac_program_info_mabs_year_round_2mo <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 0 * 7, t_end = 52 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 0, t_end = 52 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = FALSE, age_id = ind_mabs, t_start = 0, t_end = 52 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }

  make_vac_program_info_mabs_seasonal_catchup_2mo <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_HR =  list(id = TRUE, catchup = TRUE, catchupnip = FALSE, age_id_catchup = G_plus, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_LR =  list(id = TRUE, catchup = TRUE, catchupnip = FALSE, age_id_catchup = G_plus, age_id = ind_mabs, t_start = 13 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }


  make_vac_program_info_mabs_seasonal_catchup_nip_2mo <- function(seed) {
    list(
        mAB_VHR = list(id = TRUE, age_id = ind_pal, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_HR =  list(id = TRUE, catchup = FALSE, catchupnip = TRUE, age_id = ind_mabs, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9),
        mAB_LR =  list(id = TRUE, catchup = FALSE, catchupnip = TRUE, age_id = ind_mabs, t_start = 9 * 7, t_end = 34 * 7, eff = mab_eff[seed], cov = 0.9)
      )
  }


  output_season_vhr_2mo <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_vhr_seasonal, vac_par_info, 0, post)
  output_season_2mo <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_seasonal_2mo, vac_par_info, 0, post)
  output_yr_2mo <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_year_round_2mo, vac_par_info, 0, post)
  output_season_catchup_2mo <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_seasonal_catchup_2mo, vac_par_info, 0, post)
  output_season_catchup_nip_2mo <- run_sample_custom(seeds[1:S], make_vac_program_info_mabs_seasonal_catchup_nip_2mo, vac_par_info, 0, post)

  save(output_season_vhr_2mo, file = here("outputs", "nirsevimab", "impact", "admin_age", "output_season_vhr_2mov.RData"))
  save(output_season_2mo, file = here("outputs", "nirsevimab", "impact","admin_age",  "output_season_2mo.RData"))
  save(output_yr_2mo, file = here("outputs", "nirsevimab", "impact", "admin_age", "output_yr_2mo.RData"))
  save(output_season_catchup_2mo, file = here("outputs", "nirsevimab", "impact", "admin_age", "output_season_catchup_2mo.RData"))
  save(output_season_catchup_nip_2mo, file = here("outputs", "nirsevimab", "impact", "admin_age", "output_season_catchup_nip_2mo.RData"))

  list_outputs <- list(output_default_base, output_season_vhr_2mo, output_season_2mo,
    output_season_catchup_2mo, output_season_catchup_nip_2mo)
  list_names <-  list("vhr", "seasonal",
    "seasonal_and_catchup", "seasonal_and_catchup_nip")
  ppd_2mo <- cal_impact_standard(list_outputs, list_names, 20000)
  save(ppd_2mo, file = here::here("outputs", "nirsevimab", "icer", "2mo_icer.RData"))

}