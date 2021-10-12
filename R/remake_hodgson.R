#' Function which runs each of the 15 intervention programmes outlined in Hodgson et al. 2020
#' 
#' @param prog_no Programme number (1 to 15)
#' @param seed_vals A list of the seed values
#' @param all_eff A list of efficacy samples which are selected using the seed value
#' @param type Serovention times or dosing?
#' @return A list of two dataframes, the first shows the annual incidence for health outcomes, the second shows the economic metrics.
create_calendars_hodgson <- function(all_eff) {
    # age_ids
    G_base <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    G_0mo = c(1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    G_plus = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    G_2mo = c(0,0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    G_pal = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    G_2_4 =  c(0,0,0,0,0,0,0,0,0,0,0,0,0, 1.0, 1.0, 1.0,           0,0,0,0,0,0,0,0,0)
    G_5_9 = c(0,0,0,0,0,0,0,0,0,0,0,0,0, 0.0, 0.0, 0.0, 1.0,        0,0,0,0,0,0,0,0)
    G_5_14 = c(0,0,0,0,0,0,0,0,0,0,0,0,0, 0.0, 0.0, 0.0, 1.0, 1.0, 0,  0,  0,  0,0,0,0)
    G_par =  c(0,0,0,0,0,0,0,0,0,0,0,0,0, 0,   0,   0,   0,   0,   1.0,1.0,1.0,0,0,0,0)
    G_65_ = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 1.0, 1.0)
    G_75_ = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 1.0)

    # Efficacy sampels
    pal_eff <- all_eff[["pal_eff"]]
    mab_eff <- all_eff[["mab_eff"]]
    lav_eff <- all_eff[["lav_eff"]]
    mat_eff <- all_eff[["mat_eff"]]

    # uptakes rate
    up_week_o65 = c(0., 0.046, 0.177, 0.311, 0.444, 0.64, 0.737, 0.802, 0.856, 0.894, 0.920, 0.941, 0.953, 0.967, 0.978, 0.980, 0.986, 0.988, 0.988, 0.996, 1.)
    up_week_u65 = c(0., 0.031, 0.075, 0.16, 0.274, 0.387, 0.491, 0.579, 0.664, 0.726, 0.783, 0.830, 0.868, 0.906, 0.931, 0.953, 0.953, 0.965, 0.981, 0.987, 1.)
    up_week_2t3 = c(0., 0.000, 0.000, 0.00, 0.040, 0.120, 0.227, 0.369, 0.511, 0.629, 0.719, 0.809, 0.868, 0.910, 0.941, 0.965, 0.969, 0.976, 0.990, 0.993, 1.)
    up_week_preg = c(0., 0.041, 0.095, 0.192, 0.32, 0.451, 0.568, 0.664, 0.734, 0.794, 0.843, 0.881, 0.917, 0.934, 0.944, 0.967, 0.974, 0.973, 0.983, 0.990, 1.)
    i_shift <- 1

    # Defining the parameters for the 15 programmes
    P0_none_programme <- function(seed) {
        list( 
        )
    }
    P1_pal_programme <- function(seed) {
        list( 
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9)
        )
    }
    P2_mab_vhr_programme <- function(seed) {
        list( 
            mAB_VHR = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = mab_eff[seed], cov = 0.9)
        )
    }
    P3_mab_hr_programme <- function(seed) {
        list( 
            mAB_VHR = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = mab_eff[seed], cov = 0.9),
            mAB_HR = list(id = TRUE, age_id = G_0mo, catchup = FALSE, t_start = 12*7 + i_shift, t_end = (12 + 21)*7 + i_shift, eff = mab_eff[seed], cov = 0.9)
        )
    }
    P4_mab_hr_plus_programme <- function(seed) {
        list(
            mAB_VHR = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = mab_eff[seed], cov = 0.9),
            mAB_HR = list(id = TRUE, age_id = G_0mo, catchup = TRUE, age_id_catchup = G_plus, t_start = 12*7 + i_shift, t_end = (12 + 21)*7 + i_shift, eff = mab_eff[seed], cov = 0.9)
        )
    }
    P5_mab_all_programme <- function(seed) {
        list( 
            mAB_VHR = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = mab_eff[seed], cov = 0.9),
            mAB_HR = list(id = TRUE, age_id = G_0mo, catchup = FALSE, t_start = 12*7 + i_shift, t_end = (12 + 21)*7 + i_shift, eff = mab_eff[seed], cov = 0.9),
            mAB_LR = list(id = TRUE, age_id = G_0mo, catchup = FALSE, t_start = 12*7 + i_shift, t_end = (12 + 21)*7 + i_shift, eff = mab_eff[seed], cov = 0.9)
        )
    }
    P6_mab_all_plus_programme <- function(seed) {
        list(
            mAB_VHR = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = mab_eff[seed], cov = 0.9),
            mAB_HR = list(id = TRUE, age_id = G_0mo, catchup = TRUE, age_id_catchup = G_plus, t_start = 12*7 + i_shift, t_end = (12 + 21)*7 + i_shift, eff = mab_eff[seed], cov = 0.9),
            mAB_LR = list(id = TRUE, age_id = G_0mo, catchup = TRUE, age_id_catchup = G_plus, t_start = 12*7 + i_shift, t_end = (12 + 21)*7 + i_shift, eff = mab_eff[seed], cov = 0.9)
        )
    }
    P7_mat_season_programme <- function(seed) {
        list(
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            mat = list(id = TRUE, age_id = G_par, t_start = 0 + i_shift, t_end = 21*7 + i_shift, eff_inf = mat_eff[seed], eff_mat = lav_eff[seed], cov = 0.6)
        )
    }
    P8_mat_yr_programme <- function(seed) {
        list( 
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            mat = list(id = TRUE, age_id = G_par, t_start = 0 + i_shift, t_end = 52*7 + i_shift, eff_inf = mat_eff[seed], eff_mat = lav_eff[seed], cov = 0.6)
        )
    }
    P9_LAV_2mo_season_programme <- function(seed) {
        list(
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            LAV_HR = list(id = TRUE, age_id = G_2mo, uptake_type = 2, t_start = 4*7 + i_shift, t_end = (4 + 21)*7 + i_shift, eff = lav_eff[seed], cov = 0.9),
            LAV_LR = list(id = TRUE, age_id = G_2mo, uptake_type = 2, t_start = 4*7 + i_shift, t_end = (4 + 21)*7 + i_shift, eff = lav_eff[seed], cov = 0.9)
        )
    }
    P10_LAV_2mo_yr_programme <- function(seed) {
        list(
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            LAV_HR = list(id = TRUE, age_id = G_2mo, uptake_type = 2, t_start = 0 + i_shift, t_end = 52*7 + i_shift, eff = lav_eff[seed], cov = 0.9),
            LAV_LR = list(id = TRUE, age_id = G_2mo, uptake_type = 2, t_start = 0 + i_shift, t_end = 52*7 + i_shift, eff = lav_eff[seed], cov = 0.9)
        )
    }
    up_week_2t3 <- get_daily_uptake(up_week_2t3, 16)
    P11_LAV_2_4_programme <- function(seed) {
    list(
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            LAV_LR = list(id = TRUE, age_id = G_2_4, uptake_type = 1, eff = lav_eff[seed], cov = 0.45,  uptake_daily = up_week_2t3)
        )
    }
    rate_u65 <- get_daily_uptake(up_week_u65, 16)
    P12_LAV_5_9_programme <- function(seed) {
        list(
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            LAV_LR = list(id = TRUE, age_id = G_5_9, uptake_type = 1, eff = lav_eff[seed], cov = 0.6,  uptake_daily = rate_u65)
        )
    }
    P13_LAV_5_14_programme <- function(seed) {
        list(
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            LAV_LR = list(id = TRUE, age_id = G_5_14, uptake_type = 1, eff = lav_eff[seed], cov = 0.6,  uptake_daily = rate_u65)
        )
    }
    rate_o65 <- get_daily_uptake(up_week_o65, 16)
    P14_LAV_65_programme <- function(seed) {
        list(
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            LAV_LR = list(id = TRUE, age_id = G_65_, uptake_type = 1, eff = lav_eff[seed], cov = 0.7,  uptake_daily = rate_o65)
        )
    }
    P15_LAV_75_programme <- function(seed) {
        list(
            pal = list(id = TRUE, age_id = G_pal, t_start = 15*7 + i_shift, t_end = 36*7 + i_shift, eff = pal_eff[seed], cov = 0.9),
            LAV_LR = list(id = TRUE, age_id = G_75_, uptake_type = 1, eff = lav_eff[seed], cov = 0.7,  uptake_daily = rate_o65)
        )
    }

    all_programmes <- list(
        P0_none_programme = P0_none_programme,
        P1_pal_programme = P1_pal_programme,
        P2_mab_vhr_programme = P2_mab_vhr_programme,
        P3_mab_hr_programme = P3_mab_hr_programme,
        P4_mab_hr_plus_programme = P4_mab_hr_plus_programme,
        P5_mab_all_programme = P5_mab_all_programme,
        P6_mab_all_plus_programme = P6_mab_all_plus_programme,
        P7_mat_season_programme = P7_mat_season_programme,
        P8_mat_yr_programme = P8_mat_yr_programme,
        P9_LAV_2mo_season_programme = P9_LAV_2mo_season_programme,
        P10_LAV_2mo_yr_programme = P10_LAV_2mo_yr_programme,
        P11_LAV_2_4_programme = P11_LAV_2_4_programme,
        P12_LAV_5_9_programme = P12_LAV_5_9_programme,
        P13_LAV_5_14_programme = P13_LAV_5_14_programme,
        P14_LAV_65_programme = P14_LAV_65_programme,
        P15_LAV_75_programme = P15_LAV_75_programme
    )

    cov_c <- rep(0, 16)
    cov_c[8] <- 0.6
    cov_c[9] <- 0.6
    all_info <- list(programmes = all_programmes, mat_cov = cov_c)
    save(all_info, file = here::here("data", "inter_model_input", "hodgson_programmes.RData"))

    # Select correct programme
   # output <- switch(prog_no,
   #     create_calendar(P1_pal_programme, "week", 0),
    #    create_calendar(P2_mab_vhr_programme, "week", 0),
     #   create_calendar(P3_mab_hr_programme, "week", 0),
     #   create_calendar(P4_mab_hr_plus_programme, "week", 0),
     #   create_calendar(P5_mab_all_programme, "week", 0),
     #   create_calendar(P6_mab_all_plus_programme, "week", 0),
     #   create_calendar(P7_mat_season_programme, "week", 0.6),
     #   create_calendar(P8_mat_yr_programme, "week", 0.6),
     #   create_calendar(P9_LAV_2mo_season_programme, "week", 0),
     #   create_calendar(P10_LAV_2mo_yr_programme, "week", 0),
     #   create_calendar(P11_LAV_2_4_programme, "week", 0),
     #   create_calendar(P12_LAV_5_9_programme, "week", 0),
     #   create_calendar(P13_LAV_5_14_programme, "week", 0),
     #   create_calendar(P14_LAV_65_programme, "week", 0),
     #   create_calendar(P15_LAV_75_programme, "week", 0)
    #)
    

   # names <- c("pal", "mAB_VHR", "mAB_HR", "mAB_LR", "LAV_HR", "LAV_LR", "mat_LR")
   # out <- vector(mode = "list", length = length(output[[1]])) %>% setNames(names)
   # j <- 1
   # if (type == "cal") {
   #     for (i in output[[1]]) {
   #         if (is.matrix(i))
   #             out[[j]] <- i
   #         else
   #             out[[j]] <- i[[2]]
   #         j <- j + 1
   #     }
   # } else {
   #     j <- 1
   #     for (i in output[[1]]) {
    #            out[[j]] <- i
   ##         if (is.matrix(i))
    #        else 
    #            out[[j]] <- i[[1]]
    #        j <- j + 1
    #    }
    #}
    #list(out = out, cov_c = output[[2]])*/
}