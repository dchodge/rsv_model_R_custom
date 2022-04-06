
cal_impact_standard <- function(list_outputs, list_names, threshold) {
  icer_output_A <- calc_impact(list_outputs[[1]], list_outputs[[2]],  threshold, "base_mab_p1_vhr_icer", pal_base = TRUE) %>%
    mutate(names = list_names[[1]])
  icer_output_B <- calc_impact(list_outputs[[2]], list_outputs[[3]], threshold, "base_mab_p1_icer", pal_base = FALSE) %>%
    mutate(names = list_names[[2]])
  icer_output_C <- calc_impact(list_outputs[[3]], list_outputs[[4]],  threshold, "base_mab_p1_catchup_icer", pal_base = FALSE) %>%
    mutate(names = list_names[[3]])
  icer_output_D <- calc_impact(list_outputs[[4]], list_outputs[[5]],  threshold, "base_mab_p1_catchup_nip_icer", pal_base = FALSE) %>%
    mutate(names = list_names[[4]])
  icer_output <- bind_rows(icer_output_A, icer_output_B, icer_output_C, icer_output_D)
  icer_output
}

calc_impact <- function(base, programme, threshold, filename, pal_base) {
    base_sum <- base$outcomes_week_age %>% 
        group_by(outcome, seed, age_group) %>%
        summarise(incidence_base = sum(incidence))

    programme_sum <- programme$outcomes_week_age %>% 
        group_by(outcome, seed, age_group) %>%
        summarise(incidence_programme = sum(incidence))

    df_averted <- base_sum %>% left_join(programme_sum, by = c("outcome", "seed", "age_group")) %>%
        mutate(cases_averted = incidence_base - incidence_programme)

    # Total QALY loss needed
    df_qaly <- left_join(
        base$QALY %>% filter(type == "discounted", metric == "total"),
        programme$QALY %>% filter(type == "discounted", metric == "total"),
        by = c("seed", "type", "metric")
    ) %>% rename(base_qaly = value.x, inter_qaly = value.y) %>% unique

    # Need total cost different ignoring number of doses given 
    df_cost <- left_join(
        base$cost %>% filter(type == "discounted", metric == "total"),
        programme$cost %>% filter(type == "discounted", metric == "total"),
        by = c("seed", "type", "metric")
    ) %>% rename(base_cost = value.x, inter_cost = value.y) %>% unique

    if (pal_base) {
        pal_cost <- 4035.50 * sum(base$vac_cal[, 1] * exp(-0.035 / 52 * (1:521)))
        tot_doses <- sum(programme$vac_cal[, 2] * exp(-0.035 / 52 * (1:521)))
        df_icer <- left_join(df_qaly, df_cost, by = c("seed", "type", "metric")) %>%
            mutate(ppd = (threshold * (base_qaly - inter_qaly) - (inter_cost - base_cost) + pal_cost) / tot_doses)
    } else{
        tot_doses <- sum(programme$vac_cal[, 2] * exp(-0.035 / 52 * (1:521))) - sum(base$vac_cal[, 1] * exp(-0.035 / 52 * (1:521)))
        df_icer <- left_join(df_qaly, df_cost, by = c("seed", "type", "metric")) %>%
            mutate(ppd = (threshold * (base_qaly - inter_qaly) - (inter_cost - base_cost)) / tot_doses)

    }
    df_icer
}