library(brms)
library(tidybayes)
library(tidyverse)
library(readxl)
library(DescTools)
library(irr)

# Safe extraction helper functions
safe_extract <- function(obj, field) {
  if (!is.null(obj) && field %in% names(obj)) {
    return(obj[[field]])
  } else {
    return(NA_real_)
  }
}

safe_CCC <- function(x, y) {
  result <- tryCatch(DescTools::CCC(x, y), error = function(e) NULL)
  if (!is.null(result) && is.list(result) && !is.null(result$rho.c)) {
    return(list(
      est = result$rho.c[["est"]],
      lwr.ci = result$rho.c[["lwr.ci"]],
      upr.ci = result$rho.c[["upr.ci"]]
    ))
  } else {
    return(list(est = NA_real_, lwr.ci = NA_real_, upr.ci = NA_real_))
  }
}

calculate_overlap <- function(df) {
  d1 <- df %>% filter(Day == "D1") %>% pull(.epred)
  d2 <- df %>% filter(Day == "D2") %>% pull(.epred)
  
  dens1 <- density(d1)
  dens2 <- density(d2)
  
  x_common <- seq(
    from = max(min(dens1$x), min(dens2$x)),
    to = min(max(dens1$x), max(dens2$x)),
    length.out = 512
  )
  
  y1 <- approx(dens1$x, dens1$y, xout = x_common)$y
  y2 <- approx(dens2$x, dens2$y, xout = x_common)$y
  
  overlap <- sum(pmin(y1, y2)) * (x_common[2] - x_common[1])
  return(overlap)
}

# Read data
imtp <- read_xlsx("/Users/tara/Desktop/PhD/Study2/Raw data file_1.00.xlsx")
imtp <- imtp %>% 
  mutate(across(c(Method, Participant, Day), as.factor),
         Value = as.numeric(Value))

# Model settings
metrics <- c("Net Peak Vertical Force [N]", "Net Force at 100ms [N]")
methods <- c("TRAD", "SHORT")
families <- list(gaussian = gaussian(), skew_normal = skew_normal())

results_list <- list()

for (met in metrics) {
  for (method in methods) {
    for (fam_name in names(families)) {
      cat("Running:", met, method, fam_name, "\n")
      
      fam <- families[[fam_name]]
      
      df <- imtp %>% 
        filter(Metric == met, Trial %in% paste("Trial", 1:5), Method == method)
      
      formula <- bf(Value ~ (1|Participant), sigma ~ (1|Participant))
      
      m1 <- brm(formula, data = df %>% filter(Day == "D1"), family = fam, chains = 4, cores = 4, iter = 2000)
      m2 <- brm(formula, data = df %>% filter(Day == "D2"), family = fam, chains = 4, cores = 4, iter = 2000)
      
      epred_all <- bind_rows(
        df %>% filter(Day == "D1") %>% modelr::data_grid(Participant) %>% add_epred_draws(m1) %>% mutate(Day = "D1"),
        df %>% filter(Day == "D2") %>% modelr::data_grid(Participant) %>% add_epred_draws(m2) %>% mutate(Day = "D2")
      )
      
      epred_split <- epred_all %>% split(.$Participant)
      overlap_results <- map_dfr(names(epred_split), function(p) {
        tibble(Participant = p, Overlap = calculate_overlap(epred_split[[p]]))
      })
      
      ranef_m1 <- brms::ranef(m1)$Participant
      ranef_m2 <- brms::ranef(m2)$Participant
      
      df_ranef <- data.frame(
        ID = rownames(ranef_m1[, , "Intercept"]),
        Day1 = ranef_m1[, , "Intercept"][, "Estimate"],
        Day2 = ranef_m2[, , "Intercept"][, "Estimate"],
        Day1_sigma = ranef_m1[, , "sigma_Intercept"][, "Estimate"],
        Day2_sigma = ranef_m2[, , "sigma_Intercept"][, "Estimate"]
      )
      
      icc_result <- irr::icc(df_ranef[, c("Day1", "Day2")], model = "twoway", type = "agreement", unit = "single")
      icc_result_sigma <- irr::icc(df_ranef[, c("Day1_sigma", "Day2_sigma")], model = "twoway", type = "agreement", unit = "single")
      
      ccc_ranef_list <- safe_CCC(df_ranef$Day1, df_ranef$Day2)
      ccc_sigma_list <- safe_CCC(df_ranef$Day1_sigma, df_ranef$Day2_sigma)
      
      coef_m1 <- coef(m1)$Participant[, , "Intercept"]
      coef_m2 <- coef(m2)$Participant[, , "Intercept"]
      df_coefs <- data.frame(Day1 = coef_m1[, "Estimate"], Day2 = coef_m2[, "Estimate"])
      
      icc_mean <- irr::icc(df_coefs, model = "twoway", type = "agreement", unit = "single")
      ccc_mean_list <- safe_CCC(df_coefs$Day1, df_coefs$Day2)
      
      posterior_summary <- epred_all %>%
        group_by(Participant, Day) %>%
        summarise(median = median(.epred), .groups = "drop")
      
      posterior_medians <- posterior_summary %>%
        pivot_wider(names_from = Day, values_from = median) %>%
        rename(Day1 = D1, Day2 = D2)
      
      icc_median <- irr::icc(posterior_medians[, c("Day1", "Day2")], model = "twoway", type = "agreement", unit = "single")
      ccc_median_list <- safe_CCC(posterior_medians$Day1, posterior_medians$Day2)
      
      results_list[[paste(met, method, fam_name, sep = "_")]] <- tibble(
        metric = met,
        method = method,
        family = fam_name,
        mean_overlap = mean(overlap_results$Overlap),
        
        icc_ranef = icc_result$value,
        icc_ranef_lwr = safe_extract(icc_result, "lbound"),
        icc_ranef_upr = safe_extract(icc_result, "ubound"),
        ccc_ranef = safe_extract(ccc_ranef_list, "est"),
        ccc_ranef_lwr = safe_extract(ccc_ranef_list, "lwr.ci"),
        ccc_ranef_upr = safe_extract(ccc_ranef_list, "upr.ci"),
        
        icc_sigma = icc_result_sigma$value,
        icc_sigma_lwr = safe_extract(icc_result_sigma, "lbound"),
        icc_sigma_upr = safe_extract(icc_result_sigma, "ubound"),
        ccc_sigma = safe_extract(ccc_sigma_list, "est"),
        ccc_sigma_lwr = safe_extract(ccc_sigma_list, "lwr.ci"),
        ccc_sigma_upr = safe_extract(ccc_sigma_list, "upr.ci"),
        
        icc_mean_est = safe_extract(icc_mean, "value"),
        icc_mean_lwr = safe_extract(icc_mean, "lbound"),
        icc_mean_upr = safe_extract(icc_mean, "ubound"),
        ccc_mean_est = safe_extract(ccc_mean_list, "est"),
        ccc_mean_lwr = safe_extract(ccc_mean_list, "lwr.ci"),
        ccc_mean_upr = safe_extract(ccc_mean_list, "upr.ci"),
        
        icc_median_est = safe_extract(icc_median, "value"),
        icc_median_lwr = safe_extract(icc_median, "lbound"),
        icc_median_upr = safe_extract(icc_median, "ubound"),
        ccc_median_est = safe_extract(ccc_median_list, "est"),
        ccc_median_lwr = safe_extract(ccc_median_list, "lwr.ci"),
        ccc_median_upr = safe_extract(ccc_median_list, "upr.ci")
      )
    }
  }
}

final_results <- bind_rows(results_list)
write_csv(final_results, "imtp_model_fit_metrics.csv")

# Traditional metrics

traditional_results <- list()

for (met in metrics) {
  for (meth in methods) {
    df <- imtp %>%
      filter(Metric == met, Method == meth, Trial %in% paste("Trial", 1:5))
    
    summary_df <- df %>%
      group_by(Participant, Day) %>%
      summarise(
        mean_value = mean(Value, na.rm = TRUE),
        max_value = max(Value, na.rm = TRUE),
        min_value = min(Value, na.rm = TRUE), 
        trimmed_mean = mean(sort(Value, decreasing = TRUE)[1:3], na.rm = TRUE),
        median_value = median(Value, na.rm = TRUE),
        .groups = "drop"
      )
    
    wide_df <- summary_df %>%
      pivot_wider(
        names_from = Day,
        values_from = c(mean_value, max_value, min_value, trimmed_mean, median_value)
      ) %>%
      drop_na()
    
    # ICCs
    icc_mean <- icc(wide_df[, c("mean_value_D1", "mean_value_D2")],
                    model = "twoway", type = "agreement", unit = "single")
    icc_max <- icc(wide_df[, c("max_value_D1", "max_value_D2")],
                   model = "twoway", type = "agreement", unit = "single")
    icc_min <- icc(wide_df[, c("min_value_D1", "min_value_D2")],  
                   model = "twoway", type = "agreement", unit = "single")
    icc_trimmed <- icc(wide_df[, c("trimmed_mean_D1", "trimmed_mean_D2")],
                       model = "twoway", type = "agreement", unit = "single")
    icc_median <- icc(wide_df[, c("median_value_D1", "median_value_D2")],
                      model = "twoway", type = "agreement", unit = "single")
    
    # CCCs
    ccc_mean <- safe_CCC(wide_df$mean_value_D1, wide_df$mean_value_D2)
    ccc_max <- safe_CCC(wide_df$max_value_D1, wide_df$max_value_D2)
    ccc_min <- safe_CCC(wide_df$min_value_D1, wide_df$min_value_D2)  
    ccc_trimmed <- safe_CCC(wide_df$trimmed_mean_D1, wide_df$trimmed_mean_D2)
    ccc_median <- safe_CCC(wide_df$median_value_D1, wide_df$median_value_D2)
    
    traditional_results[[paste(met, meth, sep = "_")]] <- tibble(
      metric = met,
      method = meth,
      icc_mean = icc_mean$value,
      icc_max = icc_max$value,
      icc_min = icc_min$value,            
      icc_trimmed = icc_trimmed$value,
      icc_median = icc_median$value,
      ccc_mean = ccc_mean$est,
      ccc_max = ccc_max$est,
      ccc_min = ccc_min$est,              
      ccc_trimmed = ccc_trimmed$est,
      ccc_median = ccc_median$est
    )
  }
}

traditional_results_df <- bind_rows(traditional_results)

## Plotting traditional minimum and maximum ICC and CCCs

icc_long <- traditional_results_df %>%
  select(metric, method, icc_min, icc_max) %>%
  pivot_longer(cols = starts_with("icc_"),
               names_to = "stat",
               values_to = "icc") %>%
  mutate(stat = ifelse(stat == "icc_min", "Minimum", "Maximum"))

ccc_long <- traditional_results_df %>%
  select(metric, method, ccc_min, ccc_max) %>%
  pivot_longer(cols = starts_with("ccc_"),
               names_to = "stat",
               values_to = "ccc") %>%
  mutate(stat = ifelse(stat == "ccc_min", "Minimum", "Maximum"))

ggplot(icc_long, aes(x = method, y = icc, fill = stat)) +
  geom_col(position = "dodge") +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "ICC for Minimum vs Maximum",
       y = "ICC", x = "Method") +
  theme_minimal()

ggplot(ccc_long, aes(x = method, y = ccc, fill = stat)) +
  geom_col(position = "dodge") +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "CCC for Minimum vs Maximum",
       y = "CCC", x = "Method") +
  theme_minimal()
