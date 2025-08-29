library(tidyverse)
library(mclust)
library(moments)
library(brms)
library(gganimate)

data_file <- "data/all_athlete_trials_2022_2025_longjumpwomen.csv"
event <- "Long Jump"

#### DATA

tiladata <- read.csv(data_file)

#### CLEANING

tiladata <- tiladata %>%
  filter(result != "" & !is.na(result)) %>%  # remove blank results
  filter(event_name == event) # just the event of interest

# remove comps where no trial data was successfully retrieved
tiladata <- tiladata %>%  group_by(athlete, date, location, event_name) %>%
  filter(n() > 1)

# add trial numbers
tiladata <- tiladata %>%
  group_by(athlete, event_name, location, date) %>%
  mutate(trial_number = row_number()) %>%
  ungroup()

# athlete age column
tiladata <- tiladata %>% 
  mutate(
    date_has_year = str_detect(date, "\\d{4}$"),
    
    full_date_str = if_else(
      date_has_year,
      date,                            
      paste(date, year)            
    ),
    
    full_date = dmy(full_date_str),

    dob_day_month = str_remove(DOB, "\\s\\d{2}$"),
    dob_year_2digit = str_extract(DOB, "\\d{2}$"),
    dob_year_full = if_else(
      as.numeric(dob_year_2digit) > 30,
      paste0("19", dob_year_2digit),
      paste0("20", dob_year_2digit)
    ),
    dob_parsed = dmy(paste(dob_day_month, dob_year_full)),
    
    age = time_length(interval(dob_parsed, full_date), "years")
  )

tiladata$result_numeric <- as.numeric(tiladata$result)
tiladata <- tiladata %>% drop_na(result_numeric)


#### PLOTTING AGGREGATE DISTRIBUTIONS

skew_val <- skewness(tiladata$result_numeric, na.rm = TRUE)

# aggregate distributions
ggplot(tiladata, aes(x=result_numeric)) + geom_density() +
  labs(title = paste0(event, " Women's (Skewness = ", round(skew_val, 2), ")"),
       x = "Result (m)") +
  theme_minimal()


#### CLEANING WITH GAUSSIAN MIXTURE MODEL
# just necessary for long jump men's and women's

# fit a 2-component Gaussian mixture model
model <- Mclust(tiladata$result_numeric, G = 2)

tiladata$component <- model$classification
tiladata$posterior1 <- model$z[,1]  # probability of being in component 1
tiladata$posterior2 <- model$z[,2]  # probability of being in component 2

# keep the ones likely to be in component 2
tiladata_clean <- tiladata[tiladata$posterior2 > 0.5, ]

# plot cleaned against original
ggplot(tiladata, aes(x=result_numeric)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_density(data = tiladata_clean, aes(x=result_numeric), colour = "darkblue", size = 1.2) +
  theme_minimal() +
  labs(
    x="Result",
    y="Density")

#### PLOTTING INDIVIDUAL DISTRIBUTIONS

# 9 athletes with the most observations
top_athletes <- tiladata_clean %>%
  count(athlete) %>%
  arrange(desc(n)) %>%
  slice_head(n = 9)

athlete_map <- top_athletes %>%
  mutate(facet_label = paste0("Athlete ", row_number(), 
                              " (n=", n, 
                              ", skew=", round(sapply(athlete, function(a) skewness(tiladata_clean$result_numeric[tiladata_clean$athlete == a], na.rm = TRUE)), 2), ")")) %>%
  select(athlete, facet_label)

plot_data <- tiladata_clean %>%
  filter(athlete %in% top_athletes$athlete) %>%
  left_join(athlete_map, by = "athlete")

ggplot(plot_data, aes(x = result_numeric)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  facet_wrap(~facet_label, ncol = 3) +
  labs(x = "Result (m)", y = "Density") +
  theme_minimal()


#### MODELLING

m1_skew <- brm(bf(result_numeric ~ (1|athlete)), 
          data = tiladata_clean, family = skew_normal(), chains = 4, cores = 4, iter = 4000)

m1_skew_sigma <- brm(bf(result_numeric ~ (1|athlete), sigma ~ (1|athlete)), 
               data = tiladata_clean, family = skew_normal(), chains = 4, cores = 4, iter = 4000)

m1_gauss <- brm(bf(result_numeric ~ (1|athlete)), 
               data = tiladata_clean, family = gaussian(), chains = 4, cores = 4, iter = 4000)

m1_gauss_sigma <- brm(bf(result_numeric ~ (1|athlete), sigma ~ (1|athlete)), 
                data = tiladata_clean, family = gaussian(), chains = 4, cores = 4, iter = 4000)

# model evaluation
pp_check(m1_skew)
pp_check(m1_skew_sigma)
pp_check(m1_gauss)
pp_check(m1_gauss_sigma)

# other model specifications trialled
m2 <- brm(bf(result_numeric ~ (1|athlete) + (1|location), sigma ~ (1|athlete)), 
          data = tiladata_clean, family = skew_normal(), chains = 4, cores = 4, iter = 4000)

m3 <- brm(bf(result_numeric ~ age + (1|athlete), sigma ~ (1|athlete)), 
          data = tiladata_clean, family = skew_normal(), chains = 4, cores = 4, iter = 4000)

# model evaluation
pp_check(m2)
pp_check(m3)

# keeping m1_skew_sigma for following analysis
m <- m1_skew_sigma

summary(m)

#### INDIVIDUAL MODELLED DISTRIBUTION - BAYESIAN UPDATING EXAMPLE

athlete_name <- "Malaika Mihambo" 

athlete_data_all <- tiladata_clean %>%
  filter(athlete == athlete_name) %>%
  drop_na(result_numeric)

athlete_data_2025 <- athlete_data_all %>%
  filter(full_date >= as.Date("2025-01-01")) %>%
  arrange(full_date)

train_data <- tiladata_clean %>%
  drop_na(result_numeric) %>%
  filter(full_date < as.Date("2025-01-01"))

# fit model to training data
base_model <- brm(
  bf(result_numeric ~ (1 | athlete), sigma ~ (1 | athlete)),
  family = skew_normal(),
  data = train_data,
  chains = 2, iter = 2000, cores = 2
)

# group 2025 observations by date (i.e. competition)
date_groups <- athlete_data_2025 %>%
  group_by(full_date) %>%
  arrange(full_date, .by_group = TRUE) %>%
  group_split()

posterior_list <- list()
obs_points_list <- list()

# sequentially update model with data from each date
for (i in seq_along(date_groups)) {
  current_data <- bind_rows(date_groups[1:i])
  
  updated_model <- update(
    base_model,
    newdata = bind_rows(train_data, current_data),
    recompile = FALSE
  )
  
  # get posterior predictive draws for the athlete
  pred_data <- data.frame(athlete = athlete_name)
  preds <- posterior_predict(updated_model, newdata = pred_data)
  
  posterior_list[[i]] <- data.frame(
    value = as.numeric(preds),
    update_step = i,
    date = as.character(max(current_data$full_date))
  )
  
  obs_points_list[[i]] <- data.frame(
    obs_value = current_data$result_numeric,
    update_step = i,
    date = current_data$full_date
  )
}


# combine posterior draws
posterior_df <- bind_rows(posterior_list)
obs_points_df <- bind_rows(obs_points_list)

dates_for_frames <- posterior_df %>%
  distinct(update_step, date) %>%
  arrange(update_step) %>%
  pull(date)

# animated posterior plot
p <- ggplot(posterior_df, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_point(
    data = obs_points_df,
    aes(x = obs_value, y = 0, colour = as.factor(date)),
    inherit.aes = FALSE,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    title = "Posterior for Malaika Mihambo's Performance\nAfter {frame} Dates — {closest_date}",
    x = "Estimated Result (on original scale)",
    y = "Density",
    colour = "Observed Date"
  ) +
  theme_minimal(base_size = 14) +
  transition_manual(update_step, closest_date = dates_for_frames)

anim <- animate(p, width = 800, height = 600, fps = 5, renderer = gifski_renderer())
#anim_save("posterior_updates.gif", animate(p, fps = 4, width = 800, height = 600))

# static plots (a distribution per date)
p_all <- ggplot() +
  geom_density(
    data = posterior_df,
    aes(x = value),
    fill = "skyblue",
    alpha = 0.6
  ) +
  geom_point(
    data = obs_points_df,
    aes(x = obs_value, y = 0, colour = as.factor(date)),
    size = 2.5,
    alpha = 0.8
  ) +
  labs(
    title = "Posterior Predictive for Athlete 1 Across Update Steps",
    x = "Predicted Result (on original scale)",
    y = "Density",
    colour = "Observed Date"
  ) +
  facet_wrap(~ update_step, ncol = 4) +
  theme_minimal(base_size = 14)

ggsave(
  filename = "posterior_all_steps.png",
  plot = p_all,
  width = 16, height = 12
)

#### PROBABILITIES OF ACHIEVING RESULTS

# plot of final distribution (with all 2025 observations) and markings for results of interest

# just take last plot from above (so all obs)
last_step <- max(posterior_df$update_step)
plot_df <- posterior_df %>% filter(update_step == last_step)
obs_df  <- obs_points_df %>% filter(update_step == last_step)
plot_date <- unique(plot_df$date)

# Reference marks
ref_lines <- data.frame(
  label = c("Paris 2024 Bronze", "Paris 2024 Gold", "World Record"),
  value = c(6.96, 7.10, 7.52)
)

ggplot(plot_df, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_point(
    data = obs_df,
    aes(x = obs_value, y = 0, colour = as.factor(date)),
    inherit.aes = FALSE,
    size = 3,
    alpha = 0.8
  ) +
  geom_vline(
    data = ref_lines,
    aes(xintercept = value),
    linetype = "dashed",
    colour = "red"
  ) +
  geom_text(
    data = ref_lines,
    aes(x = value, y = 0.05, label = label),
    angle = 90,
    vjust = -0.5,
    hjust = 0,
    size = 4,
    colour = "red"
  ) +
  labs(
    title = paste0("Posterior Predictive for Athlete 1"),
    x = "Predicted Result (on original scale)",
    y = "Density",
    colour = "Observed Date"
  ) +
  theme_minimal(base_size = 14)

# posterior predictive draws
preds <- plot_df$value  

# probabilities
prob_bronze <- mean(preds >= 6.96)
prob_gold   <- mean(preds >= 7.10)
prob_wr     <- mean(preds >= 7.52)

# table
data.frame(
  Outcome = c("Medal", "Gold", "World Record"),
  Threshold = c(6.96, 7.10, 7.52),
  Probability = c(prob_bronze, prob_gold, prob_wr)
)

