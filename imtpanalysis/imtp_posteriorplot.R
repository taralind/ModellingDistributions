library(brms)
library(tidybayes)
library(tidyverse)
library(patchwork)

# choose metric and method to plot
met <- "Net Force at 100ms [N]"
method <- "TRAD"

imtp <- read.csv("imtpanalysis/rawdata_synthetic.csv")

imtp$Method <- as.factor(imtp$Method)
imtp$Participant <- as.factor(imtp$Participant)
imtp$Day <- as.factor(imtp$Day)
imtp$Value <- as.numeric(imtp$Value)

df <- imtp %>% 
  filter(Metric == met) %>%
  filter(Method == method)

#formula <- bf(Value ~ (1|Participant))
formula <- bf(Value ~ (1|Participant), sigma ~ (1|Participant))

m1 <- brm(formula, data = df %>% filter(Day == "D1"), family = gaussian(), chains = 4, cores = 4, iter = 4000)
m2 <- brm(formula, data = df %>% filter(Day == "D2"), family = gaussian(), chains = 4, cores = 4, iter = 4000)

# model parameters
# summary(m1)
# summary(m2)

# posteriors
epred_d1 <- df %>%
  filter(Day == "D1") %>%
  modelr::data_grid(Participant) %>%
  add_epred_draws(m1) %>%
  mutate(Day = "D1")

epred_d2 <- df %>%
  filter(Day == "D2") %>%
  modelr::data_grid(Participant) %>%
  add_epred_draws(m2) %>%
  mutate(Day = "D2")

epred_all <- bind_rows(epred_d1, epred_d2)

# summary stats
obs_summary <- df %>%
  group_by(Participant, Day) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    max_val = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(mean_val, max_val),
    names_to = "stat",
    values_to = "Value"
  )

# Left: Posterior distributions
p1 <- ggplot(epred_all, aes(x = .epred, y = Participant, fill = Day)) +
  stat_slab(alpha = 0.7) +
  labs(
    x = met,
    y = "Participant",
    title = "Posterior Distributions"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# Right: Raw summary stats
p2 <- ggplot(obs_summary, aes(x = Value, y = Participant, colour = Day, shape = stat)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(mean_val = 16, max_val = 18),
                     labels = c("Max", "Mean")) +
  labs(
    x = met,
    y = NULL,
    title = "Observed Values: Mean and Maximum",
    shape = "Statistic"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

p1 + p2 + plot_layout(ncol = 2)

