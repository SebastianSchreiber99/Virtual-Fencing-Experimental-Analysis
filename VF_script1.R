packages_to_load <- c("tidyverse", "glmmTMB", "emmeans", "ggplot2", "sf", "terra", "lubridate", "lme4", "performance", "DHARMa") # create vector of R package names that we know are needed in the rest of the code


## ----load_libraries----
lapply(packages_to_load, library, character.only = TRUE)


##-----------------------------VF Efficacy-----------------------------

#-----Escapes - GLMM----
groups <- read.csv("C:/Users/spsch/Documents/R/Virtual_Fence/Groups.csv")

esc <- read.csv(file = "C:/Users/spsch/Documents/R/Virtual_Fence/Escapes2.csv")


# We need to include zeros to give the model the information that escapes did not occur under those conditions. If you only analyze rows where an escape happened, the model thinks missing combinations are “unknown,” not “zero,” which (1) throws away most of your data, (2) inflates uncertainty, and (3) can easily lead to a non‑significant result even when the raw counts suggest a real difference.

escape_types <- c("Electric", "Virtual")

# Create full grid of combinations
full_grid_esc <- expand_grid(
  groups,
  Escape_Type = escape_types)

# Summarize escapes by animal × period
esc_summary <- esc %>%
  group_by(Animal_ID, TRT, Group, Escape_Type) %>%
  summarise(
    Total_escapes = sum(Escape_no, na.rm = TRUE),
    .groups = "drop")


# Join with full grid to retain zeroes!
esc_summary_df <- full_grid_esc %>%
  left_join(esc_summary, by = c("Animal_ID", "TRT", "Group", "Escape_Type")) %>%
  mutate(Total_escapes = replace_na(Total_escapes, 0))



# Why use GLMM rather than LMM? Data in this case are counts. Counts are not normally distributed and their variance usually depends on the mean (e.g., Poisson or negative binomial behavior).GLMMs extend LMMs by allowing different distributions (e.g., Poisson, Negative Binomial) for the response variable. Instead of modeling the raw response, GLMMs use a link function (log link for counts), ensuring predictions stay positive.

# Fit Poisson GLMM
# An offset is not needed because all animals have the same rate of exposure to potential escapes
m1 <- glmmTMB(
  Total_escapes ~ Escape_Type + (1|Group),
  family = poisson,
  data = esc_summary_df)

summary(m1)

# Rate ratio
emmeans(m1, pairwise ~ Escape_Type, type = "response")

# Model Validation
sim_res <- simulateResiduals(m1, n = 1000)
testDispersion(sim_res)


#-----Escapes - Bar plot-----

# Get estimated mean rates from your model (m_pois or m_nb)
emm <- emmeans(m1, ~ Escape_Type, type = "response")

emm_df <- as.data.frame(emm) %>%
  mutate(
    SE = SE,                # Standard Error
    ymin = rate - SE,   # Lower error bar
    ymax = rate + SE    # Upper error bar
  )


# Plot
p1 <- ggplot(emm_df, aes(x = Escape_Type, y = rate, fill = Escape_Type)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  labs(
    x = "Fence Type",
    y = "Estimated Escape Rate (escapes per animal-day)",
    title = "Escapes by Fence Type"
  ) +
  scale_fill_manual(values = c("Electric" = "firebrick", "Virtual" = "steelblue")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.4))


#ggsave(filename = "Escapes.png",
#       plot = p1,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 5,
#       height = 5,
#       units = "in")






#-----Percentage Points In/Out-----

#eshep_df1 <- read.csv(file = "C:/Users/spsch/Documents/R/Virtual_Fence/2025-07-13.csv")
#eshep_df2 <- read.csv(file = "C:/Users/spsch/Documents/R/Virtual_Fence/2025-07-20.csv")

#eshep_df <- rbind.data.frame(eshep_df1, eshep_df2) %>% # Filter to include only the three training and four testing days
#  filter(Time..UTC. > "2025-07-07 11:55:00" & Time..UTC. < "2025-07-14 12:05:00")

#eshep_full_df <- eshep_df %>%
#  left_join(groups, by = "Neckband.ID") %>%
#  filter(!is.na(Animal_ID))

#write.csv(eshep_full_df, file = "C:/Users/spsch/Documents/R/Virtual_Fence/eshep_full_df.csv")




#Read in and clean data
eshep_full_df <- read.csv(file = "C:/Users/spsch/Documents/R/Virtual_Fence/eshep_full_df.csv")

exclusion_zones <- st_read(dsn = "C:/Users/spsch/Documents/R/Virtual_Fence/Exclusion_Zones.kml")

training_VPs <- st_read(dsn = "C:/Users/spsch/Documents/R/Virtual_Fence/Training VPs.kml")



exclusion_zones$Description <- c("west", "east", "east", "east", "east", "east",
                                 "east", "west", "west", "west", "west", "west")

eshep_full_df_noNA <- eshep_full_df %>%
  na.omit(latitude) %>%
  na.omit(longitude)

# Convert df to spatial points object
eshep_full_sf <- st_as_sf(eshep_full_df_noNA, coords = c("longitude", "latitude"), crs = crs(training_VPs))


# Clean points outside pasture boundary with slight buffer zone
perimeter_buffer <- st_read(dsn = "C:/Users/spsch/Documents/R/Virtual_Fence/perimeter_buffer.kml")

eshep_full_sf <- eshep_full_sf %>%
  mutate(in_perimeter = st_within(eshep_full_sf, perimeter_buffer, sparse = FALSE)[ ,1]) %>%
  filter(in_perimeter == TRUE)


# Categorize points by experiment period
eshep_full_sf <- eshep_full_sf %>%
  mutate(period = 
           ifelse(Time..UTC. <= "2025-07-08 11:59:59", "Tr1",
                  ifelse(Time..UTC. >= "2025-07-08 12:00:00" &
                           Time..UTC. <= "2025-07-09 11:59:59", "Tr2",
                         ifelse(Time..UTC. >= "2025-07-09 12:00:00" &
                                  Time..UTC. <= "2025-07-10 11:59:59", "Tr3",
                                ifelse(Time..UTC. >= "2025-07-10 12:00:00" &
                                         Time..UTC. <= "2025-07-11 11:59:59", "E1",
                                       ifelse(Time..UTC. >= "2025-07-11 12:00:00" &
                                                Time..UTC. <= "2025-07-12 11:59:59", "W1",
                                              ifelse(Time..UTC. >= "2025-07-12 12:00:00" &
                                                       Time..UTC. <= "2025-07-13 11:59:59", "E2",
                                                     ifelse(Time..UTC. >= "2025-07-13 12:00:00" &
                                                              Time..UTC. <= "2025-07-14 12:00:00", "W2",
                                                            NA))))))))





# Training phase
training1_inside <- st_within(eshep_full_sf, training_VPs[1, ], sparse = FALSE)[, 1]
training2_inside <- st_within(eshep_full_sf, training_VPs[2, ], sparse = FALSE)[, 1]
training3_inside <- st_within(eshep_full_sf, training_VPs[3, ], sparse = FALSE)[, 1]

# Exclusion zones
west_matrix <- st_within(eshep_full_sf, filter(exclusion_zones, Description == "west"), sparse = FALSE)
west_ex_inside <- apply(west_matrix, 1, any)

east_matrix <- st_within(eshep_full_sf, filter(exclusion_zones, Description == "east"), sparse = FALSE)
east_ex_inside <- apply(east_matrix, 1, any)


# Bind the columns
eshep_full_sf_ex <- eshep_full_sf %>%
  mutate(
    training1_inside = training1_inside,
    training2_inside = training2_inside,
    training3_inside = training3_inside,
    west_ex_inside = west_ex_inside,
    east_ex_inside = east_ex_inside)


eshep_full_sf_ex <- eshep_full_sf_ex %>%
  mutate(
    inside_VF = case_when(
      period == "Tr1" ~ training1_inside,
      period == "Tr2" ~ training2_inside,
      period == "Tr3" ~ training3_inside,
      period %in% c("E1", "E2") ~ !east_ex_inside,
      period %in% c("W1", "W2") ~ !west_ex_inside,
      TRUE ~ NA))


# Calculate proportion of points inside VF for each animal and convert to df
m2_df <- eshep_full_sf_ex %>%
  group_by(Animal_ID, TRT, Group, period) %>%
  summarise(
    n_pts = n(),
    n_in = sum(inside_VF, na.rm = TRUE),
    prop_in = n_in / n_pts,
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  select(!last_col())


# Fit binomial GLMM
# Animal_ID random effect was removed because it caused model singularity. Random effect estiamte for animal_ID was 2.462e-34
m2 <- glmmTMB(
  cbind(n_in, n_pts - n_in) ~ TRT + (1|Group) + (1|period), 
  family = binomial,
  data = m2_df)

summary(m2)


# gives probability inside
emm2 <- emmeans(m2, ~ TRT, type = "response")  
print(emm2)



#-----Percent in - Bar Plot-------
emm_df <- as.data.frame(emm2) %>%
  mutate(
    SE = SE,
    ymin = asymp.LCL,
    ymax = asymp.UCL
  ) %>%
  mutate(treatment = ifelse(TRT == "CNT", "Control", "Treatment"))

p2 <- ggplot(emm_df, aes(x = treatment, y = prob, fill = treatment)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  labs(
    x = element_blank(),
    y = "% Points Inside Virtual Boundary",
    title = "Virtual Fence Efficacy by Treatment") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0.90, 1)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("darkorange", "forestgreen"))



#ggsave(filename = "Efficacy.png",
#      plot = p2,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 5,
#       height = 5,
#       units = "in")



#-----------------Rate of Learning by hour-------------------------------

# Remove any points with cues that are not near VF

# Read in shapefile for buffer zones
buffers <- st_read(dsn = "C:/Users/spsch/Documents/R/Virtual_Fence/Buffers.kml")

# Only keep points associated with cues
eshep_cues_sf_ex <- eshep_full_sf_ex %>%
  filter(No..Audios >= 1 | No..Pulses >= 1)

# Only keep points within buffer zones during phase 2
eshep_buffer_ex <- eshep_cues_sf_ex %>%
  filter(period %in% c("E1", "E2", "W1", "W2"))

eshep_buffer_ex <- st_intersection(eshep_buffer_ex, buffers)

eshep_buffer_ex <- eshep_buffer_ex %>%
  mutate(in_buffer = ifelse(Group == 1 & period %in% c("E1", "E2") & Name == "East_1_buffer", TRUE,
                            ifelse(Group == 1 & period %in% c("W1", "W2") & Name == "West_1_buffer", TRUE,
                                   ifelse(Group == 2 & period %in% c("E1", "E2") & Name == "East_2_buffer", TRUE,
                                          ifelse(Group == 2 & period %in% c("W1", "W2") & Name == "West_2_buffer", TRUE,
                                                 ifelse(Group == 3 & period %in% c("E1", "E2") & Name == "East_3_buffer", TRUE,
                                                        ifelse(Group == 3 & period %in% c("W1", "W2") & Name == "West_3_buffer", TRUE,
                                                               ifelse(Group == 4 & period %in% c("E1", "E2") & Name == "East_4_buffer", TRUE,
                                                                      ifelse(Group == 4 & period %in% c("W1", "W2") & Name == "West_4_buffer", TRUE,
                                                                             ifelse(Group == 5 & period %in% c("E1", "E2") & Name == "East_5_buffer", TRUE,
                                                                                    ifelse(Group == 5 & period %in% c("W1", "W2") & Name == "West_5_buffer", TRUE,
                                                                                           ifelse(Group == 6 & period %in% c("E1", "E2") & Name == "East_6_buffer", TRUE,
                                                                                                  ifelse(Group == 6 & period %in% c("W1", "W2") & Name == "West_1_buffer", TRUE, FALSE)))))))))))))



# Combine points from phases 1 and 2 and convert to simple df
eshep_buffer1_2 <- eshep_buffer_ex %>%
  filter(in_buffer == "TRUE") %>%
  as.data.frame() %>%
  select(Time..UTC., No..Audios, No..Pulses, Animal_ID, Group, TRT, period)

eshep_buffer2_2 <- eshep_cues_sf_ex %>%
  filter(period %in% c("Tr1", "Tr2", "Tr3")) %>%
  as.data.frame() %>%
  select(Time..UTC., No..Audios, No..Pulses, Animal_ID, Group, TRT, period)

eshep_cue_df_clean <- rbind.data.frame(eshep_buffer1_2,eshep_buffer2_2)


# Process timestamps and bin hours
eshep_cue_df_clean <- eshep_cue_df_clean %>%
  mutate(
    datetime = ymd_hms(`Time..UTC.`),
    hour = hour(datetime) + minute(datetime) / 60 + second(datetime) / 3600,
    hour = hour - 12,
    hour = if_else(hour < 0, hour + 24, hour),
    hour_bin = floor(hour)
  )


# Define all periods and hour bins
periods <- c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2")
hour_bins <- 0:23

# Create full grid of combinations
full_grid <- expand_grid(
  groups,
  period = periods,
  hour_bin = hour_bins)


# Summarize cue data by animal × period × hour
cue_summary <- eshep_cue_df_clean %>%
  group_by(Animal_ID, TRT, Group, period, hour_bin) %>%
  summarise(
    total_audios = sum(`No..Audios`, na.rm = TRUE),
    total_pulses = sum(`No..Pulses`, na.rm = TRUE),
    .groups = "drop")


# Join with full grid to retain zeroes!!!! Retaining zeros is highly important!!!!!!!!
hour_summary <- full_grid %>%
  left_join(cue_summary, by = c("Animal_ID", "TRT", "Group", "period", "hour_bin")) %>%
  mutate(
    total_audios = replace_na(total_audios, 0),
    total_pulses = replace_na(total_pulses, 0))


# Based on Shapiro-Wilk W, and visual histogram analyses, the data is clearly non-normal and strongly skewed, which confirms that a Poisson GLMM is the right modeling approach.


# Fit Poisson GLMM
m3 <- glmmTMB(
  total_audios ~ hour_bin * TRT + (1|Animal_ID) + (1|Group) + (1|period),
  family = poisson,
  data = hour_summary)

summary(m3)


# Poisson model yielded fixed effect standard errors that were all NaN: This means the model cannot reliably estimate the uncertainty of the coefficients — suggesting complete or quasi-complete separation, overfitting, or rank deficiency.
# Thus, the variance may far exceed the mean (which is likely with rare pulses), so we went with a negative binomial family.
m3.1 <- glmmTMB(
  total_pulses ~ hour_bin * TRT + (1|Animal_ID) + (1|Group) + (1|period),
  family = nbinom2,
  data = hour_summary)

summary(m3.1)


#Originally audio-shock ratio or the percentage of cues that were audio was modeled. But this is misleading because that percentage has a lower bound of 50%, since animals cannot receive an shock without an audio warning first. Thus, shocks are a subset of audio events, and a better framing of the question is: Given an audio, what’s the probability it was followed by a pulse?
#To model the probability that an audio cue leads to a pulse — i.e., how often a cue escalates from an audio-only warning to an audio+shock correction.
#This is best framed as a conditional probability:
#Given that an audio was delivered, what’s the chance it was followed by a pulse?
#This model uses a binomial response in the form cbind(successes, failures).
#In this case: Successes = total_pulses: how many audios escalated to a pulse
#Failures = total_audios - total_pulses: how many audios did not lead to a pulse
#Total trials = total_audios
#The audios are the "trials": they occur first and may or may not escalate.
#The pulses are the "successes": they occur only if the audio was ineffective.
#This model structure treats each hour (or bin) as an opportunity to observe that escalation rate.

m3.2 <- glmmTMB(
  cbind(total_pulses, total_audios - total_pulses) ~ hour_bin * TRT + 
    (1|Animal_ID) + (1|Group) + (1|period),
  family = binomial,
  data = hour_summary)

summary(m3.2)


# Model validation and over dispersion checks
check_overdispersion(m3)
check_overdispersion(m3.1)
check_overdispersion(m3.2)

res3 <- simulateResiduals(fittedModel = m3, n = 1000)
plot(res3)
testDispersion(res3)

res3.1 <- simulateResiduals(fittedModel = m3.1, n = 1000)
plot(res3.1)
testDispersion(res3.1)

res3.2 <- simulateResiduals(fittedModel = m3.2, n = 1000)
plot(res3.2)
testDispersion(res3.2)

# Validation results:
#m3 and m3.2 look solid. Their dispersion ratios are close to 1, and both performance::check_overdispersion() and DHARMa::testDispersion() agree: no overdispersion detected.
#m3.1, however, shows strong underdispersion, which likely indicates: Data artifacts (e.g., very sparse response), or overfitting (too many fixed effects or unnecessary complexity).



#-----------------Rate of Learning - Plots-------------------------------
cue_summary <- hour_summary %>%
  pivot_longer(
    cols = c(total_audios, total_pulses),
    names_to = "cue_type",
    values_to = "count"
  ) %>%
  mutate(cue_type = recode(cue_type,
                           total_audios = "Audio",
                           total_pulses = "Pulse")) %>%
  group_by(hour_bin, TRT, cue_type) %>%
  summarise(total_cues = sum(count, na.rm = TRUE), .groups = "drop")%>%
  mutate(treatment = ifelse(TRT == "TRT", "Treatment", "Control"))


# Plot of number of cues by hour
p3 <- ggplot(cue_summary, aes(x = hour_bin, y = total_cues,
                     color = treatment,
                     linetype = cue_type,
                     group = interaction(cue_type, treatment))) +
 # geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = .7) +
  labs(
    title = "Cues by Hour (smoothed)",
    x = "Hour of Period",
    y = "Cues per Animal per Hour",
    color = "Treatment",
    linetype = "Cue Type"
  ) +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_linetype_manual(values = c("Audio" = "dotted", "Pulse" = "solid")) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())


#ggsave(filename = "Cues_by_hour_np.png",
#       plot = p3,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 6,
#       height = 4,
#       units = "in")





escalation_summary <- hour_summary %>%
  mutate(
    prob_escalation = total_pulses / total_audios) %>%
  group_by(hour_bin, TRT) %>%
  summarise(
    mean_prob = mean(prob_escalation, na.rm = TRUE),
    se = sd(prob_escalation, na.rm = TRUE) / sqrt(n()),
    .groups = "drop") %>%
  mutate(
    ymin = pmax(0, mean_prob - se),
    ymax = pmin(1, mean_prob + se)) %>%
  mutate(treatment = ifelse(TRT == "TRT", "Treatment", "Control"))


# Escalation by hour
p4 <- ggplot(escalation_summary, aes(x = hour_bin, y = mean_prob, color = treatment, group = treatment)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.75) +
  labs(
    title = "Smoothed Probability of Escalation",
    x = "Hour of Period",
    y = "Probability Pulse Follows Audio",
    color = "Treatment"
  ) +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())


#ggsave(filename = "Escalation.png",
#       plot = p4,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 6,
#       height = 4,
#       units = "in")





#-----------------Rate of Learning by period-------------------------------

# Create full grid of combinations
full_grid_periods <- expand_grid(
  groups,
  period = periods)


# Summarize cue data by animal × period × period
cue_summary_by_period <- eshep_cue_df_clean %>%
  group_by(Animal_ID, TRT, Group, period) %>%
  summarise(
    total_audios = sum(`No..Audios`, na.rm = TRUE),
    total_pulses = sum(`No..Pulses`, na.rm = TRUE),
    .groups = "drop")


# Join with full grid to retain zeroes!!!! Retaining zeros is highly important!!!!!!!!
# Also convert periods to numerical days and categorize by phases
period_summary <- full_grid_periods %>%
  left_join(cue_summary_by_period, by = c("Animal_ID", "TRT", "Group", "period")) %>%
  mutate(
    total_audios = replace_na(total_audios, 0),
    total_pulses = replace_na(total_pulses, 0)) %>%
  mutate(total_cues = total_audios + total_pulses) %>%
  mutate(day = as.numeric(ifelse(period == "Tr1", 1,
                      ifelse(period == "Tr2", 2,
                             ifelse(period == "Tr3", 3,
                                    ifelse(period == "E1", 4,
                                           ifelse(period == "W1", 5,
                                                  ifelse(period == "E2", 6,
                                                         ifelse(period == "W2", 7, NA))))))))) %>%
  mutate(phase = ifelse(day <= 3, "Training", "Exclusion Zones"))



# Fit Poisson GLMM
m3.3 <- glmmTMB(
  total_audios ~ day * TRT + (1|Animal_ID) + (1|Group),
  family = poisson,
  data = period_summary)

summary(m3.3)



m3.4 <- glmmTMB(
  total_pulses ~ day * TRT + (1|Animal_ID) + (1|Group),
  family = poisson,
  data = period_summary)

summary(m3.4)



# Fit binomial GLMM
m3.5 <- glmmTMB(
  cbind(total_pulses, total_audios - total_pulses) ~ day * TRT + 
    (1|Animal_ID) + (1|Group),
  family = binomial,
  data = period_summary)

summary(m3.5)


#-----------------Rate of Learning by period - Plots-------------------------------


cue_period_long <- period_summary %>%
  pivot_longer(
    cols = c(total_audios, total_pulses),
    names_to = "cue_type",
    values_to = "count"
  ) %>%
  mutate(
    cue_type = recode(cue_type,
                      total_audios = "Audio",
                      total_pulses = "Pulse"))  %>%
  mutate(treatment = ifelse(TRT == "TRT", "Treatment", "Control"))

cue_period_long$period <- factor(cue_period_long$period,
                                 levels = c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2"))




p5 <- ggplot(cue_period_long, aes(x = period, y = count,
                            color = treatment,
                            linetype = cue_type,
                            group = interaction(cue_type, treatment))) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.4) +
  labs(
    title = "Cues by Period (smoothed)",
    x = "Period (chronological)",
    y = "Cues per Animal per Period",
    color = "Treatment",
    linetype = "Cue Type"
  ) +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_linetype_manual(values = c("Audio" = "dotted", "Pulse" = "solid")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank())



#ggsave(filename = "Cues_by_period.png",
#       plot = p5,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 6,
#       height = 4,
#       units = "in")






escalation_summary2 <- period_summary %>%
  mutate(
    prob_escalation = total_pulses / total_audios) %>%
  group_by(period, TRT) %>%
  summarise(
    mean_prob = mean(prob_escalation, na.rm = TRUE),
    se = sd(prob_escalation, na.rm = TRUE) / sqrt(n()),
    .groups = "drop") %>%
  mutate(
    ymin = pmax(0, mean_prob - se),
    ymax = pmin(1, mean_prob + se)) %>%
  mutate(treatment = ifelse(TRT == "TRT", "Treatment", "Control"))

escalation_summary2$period <- factor(escalation_summary2$period,
                                     levels = c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2"))
                                 


# Escalation by hour
p6 <- ggplot(escalation_summary2, aes(x = period, y = mean_prob, color = treatment, group = treatment)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2, span = 0.6) +
  labs(
    title = "Probability of Escalation (smoothed)",
    x = "Period (Chronological)",
    y = "Probability Pulse Follows Audio",
    color = "Treatment"
  ) +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


#ggsave(filename = "Escalation2.png",
#       plot = p6,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 6,
#       height = 4,
#       units = "in")



#-----------------Rate of Learning by phase-------------------------------

phases <- c("1: Training", "2: Exclusion Zone Testing")

# Create full grid of combinations to retain zeros
full_grid_phase <- expand_grid(
  groups,
  phase = phases)


# Summarize cue data by animal × period × phase
cue_summary_by_phase <- eshep_cue_df_clean %>%
  mutate(phase = ifelse(period %in% c("Tr1", "Tr2", "Tr3"), "1: Training",
                        ifelse(period %in% c("E1", "E2", "W1", "W2"), "2: Exclusion Zone Testing", NA))) %>%
  group_by(Animal_ID, TRT, Group, phase) %>%
  summarise(
    total_audios = sum(`No..Audios`, na.rm = TRUE),
    total_pulses = sum(`No..Pulses`, na.rm = TRUE),
    .groups = "drop")


# Join with full grid to retain zeroes!
phase_summary <- full_grid_phase %>%
  left_join(cue_summary_by_phase, by = c("Animal_ID", "TRT", "Group", "phase")) %>%
  mutate(
    total_audios = replace_na(total_audios, 0),
    total_pulses = replace_na(total_pulses, 0)) %>%
  mutate(total_cues = total_audios + total_pulses)



# Fit Poisson GLMM

#Phase 1
m3.6 <- glmmTMB(
  total_audios ~ TRT + (1|Animal_ID) + (1|Group),
  family = poisson,
  data = filter(phase_summary, phase == "1: Training"))

summary(m3.6)


#Phase 2
m3.6.1 <- glmmTMB(
  total_audios ~ TRT + (1|Animal_ID) + (1|Group),
  family = poisson,
  data = filter(phase_summary, phase == "2: Exclusion Zone Testing"))

summary(m3.6.1)


#Phase 1
m3.7 <- glmmTMB(
  total_pulses ~ TRT + (1|Animal_ID) + (1|Group),
  family = poisson,
  data = filter(phase_summary, phase == "1: Training"))

summary(m3.7)


#Phase 2
m3.7.1 <- glmmTMB(
  total_pulses ~ TRT + (1|Animal_ID) + (1|Group),
  family = poisson,
  data = filter(phase_summary, phase == "2: Exclusion Zone Testing"))

summary(m3.7.1)



# Fit binomial GLMM

#Phase 1
m3.8 <- glmmTMB(
  cbind(total_pulses, total_audios - total_pulses) ~ TRT + 
    (1|Animal_ID) + (1|Group),
  family = binomial,
  data = filter(phase_summary, phase == "1: Training"))

summary(m3.8)


#Phase 2
#Note: this result is misleading because after the first period of phase 2, no treatment animals went near the VF
# Thus, the this refelcts data from only the first exclusion zone trial when animals were initially adapting to the new phase
m3.8.1 <- glmmTMB(
  cbind(total_pulses, total_audios - total_pulses) ~ TRT + 
    (1|Animal_ID), # Group rnadom effect was removed as it caused a model convergence issue (non-positive-definite Hessian matrix)
  family = binomial,
  data = filter(phase_summary, phase == "2: Exclusion Zone Testing"))

summary(m3.8.1)


#-----------------Rate of Learning by phase - Plots-------------------------------

# Reshape to long format
phase_long <- phase_summary %>%
  pivot_longer(cols = c(total_audios, total_pulses),
               names_to = "cue_type",
               values_to = "count") %>%
  mutate(
    cue_type = recode(cue_type,
                      total_audios = "Audios",
                      total_pulses = "Pulses"),
    TRT = factor(TRT, levels = c("CNT", "TRT")),
    phase = factor(phase, levels = c("1: Training", "2: Exclusion Zone Testing")))

# Summarize with mean and 95% CI
plot_df <- phase_long %>%
  group_by(phase, TRT, cue_type) %>%
  summarise(
    mean_count = mean(count),
    se = sd(count) / sqrt(n()),
    .groups = "drop") %>%
  mutate(
    ci_lower = mean_count - 1.96 * se,
    ci_upper = mean_count + 1.96 * se,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control"))


# Plot with error bars
p7 <- ggplot(plot_df, aes(x = cue_type, y = mean_count, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ phase) +
  labs(
    title = "Cue Counts by Treatment",
    x = "Cue Type",
    y = "Mean Cues per Animal",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Control" = "forestgreen", "Treatment" = "darkorange")) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.spacing = unit(2, "lines"))


#ggsave(filename = "Cue_Counts_phase.png",
#       plot = p7,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 7,
#       height = 5,
#       units = "in")




# Extract estimated probabilities from each model
emm1.1 <- emmeans(m3.8, ~ TRT, type = "response") %>%
  as.data.frame() %>%
  mutate(phase = "1: Training")

emm2.1 <- emmeans(m3.8.1, ~ TRT, type = "response") %>%
  as.data.frame() %>%
  mutate(phase = "2: Exclusion Zone Testing")

# Combine into one data frame
emm_combined <- bind_rows(emm1.1, emm2.1) %>%
  rename(prob = prob, lower = asymp.LCL, upper = asymp.UCL) %>%
  mutate(
    phase = factor(phase, levels = c("1: Training", "2: Exclusion Zone Testing")),
    TRT = factor(TRT, levels = c("CNT", "TRT")),
    treatment = ifelse(TRT == "TRT", "Treatment", "Control"))


# Plot the estimated probabilities with error bars
# Note: this plot is misleading because after the first period of phase 2, no treatment animals went near the VF
# Thus, the right panel of the plot refelcts data from only the first exclusion zone trial when animals were initially adapting to the new phase
p8 <- ggplot(emm_combined, aes(x = phase, y = prob, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.7), width = 0.2) +
  labs(
    title = "Probability of Escalation",
    x = "Phase",
    y = "Mean Probability Pulse Follows Audio",
    fill = "Treatment"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("Control" = "forestgreen", "Treatment" = "darkorange")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank())


#ggsave(filename = "Escalation_phase.png",
#       plot = p8,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 6,
#       height = 5,
#       units = "in")



###----------------Spatial Analyses-------------------------

#---------Training Phase - Proximity to VF------------------

#This is done on training phase only to avoid spatial-visual markers given the possible association between the hay bale and the VF for control groups. Visual cues during the exclusion zone phase function less to show fence location and more to signial which fence is active.

#Read in and clean data
exclusion_zones <- st_read(dsn = "C:/Users/spsch/Documents/R/Virtual_Fence/Exclusion_Zones.kml")

fencelines <- st_read(dsn = "C:/Users/spsch/Documents/R/Virtual_Fence/Fencelines.kml")

exclusion_zones$Description <- c("west", "east", "east", "east", "east", "east",
                                 "east", "west", "west", "west", "west", "west")


# Calculate distance from fence lines for each period
eshep_distance_sf <- eshep_full_sf %>%
  mutate(distanceVF1 = st_distance(eshep_full_sf, fencelines[1, ]),
         distanceVF2 = st_distance(eshep_full_sf, fencelines[2, ]),
         distanceVF3 = st_distance(eshep_full_sf, fencelines[3, ])) %>%
  mutate(distance_curr = ifelse(period == "Tr1", distanceVF1,
                                ifelse(period == "Tr2", distanceVF2,
                                       ifelse(period == "Tr3", distanceVF3, NA)))) %>%
  filter(!is.na(distance_curr))


# Add column specifying distance of 10 m
eshep_10m_df <- eshep_distance_sf %>%
  as.data.frame() %>%
  mutate(Within_10m = ifelse(distance_curr <= 10, 1, 0))



# Redefine all periods to use as random effect later
periods_tr <- c("Tr1", "Tr2", "Tr3")

# Create full grid of combinations to retain zeros
full_grid_10m <- expand_grid(
  groups,
  period = periods_tr)


# Summarize points within 10m by animal × period
Points_within_10m <- eshep_10m_df %>%
  group_by(Animal_ID, TRT, Group, period) %>%
  summarise(
    total_points_in_10m = sum(Within_10m, na.rm = TRUE),
    .groups = "drop")


# Join with full grid to retain zeroes!
points_within_10m_summary <- full_grid_10m %>%
  left_join(Points_within_10m, by = c("Animal_ID", "TRT", "Group", "period")) %>%
  mutate(total_points_in_10m = replace_na(total_points_in_10m, 0))




# Fit GLMM Poisson
m4 <- glmmTMB(
  total_points_in_10m ~ TRT + (1|Animal_ID) + (1|Group) + (1|period),
  family = poisson,
  data = points_within_10m_summary)

summary(m4)


# Model validation and over dispersion check
check_overdispersion(m4)

res4 <- simulateResiduals(fittedModel = m4, n = 1000)
plot(res4)
testDispersion(res4)



#--------- Proximity to VF - Plot ------------------

# Summarize mean and standard error by treatment
# Summarize mean and 95% confidence interval by treatment
summary_df <- points_within_10m_summary %>%
  group_by(TRT) %>%
  summarise(
    mean_points = mean(total_points_in_10m),
    se = sd(total_points_in_10m) / sqrt(n()),
    n = n()
  ) %>%
  mutate(
    ci_lower = mean_points - 1.96 * se,
    ci_upper = mean_points + 1.96 * se,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control"))


# Plot
p9 <- ggplot(summary_df, aes(x = treatment, y = mean_points, fill = treatment)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(
    title = "Proximity to Virtual Fence during Training",
    y = "Mean Points within 10m of VF per Period",
    x = element_blank()
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.4),
    legend.position = "none") +
  scale_fill_manual(values = c("Control" = "forestgreen", "Treatment" = "darkorange"))
  

#ggsave(filename = "Points_within_10m.png",
#       plot = p9,
#       path = "C:/Users/spsch/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 5,
#       height = 5,
#       units = "in")




#-------Testing Phase - Points within "open" hay bale-------

# Add columns based on whether points occur within open hay bales
inside_open_ex_points <- eshep_full_sf_ex %>%
  mutate(inside_open_ex = ifelse(period %in% c("E1", "E2") & west_ex_inside == TRUE, TRUE,
                                 ifelse(period %in% c("W1", "W2") & east_ex_inside == TRUE, TRUE, NA))) %>%
  filter(inside_open_ex == TRUE) %>%
  as.data.frame()


# Redefine all periods to use as random effect later
periods_ex <- c("E1", "E2", "W1", "W2")

# Create full grid of combinations to retain zeros
full_grid_open_ex <- expand_grid(
  groups,
  period = periods_ex)


# Summarize points by animal × period
inside_open_ex_points <- inside_open_ex_points %>%
  group_by(Animal_ID, TRT, Group, period) %>%
  summarise(
    total_points_inside_ex = sum(inside_open_ex, na.rm = TRUE),
    .groups = "drop")


# Join with full grid to retain zeroes!
inside_open_ex_points_summary <- full_grid_open_ex %>%
  left_join(inside_open_ex_points, by = c("Animal_ID", "TRT", "Group", "period")) %>%
  mutate(total_points_inside_ex = replace_na(total_points_inside_ex, 0))



# Fit GLMM Poisson
m5 <- glmmTMB(
  total_points_inside_ex ~ TRT + (1|Animal_ID) + (1|Group) + (1|period),
  family = poisson,
  data = inside_open_ex_points_summary)

summary(m5)


# Model validation and over dispersion check
check_overdispersion(m5)

res5 <- simulateResiduals(fittedModel = m5, n = 1000)
plot(res5)
testDispersion(res5)


#########    Time interaction????? ##########



#-----------------Heatmaps for hay bale use-----------------

