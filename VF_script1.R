packages_to_load <- c("tidyverse", "glmmTMB", "emmeans", "ggplot2", "sf", "terra", "lubridate", "lme4", "performance", "DHARMa", "survival", "survminer", "coxme") # create vector of R package names that we know are needed in the rest of the code


## ----load_libraries----
lapply(packages_to_load, library, character.only = TRUE)


##-----------------------------VF Escapes and Effectiveness-----------------------------

#-----Escapes - GLMM----
groups <- read.csv("C:/Users/Sebastian/Documents/R/Virtual_Fence/Raw Data/Groups.csv")

esc <- read.csv(file = "C:/Users/Sebastian/Documents/R/Virtual_Fence/Working Data/Escapes2.csv")


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
  Total_escapes ~ Escape_Type + (1|Group) + (1|Animal_ID),
  family = poisson,
  data = esc_summary_df)

summary(m1)

# Rate ratio
emmeans(m1, pairwise ~ Escape_Type, type = "response")


# Model Validation

sim_res <- simulateResiduals(m1, n = 1000)  # parametric sims
plot(sim_res)                               # uniformity, QQ, residuals vs fitted
testUniformity(sim_res)
testDispersion(sim_res)                     # over/under-dispersion
testZeroInflation(sim_res)




#-----Escapes - Bar plot-----

# Get estimated mean rates (on response scale) and convert to data frame
emm_df <- as.data.frame(emmeans(m1, ~ Escape_Type, type = "response"))

# Plot: Dot plot with error bars and log scale
# First, ensure Escape_Type is a factor with the correct order
emm_df$Escape_Type <- factor(emm_df$Escape_Type, levels = c("Electric", "Virtual"))

# Plot with bracket and asterisk
p1 <- ggplot(emm_df, aes(x = Escape_Type, y = rate, color = Escape_Type)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  scale_y_log10(
    name = "EMM Escape Rate (escapes per animal-day)",
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  scale_color_manual(values = c("Electric" = "firebrick", "Virtual" = "steelblue")) +
  
  # Bracket-style annotation
  geom_segment(aes(x = 1, xend = 1, y = 2.2, yend = 2.5), color = "black") +  # left leg
  geom_segment(aes(x = 2, xend = 2, y = 2.2, yend = 2.5), color = "black") +  # right leg
  geom_segment(aes(x = 1, xend = 2, y = 2.5, yend = 2.5), color = "black") +  # top line
  annotate("text", x = 1.5, y = 2.6, label = "***", size = 8) +                # asterisk above
  
  labs(
    x = "Fence Type",
    title = "Escapes by Fence Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )



ggsave(filename = "EscapesLog.png",
       plot = p1,
       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
       dpi = 600,
       width = 3.5,
       height = 7,
       units = "in")






###-----Effectiveness: Percentage Points In/Out-----

#eshep_df1 <- read.csv(file = "C:/Users/Sebastian/Documents/R/Virtual_Fence/2025-07-13.csv")
#eshep_df2 <- read.csv(file = "C:/Users/Sebastian/Documents/R/Virtual_Fence/2025-07-20.csv")

#eshep_df <- rbind.data.frame(eshep_df1, eshep_df2) %>% # Filter to include only the three training and four testing days
#  filter(Time..UTC. > "2025-07-07 11:55:00" & Time..UTC. < "2025-07-14 12:05:00")

#eshep_full_df <- eshep_df %>%
#  left_join(groups, by = "Neckband.ID") %>%
#  filter(!is.na(Animal_ID))

#write.csv(eshep_full_df, file = "C:/Users/Sebastian/Documents/R/Virtual_Fence/eshep_full_df.csv")




#Read in and clean data-----------------------

# Read in data
eshep_full_df <- read.csv(
  "C:/Users/Sebastian/Documents/R/Virtual_Fence/Working Data/eshep_full_df.csv")

exclusion_zones <- st_read(
  "C:/Users/Sebastian/Documents/R/Virtual_Fence/Raw Data/KMLs/Exclusion_Zones.kml")

training_VPs <- st_read(
  "C:/Users/Sebastian/Documents/R/Virtual_Fence/Raw Data/KMLs/Training VPs.kml")

perimeter_buffer <- st_read(
  "C:/Users/Sebastian/Documents/R/Virtual_Fence/Raw Data/KMLs/perimeter_buffer.kml")

# Projected CRS in meters
target_crs <- 32612

exclusion_zones <- st_transform(exclusion_zones, target_crs)
training_VPs <- st_transform(training_VPs, target_crs)
perimeter_buffer <- st_transform(perimeter_buffer, target_crs)

exclusion_zones$Description <- c(
  "west", "east", "east", "east", "east", "east",
  "east", "west", "west", "west", "west", "west")


# Remove missing GPS fixes
eshep_full_df_noNA <- eshep_full_df %>%
  filter(!is.na(latitude), !is.na(longitude))

# Convert GPS to sf
# IMPORTANT: longitude/latitude are EPSG 4326 before transforming
eshep_full_sf <- st_as_sf(
  eshep_full_df_noNA,
  coords = c("longitude", "latitude"),
  crs = 4326) %>%
  st_transform(target_crs) %>%
  mutate(
    Time..UTC. = ymd_hms(Time..UTC.),
    period = case_when(
      Time..UTC. <= ymd_hms("2025-07-08 11:59:59") ~ "Tr1",
      Time..UTC. <= ymd_hms("2025-07-09 11:59:59") ~ "Tr2",
      Time..UTC. <= ymd_hms("2025-07-10 11:59:59") ~ "Tr3",
      Time..UTC. <= ymd_hms("2025-07-11 11:59:59") ~ "E1",
      Time..UTC. <= ymd_hms("2025-07-12 11:59:59") ~ "W1",
      Time..UTC. <= ymd_hms("2025-07-13 11:59:59") ~ "E2",
      Time..UTC. <= ymd_hms("2025-07-14 12:00:00") ~ "W2",
      TRUE ~ NA_character_))

# Keep only points inside perimeter
eshep_full_sf <- eshep_full_sf %>%
  mutate(
    in_perimeter = st_within(eshep_full_sf, perimeter_buffer, sparse = FALSE)[ ,1]) %>%
  filter(in_perimeter == TRUE)

# Buffers for GPS-error tolerant containment
# Training VP expanded: lenient for near-boundary points
# Exclusion zone contracted: lenient for near-boundary points
training_VPs_buff <- st_buffer(training_VPs, dist = 3)
exclusion_zones_buff <- st_buffer(exclusion_zones, dist = -3)


# Raw training classifications
training1_inside_raw <- st_within(eshep_full_sf, training_VPs[1, ], sparse = FALSE)[, 1]
training2_inside_raw <- st_within(eshep_full_sf, training_VPs[2, ], sparse = FALSE)[, 1]
training3_inside_raw <- st_within(eshep_full_sf, training_VPs[3, ], sparse = FALSE)[, 1]


# Raw exclusion classifications
west_matrix_raw <- st_within(eshep_full_sf, filter(exclusion_zones, Description == "west"), sparse = FALSE)
west_ex_inside_raw <- apply(west_matrix_raw, 1, any)

east_matrix_raw <- st_within(eshep_full_sf, filter(exclusion_zones, Description == "east"), sparse = FALSE)
east_ex_inside_raw <- apply(east_matrix_raw, 1, any)

# Buffered training classifications
training1_inside_buff <- st_within(eshep_full_sf, training_VPs_buff[1, ], sparse = FALSE)[, 1]
training2_inside_buff <- st_within(eshep_full_sf, training_VPs_buff[2, ], sparse = FALSE)[, 1]
training3_inside_buff <- st_within(eshep_full_sf, training_VPs_buff[3, ], sparse = FALSE)[, 1]

# Buffered exclusion classifications
west_matrix_buff <- st_within(eshep_full_sf, filter(exclusion_zones_buff, Description == "west"), sparse = FALSE)
west_ex_inside_buff <- apply(west_matrix_buff, 1, any)

east_matrix_buff <- st_within(eshep_full_sf, filter(exclusion_zones_buff, Description == "east"), sparse = FALSE)
east_ex_inside_buff <- apply(east_matrix_buff, 1, any)

# Create full classified data set
eshep_full_sf_ex <- eshep_full_sf %>%
  mutate(
    inside_VF_raw = case_when(
      period == "Tr1" ~ training1_inside_raw,
      period == "Tr2" ~ training2_inside_raw,
      period == "Tr3" ~ training3_inside_raw,
      period %in% c("E1", "E2") ~ !east_ex_inside_raw,
      period %in% c("W1", "W2") ~ !west_ex_inside_raw,
      TRUE ~ NA),
    inside_VF_buff = case_when(
      period == "Tr1" ~ training1_inside_buff,
      period == "Tr2" ~ training2_inside_buff,
      period == "Tr3" ~ training3_inside_buff,
      period %in% c("E1", "E2") ~ !east_ex_inside_buff,
      period %in% c("W1", "W2") ~ !west_ex_inside_buff,
      TRUE ~ NA))


# Buffered containment summary
m2_df_buff <- eshep_full_sf_ex %>%
  st_drop_geometry() %>%
  group_by(Animal_ID, TRT, Group, period) %>%
  summarise(
    n_pts = sum(!is.na(inside_VF_buff)),
    n_in = sum(inside_VF_buff, na.rm = TRUE),
    prop_in = n_in / n_pts,
    .groups = "drop")

# Raw containment summary
m2_df_raw <- eshep_full_sf_ex %>%
  st_drop_geometry() %>%
  group_by(Animal_ID, TRT, Group, period) %>%
  summarise(
    n_pts = sum(!is.na(inside_VF_raw)),
    n_in = sum(inside_VF_raw, na.rm = TRUE),
    prop_in = n_in / n_pts,
    .groups = "drop")



# Fit binomial GLMM
# Animal_ID random effect was removed because it caused model singularity. Random effect estimate for animal_ID was 2.462e-34
m2_buff <- glmmTMB(
  cbind(n_in, n_pts - n_in) ~ TRT + (1|Group) + (1|period), 
  family = binomial,
  data = m2_df_buff)

summary(m2_buff)


# gives probability inside
emmeans(m2_buff, pairwise ~ TRT, type = "response")



# Fit binomial GLMM
# Animal_ID random effect was removed because it caused model singularity. Random effect estimate for animal_ID was 2.462e-34
m2_raw <- glmmTMB(
  cbind(n_in, n_pts - n_in) ~ TRT + (1|Group) + (1|period), 
  family = binomial,
  data = m2_df_raw)

summary(m2_raw)


# gives probability inside
emmeans(m2_raw, pairwise ~ TRT, type = "response")







sim_res <- simulateResiduals(m2, n = 1000)  # parametric sims
plot(sim_res)                               # uniformity, QQ, residuals vs fitted
testUniformity(sim_res)
testDispersion(sim_res)                     # over/under-dispersion
testZeroInflation(sim_res)                  # usually not an issue for binomial, but quick to check










sim_res <- simulateResiduals(m2, n = 1000)  # parametric sims
plot(sim_res)                               # uniformity, QQ, residuals vs fitted
testUniformity(sim_res)
testDispersion(sim_res)                     # over/under-dispersion
testZeroInflation(sim_res)                  # usually not an issue for binomial, but quick to check







#----- Bar Plot-------

emm_df <- as.data.frame(emm2$emmeans) %>%
  mutate(
    SE = SE,
    ymin = asymp.LCL,
    ymax = asymp.UCL,
    treatment = ifelse(TRT == "CNT", "Control", "Treatment"))



p2 <- ggplot(emm_df, aes(x = treatment, y = prob, fill = treatment)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, linetype = "dashed") +
    # Bracket-style significance annotation
  geom_segment(aes(x = 1, xend = 1, y = 1.0005, yend = 1.005), color = "black") +  # left leg
  geom_segment(aes(x = 2, xend = 2, y = 1.0005, yend = 1.005), color = "black") +  # right leg
  geom_segment(aes(x = 1, xend = 2, y = 1.005, yend = 1.005), color = "black") +   # top line
  annotate("text", x = 1.5, y = 1.006, label = "**", size = 8) +                   # asterisk
  labs(
    x = element_blank(),
    y = "EMM % Points Inside Virtual Boundary",
    title = "VF Effectiveness by Treatment") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0.80, 1.004)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("forestgreen", "darkorange"))



ggsave(filename = "Effectiveness.png",
      plot = p2,
       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
       dpi = 600,
       width = 4,
       height = 7,
       units = "in")






##-----------------VF Interactions by time-------------------------------

# Remove any points with cues that are not near VF

# Read in shapefile for buffer zones
buffers <- st_read(dsn = "C:/Users/Sebastian/Documents/R/Virtual_Fence/Raw Data/KMLs/Buffers.kml")

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
    hour_bin = floor(hour),
    phase = ifelse(period %in% c("Tr1", "Tr2", "Tr3"), "Training", "EZT")
  )


# Define all periods and hour bins
periods <- c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2")
hour_bins <- 0:23

# Create full grid of combinations
full_grid <- expand_grid(
  groups,
  period = periods,
  hour_bin = hour_bins)

full_grid <- full_grid %>%
  mutate(phase = ifelse(period %in% c("Tr1", "Tr2", "Tr3"), "Training", "EZT"))


# Summarize cue data by animal × period × hour
cue_summary <- eshep_cue_df_clean %>%
  group_by(Animal_ID, TRT, Group, period, hour_bin, phase) %>%
  summarise(
    total_audios = sum(`No..Audios`, na.rm = TRUE),
    total_pulses = sum(`No..Pulses`, na.rm = TRUE),
    .groups = "drop")


# Join with full grid to retain zeroes! Retaining zeros is highly important.
VF_interaction_summary <- full_grid %>%
  left_join(cue_summary, by = c("Animal_ID", "TRT", "Group", "period", "hour_bin", "phase")) %>%
  mutate(
    total_audios = replace_na(total_audios, 0),
    total_pulses = replace_na(total_pulses, 0))

VF_interaction_summary <- VF_interaction_summary %>%
  mutate(phase = factor(phase, levels = c("Training", "EZT")),
         day_in_phase = case_when(
           period == "Tr1" ~ 1,
           period == "Tr2" ~ 2,
           period == "Tr3" ~ 3,
           period == "E1"  ~ 1,
           period == "W1"  ~ 2,
           period == "E2"  ~ 3,
           period == "W2"  ~ 4),
    day_in_phase = as.numeric(day_in_phase),
    hour_bin = as.numeric(hour_bin))



#Statistical Modeling conducted in SAS--------------------------------

#Models conducted separately for each phase

#VF interactions were modeled cumulatively within each period with a 3-parameter logistic curve (with 24 hour bins as continuous variable).

#Logistic curve parameters were modeled linearly across periods (with the 3-4 periods acting as continuous variable) with TRT as a fixed effect also.




##------------Responsiveness--------------------------------------------------


#Originally audio-shock ratio or the percentage of cues that were audio was modeled. But this is misleading because that percentage has a lower bound of 50%, since animals cannot receive an shock without an audio warning first. Thus, shocks are a subset of audio events, and a better framing of the question is: Given an audio, what’s the probability it was followed by a pulse?

#To model the probability that an audio cue leads to a pulse — i.e., how often a cue escalates from an audio-only warning to an audio+shock correction.

#This is best framed as a conditional probability:

#Given that an audio was delivered, what’s the chance it was followed by a pulse?

#This model uses a binomial response in the form cbind(successes, failures).

#In this case: Successes = total_audios - total_pulses: how many audios did not lead to a pulse

#Failures = total_pulses: how many audios lead to a pulse

#Total trials = total_audios

#The audios are the "trials": they occur first and may or may not escalate.

#The pulses are the "successes": they occur only if the audio was ineffective.

#This model structure treats each hour (or bin) as an opportunity to observe that escalation rate.

#Response = cbind(total_audios - total_pulses, total_pulses)




#Statistical modeling conducted in SAS--------------------------------

#Models conducted separately for phase 1
#Insufficient data to fit a model for phase 2

#Responsiveness was modeled linearly by hour and period with TRT as a fixed effect also.








#-----------------Audio Plot-------------------------------

# Compile observed data for 1-hr bins
obs_audio <- VF_interaction_summary %>%
  mutate(
    study_day = case_when(
      phase == "Training" & day_in_phase == 1 ~ 1,
      phase == "Training" & day_in_phase == 2 ~ 2,
      phase == "Training" & day_in_phase == 3 ~ 3,
      phase == "EZT"      & day_in_phase == 1 ~ 4,
      phase == "EZT"      & day_in_phase == 2 ~ 5,
      phase == "EZT"      & day_in_phase == 3 ~ 6,
      phase == "EZT"      & day_in_phase == 4 ~ 7),
    study_hour = (study_day - 1) * 24 + hour_bin,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control")) %>%
  group_by(treatment, study_hour) %>%
  summarise(
    obs_mean = mean(total_audios, na.rm = TRUE),
    .groups = "drop")


day_breaks <- seq(24, 144, by = 24)
phase_break <- 72
y_top <- 0.6

time_breaks <- seq(0, max(obs_audio$study_hour, na.rm = TRUE), by = 8)

time_labels <- rep(c("12PM", "8PM", "4AM"),
                  length.out = length(time_breaks))

p_audio_obs <- ggplot(
  obs_audio,
  aes(x = study_hour, y = obs_mean, color = treatment)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.2,
    span = 0.16,
    method = "loess") +
  geom_vline(xintercept = day_breaks, linetype = "dotted", alpha = 0.4) +
  geom_vline(xintercept = phase_break, linetype = "solid", linewidth = 0.8) +
  annotate("text", x = 36, y = y_top, label = "Training") +
  annotate("text", x = 120, y = y_top, label = "EZT") +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120, 144),
    labels = c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2"),
    sec.axis = dup_axis(
      breaks = time_breaks,
      labels = time_labels,
      name = element_blank())) +
  labs(
    title = "Audios by Treatment Across Time",
    y = "Audios per Animal per Hour",
    x = "Period of Study") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))


# 3-hr bins
obs_audio_3hr <- VF_interaction_summary %>%
  mutate(
    study_day = case_when(
      phase == "Training" & day_in_phase == 1 ~ 1,
      phase == "Training" & day_in_phase == 2 ~ 2,
      phase == "Training" & day_in_phase == 3 ~ 3,
      phase == "EZT"      & day_in_phase == 1 ~ 4,
      phase == "EZT"      & day_in_phase == 2 ~ 5,
      phase == "EZT"      & day_in_phase == 3 ~ 6,
      phase == "EZT"      & day_in_phase == 4 ~ 7),
    study_hour = (study_day - 1) * 24 + hour_bin,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control")) %>%
  mutate(hour_bin_3 = floor(hour_bin / 3) * 3,
         study_hour_3 = (study_day - 1) * 24 + hour_bin_3) %>%
  group_by(treatment, phase, day_in_phase, study_day, study_hour_3) %>%
  summarise(
    mean_audio = mean(total_audios, na.rm = TRUE),
    .groups = "drop")


p_audio_obs_3hr <- ggplot(
  obs_audio_3hr,
  aes(x = study_hour_3, y = mean_audio, color = treatment)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.2,
    span = 0.16,
    method = "loess") +
  geom_vline(xintercept = day_breaks, linetype = "dotted", alpha = 0.4) +
  geom_vline(xintercept = phase_break, linetype = "solid", linewidth = 0.8) +
  annotate("text", x = 36, y = y_top, label = "Training") +
  annotate("text", x = 120, y = y_top, label = "EZT") +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120, 144),
    labels = c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2"),
    sec.axis = dup_axis(
      breaks = time_breaks,
      labels = time_labels,
      name = element_blank())) +
  labs(
    title = "Audios by Treatment Across Time",
    y = "Audios per Animal per Hour",
    x = "Period of Study") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))


#ggsave(filename = "audio_3hr.png",
#       plot = p_audio_obs_3hr,
#      path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
#       dpi = 600,
#       width = 12,
#       height = 3,
#       units = "in")


#-----------------Pulse Plot-------------------------------

# Compile observed data for 1-hr bins
obs_pulse <- VF_interaction_summary %>%
  mutate(
    study_day = case_when(
      phase == "Training" & day_in_phase == 1 ~ 1,
      phase == "Training" & day_in_phase == 2 ~ 2,
      phase == "Training" & day_in_phase == 3 ~ 3,
      phase == "EZT"      & day_in_phase == 1 ~ 4,
      phase == "EZT"      & day_in_phase == 2 ~ 5,
      phase == "EZT"      & day_in_phase == 3 ~ 6,
      phase == "EZT"      & day_in_phase == 4 ~ 7),
    study_hour = (study_day - 1) * 24 + hour_bin,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control")) %>%
  group_by(treatment, study_hour) %>%
  summarise(
    obs_mean = mean(total_pulses, na.rm = TRUE),
    .groups = "drop")


day_breaks <- seq(24, 144, by = 24)
phase_break <- 72
y_top <- 0.32

time_breaks <- seq(0, max(obs_audio$study_hour, na.rm = TRUE), by = 8)

time_labels <- rep(c("12PM", "8PM", "4AM"),
                   length.out = length(time_breaks))

p_pulse_obs <- ggplot(
  obs_pulse,
  aes(x = study_hour, y = obs_mean, color = treatment)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.2,
    span = 0.16,
    method = "loess") +
  geom_vline(xintercept = day_breaks, linetype = "dotted", alpha = 0.4) +
  geom_vline(xintercept = phase_break, linetype = "solid", linewidth = 0.8) +
  annotate("text", x = 36, y = y_top, label = "Training") +
  annotate("text", x = 120, y = y_top, label = "EZT") +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120, 144),
    labels = c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2"),
    sec.axis = dup_axis(
      breaks = time_breaks,
      labels = time_labels,
      name = element_blank())) +
  labs(
    title = "Pulses by Treatment Across Time",
    y = "Pulses per Animal per Hour",
    x = "Period of Study") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))


# 3-hr bins
obs_pulses_3hr <- VF_interaction_summary %>%
  mutate(
    study_day = case_when(
      phase == "Training" & day_in_phase == 1 ~ 1,
      phase == "Training" & day_in_phase == 2 ~ 2,
      phase == "Training" & day_in_phase == 3 ~ 3,
      phase == "EZT"      & day_in_phase == 1 ~ 4,
      phase == "EZT"      & day_in_phase == 2 ~ 5,
      phase == "EZT"      & day_in_phase == 3 ~ 6,
      phase == "EZT"      & day_in_phase == 4 ~ 7),
    study_hour = (study_day - 1) * 24 + hour_bin,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control")) %>%
  mutate(hour_bin_3 = floor(hour_bin / 3) * 3,
         study_hour_3 = (study_day - 1) * 24 + hour_bin_3) %>%
  group_by(treatment, phase, day_in_phase, study_day, study_hour_3) %>%
  summarise(
    mean_pulses = mean(total_pulses, na.rm = TRUE),
    .groups = "drop")


p_pulse_obs_3hr <- ggplot(
  obs_pulses_3hr,
  aes(x = study_hour_3, y = mean_pulses, color = treatment)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.2,
    span = 0.16,
    method = "loess") +
  geom_vline(xintercept = day_breaks, linetype = "dotted", alpha = 0.4) +
  geom_vline(xintercept = phase_break, linetype = "solid", linewidth = 0.8) +
  annotate("text", x = 36, y = 0.2, label = "Training") +
  annotate("text", x = 120, y = 0.2, label = "EZT") +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120, 144),
    labels = c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2"),
    sec.axis = dup_axis(
      breaks = time_breaks,
      labels = time_labels,
      name = element_blank())) +
  labs(
    title = "Pulses by Treatment Across Time",
    y = "Pulses per Animal per Hour",
    x = "Period of Study") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))


#ggsave(filename = "pulse_3hr.png",
#       plot = p_pulse_obs_3hr,
#       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
#       dpi = 600,
#       width = 12,
#       height = 3,
#       units = "in")




#-----------------Responsiveness Plot-------------------------------
# Compile observed data for 1-hr bins
obs_resp <- VF_interaction_summary %>%
  mutate(
    study_day = case_when(
      phase == "Training" & day_in_phase == 1 ~ 1,
      phase == "Training" & day_in_phase == 2 ~ 2,
      phase == "Training" & day_in_phase == 3 ~ 3,
      phase == "EZT"      & day_in_phase == 1 ~ 4,
      phase == "EZT"      & day_in_phase == 2 ~ 5,
      phase == "EZT"      & day_in_phase == 3 ~ 6,
      phase == "EZT"      & day_in_phase == 4 ~ 7),
    study_hour = (study_day - 1) * 24 + hour_bin,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control"),
    responsiveness = if_else(
      total_audios > 0,
      (total_audios - total_pulses) / total_audios,
      NA_real_)) %>%
  group_by(treatment, study_hour) %>%
  summarise(
    obs_mean = mean(responsiveness, na.rm = TRUE),
    .groups = "drop")

day_breaks <- seq(24, 144, by = 24)
phase_break <- 72
y_top <- 0.9

time_breaks <- seq(0, max(obs_resp$study_hour, na.rm = TRUE), by = 8)

time_labels <- rep(c("12PM", "8PM", "4AM"),
                   length.out = length(time_breaks))

p_resp_obs <- ggplot(
  obs_resp,
  aes(x = study_hour, y = obs_mean, color = treatment)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.2,
    span = 0.5,
    method = "loess") +
  geom_vline(xintercept = day_breaks, linetype = "dotted", alpha = 0.4) +
  geom_vline(xintercept = phase_break, linetype = "solid", linewidth = 0.8) +
  annotate("text", x = 36, y = y_top, label = "Training") +
  annotate("text", x = 120, y = y_top, label = "EZT") +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_y_continuous(
    limits = c(0, 1.05),
    labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120, 144),
    labels = c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2"),
    sec.axis = dup_axis(
      breaks = time_breaks,
      labels = time_labels,
      name = element_blank())) +
  labs(
    title = "Responsiveness by Treatment Across Time",
    y = "Responsivness to Audio Warning",
    x = "Period of Study") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))



# Compile observed data for 3-hr bins
obs_resp_3hr <- VF_interaction_summary %>%
  mutate(
    study_day = case_when(
      phase == "Training" & day_in_phase == 1 ~ 1,
      phase == "Training" & day_in_phase == 2 ~ 2,
      phase == "Training" & day_in_phase == 3 ~ 3,
      phase == "EZT"      & day_in_phase == 1 ~ 4,
      phase == "EZT"      & day_in_phase == 2 ~ 5,
      phase == "EZT"      & day_in_phase == 3 ~ 6,
      phase == "EZT"      & day_in_phase == 4 ~ 7),
    study_hour = (study_day - 1) * 24 + hour_bin,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control"),
    responsiveness = if_else(
      total_audios > 0,
      (total_audios - total_pulses) / total_audios,
      NA_real_),
    hour_bin_3 = floor(hour_bin / 3) * 3,
    study_hour_3 = (study_day - 1) * 24 + hour_bin_3) %>%
  group_by(treatment, phase, day_in_phase, study_day, study_hour_3) %>%
  summarise(
    mean_resp = mean(responsiveness, na.rm = TRUE),
    .groups = "drop")

p_resp_obs_3hr <- ggplot(
  obs_resp_3hr,
  aes(x = study_hour_3, y = mean_resp, color = treatment)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.2,
    span = 0.5,
    method = "loess") +
  geom_vline(xintercept = day_breaks, linetype = "dotted", alpha = 0.4) +
  geom_vline(xintercept = phase_break, linetype = "solid", linewidth = 0.8) +
  annotate("text", x = 36, y = y_top, label = "Training") +
  annotate("text", x = 120, y = y_top, label = "EZT") +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_y_continuous(
    limits = c(0, 1.05),
    labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120, 144),
    labels = c("Tr1", "Tr2", "Tr3", "E1", "W1", "E2", "W2"),
    sec.axis = dup_axis(
      breaks = time_breaks,
      labels = time_labels,
      name = element_blank())) +
  labs(
    title = "Responsiveness by Treatment Across Time",
    y = "Responsivness",
    x = "Period of Study") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))


#ggsave(filename = "responsiveness_time.png",
#       plot = p_resp_obs,
#       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
#       dpi = 600,
#       width = 12.5,
#       height = 3,
#       units = "in")






###----------------Spatial Analyses-------------------------------------------------------------

#---------Training Phase - Proximity to VF------------------

#This is done on training phase only to avoid spatial-visual markers given the possible association between the hay bale and the VF for control groups. Visual cues during the exclusion zone phase function less to show fence location and more to signial which fence is active.

#Read in and clean data
exclusion_zones <- st_read(dsn = "C:/Users/Sebastian/Documents/R/Virtual_Fence/Raw Data/KMLs/Exclusion_Zones.kml")

fencelines <- st_read(dsn = "C:/Users/Sebastian/Documents/R/Virtual_Fence/Raw Data/KMLs/Fencelines.kml")

exclusion_zones$Description <- c("west", "east", "east", "east", "east", "east",
                                 "east", "west", "west", "west", "west", "west")


# Calculate distance from fence lines for each period
eshep_distance_sf <- eshep_full_sf %>%
  mutate(
    distanceVF1 = as.numeric(st_distance(geometry, fencelines[1, ])[, 1]),
    distanceVF2 = as.numeric(st_distance(geometry, fencelines[2, ])[, 1]),
    distanceVF3 = as.numeric(st_distance(geometry, fencelines[3, ])[, 1])
  ) %>%
  mutate(
    distance_curr = case_when(
      period == "Tr1" ~ distanceVF1,
      period == "Tr2" ~ distanceVF2,
      period == "Tr3" ~ distanceVF3,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(distance_curr)) %>%
  mutate(
    datetime = ymd_hms(`Time..UTC.`),
    hour = hour(datetime) +
      minute(datetime) / 60 +
      second(datetime) / 3600,
    hour = hour - 12,
    hour = if_else(hour < 0, hour + 24, hour),
    hour_bin = floor(hour))

Points_within_10m <- eshep_distance_sf %>%
  st_drop_geometry() %>%
  mutate(Within_10m = if_else(distance_curr <= 10, 1L, 0L)) %>%
  group_by(Animal_ID, TRT, Group, period, hour_bin) %>%
  summarise(
    total_points_in_10m = sum(Within_10m, na.rm = TRUE),
    .groups = "drop")


# Redefine all periods to use as random effect later
periods_tr <- c("Tr1", "Tr2", "Tr3")
hour_bin <- c(0:23)

# Create full grid of combinations to retain zeros
full_grid_10m <- expand_grid(
  groups,
  hour_bin,
  period = periods_tr)


# Join with full grid to retain zeroes!
points_within_10m_summary <- full_grid_10m %>%
  left_join(Points_within_10m, by = c("Animal_ID", "TRT", "Group", "hour_bin", "period")) %>%
  mutate(total_points_in_10m = replace_na(total_points_in_10m, 0),
         day = as.numeric(case_when(period == "Tr1" ~ 1,
                         period == "Tr2" ~ 2,
                         period == "Tr3" ~ 3)))




# Fit GLMM Poisson
m4 <- glmmTMB(
  total_points_in_10m ~ TRT * day + TRT * hour_bin + (1|Group/Animal_ID),
  family = nbinom2,
  data = points_within_10m_summary)

summary(m4)

emmeans(m4, pairwise ~ TRT, type = "response")



# Model validation and over dispersion check
check_overdispersion(m4)

res4 <- simulateResiduals(fittedModel = m4, n = 1000)
plot(res4)
testDispersion(res4)




#-------Testing Phase - Points within "open" hay bale-------



# Add columns based on whether points occur within open hay bales
inside_open_ex_points <- eshep_full_sf_ex %>%
  mutate(inside_open_ex = ifelse(period %in% c("E1", "E2") & west_ex_inside == TRUE, TRUE,
                                 ifelse(period %in% c("W1", "W2") & east_ex_inside == TRUE, TRUE, NA))) %>%
  filter(inside_open_ex == TRUE) %>%
  as.data.frame()


# Process timestamps and bin hours
inside_open_ex_points <- inside_open_ex_points %>%
  mutate(
    datetime = ymd_hms(`Time..UTC.`),
    hour = hour(datetime) + minute(datetime) / 60 + second(datetime) / 3600,
    hour = hour - 12,
    hour = if_else(hour < 0, hour + 24, hour),
    hour_bin = floor(hour))


# Redefine all periods and hour bins
periods_ex <- c("E1", "E2", "W1", "W2")
hour_bins <- 0:23


# Create full grid of combinations to retain zeros
full_grid_open_ex <- expand_grid(
  groups,
  hour_bin = hour_bins,
  period = periods_ex)


# Summarize points by animal × period
inside_open_ex_points <- inside_open_ex_points %>%
  group_by(Animal_ID, TRT, Group, period, hour_bin) %>%
  summarise(
    total_points_inside_ex = sum(inside_open_ex, na.rm = TRUE),
    .groups = "drop")


# Join with full grid to retain zeroes!
inside_open_ex_points_summary <- full_grid_open_ex %>%
  left_join(inside_open_ex_points, by = c("Animal_ID", "TRT", "Group", "period", "hour_bin")) %>%
  mutate(total_points_inside_ex = replace_na(total_points_inside_ex, 0),
         day = as.numeric(case_when(period == "E1" ~ 1,
                                    period == "W1" ~ 2,
                                    period == "E2" ~ 3,
                                    period == "W2" ~ 4)))



# Fit GLMM Poisson
m5 <- glmmTMB(
  total_points_inside_ex ~ TRT * day + TRT * hour_bin + (1|Group/Animal_ID),
  family = nbinom2,
  data = inside_open_ex_points_summary)

summary(m5)

emmeans(m5, pairwise ~ TRT, type = "response")



# Model validation and over dispersion check
check_overdispersion(m5)

res5 <- simulateResiduals(fittedModel = m5, n = 1000)
plot(res5)
testDispersion(res5)



#--------- Spatial Analyses - Plots ------------------

#Training phase proximity to VF: 3-hour plot

obs_10m_3hr <- points_within_10m_summary %>%
  mutate(
    study_day = day,                       # Tr1 = 1, Tr2 = 2, Tr3 = 3
    hour_bin_3 = floor(hour_bin / 3) * 3,
    study_hour_3 = (study_day - 1) * 24 + hour_bin_3,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control")) %>%
  group_by(treatment, study_day, study_hour_3) %>%
  summarise(
    mean_points = mean(total_points_in_10m, na.rm = TRUE),
    .groups = "drop")

day_breaks_tr <- seq(24, 48, by = 24)
y_top_tr <- max(obs_10m_3hr$mean_points, na.rm = TRUE) * 1.05

time_breaks_tr <- seq(0, max(obs_10m_3hr$study_hour_3, na.rm = TRUE), by = 8)
time_labels_tr <- rep(c("12PM", "8PM", "4AM"),
                      length.out = length(time_breaks_tr))

p_10m_obs_3hr <- ggplot(
  obs_10m_3hr,
  aes(x = study_hour_3, y = mean_points, color = treatment)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.2,
    span = 0.4,
    method = "loess") +
  geom_vline(xintercept = day_breaks_tr, linetype = "dotted", alpha = 0.4) +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_x_continuous(
    breaks = c(0, 24, 48),
    labels = c("Tr1", "Tr2", "Tr3"),
    sec.axis = dup_axis(
      breaks = time_breaks_tr,
      labels = time_labels_tr,
      name = NULL)) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    title = "Training Phase: Points in Proximity to VF",
    y = "Points within 10m of VF per Animal per Hour",
    x = "Period of Study") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))



#ggsave(filename = "Points_in_10m.png",
#       plot = p_10m_obs_3hr,
#       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
#       dpi = 600,
#       width = 6.5,
#       height = 4,
#       units = "in")


# EZT phase open hay bale use: 3-hour plot

obs_open_3hr <- inside_open_ex_points_summary %>%
  mutate(
    study_day = day,                       # E1 = 1, W1 = 2, E2 = 3, W2 = 4
    hour_bin_3 = floor(hour_bin / 3) * 3,
    study_hour_3 = (study_day - 1) * 24 + hour_bin_3,
    treatment = ifelse(TRT == "TRT", "Treatment", "Control")) %>%
  group_by(treatment, study_day, study_hour_3) %>%
  summarise(
    mean_points = mean(total_points_inside_ex, na.rm = TRUE),
    .groups = "drop")

day_breaks_ex <- seq(24, 72, by = 24)
y_top_ex <- max(obs_open_3hr$mean_points, na.rm = TRUE) * 1.05

time_breaks_ex <- seq(0, max(obs_open_3hr$study_hour_3, na.rm = TRUE), by = 8)
time_labels_ex <- rep(c("12PM", "8PM", "4AM"),
                      length.out = length(time_breaks_ex))

p_open_obs_3hr <- ggplot(
  obs_open_3hr,
  aes(x = study_hour_3, y = mean_points, color = treatment)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_smooth(
    se = FALSE,
    linewidth = 1.2,
    span = 0.4,
    method = "loess") +
  geom_vline(xintercept = day_breaks_ex, linetype = "dotted", alpha = 0.4) +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72),
    labels = c("E1", "W1", "E2", "W2"),
    sec.axis = dup_axis(
      breaks = time_breaks_ex,
      labels = time_labels_ex,
      name = NULL)) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    title = "EZT Phase: Use of Open Hay Bale",
    y = "Points inside Open E.Z. per Animal per Hour",
    x = "Period of Study") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))


#ggsave(filename = "Points_in_open_bale.png",
#       plot = p_open_obs_3hr,
#       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
#      dpi = 600,
#       width = 8,
#       height = 4,
#       units = "in")



#-----------------Heatmaps for hay bale use-----------------

# Read and transform spatial layers
heatmap_shapes <- st_read("C:/Users/Sebastian/Documents/R/Virtual_Fence/heatmap_shapes.kml") %>%
  st_transform(32614) %>%
  mutate(pasture = c(1, 2, 3, 4, 5, 6))

eshep_full_proj <- st_transform(eshep_full_sf, crs = 32614)
exclusion_zones_proj <- st_transform(exclusion_zones, crs = 32614)


# Create hay bale centroids (retain geometry)
hay_bales_east <- exclusion_zones_proj %>%
  filter(Description == "east") %>%
  mutate(pasture = c(6, 5, 4, 3, 2, 1)) %>%
  st_centroid()

hay_bales_west <- exclusion_zones_proj %>%
  filter(Description == "west") %>%
  mutate(pasture = c(6, 1, 2, 3, 4, 5)) %>%
  st_centroid()


# Label pasture membership for GPS points
eshep_pts_west <- eshep_full_proj %>%
  filter(period %in% c("E1", "E2")) %>%
  st_join(heatmap_shapes, join = st_intersects)

eshep_pts_east <- eshep_full_proj %>%
  filter(period %in% c("W1", "W2")) %>%
  st_join(heatmap_shapes, join = st_intersects)


# Extract original GPS coordinates before changing geometry
eshep_pts_east <- eshep_pts_east %>%
  mutate(point_coords = st_coordinates(geometry))

eshep_pts_west <- eshep_pts_west %>%
  mutate(point_coords = st_coordinates(geometry))


# Join hay bale geometry to each GPS point by pasture ID
eshep_pts_east <- eshep_pts_east %>%
  left_join(st_drop_geometry(hay_bales_east), by = "pasture") %>%
  mutate(hay_coords = st_coordinates(st_geometry(hay_bales_east)[match(pasture, hay_bales_east$pasture)])) %>%
  filter(!is.na(pasture))

eshep_pts_west <- eshep_pts_west %>%
  left_join(st_drop_geometry(hay_bales_west), by = "pasture") %>%
  mutate(hay_coords = st_coordinates(st_geometry(hay_bales_west)[match(pasture, hay_bales_west$pasture)])) %>%
  filter(!is.na(pasture))


# Calculate shifted coordinates
eshep_pts_east <- eshep_pts_east %>%
  mutate(
    x_shifted = point_coords[,1] - hay_coords[,1],
    y_shifted = point_coords[,2] - hay_coords[,2])

eshep_pts_west <- eshep_pts_west %>%
  mutate(
    x_shifted = point_coords[,1] - hay_coords[,1],
    y_shifted = point_coords[,2] - hay_coords[,2])



# Plot the overlaid heatmap - East
heatmap_east_trt <- ggplot(filter(eshep_pts_east, TRT == "TRT"),
       aes(x = x_shifted, y = y_shifted)) +
  stat_density_2d_filled(
    contour_var = "density",
    adjust = 1.2,
    alpha = 0.9) +
  scale_fill_viridis_d(option = "C", direction = -1, name = "Density") +
  coord_equal() +
  labs(
    title = "East (Treatment)") +
  theme_minimal(base_size = 14)

# Plot the overlaid heatmap - West
heatmap_west_trt <- ggplot(filter(eshep_pts_west, TRT == "TRT"),
       aes(x = x_shifted, y = y_shifted)) +
  stat_density_2d_filled(
    contour_var = "density",
    adjust = 1.2,
    alpha = 0.9) +
  scale_fill_viridis_d(option = "C", direction = -1, name = "Density") +
  coord_equal() +
  labs(
    title = "West (Treatment)") +
  theme_minimal(base_size = 14)



# Plot the overlaid heatmap - East
heatmap_east_cnt <- ggplot(filter(eshep_pts_east, TRT == "CNT"),
       aes(x = x_shifted, y = y_shifted)) +
  stat_density_2d_filled(
    contour_var = "density",
    adjust = 1.2,
    alpha = 0.9) +
  scale_fill_viridis_d(option = "C", direction = -1, name = "Density") +
  coord_equal() +
  labs(
    title = "East (Control)") +
  theme_minimal(base_size = 14)

# Plot the overlaid heatmap - West
heatmap_west_cnt <- ggplot(filter(eshep_pts_west, TRT == "CNT"),
       aes(x = x_shifted, y = y_shifted)) +
  stat_density_2d_filled(
    contour_var = "density",
    adjust = 1.2,
    alpha = 0.9) +
  scale_fill_viridis_d(option = "C", direction = -1, name = "Density") +
  coord_equal() +
  labs(
    title = "West (Control)") +
  theme_minimal(base_size = 14)



#ggsave(filename = "heatmap_east_cnt.png",
#       plot = heatmap_east_cnt,
#       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
#       dpi = 300,
#       width = 8,
#       height = 6,
#       units = "in")




#-----------Conditioned Response Extinction--------------------------

extinction <- read.csv(file = "C:/Users/Sebastian/Documents/R/Virtual_Fence/Raw Data/Extinction.csv")[, 1:5]

groups <- groups %>%
  mutate(TRT = as.character(TRT))


# Define extinction start time
extinction_start <- mdy_hm("07/14/2025 12:00")


# Clean extinction data: assign reentry time
extinction_clean <- extinction %>%
  select(Animal_ID, Group, Day, Time) %>%
  mutate(
    DateTime = mdy_hm(paste(Day, Time)),
    reentered = 1) %>%
  select(Animal_ID, Group, DateTime, reentered)


# Join extinction records to full animal list, pulling TRT from groups only
surv_data <- groups %>%
  select(Animal_ID, Group, TRT) %>%
  left_join(extinction_clean, by = c("Animal_ID", "Group")) %>%
  mutate(
    # Calculate time to event or assign max time if no reentry
    time_to_reentry = as.numeric(difftime(DateTime, extinction_start, units = "mins")),
    reentered = ifelse(is.na(reentered), 0, reentered),
    time_to_reentry = ifelse(is.na(time_to_reentry), 1229, time_to_reentry))

surv_data <- surv_data %>%
  group_by(Animal_ID) %>%
  slice_min(time_to_reentry, n = 1, with_ties = FALSE) %>%
  ungroup()



## Survivor model is appropriate since there is not definitive time cut-off

# Fit survival model
surv_obj <- Surv(surv_data$time_to_reentry, surv_data$reentered)


# Fit Kaplan-Meier survival curves
km_fit <- survfit(surv_obj ~ TRT, data = surv_data)


# Fit a Cox Model with random effects to Estimate Effects
coxme_model <- coxme(Surv(time_to_reentry, reentered) ~ TRT + (1 | Group), data = surv_data)

summary(coxme_model)




#--- Survivor Plot -----------------------------------

surv_data <- surv_data %>%
  mutate(
    hours = time_to_reentry / 60,
    time_of_day = format(mdy_hm("07/14/2025 12:00") + minutes(time_to_reentry), "%I:%M %p"))


km_fit <- survfit(Surv(time_to_reentry, reentered) ~ TRT, data = surv_data)


# Convert KM model to data frame
km_df <- survminer::ggsurvplot(km_fit, data = surv_data, conf.int = FALSE)

km_df <- km_df$data.survplot


# Adjust labels and extract treatment
km_df <- km_df %>%
  mutate(
    hours = time / 60,
    time_of_day = format(mdy_hm("07/14/2025 12:00") + minutes(time), "%I:%M %p"),
    TRT = gsub("TRT=", "", strata))

# Update treatment labels
km_df <- km_df %>%
  mutate(Treatment = ifelse(TRT == "TRT", "Treatment", "Control"))



# Plot with readable time axis and proper coloring
survivor_plot <- ggplot(km_df, aes(x = hours, y = 1 - surv, color = Treatment)) +
  geom_line(size = 1.2, alpha = 0.8) +
  scale_y_continuous(
    name = "Proportion of Animals Re-entered",
    limits = c(0, .5),
    expand = c(0, 0)) +
  scale_x_continuous(
    name = "Hours since VF deactivation",
    breaks = seq(0, 24, by = 2),
    sec.axis = dup_axis(
      breaks = seq(0, 24, by = 4),  # Choose sensible numeric breaks
      labels = function(x) format(lubridate::mdy_hm("07/14/2025 12:00") + lubridate::hours(x), "%I %p"),
      name = "Time of Day")) +
  scale_color_manual(values = c("Treatment" = "darkorange", "Control" = "forestgreen")) +
  theme_minimal(base_size = 14) +
  labs(title = "Re-entry into Exclusion Zone by Treatment") +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))


ggsave(filename = "survivor_plot.png",
       plot = survivor_plot,
       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence",
       dpi = 300,
       width = 7,
       height = 4,
       units = "in")








###-----------------------------SAS Results and Summary Plots---------------------------------------------

Full_results <- readxl::read_xlsx("C:/Users/Sebastian/Documents/R/Virtual_Fence/Results/Combined_Results.xlsx",
                                  sheet = 1)

#------------Plot 1--------------

sig_lookup <- Full_results %>%
  dplyr::filter(
    Phase == "Training",
    !Period %in% c("1", "2", "3"),
    `Statistic Type` %in% c("LSMean Contrast", "Slope contrast") 
  ) %>%
  mutate(
    p_clean = str_replace(as.character(`P-value`), "<", ""),
    p_num = as.numeric(p_clean),
    sig = case_when(
      is.na(p_num)  ~ "",
      p_num < 0.001 ~ "***",
      p_num < 0.01  ~ "**",
      p_num < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(
    `Response Variable`,
    `Model Parameter`,
    `P-value`,
    sig)


Training_parms <- Full_results %>%
  filter(
    Phase == "Training",
    is.na(Period) | !as.character(Period) %in% c("1", "2", "3"),
    `Statistic Type` %in% c("LSMean", "Slope")
  ) %>%
  left_join(sig_lookup, by = c("Response Variable", "Model Parameter")) %>%
  mutate(
    `Response Variable` = factor(
      `Response Variable`,
      levels = c("Audios", "Pulses", "Points within 10m"),
      labels = c("Audio", "Pulse", "Within 10m")
    ),
    `Model Parameter` = factor(
      `Model Parameter`,
      levels = c("A", "r", "t50",
                 "Period slope for A", "Period slope for r", "Period slope for t50"),
      labels = c("A EMM", "r EMM", "t50 EMM",
                 "Period slope of A", "Period slope of r", "Period slope of t50")
    ),
    `Treatment Group` = factor(`Treatment Group`, levels = c("CNT", "TRT")),
    lower = Estimate - 1.96 * `Standard Error`,
    upper = Estimate + 1.96 * `Standard Error`
  )

sig_labels <- Training_parms %>%
  group_by(`Response Variable`, `Model Parameter`) %>%
  summarise(
    sig = first(sig),
    x = mean(Estimate),
    .groups = "drop"
  ) %>%
  filter(sig != "")



training_plot1 <- ggplot(
  Training_parms,
  aes(
    x = Estimate,
    y = `Treatment Group`,
    color = `Treatment Group`
  )
) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    height = 0.18,
    linewidth = 1
  ) +
  geom_point(size = 4) +
  facet_grid(
    `Response Variable` ~ `Model Parameter`,
    scales = "free",
    switch = "y"
  ) +
  scale_color_manual(
    values = c(
      "TRT" = "darkorange",
      "CNT" = "forestgreen"
    )
  ) +
  geom_text(
    data = sig_labels,
    aes(
      x = x,
      y = 1.5,
      label = sig
    ),
    inherit.aes = FALSE,
    size = 6,
    color = "black"
  ) +
  labs(
    x = "Slope estimate ± 95% CI",
    y = NULL,
    color = "Treatment",
    title = "Training Phase Slope Estimates"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 


training_plot1

ggsave(filename = "training_plot1_sig.png",
       plot = training_plot1,
       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence/Results",
       dpi = 600,
       width = 9,
       height = 4.5,
       units = "in")


#-----------Plot 2-----------

sig_lookup2 <- Full_results %>%
  dplyr::filter(
    `Response Variable` == "Responsiveness",
    `Statistic Type` %in% c("LSMean Contrast", "Period slope contrast")) %>%
  mutate(
    p_clean = str_replace(as.character(`P-value`), "<", ""),
    p_num = as.numeric(p_clean),
    sig = case_when(
      is.na(p_num)  ~ "",
      p_num < 0.001 ~ "***",
      p_num < 0.01  ~ "**",
      p_num < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(
    `Statistic Type`,
    Period,
    `P-value`,
    sig)


Training_resp_parms <- Full_results %>%
  dplyr::filter(
      `Response Variable` == "Responsiveness",
      `Statistic Type` %in% c("Period slope", "LSMean")) %>%
  left_join(sig_lookup2, by = c("Period")) %>%
  mutate(
    Label = case_when(Period == "NA" ~ "Period Slope",
                      Period == "1" ~ "Period 1 EMM",
                      Period == "2" ~ "Period 2 EMM",
                      Period == "3" ~ "Period 3 EMM"),
    `Treatment Group` = factor(`Treatment Group`, levels = c("CNT", "TRT")),
    lower = Estimate - 1.96 * `Standard Error`,
    upper = Estimate + 1.96 * `Standard Error`)



sig_labels2 <- Training_resp_parms %>%
  group_by(Label) %>%
  summarise(
    sig = first(sig),
    x = mean(Estimate),
    .groups = "drop"
  ) %>%
  filter(sig != "")



training_plot2 <- ggplot(
  Training_resp_parms,
  aes(
    x = Estimate,
    y = `Treatment Group`,
    color = `Treatment Group`
  )
) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    height = 0.18,
    linewidth = 1
  ) +
  geom_point(size = 4) +
  facet_grid(`Response Variable` ~ Label,
    scales = "free",
    switch = "y"
  ) +
  scale_color_manual(
    values = c(
      "TRT" = "darkorange",
      "CNT" = "forestgreen"
    )
  ) +
  geom_text(
    data = sig_labels2,
    aes(
      x = x,
      y = 1.5,
      label = sig
    ),
    inherit.aes = FALSE,
    size = 6,
    color = "black"
  ) +
  labs(
    x = "Slope estimate ± 95% CI",
    y = NULL,
    color = "Treatment",
    title = "Responsiveness (Training Phase Only)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 


training_plot2

ggsave(filename = "training_plot2_sig.png",
       plot = training_plot2,
       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence/Results",
       dpi = 600,
       width = 7,
       height = 2.7,
       units = "in")





#------------Plot 3--------------

sig_lookup3 <- Full_results %>%
  dplyr::filter(
    Phase == "EZT",
    !Period %in% c("1", "2", "3", "4"),
    `Statistic Type` %in% c("LSMean Contrast", "Slope contrast") 
  ) %>%
  mutate(
    p_clean = str_replace(as.character(`P-value`), "<", ""),
    p_num = as.numeric(p_clean),
    sig = case_when(
      is.na(p_num)  ~ "",
      p_num < 0.001 ~ "***",
      p_num < 0.01  ~ "**",
      p_num < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(
    `Response Variable`,
    `Model Parameter`,
    `P-value`,
    sig)


EZT_parms <- Full_results %>%
  filter(
    Phase == "EZT",
    is.na(Period) | !as.character(Period) %in% c("1", "2", "3", "4"),
    `Statistic Type` %in% c("LSMean", "Slope")
  ) %>%
  left_join(sig_lookup3, by = c("Response Variable", "Model Parameter")) %>%
  mutate(
    `Response Variable` = factor(
      `Response Variable`,
      levels = c("Audios", "Pulses", "Points in Open EZ"),
      labels = c("Audio", "Pulse", "Open Bale")
    ),
    `Model Parameter` = factor(
      `Model Parameter`,
      levels = c("A", "r", "t50",
                 "Period slope for A", "Period slope for r", "Period slope for t50"),
      labels = c("A EMM", "r EMM", "t50 EMM",
                 "Period slope of A", "Period slope of r", "Period slope of t50")
    ),
    `Treatment Group` = factor(`Treatment Group`, levels = c("CNT", "TRT")),
    lower = Estimate - 1.96 * `Standard Error`,
    upper = Estimate + 1.96 * `Standard Error`
  )

sig_labels3 <- EZT_parms %>%
  group_by(`Response Variable`, `Model Parameter`) %>%
  summarise(
    sig = first(sig),
    x = mean(Estimate),
    .groups = "drop"
  ) %>%
  filter(sig != "")



EZT_plot <- ggplot(
  EZT_parms,
  aes(
    x = Estimate,
    y = `Treatment Group`,
    color = `Treatment Group`
  )
) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    height = 0.18,
    linewidth = 1
  ) +
  geom_point(size = 4) +
  facet_grid(
    `Response Variable` ~ `Model Parameter`,
    scales = "free",
    switch = "y"
  ) +
  scale_color_manual(
    values = c(
      "TRT" = "darkorange",
      "CNT" = "forestgreen"
    )
  ) +
  geom_text(
    data = sig_labels3,
    aes(
      x = x,
      y = 1.5,
      label = sig
    ),
    inherit.aes = FALSE,
    size = 6,
    color = "black"
  ) +
  labs(
    x = "Slope estimate ± 95% CI",
    y = NULL,
    color = "Treatment",
    title = "EZT Phase Estimates"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 


EZT_plot

ggsave(filename = "EZT_plot_sig.png",
       plot = EZT_plot,
       path = "C:/Users/Sebastian/Documents/R/Virtual_Fence/Results",
       dpi = 600,
       width = 9,
       height = 4.5,
       units = "in")



