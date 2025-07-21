packages_to_load <- c("tidyverse", "glmmTMB", "emmeans", "ggplot2", "sf", "terra") # create vector of R package names that we know are needed in the rest of the code


## ----load_libraries----
lapply(packages_to_load, library, character.only = TRUE)


##-----------------------------VF Efficacy-----------------------------

#-----Escapes - GLMM----
esc <- read.csv(file = "C:/Users/spsch/Documents/R/Virtual_Fence/Escapes.csv")

esc <- esc %>%
  mutate(Escape_Type = as.factor(Escape_Type)) %>%
  mutate(Group = as.factor(Group))

# This data includes zeros to give the model the information that escapes did not occur under those conditions. If you only analyze rows where an escape happened, the model thinks missing combinations are “unknown,” not “zero,” which (1) throws away most of your data, (2) inflates uncertainty, and (3) can easily lead to a non‑significant result even when the raw counts suggest a real difference.
esc_summ <- read.csv(file = "C:/Users/spsch/Documents/R/Virtual_Fence/Escapes_summary.csv")


groups <- read.csv("C:/Users/spsch/Documents/R/Virtual_Fence/Groups.csv")


# Why use GLMM rather than LMM? Data in this case are counts. Counts are not normally distributed and their variance usually depends on the mean (e.g., Poisson or negative binomial behavior).GLMMs extend LMMs by allowing different distributions (e.g., Poisson, Negative Binomial) for the response variable. Instead of modeling the raw response, GLMMs use a link function (log link for counts), ensuring predictions stay positive.

# Fit Poisson GLMM
m1 <- glmmTMB(
  total_esc ~ esc_type + (1|group),
  offset = log(animal_days), # An offset in a Poisson model adjusts for exposure time (or n) to estimate rate instead of raw counts. 
  family = poisson,
  data = esc_summ
)

summary(m1)
# Rate ratio
emmeans(m1, pairwise ~ esc_type, offset = log(1), type = "response")# offset = log(1) sets the rate to escapes per 1 animal-day

# Model Validation
sim_res <- simulateResiduals(m1, n = 1000)
testDispersion(sim_res)


#-----Escapes - Bar plot-----

# Get estimated mean rates from your model (m_pois or m_nb)
emm <- emmeans(m1, ~ esc_type, offset = log(1), type = "response")

emm_df <- as.data.frame(emm) %>%
  mutate(
    SE = SE,                # Standard Error
    ymin = rate - SE,   # Lower error bar
    ymax = rate + SE    # Upper error bar
  )


# Plot
ggplot(emm_df, aes(x = esc_type, y = rate, fill = esc_type)) +
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









#-----Percentage Points In/Out-----

#eshep_df1 <- read.csv(file = "C:/Users/spsch/Documents/R/Virtual_Fence/2025-07-13.csv")
#eshep_df2 <- read.csv(file = "C:/Users/spsch/Documents/R/Virtual_Fence/2025-07-20.csv")

#eshep_df <- rbind.data.frame(eshep_df1, eshep_df2) %>% # Filter to include only the three training and four testing days
#  filter(Time..UTC. > "2025-07-07 11:55:00" & Time..UTC. < "2025-07-14 12:05:00")

#eshep_full_df <- eshep_df %>%
#  left_join(groups, by = "Neckband.ID") %>%
#  filter(!is.na(Animal_ID))

#write.csv(eshep_full_df, file = "C:/Users/spsch/Documents/R/Virtual_Fence/eshep_full_df.csv")



#----Read in and clean data-------
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

# Add columns based on whether points occur within polygons of interest
training1_inside <- st_within(eshep_full_sf, training_VPs[1, ], sparse = FALSE)[ , 1]
training2_inside <- st_within(eshep_full_sf, training_VPs[2, ], sparse = FALSE)[ , 1]
training3_inside <- st_within(eshep_full_sf, training_VPs[3, ], sparse = FALSE)[ , 1]
west_ex_inside <- st_within(eshep_full_sf,
                            filter(exclusion_zones, Description == "west"),
                                   sparse = FALSE)[ , 1]
east_ex_inside <- st_within(eshep_full_sf,
                            filter(exclusion_zones, Description == "east"),
                            sparse = FALSE)[ , 1]


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


# Match point occurrence to experiment period
eshep_full_sf <- eshep_full_sf %>%
  mutate(training1_inside = ifelse(period == "Tr1", training1_inside, NA)) %>%
  mutate(training2_inside = ifelse(period == "Tr2", training2_inside, NA)) %>%
  mutate(training3_inside = ifelse(period == "Tr3", training3_inside, NA)) %>%
  mutate(west_ex_inside = ifelse(period %in% c("W1", "W2"), west_ex_inside, NA)) %>%
  mutate(east_ex_inside = ifelse(period %in% c("E1", "E2"), east_ex_inside, NA))

eshep_full_sf <- eshep_full_sf %>%
  mutate(
    inside_VF = case_when(
      period == "Tr1" ~ training1_inside,
      period == "Tr2" ~ training2_inside,
      period == "Tr3" ~ training3_inside,
      period %in% c("E1","E2") ~ !east_ex_inside,
      period %in% c("W1","W2") ~ !west_ex_inside))


# Calculate proportion of points inside VF for each animal and convert to df
m2_df <- eshep_full_sf %>%
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
emm <- emmeans(m2, ~ TRT, type = "response")  
print(emm)



#-----Bar Plot-------
emm_df <- as.data.frame(emm) %>%
  mutate(
    SE = SE,
    ymin = asymp.LCL,
    ymax = asymp.UCL
  ) %>%
  mutate(treatment = ifelse(TRT == "CNT", "Control", "Treatment"))

ggplot(emm_df, aes(x = treatment, y = prob, fill = treatment)) +
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
  scale_fill_manual(values = c("darkorange", "lightgreen"))



