# Purpose
# This is an attempt to analyze the community structure of the FMWT datasets,
# for the random and fixed stations. If successful, the two community structures
# will be compared to assess if there are differences between the two designs
# on the species community level.

# Workflow ----------------------------------------------------------------

# 1. Read in the data
# 2. Clean the data
# 3. Define functions required for the analysis
# 4. Run the models
# 5. Optimize the models
# 6. Validate the models
# 7. Compare the models
# 8. Conclusion

# Required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(ggdark)
library(PTAk)
library(cluster)
library(dendextend)
library(gifski)
library(ggridges)

options(scipen = 99999)

myTheme <- dark_theme_bw(base_size = 24) +
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#2c2828", color = NA),
        panel.background = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.major = element_line(color = "#646464", size = 0.2),
        legend.background = element_blank(),
        legend.key = element_blank())

theme_set(myTheme)

# Reading in the data -----------------------------------------------------

# # This data has both the FMWT and special study in it. Just catch matrix of 2021-2022
# data <- read_xlsx(list.files(file.path("data"), "*.xlsx", full.names = T)) %>% 
#   # Change catch into CPUE
#   # The provided CPUE is already correct
#   mutate(same = ifelse(all.equal(Catch / Volume * 1000, CPUE), T, F),
#          # Creating a month column to define the time component as month-year
#          month = case_when(SurveyNumber == 3 ~ 9,
#                            SurveyNumber == 4 ~ 10,
#                            SurveyNumber == 5 ~ 11,
#                            SurveyNumber == 6 ~ 12)
#          # # Combine all striped bass into a singular
#          # Species = ifelse(Species %in% c("Striped Bass age-0", "Striped Bass age-1", 
#          #                                 "Striped Bass age-2", "Striped Bass age-3"), 
#          #                  "Striped Bass", Species)
#          ) %>% 
#   # Removing "Goby (unid)" from the analysis
#   filter(!Species %in% c("Goby (unid)", "Jellyfish", "Unid")) %>% 
#   rename(Station = StationCode) %>% 
#   # Split the data into a list of two elements
#   split(., .$Study)

data <- readRDS(file.path("data", "SSdata.Rds")) %>% 
  mutate(same = ifelse(all.equal(Catch / Volume * 1000, CPUE), T, F),
         # Creating a month column to define the time component as month-year
         month = case_when(SurveyNumber == 3 ~ 9,
                           SurveyNumber == 4 ~ 10,
                           SurveyNumber == 5 ~ 11,
                           SurveyNumber == 6 ~ 12)
         # # Combine all striped bass into a singular
         # Species = ifelse(Species %in% c("Striped Bass age-0", "Striped Bass age-1", 
         #                                 "Striped Bass age-2", "Striped Bass age-3"), 
         #                  "Striped Bass", Species)
  ) %>% 
  # Removing "Goby (unid)" from the analysis
  filter(!Species %in% c("Goby (unid)", "Jellyfish", "Unid")) %>% 
  rename(Station = StationCode) %>% 
  # Split the data into a list of two elements
  split(., .$Study)

# Are there important stations that should not be removed?
importantStations <- NA


# Whole dataset to run the prop 0 analysis --------------------------------

data$all <- read.csv(file.path("data", "FMWT 1967-2022 Catch Matrix_updated_tidy.csv")) %>% 
  mutate(Region = case_when(StationCode %in% c(305:339, 401:407) ~ "San Pablo Bay and Carquinez Strait",
                            StationCode %in% c(340:341) ~ "Napa River",
                            StationCode %in% c(408:418, 501:509, 515:519, 601:604) ~ "Suisun and Honker Bays",
                            StationCode %in% c(605:608) ~ "Suisun Marsh",
                            StationCode %in% c(510:513, 701:711, 802:813) ~ "Confluence",
                            StationCode %in% c(713:716, 719, 722:723) ~ "Cache Slough",
                            StationCode %in% c(795:797) ~ "Sacramento Ship Channel"),
         CPUE = (Catch/Volume) * 10000)

# Checking data distributions ---------------------------------------------

# Finding stations that were sampled per time point
stationFilterYear <- function(data, yearStart, perSurvey = F) {
  
  if (perSurvey) {
    elapsedTime <- (max(data$Year) - yearStart + 1) * 4 # 4 surveys
    
    dataDistinct <- data %>% 
      distinct(Year, month, Station)
  } else {
    elapsedTime <- (max(data$Year) - yearStart + 1)
    
    dataDistinct <- data %>% 
      distinct(Year, Station)
  }
  
  dataDistinct %>% 
    filter(Year >= yearStart) %>% 
    group_by(Station) %>% 
    count() %>% 
    filter(n == elapsedTime) %>% 
    pull(Station)
}

# Plotting the number of stations per time step....
data$FMWT %>% 
  distinct(Station, month, Year) %>% 
  mutate(dateDummy = as.Date(paste(Year, month, "01", sep = "-"))) %>% 
  arrange(dateDummy) %>% 
  mutate(dateDummy = factor(as.character(dateDummy), levels = unique(as.character(dateDummy))),
         stationFactor = as.numeric(factor(Station, levels = unique(Station)))) %>% 
  {
    ggplot(., aes(dateDummy, stationFactor)) +
      geom_tile(color = "black") +
      theme(panel.grid.major = element_blank()) +
      scale_y_continuous(sec.axis = sec_axis(~.*1,
                                             breaks = seq(2, max(.$stationFactor), 2),
                                             labels = filter(., stationFactor %in%
                                                               seq(2, max(.$stationFactor), 2)) %>%
                                               pull(Station) %>%
                                               unique() %>%
                                               as.character()),
                         breaks = seq(1, max(.$stationFactor), 2),
                         labels = filter(., stationFactor %in% seq(1, max(.$stationFactor), 2)) %>%
                           pull(Station) %>%
                           unique() %>%
                           as.character(),
                         expand = expansion(add = c(0.6, 0.6)))
  } +
  theme(axis.text.y = element_text(size = 15))

proportionCPUE <- data$all %>% 
  filter(Year >= 2013) %>% 
  group_by(Region = factor(Region, 
                           levels = c("Sacramento Ship Channel", "Cache Slough", "Confluence", "Suisun Marsh",
                                      "Suisun and Honker Bays", "San Pablo Bay and Carquinez Strait", "Napa River")),
           Species) %>%
  summarise(totalCPUE = sum(CPUE, na.rm = T)) %>% 
  # Remove stations not assigned to a region
  filter(!is.na(Region)) %>% 
  group_by(Species) %>% 
  mutate(proportionCPUE = totalCPUE / sum(totalCPUE),
         highestRegion = ifelse(totalCPUE == max(totalCPUE) & totalCPUE > 0, as.character(Region), NA)) %>% 
  # Highest region is the region in which the species was found the most often
  fill(highestRegion) %>% 
  ungroup() %>% 
  # I want to order the x-axis based on prop per region, starting with Sacramento Ship Channel
  mutate(xIndex = case_when(proportionCPUE == 1 & Region == "Sacramento Ship Channel" ~ 1,
                            proportionCPUE == 1 & Region == "Cache Slough" ~ 2,
                            proportionCPUE == 1 & Region == "Confluence" ~ 3,
                            proportionCPUE == 1 & Region == "Suisun Marsh" ~ 4,
                            proportionCPUE == 1 & Region == "Suisun and Honker Bays" ~ 5,
                            proportionCPUE == 1 & Region == "San Pablo Bay and Carquinez Strait" ~ 6,
                            proportionCPUE == 1 & Region == "Napa River" ~ 7,
                            # Now for species that are not 1
                            proportionCPUE > 0 & highestRegion == "Sacramento Ship Channel" ~ 8,
                            proportionCPUE > 0 & highestRegion == "Cache Slough" ~ 9,
                            proportionCPUE > 0 & highestRegion == "Confluence" ~ 10,
                            proportionCPUE > 0 & highestRegion == "Suisun Marsh" ~ 11,
                            proportionCPUE > 0 & highestRegion == "Suisun and Honker Bays" ~ 12,
                            proportionCPUE > 0 & highestRegion == "San Pablo Bay and Carquinez Strait" ~ 13,
                            proportionCPUE > 0 & highestRegion == "Napa River" ~ 14,
                            # Now for all other species with 0 catch
                            TRUE ~ 15)) %>% 
  group_by(Species) %>%
  # I want to order by highest occurring species first per index
  mutate(xIndex = min(xIndex),
         maxCPUE = max(totalCPUE)) %>% 
  group_by(xIndex, reorder(Species, -maxCPUE)) %>% 
  mutate(orderedIndex = cur_group_id()) %>% 
  ungroup()

# Plotting all species
ggplot(proportionCPUE, aes(orderedIndex, Region, fill = proportionCPUE)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis_c(option = "turbo")

# Plotting only species with catch
proportionCPUE %>% 
  filter(xIndex < 15) %>% 
  ggplot(aes(orderedIndex, Region, fill = proportionCPUE)) +
  geom_tile(color = "#646464") +
  scale_y_discrete(limits = rev,
                   expand = expansion(mult = c(0, 0))) +
  scale_fill_viridis_c(option = "turbo") +
  scale_x_continuous(breaks = seq(1, length((proportionCPUE %>% 
                                               filter(xIndex < 15) %>% 
                                               distinct(Species, orderedIndex) %>% 
                                               arrange(orderedIndex) %>% 
                                               pull(Species))), 1),
                     labels = (proportionCPUE %>% 
                                 filter(xIndex < 15) %>% 
                                 distinct(Species, orderedIndex) %>% 
                                 arrange(orderedIndex) %>% 
                                 pull(Species)),
                     expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(fill = "Proportional CPUE",
       x = "Species")
# Plotting species that are dominant in ship channel, cache slough, and confluence only (> 0.5)
speciesCatchRateAll <- proportionCPUE %>% 
  mutate(regionOfInterest = ifelse(Region %in% c("Sacramento Ship Channel", "Cache Slough", "Confluence"), T, F)) %>%
  group_by(Species, regionOfInterest) %>% 
  mutate(sumCPUE = sum(proportionCPUE, na.rm = T)) %>% 
  ungroup() %>% 
  filter(regionOfInterest == T, sumCPUE > 0.5) %>% 
  pull(Species) %>% 
  unique()

proportionCPUE %>% 
      filter(Species %in% speciesCatchRateAll) %>% 
      group_by(orderedIndex) %>% 
      mutate(orderedIndex = cur_group_id()) %>% 
      ggplot(aes(reorder(Species, orderedIndex), Region, fill = proportionCPUE)) +
      geom_tile(color = "#646464") +
      scale_y_discrete(limits = rev,
                       expand = expansion(mult = c(0, 0))) +
      scale_fill_viridis_c(option = "turbo") +
      # scale_x_continuous(breaks = seq(1, length(speciesCatchRateAll), 1),
      #                    labels = speciesCatchRateAll,
      #                    expand = expansion(mult = c(0, 0))) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(fill = "Proportional CPUE",
           x = "Species",
           title = "Species dominant in Sacramento Ship Channel, Cache Slough, and Confluence")

# Plotting species that are dominant in up to Suisun Bay (> 0.5)
proportionCPUE %>% 
  filter(Region %in% c("Sacramento Ship Channel", "Cache Slough", "Confluence", 
                       "Suisun Marsh", "Suisun and Honker Bays"),
         proportionCPUE > 0.5) %>% 
  pull(Species) %>% 
  {
    speciesOfInterest <- .
    proportionCPUE %>% 
      filter(Species %in% speciesOfInterest) %>% 
      group_by(orderedIndex) %>% 
      mutate(orderedIndex = cur_group_id()) %>% 
      ggplot(aes(orderedIndex, Region, fill = proportionCPUE)) +
      geom_tile(color = "#646464") +
      scale_y_discrete(limits = rev,
                       expand = expansion(mult = c(0, 0))) +
      scale_fill_viridis_c(option = "turbo") +
      scale_x_continuous(breaks = seq(1, length(speciesOfInterest), 1),
                         labels = speciesOfInterest,
                         expand = expansion(mult = c(0, 0))) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(fill = "Proportional CPUE",
           x = "Species",
           title = "Species dominant in Sacramento Ship Channel, Cache Slough, \nConfluence, Suisun Marsh, and Suisun and Honker Bays")
  }

# More stations in 2022 than 2021
# This calculates the number of unique stations that exists for the entire dataset from year X to present
# E.g., Station 501 exists in 2021-2022 = makes it to 2021 year but if 501 only existed in 2022, would exists
# only in 2022
data$FMWTRegion <- data$FMWT %>% 
  mutate(Station = Region) %>% 
  group_by(Study, Year, SurveyNumber, month, Station, Region, Species) %>% 
  summarise(CPUE = mean(CPUE, na.rm = T), .groups = "drop")

stationConsistentYear <- lapply(c(data$FMWTRegion$Year %>% unique()),
                                function(x) {
                                  data$FMWTRegion %>% 
                                    stationFilterYear(yearStart = x,
                                                      perSurvey = T) %>% 
                                    data.frame(Station = .)
                                }) %>% 
  setNames(data$FMWTRegion$Year %>% unique()) %>% 
  bind_rows(.id = "Year") %>% 
  group_by(Year) %>% 
  mutate(cumulativeCount = 1:n()) %>% 
  arrange(Year, Station)

fishOccurencesYearly <- data$FMWTRegion %>% 
  # Must take into account the starting year to ensure only stations sampled in all years are used
  filter(Station %in% (filter(stationConsistentYear, Year %in% 2022) %>% pull(Station))) %>% 
  distinct(Year, Station, Species) %>% 
  group_by(Station, Species) %>% 
  count() %>% 
  ungroup() %>% 
  distinct(Species, n) %>% 
  group_by(Species) %>% 
  slice_max(n) %>% 
  arrange(-n) %>% 
  ungroup()

# Find the proportion of zeroes per species
# Only want from the top 5 stations that each species occur at
speciesPropZero <- data$FMWTRegion %>% 
  # Only want stations that have sampling across all time steps
  filter(Station %in% pull(filter(stationConsistentYear, Year == 2022), Station)) %>% 
  group_by(Station, Region, Species) %>% 
  summarise(sumCPUE = sum(CPUE), .groups = "drop") %>% 
  # Sort
  arrange(-sumCPUE) %>% 
  # Per species, pick highest CPUE stations, 5 here
  group_by(Species) %>% 
  slice(1:2) %>% 
  # Change names to bind
  transmute(topStations = Station, Species) %>% 
  ungroup() %>% 
  # Now join to beginning dataset again
  full_join(data$FMWTRegion,
            by = "Species") %>% 
  # Filter for only the top stations
  filter(Station == topStations) %>% 
  # Determine the proportion of zeroes
  group_by(Species) %>% 
  mutate(isZero = ifelse(CPUE == 0, 1, 0),
         noZero = ifelse(CPUE > 0, 1, 0)) %>% 
  summarise(sumNotZero = sum(noZero),
            propZero = sum(isZero)/max(n()), 
            .groups = "drop") %>% 
  mutate(Species = as.character(Species)) %>% 
  arrange(propZero) %>% 
  filter(propZero < 1,
         # Need to be caught more than just once
         sumNotZero > 1)

# It appears as 0.27 may work for this dataset as well
speciesPropZero %>% 
  mutate(remove = ifelse(propZero <= quantile(propZero, 0.27), T, F)) %>% 
  {
    # Plotting it
    ggplot(., aes(reorder(Species, -propZero), propZero, color = remove)) +
      geom_point(size = 3, color = "#FFFFB8") +
      # geom_vline(xintercept = c("Prickly Sculpin", "Speckled Sanddab")) +
      # scale_color_manual(values = c("TRUE" = "#FF4D00", "FALSE" = "#FFFFB8")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            panel.grid.major = element_blank(),
            legend.position = "none") +
      labs(y = "Proportion of Zeroes",
           x = "Common Name")
  }

# Proportion of zero as done by Jereme
# All species found as a proportion of total number of sampling events
# How many distinct sampling events in 2022
speciesCatchRateFMWT <- data$FMWT %>% 
  filter(Year %in% 2022, Region %in% stationConsistentYear$Station) %>% 
  group_by(SurveyNumber, Station) %>% 
  mutate(sampleEvent = cur_group_id(),
         catch = ifelse(CPUE > 0, 1, 0)) %>% 
  group_by(Species, sampleEvent) %>% 
  mutate(catch = sum(ifelse(sum(CPUE > 0, na.rm = T) > 0, T, F))) %>% 
  group_by(Species) %>% 
  summarise(totalSample = max(sampleEvent),
            overallCatch = sum(catch), .groups = "drop") %>%
  filter(overallCatch > 0)
         # Species %in% speciesCatchRateAll)

speciesCatchRateSS <- data$SS %>% 
  filter(Year %in% 2022, Region %in% stationConsistentYear$Station) %>% 
  group_by(SurveyNumber, Station) %>% 
  mutate(sampleEvent = cur_group_id(),
         catch = ifelse(CPUE > 0, 1, 0)) %>% 
  group_by(Species, sampleEvent) %>% 
  mutate(catch = sum(ifelse(sum(CPUE > 0, na.rm = T) > 0, T, F))) %>% 
  group_by(Species) %>% 
  summarise(totalSample = max(sampleEvent),
            overallCatch = sum(catch), .groups = "drop") %>%
  filter(overallCatch > 0)
         # Species %in% speciesCatchRateAll)

speciesCatchRate <- speciesCatchRateFMWT %>% 
  transmute(Species, overallCatchFMWT = overallCatch) %>% 
  full_join(speciesCatchRateSS %>% 
              transmute(Species, overallCatchSS = overallCatch),
            by = "Species") %>% 
  rename(FMWT = overallCatchFMWT, SS = overallCatchSS) %>% 
  pivot_longer(-Species, names_to = "study", values_to = "catch")

speciesCatchRate %>% 
  ggplot(aes(reorder(Species, catch), catch)) +
  geom_col() +
  geom_text(aes(label = catch), vjust = -1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap(~study) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.07))) +
  labs(x = "Species", y = "# of Catch Events", title = "Filtered for only Cache Slough, Confluence, and Sacramento Ship Channel")

speciesCatchRate <- speciesCatchRate %>% 
  pivot_wider(names_from = study, values_from = catch) %>% 
  filter(!is.na(FMWT), !is.na(SS)) %>% 
  mutate(overallCatch = FMWT + SS)

speciesCatchRate %>% 
  ggplot(aes(reorder(Species, overallCatch), overallCatch)) +
  geom_col() +
  geom_text(aes(label = overallCatch), vjust = -1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.07))) +
  labs(x = "Species", y = "# of Catch Events", title = "Filtered for only Cache Slough, Confluence, and Sacramento Ship Channel")

# Transformations
transformation <- function(x, transformation = c("log", "squareRoot", "cubicRoot", "fourthRoot"), zscore = F) {
  transformation <- match.arg(transformation)
  y <- switch(transformation,
              log = log(x + 1),
              squareRoot = sqrt(x),
              cubicRoot = x^(1/3),
              fourthRoot = x^(1/4))
  
  if (isTRUE(zscore)) {
    y = (y - mean(y))/sd(y)
  }
  y
}

# lapply(c("log", "squareRoot", "cubicRoot", "fourthRoot"),
#        function(y) {
#          data$FMWTRegion %>% 
#            filter(Species %in% speciesPropZero$Species,
#                   CPUE > 0) %>% 
#            mutate("transformedCPUE" := transformation(CPUE, y, zscore = T),
#                   transformation = y)
#        }) %>% 
#   bind_rows() %>% 
#   ggplot(aes(x = transformedCPUE, y = Species)) +
#   geom_density_ridges(panel_scaling = F) +
#   # scale_x_continuous(limits = c(-1, 5)) +
#   facet_wrap(~factor(transformation, levels = c("log", "squareRoot", "cubicRoot", "fourthRoot"))) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# # Toss up between log and cubicRoot

# Calculating overlapping coefficient here
overlapping_coefficient <- function(x, y, 
                                    transformation = c("raw", "log", "squareRoot", "cubicRoot", "fourthRoot"),
                                    centered = T, figure = F){
  
  # First, calculate skewness
  skew <- lapply(list(x, y), function(var) round(e1071::skewness(var), 4))
  
  transformation <- match.arg(transformation)
  
  xTransformed <- switch(transformation,
                         raw = x,
                         log = log(x + 1),
                         squareRoot = sqrt(x),
                         cubicRoot = x^(1/3),
                         fourthRoot = x^(1/4))
  
  yTransformed <- switch(transformation,
                         raw = y,
                         log = log(y + 1),
                         squareRoot = sqrt(y),
                         cubicRoot = y^(1/3),
                         fourthRoot = y^(1/4))
  
  if (isTRUE(centered)) {
    xTransformed <- (xTransformed - mean(xTransformed))/sd(xTransformed)
    yTransformed <- (yTransformed - mean(yTransformed))/sd(yTransformed)
  }
  
  if (transformation != "raw") {
    skewTransformed <- lapply(list(xTransformed, yTransformed), function(var) round(e1071::skewness(var), 4))
  } else {
    skewTransformed <- "raw"
  }
  
  # jointDistribution <- function(xDistribution, yDistribution) {
  #   combined <- c(xDistribution, yDistribution)
  #   
  #   densityX <- density(xDistribution, from = min(combined), to = max(combined))
  #   densityY <- density(yDistribution, from = min(combined), to = max(combined))
  #   
  #   joint <- pmin(densityX$y, densityY$y)
  #   
  #   lengthCombined <- seq(min(combined), max(combined), length.out = length(joint))
  #   
  #   xArea <- integrate(approxfun(densityX), min(combined), max(combined), stop.on.error = F)
  #   yArea <- integrate(approxfun(densityY), min(combined), max(combined), stop.on.error = F)
  #   overlapArea <- integrate(approxfun(lengthCombined, joint),
  #                            min(combined), max(combined), stop.on.error = F)
  #   
  #   xArea <- ifelse(xArea$message == "OK", xArea$value, NA)
  #   yArea <- ifelse(yArea$message == "OK", yArea$value, NA)
  #   overlapArea <- ifelse(overlapArea$message == "OK", overlapArea$value, NA)
  #   
  #   
  #   val <- mean(c(overlapArea/xArea, overlapArea/yArea))
  #   
  #   list(densityX = densityX,
  #        densityY = densityY,
  #        joint = joint,
  #        val = val)
  # }
  
  combined <- c(xTransformed, yTransformed)
  
  # integrateCalc <- jointDistribution(xTransformed, yTransformed)
  densityX <- density(xTransformed, from = min(combined), to = max(combined), n = 2^10)
  densityY <- density(yTransformed, from = min(combined), to = max(combined), n = 2^10)
  
  joint <- pmin(densityX$y, densityY$y)
  overlap <- mean(sum(joint)/sum(densityX$y),
                  sum(joint)/sum(densityY$Y))
  
  # p <- data.frame(x = rep(integrateCalc$densityX$x, 3),
  #                 y = c(integrateCalc$densityX$y, integrateCalc$densityY$y, integrateCalc$joint),
  #                 data = factor(rep(c("x", "y", "overlap"), each = length(integrateCalc$densityX$x)),
  #                               levels = c("x", "y", "overlap"))) %>% 
  #   ggplot(aes(x, y, fill = data)) +
  #   geom_area(position = position_identity(), color = "black") +
  #   scale_fill_brewer(palette = "Dark2") +
  #   labs(subtitle = paste0("Overlap: ", round(integrateCalc$val, 4), "\nSkewness Old (x, y): ",
  #                          paste(skew, collapse = ", "), "\nSkewness New (x, y): ",
  #                          ifelse(is.character(skewTransformed), "No change", 
  #                                 paste(skewTransformed, collapse = ", "))))
  
  p <- data.frame(x = rep(densityX$x, 3),
                  y = c(densityX$y, densityY$y, joint),
                  data = factor(rep(c("x", "y", "overlap"), each = length(densityX$x)),
                                levels = c("x", "y", "overlap"))) %>% 
    ggplot(aes(x, y, fill = data)) +
    geom_area(position = position_identity(), color = "black") +
    scale_fill_brewer(palette = "Dark2") +
    labs(subtitle = paste0("Overlap: ", round(overlap, 4), "\nSkewness Old (x, y): ",
                           paste(skew, collapse = ", "), "\nSkewness New (x, y): ",
                           ifelse(is.character(skewTransformed), "No change", 
                                  paste(skewTransformed, collapse = ", "))))
  if (isTRUE(figure)) {
    return(p)
  }
  
  data.frame(transformation = transformation,
             skewOld = unlist(skew),
             skewNew = ifelse(is.character(skewTransformed), NA, unlist(skewTransformed)),
             overlap = overlap)
}

# overlapping_coefficient(with(data$FMWT, CPUE[Species %in% "Longfin Smelt"]),
#                         with(data$FMWT, CPUE[Species %in% "American Shad"]), 
#                         transformation = "raw")

transformations <- lapply(combn(speciesCatchRate$Species, 2, simplify = F), 
                          function(x) {
                            
                            raw <- overlapping_coefficient(with(data$FMWTRegion, CPUE[Species %in% x[1]]),
                                                           with(data$FMWTRegion, CPUE[Species %in% x[2]]),
                                                           transformation = "raw")
                            log <- overlapping_coefficient(with(data$FMWTRegion, CPUE[Species %in% x[1]]),
                                                           with(data$FMWTRegion, CPUE[Species %in% x[2]]),
                                                           transformation = "log")
                            squareRoot <- overlapping_coefficient(with(data$FMWTRegion, CPUE[Species %in% x[1]]),
                                                                  with(data$FMWTRegion, CPUE[Species %in% x[2]]),
                                                                  transformation = "squareRoot")
                            cubicRoot <- overlapping_coefficient(with(data$FMWTRegion, CPUE[Species %in% x[1]]),
                                                                 with(data$FMWTRegion, CPUE[Species %in% x[2]]),
                                                                 transformation = "cubicRoot")
                            fourthRoot <- overlapping_coefficient(with(data$FMWTRegion, CPUE[Species %in% x[1]]),
                                                                  with(data$FMWTRegion, CPUE[Species %in% x[2]]),
                                                                  transformation = "fourthRoot")
                            
                            data.frame(species = rep(c(x), 5),
                                       speciesCombination = 
                                         rep(paste0(janitor::make_clean_names(x, case = "lower_camel"), 
                                                    collapse = "."), 5)) %>% 
                              bind_cols(bind_rows(raw, log, squareRoot, cubicRoot, fourthRoot))
                          }) %>% 
  bind_rows() %>% 
  mutate(transformation = factor(transformation, levels = c("raw", "log", "squareRoot", "cubicRoot", "fourthRoot")))

library(patchwork)
{transformations %>% 
    pivot_longer(c(skewOld, skewNew), names_to = "skewType", values_to = "skew") %>% 
    mutate(type = factor(skewType, levels = c("skewOld", "skewNew"))) %>% 
    ggplot(aes(transformation, skew, color = type)) +
    geom_boxplot(position = position_dodge(preserve = "single"), fill = "#2C2828FF", size = 1) +
    scale_color_manual(values = c("#FF4D00", "#FFFFB8")) +
    theme(axis.text.x.bottom = element_blank(),
          axis.title.x.bottom = element_blank())}/
  transformations %>% 
  ggplot(aes(transformation, overlap)) +
  geom_boxplot(color = "white", fill = "#2C2828FF") +
  geom_text(data = transformations %>% group_by(transformation) %>% summarise(n = sum(is.na(overlap))),
            aes(x = transformation, y = 0.05, label = paste0("NAs = ", n)), size = 6)
# will simply cube root the results from these results


# Defining helper functions -----------------------------------------------

filterDataset <- function(data, 
                          stationYear, 
                          speciesThreshold, propThres,
                          # These two tables should be created beforehand
                          stationAdd = NULL, stationRemove = NULL,
                          speciesAdd = NULL, speciesRemove = NULL,
                          # If you wanted a static list of stations or species
                          stationFixed = NULL,
                          speciesFixed = NULL,
                          stationDF = stationConsistentYear,
                          fishDF = fishOccurencesYearly) {
  # stationYear will depend on the stationConsistentYear data frame; detects which year to use
  # The list of stations will come from the data frame
  # speciesThreshold will depend on the fishOccurencesYearly table, in which total number
  # of catch per year will be used as the filter. Species with values >= provided number
  # will be returned
  
  if (is.null(stationFixed)) {
    stationsToKeep <- stationDF %>% 
      filter(Year == stationYear) %>% 
      pull(Station) %>% 
      c(stationAdd, .) %>% 
      unique()
    
    if (!is.null(stationRemove)) {
      stationsToKeep <- stationsToKeep[-which(stationsToKeep %in% stationRemove)]
    }
  } else {
    stationsToKeep <- stationFixed
  }

  if (is.null(speciesFixed)) {
    # speciesToKeep <- fishDF %>% 
    #   filter(n >= speciesThreshold) %>% 
    #   pull(Species) %>% 
    #   c(speciesAdd, .) %>% 
    #   unique()
    
    speciesToKeep <- speciesPropZero %>% 
      filter(propZero <= quantile(propZero, propThres)) %>% 
      pull(Species) %>% 
      c(speciesAdd, .) %>% 
      unique()
    
  } else {
    speciesToKeep <- speciesFixed
  }
  
  if (!is.null(speciesRemove)) {
    speciesToRemove <- which(speciesToKeep %in% speciesRemove)
    
    if (length(speciesToRemove) > 0) {
      speciesToKeep <- speciesToKeep[-speciesToRemove]
    }
  }
  
  data %>% 
    filter(Station %in% stationsToKeep,
           Species %in% speciesToKeep,
           Year >= stationYear)
}

buildArray <- function(data) {

  fishListSurvey <- data %>% 
    mutate(Species = as.character(Species)) %>% 
    pivot_wider(names_from = timeStep,
                values_from = CPUE)
  
  # Are there stations that are NAs...More a problem as you get to less and less data
  # also a problem as the station filter (filterDataset) is on a YEARLY time step, 
  # but if there are entire surveys missing, the filter will not pick it up
  missingData <- fishListSurvey %>% 
    filter(if_any(everything(), is.na)) %>% 
    pivot_longer(-c(Station, Species), names_to = "Year", values_to = "CPUE") %>% 
    filter(is.na(CPUE)) %>% 
    distinct(Station, Year)
  
  if (nrow(missingData) > 0) {
    warning("NAs found for station(s) ", paste(unique(missingData$Station), collapse = ", "), 
            " for years(s) ", paste(unique(missingData$Year), collapse = ", "), ". Removing these stations",
            call. = F)
  }
  
  fishListSurvey <- fishListSurvey %>% 
    filter(!Station %in% missingData$Station) %>% 
    {split(.[, !names(.) %in% "Species"], .$Species)} %>% 
    {lapply(., function(x) {
      tibble::column_to_rownames(.data = x, var = "Station")
    })}
  
  fishArraySurvey <- fishListSurvey %>% 
    # Number of station, Number of Survey, Number of species
    {array(unlist(.), dim = c((data$Station %>% unique() %>% length()) - length(unique(missingData$Station)),
                              data$timeStep %>% unique() %>% length(),
                              data$Species %>% unique() %>% length()))}
  
  dimnames(fishArraySurvey) <- list(rownames(fishListSurvey[[1]]), # Any of the array element works here
                                    names(fishListSurvey[[1]]),
                                    names(fishListSurvey))
  
  fishArraySurvey
}


# Running the models ------------------------------------------------------

set.seed(135)
seed <- sample(1:1000000, 100)

# inputsData <- lapply(rev(seq(0.02, 0.50, 0.01)), function(x) {
#   data.frame(iteration = which(rev(seq(0.02, 0.50, 0.01)) %in% x),
#              scenario = x,
#              removeStations = NA[1:9],
#              yearStart = 2021[1:9],
#              # speciesThres = 6[1:9],
#              propThres = x[1:9],
#              speciesRemove = NA[1:9],
#              modeNames = c("Station", "Survey", "Species")[1:9])  
# })

# lengthSpecies <- 1:unique(length(speciesCatchRate$overallCatch))

selectSpecies <- function(data, index) {
  speciesCatchRate %>% 
    group_by(overallCatch) %>% 
    mutate(groupIndex = cur_group_id()) %>% 
    filter(groupIndex >= index) %>% 
    pull(Species)
}

maxRows <- length(selectSpecies(speciesCatchRate, 1))

inputsData <- lapply(1:maxRows, function(x) {
  data.frame(iteration = which(1:maxRows %in% x),
             scenario = which(1:maxRows %in% x),
             removeStations = (NA)[1:maxRows],
             yearStart = (2022)[1:maxRows],
             # speciesThres = 6[1:9],
             speciesFixed = selectSpecies(speciesCatchRate, x)[1:maxRows],
             propThres = (NA)[1:maxRows],
             speciesRemove = (NA)[1:maxRows],
             modeNames = c("Station", "Survey", "Species")[1:maxRows])  
})

names(inputsData) <- sapply(inputsData, function(x) unique(x$scenario))

simulateRemoval <- function(data, removeStations = NULL,
                            timeStep = c("Year", "Survey"),
                            yearStart, 
                            speciesThres = NULL,
                            propThres = NULL,
                            removeSpecies = NULL,
                            modeNames = c("Station", "Survey", "Species"),
                            iteration = NULL,
                            ...) {

  # Determine stations to remove, checking with those that cannot be removed
  if (any(removeStations %in% importantStations)) {
    cantRemove <- removeStations[which(removeStations %in% importantStations)]
    cat("Station(s)", paste(c(cantRemove), collapse = ", "), "will not be removed. \n")
    removeStations <- c(removeStations[which(!removeStations %in% cantRemove)])
  } else {
    if (length(removeStations) == 0) {
      removeStations = NULL
    }
    cantRemove = NULL
  }
  
  if (length(removeSpecies) == 0) removeSpecies = NULL

  # Design the simulation array with the appropriate filters
  preArray <- data %>% 
    filterDataset(stationYear = yearStart,
                  stationRemove = removeStations,
                  speciesThreshold = speciesThres,
                  propThres = propThres,
                  speciesRemove = removeSpecies, ...) %>% 
    mutate(CPUE = transformation(CPUE, "log", zscore = F))

  # Creating the relevant time step:
  if (timeStep == "Year") {
    preArray <- preArray %>% 
      mutate(timeStep = Year) %>% 
      group_by(timeStep, Station, Species) %>% 
      summarise(CPUE = sqrt(ceiling(sum(CPUE))), .groups = "drop") 
  } else {
    if (timeStep == "Survey") {
      preArray <- preArray %>%
        mutate(timeStep = SurveyNumber) %>% 
        select(timeStep, Station, Species, CPUE)
    } else {
      stop("Check your timeStep value.", call. = F)
    }
  }
  
  simulation <- preArray %>% 
    buildArray()
  
  cat("Model parameters: ", paste(c(unique(iteration), yearStart, removeStations, 
                                    speciesThres, removeSpecies), 
                                  collapse = ", "), "\n\n")
  
  # Scale the simulation array and then run the PTA
  set.seed(135)
  # This Multcent function scales and centers the data...
  # Starts by finding the mean catch of each species across all space/time
  # Then subtracts that mean from each value
  # Then, find sd per species across all space/time and scale each entry
  ptaModel <- Multcent(dat = simulation, 
                       bi = NULL, by = 3, centre = mean,
                       centrebyBA = c(TRUE, FALSE), 
                       scalebyBA = c(TRUE, FALSE)) %>% 
    PTA3(nbPT = 3, nbPT2 = 3,
         minpct = 0.1)
  
  # There are 3 modes here:
  # 1 = station (spatial mode)
  # 2 = surveys
  # 3 = species 
  # JUST NOTE THOUGH, the last "mode" will also contain other data saved by the object
  names(ptaModel) <- modeNames
  
  # Summarize the PTA
  modelScaleSummary <- summary(ptaModel, testvar = 0)
  
  # The scree plot to determine # of tensors to keep
  plotScree <- function() {
    par(mfrow=c(1,1), mar=c(4.5,4.5,1.5,4.5)+0.1)
    plot(ptaModel, scree = TRUE)
  }
  
  list(cantRemove = cantRemove,
       removeStations = removeStations,
       preArray = preArray,
       ptaModel = ptaModel,
       modelScaleSummary = modelScaleSummary,
       plotScree = plotScree)
}

models <- Map(simulateRemoval, 
              removeStations = lapply(inputsData[1:12], function(x) x$removeStations %>% na.omit),
              yearStart = lapply(inputsData[1:12], function(x) x$yearStart %>% na.omit),
              # speciesThres = lapply(inputsData, function(x) x$speciesThres %>% na.omit),
              propThres = lapply(inputsData[1:12], function(x) x$propThres %>% na.omit()),
              speciesFixed = lapply(inputsData[1:12], function(x) x$speciesFixed), 
              removeSpecies = lapply(inputsData[1:12], function(x) x$speciesRemove %>% na.omit),
              iteration = lapply(inputsData[1:12], function(x) x$iteration %>% na.omit),
              MoreArgs = list(data = data$FMWTRegion,
                              timeStep = "Survey"))

sapply(models, function(x) {
  evalMetricSum <- x$modelScaleSummary %>% 
    data.frame() %>% 
    pull(`X..Global.Pct..`) %>% 
    sum()
  
  numSpecies <- x$preArray$Species %>% 
    unique() %>% 
    length()
  
  data.frame(evalMetric = evalMetricSum,
             numSpecies = numSpecies)
}) %>% 
  t() %>% 
  data.frame() %>% 
  tibble::rownames_to_column(var = "threshold") %>% 
  mutate(across(everything(), ~as.numeric(.x))) %>% 
  pivot_longer(-threshold, values_to = "values", names_to = "metric") %>% 
  ggplot(aes(threshold, values)) +
  # geom_vline(xintercept = 0.36) +
  geom_point(size = 4) +
  facet_wrap(~metric, scales = "free_y")

# Only top X tensors:
sapply(models, function(x, n = 3) {
  evalMetricSum <- x$modelScaleSummary %>% 
    data.frame() %>% 
    arrange(-X..Global.Pct..) %>% 
    slice(1:n) %>% 
    pull(`X..Global.Pct..`) %>% 
    sum()
  
  numSpecies <- x$preArray$Species %>% 
    unique() %>% 
    length()
  
  data.frame(evalMetric = evalMetricSum,
             numSpecies = numSpecies)
}) %>% 
  t() %>% 
  data.frame() %>% 
  tibble::rownames_to_column(var = "threshold") %>% 
  mutate(across(everything(), ~as.numeric(.x))) %>% 
  pivot_longer(-threshold, values_to = "values", names_to = "metric") %>% 
  ggplot(aes(threshold, values)) +
  # geom_vline(xintercept = 0.36) +
  geom_point(size = 4) +
  facet_wrap(~metric, scales = "free_y")
# 4 or 5

# Finding tensors you want to keep
keepTensors <- function(data, start = 1, end) {
  data %>% 
    data.frame() %>% 
    arrange(-`X..Global.Pct..`) %>% 
    slice(start:end) %>% 
    pull(X.no.)
}

keep <- Map(keepTensors, 
            data = lapply(models, function(x) x$modelScaleSummary),
            end = 4)

# models$`135`$modelScaleSummary %>% .[, 5] %>% sort() %>% rev() %>% .[1:3] %>% sum()
# lapply(1:49, function(x) {
#   cat(names(models)[x], "\n")
#   
#   top3 <- models[[x]]$modelScaleSummary %>% 
#     .[, 5] %>% 
#     sort() %>% 
#     rev() %>% 
#     .[1:3] %>% 
#     sum()
#   
#   total <- models[[x]]$modelScaleSummary %>% 
#     .[, 5] %>% 
#     sort() %>% 
#     rev() %>%
#     sum()
#   
#   data.frame(name = names(models)[x],
#              top3 = top3,
#              all = total,
#              top3per = top3/100,
#              prop = top3/total)
#   
# }) %>% 
#   bind_rows() %>% 
#   pivot_longer(-c(name, top3, all), names_to = "eval", values_to = "values") %>% 
#   ggplot(aes(name, values, color = eval)) +
#   geom_point() +
#   geom_hline(yintercept = 0.5) +
#   theme(axis.text.x = element_text(size = 11))

# Creating the projection DF
createProjectionDF <- function(tensorsToKeep, model, mode = "Species") {
  
  if (is.null(names(model))) 
    stop("Model modes are not named. Either name them or provide index value.", call. = F)
  
  projectionDF <- t(model[[mode]]$v[c(tensorsToKeep),])
  rownames(projectionDF) <- model[[mode]]$n
  
  projectionDF
}

projectionDF <- Map(createProjectionDF,
                    tensorsToKeep = keep,
                    model = lapply(models, function(x) x$ptaModel))

# Begin clustering --------------------------------------------------------
optimizeCluster <- function(data, distMethod = "euclidean") {
  # data = projection DF
  
  # Calculate evaluation metrics for the various clustering methods:
  # ward.D2, single, complete, and average
  # The evaluation metrics are cophenetic correlation and Gower's distance
  evalTable <- function(data, distMethod = "euclidean") {
    
    dist1 = dist(data, method = distMethod)
    
    ##  Ward_2 clustering
    hclust.ward2<- hclust(dist1, method = "ward.D2")
    ward2.coph <- cophenetic(hclust.ward2)
    cophcorr.ward2 <- cor(dist1, ward2.coph)
    
    ##  Single linkage
    hclust.single <- hclust(dist1, method="single")
    single.coph <- cophenetic(hclust.single)
    cophcorr.single <- cor(dist1, single.coph)
    
    ##  Complete linkage
    hclust.complete <- hclust(dist1, method="complete")
    complete.coph <- cophenetic(hclust.complete)
    cophcorr.complete <- cor(dist1, complete.coph)
    
    ##  Average clustering
    hclust.avg <- hclust(dist1, method="average")
    avg.coph <- cophenetic(hclust.avg)
    cophcorr.avg <- cor(dist1, avg.coph)
    
    #####
    
    gow.dist.single <- sum((dist1 - single.coph)^2)
    gow.dist.complete <- sum((dist1 - complete.coph)^2)
    gow.dist.average <- sum((dist1 - avg.coph)^2)
    gow.dist.ward2 <- sum((dist1 - ward2.coph)^2)
    
    ###
    #####  Combine coph corr and Gower dist in a dataframe   #####
    ###
    
    fit.metrics = data.frame(clustmethod = c("single", "complete", "avg", "ward2"),
                             cophcorr = c(cophcorr.single, cophcorr.complete,
                                          cophcorr.avg, cophcorr.ward2),
                             gowdist = c(gow.dist.single, gow.dist.complete,
                                         gow.dist.average, gow.dist.ward2))
    
    bestClusterMethod <- filter(fit.metrics, cophcorr == max(cophcorr), gowdist == min(gowdist)) %>% 
      pull(clustmethod)
    
    if (length(bestClusterMethod) == 0) {
      warning("Selection based on cophenetic correlation and Gower's distance are different. Defaulting to avg", call. = F)
      bestClusterMethod <- "manual"
    }
    
    print(bestClusterMethod)
    
    if (bestClusterMethod == "single") hclustBest = hclust.single
    else if (bestClusterMethod == "complete") hclustBest = hclust.complete
    else if (bestClusterMethod == "avg") hclustBest = hclust.avg
    else if (bestClusterMethod == "ward2") hclustBest = hclust.ward2
    else if (bestClusterMethod == "manual") hclustBest = hclust.avg
    
    list(fitTable = fit.metrics,
         bestClusterMethod = bestClusterMethod,
         hclustBest = hclustBest,
         dist = dist1)
  }
  
  evalMetric <- evalTable(data) 
  
  # Calculate optimal cluster via silhouette length
  asw = numeric(nrow(data))
  # write values
  for(k in 2:(nrow(data)-1)) {
    sil = silhouette(cutree(evalMetric$hclustBest, k = k), evalMetric$dist)
    asw[k] = summary(sil)$avg.width
  }
  
  # best (largest Silhouette width)
  k.best.silhouette = which.max(asw)
  
  # Plotting the best cluster
  p <- data.frame(numClust = factor(1:nrow(data)),
                  avgSilDist = asw) %>% 
    mutate(bestSil = ifelse(avgSilDist == max(avgSilDist), T, NA)) %>% 
    {ggplot(., aes(numClust, avgSilDist, fill = bestSil)) +
        geom_col(width = 0.5, show.legend = F) +
        geom_label(data = filter(., bestSil), aes(x = 0.9 * nrow(data), 0.9 * max(avgSilDist), 
                                                  label = paste0("Best = ", numClust)), 
                   size = 10, label.padding = unit(1, "lines"), 
                   fill = "#2c2828", inherit.aes = F)}
  
  # Resulting cluster assignments to each species
  clusterAssignments <- cutree(evalMetric$hclustBest, k.best.silhouette) %>% 
    sort()
  
  taxaDendro <- as.dendrogram(evalMetric$hclustBest)
  
  list(fitTable = evalMetric$fitTable,
       bestClusterMethod = evalMetric$bestClusterMethod,
       hclustBest = evalMetric$hclustBest,
       dist = evalMetric$dist,
       k.best.silhouette = k.best.silhouette,
       p = p,
       clusterAssignments = clusterAssignments,
       taxaDendro = taxaDendro)
}

clustered <- Map(optimizeCluster,
                 data = projectionDF)

data.frame(cluster = sapply(clustered, function(x) x$k.best.silhouette),
           species = as.numeric(names(clustered))) %>% 
  ggplot(aes(x = species, cluster)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, 40, 1)) +
  scale_y_continuous(breaks = seq(1, 10, 1))

# Visualize the dendogram -------------------------------------------------
cols <- c('#01665e', '#8c510a', '#b2182b', '#5ab4ac', '#d8b365', '#d6604d', '#c7eae5', '#f6e8c3', '#f4a582','white')

evaluateCommunity <- function(dendo, k, dist, preArray, colors,
                              fullModelDendro, untangleMethod,
                              iteration = NULL,
                              plot = F,
                              ...) {
  
  if (plot) {
    
    # # old code
    # p <- dendo %>% 
    #   set("branches_col", "#878787") %>% 
    #   set("branches_lwd", 1.1) %>% 
    #   color_branches(k = k, col = colors, groupLabels = T) %>% 
    #   color_labels(col = "white") %>% 
    #   # raise.dendrogram(0.1) %>%
    #   as.ggdend() %>% 
    #   ggplot(horiz = T, offset_labels = -0.01) +
    #   # geom_point(data = what$segments, aes(x, y, color = col), show.legend = T) +
    #   theme(panel.border = element_blank(),
    #         plot.margin = margin(l = -65,
    #                              b = -15))
    
    ggdendObject <- dendo %>% 
      set("branches_col", "#878787") %>% 
      set("branches_lwd", 1.1) %>% 
      color_branches(k = k, col = colors, groupLabels = T) %>% 
      color_labels(k = k, col = colors) %>% 
      # raise.dendrogram(0.1) %>%
      as.ggdend()
    
    dendoData <- ggdendObject$labels %>% 
      left_join(data.frame(cluster = cutree(dendo, k)) %>% 
                  mutate(label = rownames(.)), by = "label") %>% 
      distinct(col, cluster) %>% 
      right_join(ggdendObject$segments, by = 'col')
    
    p <- ggplot() +
      geom_segment(data = dendoData %>% 
                     mutate(cluster = factor(cluster, levels = na.omit(sort(unique(cluster))))),
                   aes(x, y, xend = xend, yend = yend, 
                       color = cluster), 
                   size = 1.1) +
      geom_text(data = ggdendObject$labels,
                aes(x, y, label = label), color = "white", hjust = -0.05, size = 6) +
      scale_color_manual(values = unique(dendoData$col), 
                         breaks = sort(na.omit(unique(dendoData$cluster)))) +
      coord_flip() +
      scale_y_reverse(expand = c(0.2, 0)) +
      guides(color = guide_legend(nrow = 1, title.vjust = 1.3, title = "Clusters: ")) +
      theme_dendro() +
      theme(legend.position = "bottom",
            legend.margin = margin(t=-25),
            panel.border = element_blank(),
            plot.margin = margin(l = -125,
                                 b = 35, t = 15))
    
    spatialPlot <- dendo %>% 
      cutree(k = k) %>% 
      data.frame(cluster = .) %>% 
      mutate(Species = rownames(.), .before = cluster) %>% 
      right_join(preArray, by = "Species") %>% 
      # Change this to propotion across all stations...
      group_by(cluster, Species) %>% 
      mutate(propCPUE = CPUE/sum(CPUE)) %>% 
      ungroup() %>% 
      mutate(stationIndex = as.numeric(as.factor(Station))) %>% 
      # floorStation = ifelse(as.numeric(Station) <= 100, 
      #                       ceiling(as.numeric(Station)/100) * 100,
      #                       floor(as.numeric(Station)/100) * 100)) %>%
      # group_by(Station) %>% 
      # mutate(occurenceIndex = 1:n(),
      #        stationLabel = ifelse(occurenceIndex == 1, Station, NA),
      #        stationBreak = ifelse(occurenceIndex == 1, stationIndex, NA)) %>% 
      # ungroup() %>% 
      {
        ggplot(., aes(Station, propCPUE, fill = Species)) +
          geom_col(show.legend = F) +
          facet_wrap(~cluster, scales = "free_y") +
          scale_fill_viridis_d(option = "plasma") +
          # scale_x_continuous(labels = na.omit(pull(., Station)),
          #                    expand = expansion(add = (0.6))) +
          theme(panel.grid.major = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      }
    
    temporalPlot <- dendo %>% 
      cutree(k = k) %>% 
      data.frame(cluster = .) %>% 
      mutate(Species = rownames(.), .before = cluster) %>% 
      right_join(preArray, by = "Species") %>% 
      group_by(cluster, Species) %>% 
      mutate(propCPUE = CPUE/sum(CPUE)) %>% 
      ggplot(aes(timeStep, propCPUE, fill = Species)) +
      geom_col(show.legend = F) +
      facet_wrap(~cluster, scales = "free_y") +
      scale_fill_viridis_d(option = "plasma") +
      scale_x_continuous(expand = expansion(mult = 0.01), breaks = c(3:6, 1)) +
      theme(panel.grid.major = element_blank())
    
    treeCut <- cutree(dendo, k = k)
    
    sil <- silhouette(treeCut, dist)
    
    silPlot <- sil %>% 
      .[, 1:3] %>% 
      data.frame() %>% 
      mutate(Species = names(treeCut)) %>% 
      group_by(cluster) %>% 
      arrange(-sil_width, .by_group = T) %>% 
      mutate(Species = factor(Species, levels = Species)) %>% 
      ggplot(aes(cluster, sil_width, group = Species)) +
      geom_col(position = position_dodge2(preserve = "single")) +
      geom_text(aes(y = -0.01, label = Species), vjust = 0.2, hjust = 1, position = position_dodge2(width = 0.9, preserve = "single")) +
      annotate("segment", y = 0, yend = 0, x = 0.5, xend = Inf, size = rel(1), color = "#CCCCCC") +
      annotate("segment", y = 0, yend = Inf, x = 0.5, xend = 0.5, size = rel(1), color = "#CCCCCC") +
      scale_y_continuous(expand = expansion(mult = c(0.1))) +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
      coord_flip() +
      theme(axis.text.y = element_blank(), 
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.border = element_blank())
    
    # Finally to create tanglegram to full model
    plotTanglegram <- function() {
      dendlist(fullModelDendro = fullModelDendro, reducedModelDendro = dendo) %>% 
        untangle(method = untangleMethod) %>% 
        tanglegram(...)
    }
  } else {
    p <- NULL
    spatialPlot <- NULL
    temporalPlot <- NULL
    sil <- NULL
    silPlot <- NULL
    plotTanglegram <- NULL
  }
  
  # Now to pull eval metric values 
  compareCommunity <- function(fullModelDendro,
                               reducedModelDendro,
                               k) {
    
    # Calculating coph cor
    compareDendograms <- intersect_trees(fullModelDendro, 
                                         reducedModelDendro)
    
    copheneticCorrelation <- cor_cophenetic(compareDendograms, method = "kendall")
    
    # calculating prop correct clustered
    fullClustered <- cutree(fullModelDendro, k = k) %>% 
      data.frame(clusterFull = .) %>% 
      tibble::rownames_to_column(var = "Species")
    
    reducedClustered <- cutree(reducedModelDendro, k = k) %>% 
      data.frame(clusterReduced = .) %>% 
      tibble::rownames_to_column(var = "Species")
    
    
    proportionCorrect <- fullClustered %>% 
      full_join(reducedClustered,
                by = "Species") %>% 
      mutate(clusterCorrect = ifelse(clusterFull == clusterReduced, T, F)) %>% 
      group_by(clusterCorrect) %>% 
      count() %>% 
      ungroup() %>% 
      mutate(propCorrect = n/sum(n)) %>% 
      filter(clusterCorrect) %>% 
      pull(propCorrect)
    
    data.frame(copheneticCorrelation = copheneticCorrelation,
               proportionCorrect = proportionCorrect)
  }
  
  evalCompare <- compareCommunity(fullModelDendro = fullModelDendro,
                                  reducedModelDendro = dendo,
                                  k = k)
  # evalCompare <- NA
  if (!is.null(iteration)) cat(unique(iteration), "\n")
  
  list(dendogram = p,
       spatialPlot = spatialPlot,
       temporalPlot = temporalPlot,
       sil = sil,
       silPlot = silPlot,
       evalCompare = evalCompare,
       tanglegram = plotTanglegram)
}

assessmentPlots <- lapply(1:length(clustered), function(x) {
  evaluateCommunity(dendo = clustered[[x]]$taxaDendro,
                    preArray = models[[x]]$preArray,
                    fullModelDendro = clustered[[x]]$taxaDendro,
                    k = clustered[[x]]$k.best.silhouette,
                    dist = clustered[[x]]$dist,
                    plot = T,
                    untangleMethod = "labels",
                    columns_width = c(5, 2, 5), 
                    margin_inner = 12, 
                    k_labels = clustered[[x]]$k.best.silhouette,
                    colors = cols)
}) %>% 
  setNames(names(clustered))
stop()
save.image(file = "FMWTRegions.RData")
save(models, clustered, assessmentPlots, file = "FMWTRegionMain.RData")

assessmentPlots$`18`$dendogram
assessmentPlots$`18`$spatialPlot
assessmentPlots$`18`$temporalPlot
assessmentPlots$`18`$silPlot
clustered$`18`$clusterAssignments

assessmentPlots$`18`$spatialPlot$data %>% 
  filter(cluster == 3) %>% 
  ggplot(aes(stationIndex, propCPUE)) +
  geom_col() +
  facet_wrap(~Species)

# Gif of space and time timing per species --------------------------------
plotTimingSpecies <- function(dendo, preArray, numClust, species) {
  
  data <- dendo %>% 
    cutree(k = numClust) %>% 
    data.frame(cluster = .) %>% 
    mutate(Species = rownames(.), .before = cluster) %>% 
    right_join(preArray, by = "Species") %>% 
    group_by(cluster, Species) %>% 
    mutate(propCPUE = CPUE/sum(CPUE)) %>% 
    group_by(Species, cluster, timeStep) %>% 
    summarise(propCPUE = sum(propCPUE), .groups = "drop") %>% 
    mutate(speciesTarget = ifelse(Species %in% species, Species, NA))
  
  if (all(unique(data$speciesTarget) %in% NA)) {
    stop("No species was matched. Check your spelling.", call. = F)
  }
  
  dataWithoutSpecies <- data %>% 
    filter(!Species %in% species)
  dataWithSpecies <- data %>% 
    filter(Species %in% species)
  
  
  p <- ggplot(dataWithoutSpecies, aes(timeStep, propCPUE, 
                                      group = Species, 
                                      color = speciesTarget)) +
    geom_line(show.legend = F) +
    geom_line(data = dataWithSpecies, 
              aes(timeStep, propCPUE, color = Species), size = 1.3,
              show.legend = F) +
    facet_wrap(~cluster, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    scale_color_manual(values = c("#E64B35FF", "")) +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = species)
  
  print(p)
}

plotSpaceSpecies <- function(dendo, preArray, numClust, species) {
  
  data <- dendo %>% 
    cutree(k = numClust) %>% 
    data.frame(cluster = .) %>% 
    mutate(Species = rownames(.), .before = cluster) %>% 
    right_join(preArray, by = "Species") %>% 
    # Change this to propotion across all stations...
    group_by(cluster, Species) %>% 
    mutate(propCPUE = CPUE/sum(CPUE)) %>% 
    # Removing survey now
    group_by(Species, cluster, Station) %>% 
    summarise(propCPUE = sum(propCPUE), .groups = "drop") %>% 
    mutate(stationIndex = as.numeric(as.factor(Station)),
           floorStation = ifelse(as.numeric(Station) <= 100, 
                                 ceiling(as.numeric(Station)/100) * 100,
                                 floor(as.numeric(Station)/100) * 100)) %>%
    group_by(floorStation) %>% 
    mutate(occurenceIndex = 1:n(),
           stationLabel = ifelse(occurenceIndex == 1, Station, NA),
           stationBreak = ifelse(occurenceIndex == 1, stationIndex, NA),
           speciesTarget = ifelse(Species %in% species, Species, NA)) %>% 
    ungroup()
  
  if (all(unique(data$speciesTarget) %in% NA)) {
    stop("No species was matched. Check your spelling.", call. = F)
  }
  
  dataWithoutSpecies <- data %>% 
    filter(!Species %in% species)
  dataWithSpecies <- data %>% 
    filter(Species %in% species)
  
  p <- ggplot(dataWithoutSpecies, aes(stationIndex, propCPUE, 
                                      group = Species, 
                                      color = speciesTarget)) +
    geom_line(show.legend = F) +
    geom_line(data = dataWithSpecies, 
              aes(stationIndex, propCPUE, color = Species), size = 1.3,
              show.legend = F) +
    facet_wrap(~cluster, scales = "free_y") +
    scale_x_continuous(breaks = na.omit(pull(data, stationBreak)),
                       labels = na.omit(pull(data, stationLabel)),
                       expand = expansion(add = (0.6))) +
    scale_color_manual(values = c("#E64B35FF", "")) +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = species)
  
  print(p)
}

unlink(file.path(tempdir(), list.files(tempdir(), pattern = "frame.*.png$")), recursive = T)
# Timing gif
Cairo::CairoPNG(filename = file.path(tempdir(), "frame%03d.png"),
                width = 1920, height = 1080, res = 144)
for (i in cutree(clustered$`18`$hclustBest, k = clustered$`18`$k.best.silhouette) %>% sort() %>% names()) {
  plotTimingSpecies(clustered$`18`$taxaDendro,
                    models$`18`$preArray,
                    clustered$`18`$k.best.silhouette,
                    i)
}
dev.off()

gifski(paste0(tempdir(), "\\", list.files(tempdir(), pattern = "^frame.*png$")), 
       gif_file = "surveyDistributionFMWT2022.gif", width = 1920, height = 1080)

# Space gif
Cairo::CairoPNG(filename = file.path(tempdir(), "frameSpace%03d.png"),
                width = 1920, height = 1080, res = 144)
for (i in cutree(clustered$`18`$hclustBest, k = 5) %>% sort() %>% names()) {
  plotSpaceSpecies(clustered$`18`$taxaDendro,
                   models$`18`$preArray,
                   5,
                   i)
}
dev.off()

gifski(paste0(tempdir(), "\\", list.files(tempdir(), pattern = "^frameSpace.*png$")), 
       gif_file = "surveyDistributionFMWT2022_space.gif", width = 1920, height = 1080)

# Combining the gifs
library(magick)
combineGif <- function(gif1Path, gif2Path, ...) {
  gif1 <- image_read(gif1Path)
  gif2 <- image_read(gif2Path)
  
  gif1Info <- image_info(gif1)
  gif2Info <- image_info(gif2)
  
  if (nrow(gif1Info) != nrow(gif2Info)) {
    stop("Images do not have the same number of frames", call. = F)
  }
  
  newGif <- image_append(c(gif1[1], gif2[1]), ...)
  for (i in 2:nrow(gif1Info)) {
    combined <- image_append(c(gif1[i], gif2[i]), ...)
    newGif <- c(newGif, combined)
  }
  newGif
}

combinedGif <- combineGif("surveyDistributionFMWT2022.gif", "surveyDistributionFMWT2022_space.gif",
                          stack = T)
image_write_gif(combinedGif, "combinedFMWT.gif", delay = 1)
