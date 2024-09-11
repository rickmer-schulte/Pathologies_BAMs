# Data Reading and Preparation
load_covid_dat <- function(data_folder_path, sf=FALSE){
df <- read.csv(paste0(data_folder_path, "/aggregated_us_lev2_subset.csv"))
df <- df %>% filter(new_confirmed >= 0 & date <= 863)
df <- df %>% 
  mutate(
    prevalence = log(prevalence*1000 + 1),
    population = log(population)) %>%
  rename(
    temp = average_temperature_celsius,
    humid = relative_humidity
  )
df <- df %>% filter(date <= 726)

# Convert days since 2020-01-01 to Date
df$dates <- as.Date("2020-01-01") + df$date

# Select only data for San Francisco if desired
if(sf==TRUE){df <- df %>% filter(latitude == 37.44000 & longitude == -122.36)}

return(df)
}