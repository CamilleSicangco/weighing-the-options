# Main function to clean soil water content data ##############################
clean_SWC = function(input_folder)
{
  # Clean individual data frames and combine into a single one
  ls = list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
  clean_ls = lapply(ls, clean_SWC0)
  clean_df = bind_rows(clean_ls)
  
  # Calculate the 30-minute mean
  sw = .colMeans(clean_df$sw, 2, length(clean_df$sw) / 2)
  sw0 = .colMeans(clean_df$sw0, 2, length(clean_df$sw0) / 2)
  SWC_df = clean_df[seq(1, nrow(clean_df), 2),] %>% select(-sw, -sw0)
  SWC_df = SWC_df %>% mutate(sw = sw, sw0 = sw0)
  
  return(SWC_df)
}

# Internal function for cleaning SWC data ######################################
clean_SWC0 = function(path) 
{
  # Load raw data
  headers = read.csv(path, skip = 1, header = F, nrows = 1, as.is = T)
  raw = read.csv(path, skip = 4, header = F)
  colnames(raw)= headers
  
  # Reformat date and time
  raw1 = raw %>% select(TIMESTAMP, `VW_Avg(1)`, `VW_Avg(3)`) %>%
    rename(sw0 = `VW_Avg(1)`, sw = `VW_Avg(3)`) %>%
    mutate(DateTime_hr = ymd_hms(TIMESTAMP))
  #datetime = extract_date_time(raw1)
  #raw2 = cbind(datetime, raw1)
  
  # Select data of interest
  #clean = raw2 %>% select(year, doy, hod, sw, sw0) %>%
  #  filter(doy >= 303 & doy <= 320) 
  clean = raw1 %>% select(DateTime_hr, sw, sw0)
  
  # Add variable for chamber number
  chamber_no = str_sub(path, 26, 28)
  clean = clean %>%
    mutate(chamber = rep(chamber_no, times = nrow(clean)))
  
  return(clean)
}

# Match a variable from one data frame to another data frame based on chamber and DateTime_hr
match_data = function(df_orig, df_match, var)
{
  match = sapply(1:nrow(df_orig), function(i)
  {
    j = which(df_match$chamber == df_orig$chamber[i] & 
              df_match$DateTime_hr == df_orig$DateTime_hr[i])
    var = df_match[[var]][j]
  })
  
  # Unlist if necessary and replace 0 values with NA's
  if (is.list(match)) {
    k <- !(sapply(match, length))
    # replace these values with NA
    match[k] <- NA
    # transform list to vector
    match <- unlist(match)
  } else {
    match = match
  }
  
  return(match)
}
