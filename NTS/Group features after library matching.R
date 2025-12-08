# Set up directory
setwd("C:/Users/wangha/Desktop/vEDA_AhR_version_9/Zenodo/Figures/Figure S10-11/03_Library matching output after filter for grouping")
library(dplyr)

# The input data are the curated data after appending InChiKey and CID, as well as AhR classification and EC10 prediction
# List CSV files in the current directory
csv_files <- list.files(pattern = "\\.csv$")
# show the files in the current directory
#list.files()

df <- read.csv("20241202_19_Schmilka.csv", sep = ",", fileEncoding = "UTF-8")

df <- df %>%
  filter(prob >= 4)

# To reduce false positives and ensure a definitive library matching. Prob value in the library matching results can be used
# According to a literature, https://doi.org/10.1021/jasms.4c00441
# The identification probability calculation is more complicated, as it also considers the matching factors of other competing compounds, but it
# produces much better identification results, particularly if the identification probability of #1 is much higher than of #2.
# Moreover, based on the first round of identification, if the prob value is too low, the deconvoluted spectra is most likely messy.
# Therefore, a threshold of 5-10 is used here to achieve a more definitive matching result.

Unknown_ID = unique(df$Unknown_ID)
#-------------------------------------------------------------------------------
# Group by mw of candidates from library matching for the whole hitlist
# split into a list of dataframes
# Each element in 'grouped_list' is a dataframe of rows with the same mw
grouped_list <- df %>%
  group_by(mw) %>%
  group_split()

# Get modelion values for each row
group_info_list <- list()

for (i in seq_along(grouped_list)) {
  temp_df <- grouped_list[[i]]
  
  # Get unique 'Unknown_ID's and associated data
  unique_ids <- unique(temp_df$Unknown_ID)
  RTs <- RI_vals <- modelions <- modelionareas <- rep(NA, length(unique_ids))
  for (j in seq_along(unique_ids)) {
    id_val <- unique_ids[j]
    match_row <- temp_df %>% filter(Unknown_ID == id_val) %>% slice(1)
    RTs[j] <- match_row$RT
    RI_vals[j] <- match_row$RI
    modelions[j] <- match_row$modelion
    modelionareas[j] <- match_row$modelionarea
  }
  # Create combined data frame
  combined_df <- data.frame(
    Unknown_ID = unique_ids,
    RT = RTs,
    RI = RI_vals,
    modelion = modelions,
    modelionarea = modelionareas,
    stringsAsFactors = FALSE
  )
  
  # Save to list
  group_info_list[[i]] <- combined_df
}

# Loop through each sublist
for (i in seq_along(group_info_list)) {
    # Format the 'modelionarea' column to scientific notation with 2 digits
    group_info_list[[i]]$modelionarea <- formatC(group_info_list[[i]]$modelionarea, format = "e", digits = 2)
  }

#--------------------------------------------------------------------------------
# re-group the features based on modelion.
separated_sublists <- list()
# Loop over each sublist in group_info_list
for (i in seq_along(group_info_list)) {
  sub_df <- group_info_list[[i]]
  
  # Create a new column with rounded modelion values
  sub_df$rounded_modelion <- round(sub_df$modelion)
  
  # Check how many unique rounded modelion values
  unique_rounded_modelions <- unique(sub_df$rounded_modelion)
  num_unique_ions <- length(unique_rounded_modelions)
  
  if (num_unique_ions > 1) {
    # Separate rows where rounded_modelion is the same
    for (m in unique_rounded_modelions) {
      group_df <- sub_df %>% filter(rounded_modelion == m)
      separated_sublists <- append(separated_sublists, list(group_df))
    }
  } else {
    # Only one rounded modelion group
    separated_sublists <- append(separated_sublists, list(sub_df))
  }
}
#-------------------------------------------------------------------------------
# The library matching results have been clustered by molecular weight of candidates
# and 'quant mass' of features
summary_df <- do.call(rbind, lapply(seq_along(separated_sublists), function(i) {
  data.frame(
    index = i,
    row_num = nrow(separated_sublists[[i]]),  # number of rows in each sublist
    stringsAsFactors = FALSE
  )
}))

#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# If there are more than one features sharing the same quant mass, it suggests there are isomers
filtered_df <- filter(summary_df, row_num >= 2)
index = filtered_df$index
filtered_sublists <- separated_sublists[index]

# append library matching results (name, RI_extracted, RI_type_char, AhR classification) for each sublist
# in the filtered_sublists
group_info_list_combine <- list()

for (i in seq_along(filtered_sublists)) {
  temp_df <- filtered_sublists[[i]]
  unique_ids <- unique(temp_df$Unknown_ID)
  
  all_rows_list <- list()
  
  for (uid in unique_ids) {
    rows <- df %>%
      filter(Unknown_ID == uid) %>%
      select(name, RI_extracted, RI_type_char, mw, AhR_classification, EC10_prediction) %>%
      distinct(name, .keep_all = TRUE)
    all_rows_list[[as.character(uid)]] <- rows
  }
  
  combined_df <- bind_rows(all_rows_list) %>%
    distinct(name, .keep_all = TRUE) # Remove duplicates across combined data
  
  # Save unsorted combined data
  group_info_list_combine[[i]] <- combined_df
}

# sort all data frames by 'RI_extracted' ascending
group_info_list_combine <- lapply(group_info_list_combine, function(df) {
  df %>% arrange(RI_extracted)
})

#-----------------------------------------------------------------------------------------------------
# concatenate each pair of sublists at the same index in group_info_list_combine and filtered_sublists
# combine them by columns and fill empty spaces in columns with shorter length with NA
# Function to pad a data frame with NAs to reach target_rows
pad_dataframe <- function(df, target_rows) {
  n <- nrow(df)
  if (n < target_rows) {
    na_rows <- as.data.frame(matrix(NA, nrow = target_rows - n, ncol = ncol(df)))
    colnames(na_rows) <- colnames(df)
    df <- bind_rows(df, na_rows)
  }
  return(df)
}

n <- length(filtered_sublists)
summary_list <- vector("list", n)

for (i in seq_len(n)){
  df1 <- filtered_sublists[[i]]
  df2 <- group_info_list_combine[[i]]
  
  max_rows <- max(nrow(df1), nrow(df2))
  df1 <- pad_dataframe(df1, max_rows)
  df2 <- pad_dataframe(df2, max_rows)
  
  combined_df <- bind_cols(df1, df2)
  summary_list[[i]] <- combined_df
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# filter the sublists with only one feature in each sublist
filtered_df_one <- filter(summary_df, row_num == 1)
# Extract indices of sublists with exactly one row
index_one <- filtered_df_one$index

# Use these indices to filter the original list
filtered_sublists_one <- separated_sublists[index_one]

# Append library matching results to each sublist in filtered_sublists_one from df
appended_sublist <- list()

for (i in seq_along(filtered_sublists_one)) {
  temp_df <- filtered_sublists_one[[i]]
  unique_ids <- temp_df$Unknown_ID[1]
  
  all_rows_list <- list()
  
  for (uid in unique_ids) {
    rows <- df %>%
      filter(Unknown_ID == uid) %>%
      select(name, RI_extracted, RI_type_char, mw, prob, rank, RI_tolerance, AhR_classification, EC10_prediction) 
    all_rows_list[[as.character(uid)]] <- rows
  }
  
  combined_df <- bind_rows(all_rows_list)
   
  # Save unsorted combined data
  appended_sublist[[i]] <- combined_df
}


m <- length(filtered_sublists_one)
summary_list_one <- vector("list", m)

for (i in seq_len(m)){
  df1 <- filtered_sublists_one[[i]]
  df2 <- appended_sublist[[i]]
  
  max_rows <- max(nrow(df1), nrow(df2))
  df1 <- pad_dataframe(df1, max_rows)
  df2 <- pad_dataframe(df2, max_rows)
  
  combined_df <- bind_cols(df1, df2)
  summary_list_one[[i]] <- combined_df
}


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Filter out sublists where all AhR_prediction are zero for summary_list containing isomers
summary_list_filtered <- list()

for (i in seq_along(summary_list)) {
  df <- summary_list[[i]]
  # Check if all AhR_prediction are zero
  if (!all(df$AhR_classification == 0, na.rm = TRUE)) {
    summary_list_filtered[[length(summary_list_filtered) + 1]] <- df
  }
}

# Convert list of lists to a data frame
DF1 <- bind_rows(lapply(summary_list_filtered, as.data.frame))
# Write to CSV
#write.csv(DF1, "summary_list_filtered.csv", row.names = FALSE)


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Filter out sublists where all AhR_prediction are zero for summary_list containing only one feature
summary_list_one_filtered <- list()

for (i in seq_along(summary_list_one)) {
  df <- summary_list_one[[i]]
  # Check if all AhR_prediction are zero
  if (!all(df$AhR_classification == 0, na.rm = TRUE)) {
    summary_list_one_filtered[[length(summary_list_one_filtered) + 1]] <- df
  }
}

# remove duplicated sublists in summary_list_one_filtered
ids <- sapply(summary_list_one_filtered, function(df) {
  paste(capture.output(str(df)), collapse = "")
})

# Keep only unique data frames
unique_indices <- !duplicated(ids)
summary_list_unique <- summary_list_one_filtered[unique_indices]


# Convert list of lists to a data frame
DF2 <- bind_rows(lapply(summary_list_unique, as.data.frame))
# Write to CSV
write.csv(DF2, "summary_list_unique.csv", row.names = FALSE)

##--------------------------------------------------------------------------------------------
## Group statistics
# Number of elements in each list
num_separated_sublists <- length(separated_sublists)
num_summary_list_one <- length(summary_list_one)
num_summary_list <- length(summary_list)
num_summary_list_unique <- length(summary_list_unique)
num_summary_list_filtered <- length(summary_list_filtered)

# Create a data frame with one row of these values
counts_df <- data.frame(
  Total_groups = num_separated_sublists,
  One_feature_group_count = num_summary_list_one,
  Isomer_feature_group_count = num_summary_list,
  One_feature_group_active_count = num_summary_list_unique,
  Isomer_feature_group_active_count = num_summary_list_filtered
)
# Export to CSV
write.csv(counts_df, "element_counts_summary.csv", row.names = FALSE)

