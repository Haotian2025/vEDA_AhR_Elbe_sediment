# Set up directory
setwd("Zenodo/Related files for suspect screening analysis/07_Quantification&iceberg_modelling")
library(dplyr)
library(tidyr)
library(ggplot2)
library(UpSetR)
library(gridExtra)
library(ggforce)
library(cowplot)
library(magick)
library(grid)
#--------------------------------------------------------------------------------
# After alignment using Schmilka and Wittenberg, the input needs to be edited a little bit.
# Dibenz[a,h]anthracene needs to be mannually added into the ID & Area list as well as Target_list_slope
# B[b]F and B[j]F elutes in the same chromatography, therefore they share the same unknonwn ID and splite the Area equally,
# B[j]F needs to be mannually added into the lists.


df_Align <- read.csv("Alignment_Group_REP.csv", sep = ",", fileEncoding = "UTF-8")
df_target_list <- read.csv("Target_list_slope.csv", sep = ",", fileEncoding = "UTF-8")

df_align_sch_ID <- read.csv("Summary_ID_schmilka.csv", sep = ",", fileEncoding = "UTF-8")
df_align_wit_ID <- read.csv("Summary_ID_witt.csv", sep = ",", fileEncoding = "UTF-8")
df_align_sch_Area <- read.csv("Summary_Area_schmilka.csv", sep = ",", fileEncoding = "UTF-8")
df_align_wit_Area <- read.csv("Summary_Area_witt.csv", sep = ",", fileEncoding = "UTF-8")

colnames(df_align_sch_ID)[colnames(df_align_sch_ID) == "Schmilka"] <- "Unknown ID"
colnames(df_align_wit_ID)[colnames(df_align_wit_ID) == "Wittenberg"] <- "Unknown ID"
colnames(df_align_sch_Area)[colnames(df_align_sch_Area) == "Schmilka"] <- "Unknown ID"
colnames(df_align_wit_Area)[colnames(df_align_wit_Area) == "Wittenberg"] <- "Unknown ID"
# Stack data frames vertically
df_ID <- rbind(df_align_sch_ID, df_align_wit_ID)
df_Area <- rbind(df_align_sch_Area, df_align_wit_Area)


column_names <- colnames(df_ID)
column_names <- column_names[-c(1, 2)]

#create a separate data frame (df) for each col_name in column_names that contains the Chemicals, 
#and the corresponding columns from df_ID and df_Area.
#-----------------------------------------------------------------------------------------------

# Append the values from df_ID[[col_name]] to df[["UnkonwnID"]] using Chemicals as a key.
# Append the values from df_Area[[col_name]] to df[["UnkonwnID"]] using Chemicals as a key.
list_of_dfs <- list()
for (col_name in column_names) {
  # Create initial df with Chemicals
  df <- data.frame(Chemicals = df_ID$Chemicals)
  
  # Initialize columns
  df[["UnknownID"]] <- NA
  df[["Area"]] <- NA
  
  # Loop through each row
  for (i in seq_len(nrow(df))) {
    chem <- df$Chemicals[i]
    
    # Find index in df_ID
    idx_id <- match(chem, df_ID$Chemicals)
    if (!is.na(idx_id)) {
      df$UnknownID[i] <- df_ID[[col_name]][idx_id]
    }
    
    # Find index in df_Area
    idx_area <- match(chem, df_Area$Chemicals)
    if (!is.na(idx_area)) {
      df$Area[i] <- df_Area[[col_name]][idx_area]
    }
  }
  
  # Store in list
  list_of_dfs[[col_name]] <- df
}

##---------------------------------------------------------------------------------
# Loop through each data frame in the list and merge and append REP and Group info
for (name in names(list_of_dfs)) {
  df <- list_of_dfs[[name]]
  
  # Merge with selected columns from df_align
  merged_df <- merge(df, df_Align, by = "Chemicals", all.x = TRUE)
  
  # Save back to the list
  list_of_dfs[[name]] <- merged_df
}


##-------------------------------------------------------------------------------------
##Clean up each df by removing rows with Not Found in the second column

for (name in names(list_of_dfs)) {
  df <- list_of_dfs[[name]]

  # Remove rows where second column equals 'Not Found'
  df <- df[df[[2]] != 'Not Found', ]
  
  # Update the list with cleaned df
  list_of_dfs[[name]] <- df
}
#--------------------------------------------------------------------------------------
#Normalize the area using IS

df_IS <- read.csv("IS.csv", sep = ",", fileEncoding = "UTF-8")

for (name in names(list_of_dfs)) {
  df <- list_of_dfs[[name]]
  df[["Normalized_area"]] <- NA
  
  # Loop through each row
  for (i in seq_len(nrow(df))) {
    reference <- df$Reference[i]
    
    # Check if this reference chemical exists in df_IS
    match_idx <- which(df_IS$Chemicals == reference)
    
    if (length(match_idx) == 1) {
      # Extract column name from df_IS corresponding to the sample (assuming 'name' matches)
      reference_column <- match(name, colnames(df_IS))
      
      if (!is.na(reference_column)) {
        reference_value <- df_IS[match_idx, reference_column]
        
        if (!is.na(reference_value) && reference_value != 0) {  # Ensure reference_value is valid
          # Normalize the Area
          df$Normalized_area[i] <- df$Area[i] / reference_value
        } else {
          # Handle case where reference value is NA or 0
          df$Normalized_area[i] <- NA  # or some other handling like an error message
        }
      }
    } else {
      # Handle case where no unique match is found
      df$Normalized_area[i] <- NA
    }
  }
  
  # Save the updated df back to the list
  list_of_dfs[[name]] <- df
}
    
#--------------------------------------------------------------------------------------
## Calculate concentration in ACN and provide a label in a new column
# Loop through each data frame in the list
for (name in names(list_of_dfs)) {
  df <- list_of_dfs[[name]]
  
  # Initialize the new column with proper name syntax
  df[["Concentration(ng/mL ACN)"]] <- NA
  df[["Label"]] <- NA
  
  # Loop through each row
  for (i in seq_len(nrow(df))) {
    chem <- df$Chemicals[i]
    
    # Check if this chemical exists in df_target_list
    match_idx <- which(df_target_list$Chemicals == chem)
    
    if (length(match_idx) > 0) {
      area_value <- df[i, "Normalized_area"]  
      intercept_value <- df_target_list$Intercept[match_idx]
      slope_value <- df_target_list$Slope[match_idx]
      
      df[["Concentration(ng/mL ACN)"]][i] <- (area_value - intercept_value) / slope_value
      # Set label
      df[["Label"]][i] <- "in target list"
    } else {
      # No matching chemical, fallback to Area / Slope in df
      area_value <- df[i, "Normalized_area"]
      slope_value <- df$Slope[i]
      df[["Concentration(ng/mL ACN)"]][i] <- area_value / slope_value
      
      # Set label
      df[["Label"]][i] <- "not in target list"
    }
  }
  
  # Save updated df back into the list
  list_of_dfs[[name]] <- df
}


# export each df in the list as a csv file
#for (name in names(list_of_dfs)) {
#  df <- list_of_dfs[[name]]
  
  # Create filename (you can customize the path if needed)
#  filename <- paste0(name, ".csv")
  
  # Write to CSV without row names
#  write.csv(df, file = filename, row.names = FALSE)
#}




#-------------------------------------------------------------------------------------
# Calculate the concentration in the sediment using the enrichment factor
# Calculate the BEQchem by REP*Concentration in the sediment

df_EF <- read.csv("Enrichment factor.csv", sep = ",", fileEncoding = "UTF-8")

# Remove the first and last elements
list_of_dfs_sample <- list_of_dfs[-c(1, length(list_of_dfs))]


# Initialize a new list for renamed data frames (to preserve original if needed)
renamed_list <- list()

for (name in names(list_of_dfs_sample)) {
  df <- list_of_dfs_sample[[name]]
  
  # Remove first 13 characters
  temp_name <- substring(name, 14)
  
  # Remove last 4 characters if long enough
  if (nchar(temp_name) > 4) {
    filename_name <- substr(temp_name, 1, nchar(temp_name) - 4)
  } else {
    filename_name <- temp_name
  }
  
  # Assign the new name to the list
  renamed_list[[filename_name]] <- df
}


for (name in names(renamed_list)) {
  df <- renamed_list[[name]]

  df[["Concentration(ng/gsed,dw)"]] <- NA
  df[["BEQchem(ng TCDD/gsed,dw)"]] <- NA
 
  
   # Retrieve EF value from df_EF based on the current dataframe's name
  EF_value <- df_EF$EF[df_EF$Sample == name]
  
  df[["Concentration(ng/gsed,dw)"]] <- df$`Concentration(ng/mL ACN)` / EF_value
  df[["BEQchem(ng TCDD/gsed,dw)"]] <- df$`Concentration(ng/gsed,dw)` * df$`REP`
  
  # Save changes back to list
  renamed_list[[name]] <- df
}
  
  
#-----------------------------------------------------------------------------------------------
# Calculate REP with uncertainty, if chemicals are not in the target list, and REP > 2.95937E-07
# Then multiply 10.
# Recalculate BEQchem with uncertainty
# Initialize REP_uncertainty based on label.
# If label is 'in target list', assign REF.
#Else:
#  Check if REF > threshold, then multiply REF by 10.
#Else, assign REF as-is.

for (name in names(renamed_list)) {
  df <- renamed_list[[name]]
  
  
  # Remove rows with negative Concentration(ng/mL ACN)
  df <- df[df$`Concentration(ng/mL ACN)` >= 0, ]
  
  # Initialize columns
  df[["REP_uncertainty"]] <- NA
  df[["BEQchem_uncertainty(ng TCDD/gsed,dw)"]] <- NA
  df[["BEQchem(ng BaP/gsed,dw)"]] <- NA
  df[["BEQchem_uncertainty(ng BaP/gsed,dw)"]] <- NA
  
  
  
  for (i in seq_len(nrow(df))) {
    label <- df$Label[i]
    ref_value <- df$REP[i]
    
    if (label == 'in target list') {
      # Assign REF directly
      df$REP_uncertainty[i] <- ref_value
    } else {
      # Only proceed if ref_value is present
      if (!is.na(ref_value)) {
        if (ref_value > 2.95937E-07) {
          df$REP_uncertainty[i] <- ref_value * 10
        } else {
          df$REP_uncertainty[i] <- ref_value
        }
      } else {
        # If REF is NA, then set as NA or default
        df$REP_uncertainty[i] <- NA
      }
    }
  }
  
  # Compute BEQchem_uncertainty
  df[["BEQchem_uncertainty(ng TCDD/gsed,dw)"]] <- df$`Concentration(ng/gsed,dw)` * df$`REP_uncertainty`
  df[["BEQchem(ng BaP/gsed,dw)"]] <- df$`BEQchem(ng TCDD/gsed,dw)` * 35790
  df[["BEQchem_uncertainty(ng BaP/gsed,dw)"]] <- df$`BEQchem_uncertainty(ng TCDD/gsed,dw)` * 35790
  
  # Save back
  renamed_list[[name]] <- df
}


#Summary
df_summary <- data.frame(
  Sample = character(),
  BEQchem_TCDD = numeric(),
  BEQchem_TCDD_uncertainty = numeric(),
  BEQchem_BaP = numeric(),
  BEQchem_BaP_uncertainty = numeric(),
  BEQbio_TCDD = numeric(),
  BEQbio_BaP = numeric(),
  Explained_percent_TCDD = numeric(),
  Explained_percent_TCDD_uncertainty = numeric(),
  Explained_percent_BaP = numeric(),
  Explained_percent_BaP_uncertainty = numeric(),
  Uncertainty = numeric(),
  Unexplained = numeric(),
  stringsAsFactors = FALSE
)

for (name in names(renamed_list)) {
  df <- renamed_list[[name]]
  numeric_cols <- sapply(df, is.numeric)
  df[, numeric_cols] <- lapply(df[, numeric_cols], function(col) {
    col[col < 0] <- NA
    return(col)
  })
  # Calculate sums
  BEQchem_TCDD <- sum(df$`BEQchem(ng TCDD/gsed,dw)`, na.rm = TRUE)
  BEQchem_TCDD_uncertainty <- sum(df$`BEQchem_uncertainty(ng TCDD/gsed,dw)`, na.rm = TRUE)
  BEQchem_BaP <- sum(df$`BEQchem(ng BaP/gsed,dw)`, na.rm = TRUE)
  BEQchem_BaP_uncertainty <- sum(df$`BEQchem_uncertainty(ng BaP/gsed,dw)`, na.rm = TRUE)
  
  # Fetch bio estimates
  BEQbio_TCDD <- df_EF$TCDD.EQ[df_EF$Sample == name]
  BEQbio_BaP <- df_EF$BaP.EQ[df_EF$Sample == name]
  
  # Calculate explained percent
  explained_TCDD <- ifelse(length(BEQbio_TCDD) == 1 && !is.na(BEQbio_TCDD),
                           BEQchem_TCDD / BEQbio_TCDD * 100,
                           NA)
  
  explained_TCDD_uncertainty <- ifelse(length(BEQbio_TCDD) == 1 && !is.na(BEQbio_TCDD),
                           BEQchem_TCDD_uncertainty / BEQbio_TCDD *100,
                           NA)
  
  explained_BaP <- ifelse(length(BEQbio_BaP) == 1 && !is.na(BEQbio_BaP),
                          BEQchem_BaP / (BEQbio_BaP * 1000) * 100,
                          NA)
  
  explained_BaP_uncertainty <- ifelse(length(BEQbio_BaP) == 1 && !is.na(BEQbio_BaP),
                          BEQchem_BaP_uncertainty / (BEQbio_BaP * 1000) * 100,
                          NA)
  
  uncertainty <- explained_BaP_uncertainty - explained_BaP
  unexplained <- 100 - explained_BaP_uncertainty
  
  # Append to results data frame
  df_summary <- rbind(df_summary, data.frame(
    Sample = name,
    BEQchem_TCDD = BEQchem_TCDD,
    BEQchem_TCDD_uncertainty = BEQchem_TCDD_uncertainty,
    BEQchem_BaP = BEQchem_BaP,
    BEQchem_BaP_uncertainty = BEQchem_BaP_uncertainty,
    BEQbio_TCDD = BEQbio_TCDD,
    BEQbio_BaP = BEQbio_BaP,
    Explained_percent_TCDD = explained_TCDD,
    Explained_percent_TCDD_uncertainty = explained_TCDD_uncertainty,
    Explained_percent_BaP = explained_BaP,
    Explained_percent_BaP_uncertainty = explained_BaP_uncertainty,
    Uncertainty = uncertainty,
    Unexplained = unexplained,
    stringsAsFactors = FALSE
  ))
}
#-------------------------------------------------------------------------------
# Plot Figure S17
# Prioritize the top 10 AhR agonists and calculate thier individual explained percent
top10_list <- list()

for (name in names(renamed_list)) {
  df <- renamed_list[[name]]
  
  # Sum over all data for total
  BEQchem_BaP_total <- sum(df$`BEQchem(ng BaP/gsed,dw)`, na.rm = TRUE)
  
  # Select the top 10 rows with the highest 'BEQchem(ng BaP/gsed,dw)'
  top10_df <- df[order(-df$`BEQchem(ng BaP/gsed,dw)`), ]
  top10_df <- head(top10_df, 10)
  top10_df[["explained_percent"]] <- NA
  
  #top1_df <- df[which.max(df$`BEQchem(ng BaP/gsed,dw)`), ]
  
  # Sum of top 10
  BEQchem_BaP_top_10 <- sum(top10_df$`BEQchem(ng BaP/gsed,dw)`, na.rm = TRUE)
  
  top10_df[["explained_percent"]] <- top10_df$`BEQchem(ng BaP/gsed,dw)` / BEQchem_BaP_total * 100

  # Save the top 10 data frame
  top10_list[[name]] <- top10_df
}


library(gridExtra)
patch_colors <- c(
  'red',
  "dodgerblue3",
  "#4DAF4A",
  "lightpink",
  "lightcoral",
  "aquamarine",
  "lightcyan",
  "#69b3a2",
  "beige",
  "slategrey"
)
plots <- list()
names(top10_list) <- c('Prossen',
                       'Riesa',
                       'Werben',
                       'Magdeburg',
                       'Lauenburg',
                       'Wittenberg',
                       'Torgau',
                       'Dresden',
                       'Dessau',
                       'Domitz')


ordered_plot_names <- c('Prossen',
                        'Dresden',
                        'Riesa',
                        'Torgau',
                        'Wittenberg',
                        'Dessau',
                        'Magdeburg',
                        'Werben',
                        'Domitz',
                        'Lauenburg')



for (name in ordered_plot_names) {
  df <- top10_list[[name]]
  # Preserve order in legend
  df$label <- factor(df$`Chemicals`, levels = df$`Chemicals`)
  
  p <- ggplot(df, aes(x = "", y = explained_percent, fill = label)) +
    geom_bar(stat = "identity", width = 1, color = 'black', size = 0.5) +
    coord_polar(theta = "y") +
    theme_void() +
    ggtitle(name) +
    theme(legend.position = "right") +
    scale_fill_manual(values = patch_colors, name = "Chemicals")
  plots[[name]] <- p
}
grid.arrange(grobs = plots, ncol = 3)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Calculate detection frequency

# Define the start and end column names
start_col <- "X20241202_19_Schmilka.msp"
end_col <- "X20241202_28_Domitz.msp"

# Find the positions (indices) of these columns
start_idx <- which(names(df_ID) == start_col)
end_idx <- which(names(df_ID) == end_col)



# Get the sequence of column indices between start and end (inclusive)
cols_range <- start_idx:end_idx

# Get the names of these columns
subset_cols <- names(df_ID)[cols_range]

# Calculate detection frequency for each row
df_ID$Detection_Freq <- apply(df_ID[, subset_cols], 1, function(row) {
  total_cols <- length(subset_cols)  # number of columns in subset
  not_found_count <- sum(row == "Not Found")
  (total_cols - not_found_count) / total_cols
})


#Get chemicals with detection frequency > 0.9
high_freq_chemicals <- df_ID$Chemicals[df_ID$Detection_Freq >= 0.9]

#Function to extract concentration for a chemical in a dataframe
get_concentration <- function(chem_name, df) {
  row_idx <- which(df$Chemicals == chem_name)
  if(length(row_idx) == 0) {
    return(NA)
  }
  value <- df$`Concentration(ng/gsed,dw)`[row_idx]
  if(length(value) == 0 || is.na(value) || value == "") {
    return(NA)
  } else {
    return(value)
  }
}

# Create a new df to restore the results
result_df <- data.frame(Chemical = high_freq_chemicals, stringsAsFactors = FALSE)

# Loop over each dataframe to fill in concentration values
for (df_name in names(renamed_list)) {
  # Extract concentration for each chemical in the list
  concentrations <- sapply(high_freq_chemicals, function(chem) get_concentration(chem, renamed_list[[df_name]]))
  # Add as a new column
  result_df[[df_name]] <- concentrations
}


#-------------------------------------------------------------------------------
# seperate results_df by labels 'in target list' and 'not in target list'

result_in_target_list <- result_df[result_df$Chemical %in% df_target_list$Chemicals, ]
result_not_in_target_list <- result_df[!(result_df$Chemical %in% df_target_list$Chemicals), ]
# Prepare the plot for `result_in_target_list`
Quant_target <- result_in_target_list
long_data_target <- Quant_target %>%
  pivot_longer(cols = -Chemical, names_to = "Location", values_to = "Value") %>%
  left_join(df_Align[, c("Chemicals", "Group")], by = c("Chemical" = "Chemicals")) %>%
  mutate(Group = as.character(Group))

# Reorder chemicals based on median value within each group
median_by_group_target <- long_data_target %>%
  group_by(Group, Chemical) %>%
  summarize(median_value = median(Value, na.rm = TRUE)) %>%
  arrange(Group, desc(median_value))
ordered_chemicals_target <- median_by_group_target %>%
  arrange(Group, desc(median_value)) %>%
  pull(Chemical)

long_data_target$Chemical <- factor(long_data_target$Chemical, levels = unique(ordered_chemicals_target))

# Prepare the plot for `result_not_in_target_list`
Quant_notarget <- result_not_in_target_list
long_data_notarget <- Quant_notarget %>%
  pivot_longer(cols = -Chemical, names_to = "Location", values_to = "Value") %>%
  left_join(df_Align[, c("Chemicals", "Group")], by = c("Chemical" = "Chemicals")) %>%
  mutate(Group = as.character(Group))

# Reorder chemicals based on median value within each group
median_by_group_notarget <- long_data_notarget %>%
  group_by(Group, Chemical) %>%
  summarize(median_value = median(Value, na.rm = TRUE)) %>%
  arrange(Group, desc(median_value))
ordered_chemicals_notarget <- median_by_group_notarget %>%
  arrange(Group, desc(median_value)) %>%
  pull(Chemical)

long_data_notarget$Chemical <- factor(long_data_notarget$Chemical, levels = unique(ordered_chemicals_notarget))
#-------------------------------------------------------------------------------
result_df[result_df < 0] <- NA
Quant <- result_df
# Reshape data with pivot_longer
long_data <- Quant %>%
  pivot_longer(cols = -Chemical, names_to = "Location", values_to = "Value")

# Look up Group values from df_Align using Chemical names
long_data <- long_data %>%
  left_join(df_Align[, c("Chemicals", "Group")], by = c("Chemical" = "Chemicals"))

long_data$Group <- as.character(long_data$Group)

median_by_group <- long_data %>%
  group_by(Group, Chemical) %>%
  summarize(median_value = median(Value, na.rm = TRUE)) %>%
  arrange(Group, desc(median_value))

# Order chemicals within each group by descending median
ordered_chemicals <- median_by_group %>%
  arrange(Group, desc(median_value)) %>%
  pull(Chemical)

# Relevel the 'Chemical' factor based on this ordered list
long_data$Chemical <- factor(long_data$Chemical, levels = unique(ordered_chemicals))
long_data$Group <- as.character(long_data$Group)


group_labels <- c(
  "1" = "Unsubstituted PAHs",
  "2" = "Alkylated PAHs",
  "3" = "Oxygen-containing PAHs",
  "4" = "Sulfur-containing PAHs",
  "5" = "Nitrogen-containing PAHs",
  "6" = "Others"
)

group_colors <- c(
  "1" = "#69b3a2",
  "2" = "dodgerblue3",
  "3" = "#4DAF4A",
  "4" = "lightpink",
  "5" = "lightcoral",
  "6" = "aquamarine"
)
# Plot
ggplot(long_data, aes(x = Chemical, y = Value, fill = Group)) +
  geom_boxplot(color = "black") +  # border color set to black
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.5) +  # Add individual points
  scale_y_log10() +
  labs(x = "Chemical", y = "Concentration in the sediment (ng/gsed,dw)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = c(0.12, 0.2),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  scale_fill_manual(values = group_colors, name = "Chemical group", labels = group_labels)

ggsave("Quantification.jpg", dpi = 600, width = 16, height = 12, units = "in")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot Figure 5C-D
# Plot quantification and semiquantification seperately
group_labels <- c(
  "1" = "Unsubstituted PAHs",
  "2" = "Alkylated PAHs",
  "3" = "Oxygen-containing PAHs",
  "4" = "Sulfur-containing PAHs",
  "5" = "Nitrogen-containing PAHs",
  "6" = "Others"
)

group_colors <- c(
  "1" = "#FF8081",
  "2" = "#066943",
  "3" = "#B5EAD9",
  "4" = "#FFD0A0",
  "5" = "#80B7FB",
  "6" = "aquamarine"
)
# Plot
p1<-ggplot(long_data_target, aes(x = Chemical, y = Value, fill = Group)) +
  geom_boxplot(color = "black", width = 0.6) +  # border color set to black
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.5) +  # Add individual points
  scale_y_log10(limits = c(1, 8000)) +
  labs(x = "Chemical", y = "Concentration in the sediment\n(ng/gsed,dw)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "none"  # Remove here too
  ) +
  scale_fill_manual(values = group_colors, name = "Chemical group", labels = group_labels)

#ggsave("Quantification_target.jpg", dpi = 600, width = 10, height = 12, units = "in")

#---------------------------------------------------------------------------------------------
# Plot
p2<-ggplot(long_data_notarget, aes(x = Chemical, y = Value, fill = Group)) +
  geom_boxplot(color = "black", width = 0.6) +  # border color set to black
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.5) +  # Add individual points
  scale_y_log10(limits = c(1, 10000)) +
  labs(x = "Chemical", y = "Concentration in the sediment (ng/gsed,dw)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 30),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    legend.position = "none"  # Remove here too
  ) +
  scale_fill_manual(values = group_colors, name = "Chemical group", labels = group_labels)

#ggsave("Quantification_notarget.jpg", dpi = 600, width = 16, height = 12, units = "in")


plot_grid(
  p1, p2,
  ncol = 2,
  align = "h",          # Align the plots horizontally
  rel_widths = c(1, 1.4) # Adjust the relative widths
)

ggsave("Combine.jpg", dpi = 600, width = 32, height = 15, units = "in")
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
# Plot Figure S15_Intersection of AhR agonists across samples
# Replace NAs with 0
df_ID_plot <- df_ID %>% select(ends_with('.msp'))
df[is.na(df)] <- 0

# Replace "Not Found" with 0 in the entire data frame
df_ID_plot[] <- lapply(df_ID_plot, function(x) {
  if (is.character(x)) {
    x[x == "Not Found"] <- 0
  }
  return(x)
})

#Replace all non-zero values with 1
df_ID_plot[] <- lapply(df_ID_plot, function(x) {
  if (is.character(x)) as.numeric(x) else x
})

df_ID_plot[] <- lapply(df_ID_plot, function(x) {
  if (is.numeric(x)) {
    x[x != 0] <- 1
  }
  return(x)
})

# Convert all numeric columns in df to integers
df_ID_plot[] <- lapply(df_ID_plot[], function(x) {
  if (is.numeric(x)) {
    as.integer(x)
  } else {
    x
  }
})
#append chemical column
df_ID_plot$Chemicals <- df_ID$Chemicals
#reorder last first
df_ID_plot <- df_ID_plot[, c(ncol(df_ID_plot), 1:(ncol(df_ID_plot)-1))]
#rename the columns
names(df_ID_plot) <- c('Chemicals',
                       'Blank1',
                       "Prossen",
                       "Riesa",
                       "Werben",
                       "Magdeburg",
                       "Lauenburg",
                       'Wittenberg',
                       'Torgau',
                       'Dresden',
                       'Dessau',
                       'Domitz',
                       'Blank2')

# remove PBs in the second and last columns
df_ID_plot <- df_ID_plot[, -c(2, ncol(df_ID_plot))]
png("up-set_diagram.png", width = 10*300, height = 12*300, res = 300)
upset(df_ID_plot, 
      sets = c("Prossen", "Riesa", "Werben", "Magdeburg", "Lauenburg", 'Wittenberg', 'Torgau', 'Dresden', 'Dessau', 'Domitz'),
      sets.bar.color = "black",
      nintersects = 10,
      sets.x.label = 150,
      text.scale = 2.5,
      order.by = "freq",
      empty.intersections = "on")
dev.off()
#ggsave("up-set_digram.tiff", dpi = 600, width = 16, height = 12, units = "in")

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# export each df in the list as a csv file
for (name in names(renamed_list)) {
  df <- renamed_list[[name]]
  
  # Create filename (you can customize the path if needed)
  filename <- paste0(name, ".csv")
  
  # Write to CSV without row names
  write.csv(df, file = filename, row.names = FALSE)
}

write.csv(df_ID, "df_ID_output.csv", row.names = FALSE)
write.csv(df_summary, "Summary_iceberg_model.csv", row.names = FALSE)

