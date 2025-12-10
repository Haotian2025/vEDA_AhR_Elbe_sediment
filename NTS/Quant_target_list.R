## Qantification and Semiquantification by Haotian, 20250729
##------------------------------------------------------------------------------
## Inputs are peak lists in msp files exported from MS-DAIL (4.9) 
##------------------------------------------------------------------------------
## This is used to quantify the chemicals in a target list, which need to be
## identified first in CAL1000 ppb (actually, 950 ppb). The Unknown ID needs to
## be extracted and used for alignment.
## 16_EPA_PAHs are used as a case study here.
##------------------------------------------------------------------------------
# Set up library
# Needs R version > 4.3
# Install 'mssearchr' from CRAN:
# https://cran.r-project.org/web/packages/mssearchr/index.html
install.packages("mssearchr")
library("mssearchr")
library(tidyr)
library(dplyr)
library(readr)
library(purrr)
# Set up directory here, and put all the relevant files into this directory
setwd("Zenodo/Related files for suspect screening analysis/06_Alignment_16PAHs")

## Import data
#-------------------------------------------------------------------------------
df <- read_csv("16_EPA_PAHs.csv")

# show the files in the current directory
list.files()

#Put the file name into a vector
file_list <- c(list.files())
# remove the last element in the vector as it is the R script
#msp_list <- file_list[-length(file_list)]
# Remove the first element
msp_list <- file_list[-1]

# Remove the last two elements
#msp_list <- msp_list[-(length(msp_list)-1:0)]
#msp <- '20241202_12_CAL_1000_ppb.msp'
#msp_objs <- ReadMsp(msp)
#-------------------------------------------------------------------------------
## Read msp files and target list
msp <- list()
for (i in seq_along(msp_list)) {
  file_name <- msp_list[i]
  
  #Read the MSP file
  msp_objs <- ReadMsp(file_name)
  
  #Extract UnknownID
  for (j in 1:length(msp_objs)) {
    name_value <- msp_objs[[j]]$name
    unknown_id <- sub("UnknownID=(\\d+)\\|.*", "\\1", name_value)
    msp_objs[[j]]$UnknownID <- as.numeric(unknown_id)
  }

  element_name <- basename(file_name)
  msp[[element_name]] <- list(msp_objs = msp_objs)
}

##------------------------------------------------------------------------------
##Read target list in calibration curve  
##Put unknown id in 1000 pb in a vector for alignment
targetlist <- df$`Unknown Id in CAL1000ppb`
CAL1000 <- msp[[11]]
result <- Filter(function(x) x$UnknownID %in% targetlist, CAL1000$msp_objs)
#-------------------------------------------------------------------------------
# compare each element in result with each list in msp, 
# storing all matches within the specified tolerances
# Initialize a named list to store matched elements for each list in msp
# The mass error is set up at 0.005 Da and retention time tolerance is set up at 0.03 min
matched_elements <- list()
# Loop through each named list in msp
for (i in seq_along(msp)) {
  current_list <- msp[[i]]
  all_elements <- current_list$msp_objs
  
  # Get the name of current list
  list_name <- names(msp)[i]
  matched_elements[[list_name]] <- list()
  
  # Loop over each element in result (targets)
  for (j in seq_along(result)) {
    target_model_ion <- as.numeric(result[[j]]$modelion)
    target_retention_time <- as.numeric(result[[j]]$retentiontime)
    
    # Extract UnknownID from result for naming
    unknown_id_name <- as.character(result[[j]]$UnknownID)
    
    # Find all matching elements within tolerance
    matches <- Filter(function(x) {
      # Convert x$modelion from chr to numeric
      modelion_num <- as.numeric(x$modelion)
      retentiontime <- as.numeric(x$retentiontime)
      abs(modelion_num - target_model_ion) <= 0.005 &&
        abs(retentiontime - target_retention_time) <= 0.03
    }, all_elements)
    
    
    if (length(matches) > 0) {
      # For each match, extract details into a sublist with the UnknownID as the name
      details_list <- list()
      for (m in seq_along(matches)) {
        match_obj <- matches[[m]]
        details_list[[as.character(match_obj$UnknownID)]] <- list(
          UnknownID = match_obj$UnknownID,
          retentiontime = match_obj$retentiontime,
          modelion = match_obj$modelion,
          modelionarea = match_obj$modelionarea
        )
      }
      # Store in main list under the name of the current UnknownID
      matched_elements[[list_name]][[unknown_id_name]] <- details_list
    } else {
      # If no matches, store 'Not Found'
      matched_elements[[list_name]][[unknown_id_name]] <- 'Not Found'
    }
  }
}
#-------------------------------------------------------------------------------
# Initialize a vector to store innermost list names
innermost_names <- c()

# Loop over each top-level element in your main list
for (main_str in names(matched_elements)) {
  inner_list <- matched_elements[[main_str]]
  
  # Loop over each sublist within the main list
  for (sub_str in names(inner_list)) {
    # Store the innermost list name, which is 'sub_str'
    innermost_names <- c(innermost_names, sub_str)
  }
}

# Remove duplicates if any
innermost_names <- unique(innermost_names)


# Initialize an empty list to store the results
extracted_list <- list()

for (list_name in names(matched_elements)) {
  sublist <- matched_elements[[list_name]]
  
  # Initialize a list to store the first element of innermost sublist
  first_elements <- list()
  
  # Loop over all innermost names
  for (innermost_name in innermost_names) {
    value <- NA  # Default value
    
    if (innermost_name %in% names(sublist)) {
      item <- sublist[[innermost_name]]
      
      # Handle 'Not Found'
      if (is.character(item) && item == "Not Found") {
        value <- "Not Found"
      } else {
        # Extract the first element from the innermost list or vector
        if (is.list(item)) {
          value <- item[[1]]
        } else if (is.vector(item)) {
          value <- item[1]
        }
      }
    }
    # Store the value for this innermost name
    first_elements[[innermost_name]] <- value
  }
  
  # Assign the list of first elements to the main result list
  extracted_list[[list_name]] <- first_elements
}
#-------------------------------------------------------------------------------


## Summary for UnknownID each sample
#-------------------------------------------------------------------------------
Align_ID <- list()
for (list_name in names(extracted_list)) {
  sublist <- extracted_list[[list_name]]
  
  # Initialize a vector to hold IDs for this sample
  IDs <- c()
  
  # Loop over all innermost names
  for (innermost_name in innermost_names) {
    value <- NA  # Default
    
    if (innermost_name %in% names(sublist)) {
      item <- sublist[[innermost_name]]
      
      if (is.character(item) && item == "Not Found") {
        value <- "Not Found"
      } else if (is.list(item) && "UnknownID" %in% names(item)) {
        value <- item$UnknownID
      } else if (is.list(item) && length(item) >= 1) {
        value <- item[[1]]
      } else {
        value <- NA
      }
    } else {
      value <- NA
    }
    IDs <- c(IDs, value)
  }
  
  # Save IDs for this sample, with names
  names(IDs) <- innermost_names
  Align_ID[[list_name]] <- IDs
}

df_list_ID <- list()

# Loop through each sample in Align_ID
for (sample_name in names(Align_ID)) {
  # Retrieve the vector of IDs, which is named by innermost_names
  id_vector <- Align_ID[[sample_name]]
  
  # Convert to a data frame
  df_temp <- data.frame(
    innermost_name = names(id_vector),
    stringsAsFactors = FALSE
  )
  
  # Add each sample's IDs as a column
  df_temp[[sample_name]] <- id_vector
  
  df_list_ID[[sample_name]] <- df_temp
}

# Merge all data frames by 'innermost_name' to get one combined data frame
Summary_ID <- reduce(df_list_ID, full_join, by = "innermost_name")

# Rename 'innermost_name' to just 'innermost_name' or keep it
colnames(Summary_ID)[colnames(Summary_ID) == "innermost_name"] <- "CAL1000ppb"

# Create a mapping from 'Unknown Id in 1000ppb' to 'name' in df
name_mapping <- setNames(df$name, df$`Unknown Id in CAL1000ppb`)

# Map the IDs in Summary_ID to get the corresponding names
Summary_ID$Chemicals <- name_mapping[as.character(Summary_ID$CAL1000ppb)]

# Ensure that any unmatched IDs are handled (e.g., set to NA or keep as is)
Summary_ID$Chemicals[is.na(Summary_ID$Chemicals)] <- NA
Summary_ID <- Summary_ID %>% select(Chemicals, everything())
#-------------------------------------------------------------------------------


## Summary for modelion area in each sample
#-------------------------------------------------------------------------------
Align_area <- list()
for (list_name in names(extracted_list)) {
  sublist <- extracted_list[[list_name]]
  
  # Initialize a vector to hold IDs for this sample
  IDs <- c()
  
  # Loop over all innermost names
  for (innermost_name in innermost_names) {
    value <- NA  # Default
    
    if (innermost_name %in% names(sublist)) {
      item <- sublist[[innermost_name]]
      
      if (is.character(item) && item == "Not Found") {
        value <- "Not Found"
      } else if (is.list(item) && "UnknownID" %in% names(item)) {
        value <- item$modelionarea
      } else if (is.list(item) && length(item) >= 1) {
        value <- item[[1]]
      } else {
        value <- NA
      }
    } else {
      value <- NA
    }
    IDs <- c(IDs, value)
  }
  
  # Save IDs for this sample, with names
  names(IDs) <- innermost_names
  Align_ID[[list_name]] <- IDs
}

df_list_area <- list()

# Loop through each sample in Align_ID
for (sample_name in names(Align_ID)) {
  # Retrieve the vector of IDs, which is named by innermost_names
  id_vector <- Align_ID[[sample_name]]
  
  # Convert to a data frame
  df_temp <- data.frame(
    innermost_name = names(id_vector),
    stringsAsFactors = FALSE
  )
  
  # Add each sample's IDs as a column
  df_temp[[sample_name]] <- id_vector
  
  df_list_area[[sample_name]] <- df_temp
}

# Merge all data frames by 'innermost_name' to get one combined data frame
Summary_area <- reduce(df_list_area, full_join, by = "innermost_name")

# Rename 'innermost_name' to just 'innermost_name' or keep it
colnames(Summary_area)[colnames(Summary_area) == "innermost_name"] <- "CAL1000ppb"
# Map the IDs in Summary_ID to get the corresponding names
Summary_area$Chemicals <- name_mapping[as.character(Summary_ID$CAL1000ppb)]
Summary_area <- Summary_area %>% select(Chemicals, everything())
Summary_area[Summary_area == "Not Found"] <- NA


## check the Summary_area and Summary_ID
write.csv(Summary_ID, file = "Summary_ID.csv", row.names = FALSE)
write.csv(Summary_area, file = "Summary_Area.csv", row.names = FALSE)
##------------------------------------------------------------------------------



#Still need to mannually inspect several missing values as they might have a different quant mass
#E.g., anthracene might have a quant mass of 177 instead of 178


##Linear regression for calibration curve and do quantification in samples
#-------------------------------------------------------------------------------
Summary_area[, 3:ncol(Summary_area)] <- lapply(Summary_area[, 3:ncol(Summary_area)], as.numeric)

cols_to_check <- 3:ncol(Summary_area)
rows_to_check <- c('Acenaphthene-D10', 'PCB118-13C12', 'Benzo[a]pyrene-D12')

# apply the NA imputation only to specific rows of Summary_area where the 'Chemicals' column matches entries in rows_to_check. 
# Find the indices of the matching rows
rows_indices <- which(Summary_area$Chemicals %in% rows_to_check)

# Loop through only these rows
for (i in rows_indices) {
  row_values <- as.numeric(as.character(Summary_area[i, cols_to_check]))
  na_indices <- which(is.na(row_values))
  # If there are NA values in the row
  if (length(na_indices) > 0) {
    # Compute the mean of non-NA values
    mean_value <- mean(row_values, na.rm = TRUE)
    # Replace NA with the mean
    Summary_area[i, cols_to_check][na_indices] <- mean_value
  }
}

Chemical_names <- Summary_area[[1]]
Chemical_names_filtered <- Chemical_names[!Chemical_names %in% rows_to_check]

#normalize each row in Summary_area from the 3rd column onwards based on specific reference rows, 
#determined by whether the 'Chemical' matches certain groups of elements in Chemical_names_filtered.
#If 'Chemical' matches any of the first 4 elements: divide by the Acenaphthene-D10 reference row.
#If 'Chemical' matches any of the last 6 elements: divide by the Benzo[a]pyrene-D12 reference row.


# Define the reference rows
ref_row_Acenaphthene <- which(Summary_area$Chemicals == "Acenaphthene-D10")
ref_row_BenzoPA <- which(Summary_area$Chemicals == "Benzo[a]pyrene-D12")
ref_row_PCB118 <- which(Summary_area$Chemicals == "PCB118-13C12")

# Extract the reference rows
ref_Acenaphthene <- Summary_area[ref_row_Acenaphthene, ]
ref_BenzoPA <- Summary_area[ref_row_BenzoPA, ]
ref_PCB118 <- Summary_area[ref_row_PCB118, ]
reference_df <- rbind(ref_Acenaphthene, ref_BenzoPA, ref_PCB118)
write.csv(reference_df, "reference_rows.csv", row.names = FALSE)
# Define your lists of chemicals
first_group <- head(Chemical_names_filtered, 4)
last_group  <- tail(Chemical_names_filtered, 6)

# Loop through each row for normalization
for (i in 1:nrow(Summary_area)) {
  chem_name <- Summary_area$Chemicals[i]
  
  if (chem_name %in% first_group) {
    ref_values <- ref_Acenaphthene[, 3:ncol(Summary_area)]
  } else if (chem_name %in% last_group) {
    ref_values <- ref_BenzoPA[, 3:ncol(Summary_area)]
  } else {
    # For all other chemicals, use PCB118 reference
    ref_values <- ref_PCB118[, 3:ncol(Summary_area)]
  }
  
  # Normalize only numeric columns (columns 3 onward)
  Summary_area[i, 3:ncol(Summary_area)] <- Summary_area[i, 3:ncol(Summary_area)] / as.numeric(ref_values)
}







# Initialize list to store models
models_list <- list()

for (i in 1:nrow(Summary_area)) {
  # Extract response data for the current row
  data_row <- Summary_area[i, 3:13]
  data_row <- as.numeric(data_row)
  
  # Create dataframe for regression
  df <- data.frame(
    Concentration = c(0.1, 0.5, 1, 5, 10, 20, 50, 100, 200, 500, 950),
    Response = data_row
  )
  
  # Fit linear model
  lm_fit <- lm(Response ~ Concentration, data = df)
  
  # Get model name (assumes column 1 is 'Chemicals')
  model_name <- Summary_area[i, 1]
  
  # Store the model
  models_list[[model_name]] <- lm_fit
}

# Create a list to store coefficients and R2
coefficients_list <- list()

for (name in names(models_list)) {
  model <- models_list[[name]]
  # Extract coefficients
  coeffs <- coef(model)
  # Extract R-squared
  r2_value <- summary(model)$r.squared
  # Store in list
  coefficients_list[[name]] <- list(
    coeffs = coeffs,
    r2 = r2_value
  )
}

# Transform list into a dataframe
df_coeffs <- do.call(rbind, lapply(names(coefficients_list), function(name) {
  item <- coefficients_list[[name]]
  coeffs <- item$coeffs
  r2_value <- item$r2
  # Convert to data frame row
  data.frame(
    Chemicals = name,
    Slope = coeffs["Concentration"],
    Intercept = coeffs["(Intercept)"],
    R2 = r2_value,
    stringsAsFactors = FALSE
  )
}))
colnames(df_coeffs)[colnames(df_coeffs) == "Concentration"] <- "Slope"
colnames(df_coeffs)[colnames(df_coeffs) == "Model"] <- "Chemicals"
#-------------------------------------------------------------------------------

Summary_area_unknown <- Summary_area[, c(1, 14:ncol(Summary_area))]
# Merge the dataframes on 'Chemicals'
df_summary_unknown <- cbind(df_coeffs, Summary_area_unknown[ , !(names(Summary_area_unknown) %in% names(df_coeffs))])
for (col in 5:ncol( df_summary_unknown)) {
  # Extract original column name
  orig_name <- names(df_summary_unknown)[col]
  
  # Remove first 12 characters and last 4 characters
  new_name <- substr(orig_name, 10, nchar(orig_name) - 4)
  
  
  #Calculate the concentrations in real samples based on the CAL curve
  new_col_name <- paste0(new_name, "_calc")
  df_summary_unknown[[new_col_name]] <- (df_summary_unknown[[col]] - df_summary_unknown[[3]]) / df_summary_unknown[[2]]
}

write.csv(df_summary_unknown, file = "Quantification.csv", row.names = FALSE)




