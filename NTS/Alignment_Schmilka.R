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
setwd("Zenodo/Related files for suspect screening analysis/06_Alignment_Schmilka")

## Import data
#-------------------------------------------------------------------------------
df <- read_csv("Alignment.csv")

# show the files in the current directory
list.files()

#Put the file name into a vector
file_list <- c(list.files())

# remove the last element in the vector as it is the R script
msp_list <- file_list[-length(file_list)]


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
targetlist <- df$`Unknown Id in Schmilka`
Schmilka <- msp[[2]]
result <- Filter(function(x) x$UnknownID %in% targetlist, Schmilka$msp_objs)
# Access the first list in msp
#fluoranthene <- msp[[1]]
# Filter the sublist items with UnknownID == 1268
#result <- Filter(function(x) x$UnknownID == 1268, fluoranthene$msp_objs)
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
colnames(Summary_ID)[colnames(Summary_ID) == "innermost_name"] <- "Schmilka"



# Create a mapping from 'Unknown Id in Schmilka' to 'name' in df
name_mapping <- setNames(df$name, df$`Unknown Id in Schmilka`)

# Map the IDs in Summary_ID to get the corresponding names
Summary_ID$Chemicals <- name_mapping[as.character(Summary_ID$Schmilka)]

# Ensure that any unmatched IDs are handled (e.g., set to NA or keep as is)
Summary_ID$Chemicals[is.na(Summary_ID$Chemicals)] <- NA

# Reorder columns if necessary
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
colnames(Summary_area)[colnames(Summary_area) == "innermost_name"] <- "Schmilka"

# Map the IDs in Summary_ID to get the corresponding names
Summary_area$Chemicals <- name_mapping[as.character(Summary_ID$Schmilka)]
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
Concentration <- c(0.1,
       0.5,
       1,
       5,
       10,
       20,
       50,
       100,
       200,
       500,
       950)

#first_row <- Summary_area[1, ]
#subset_elements <- first_row[ , 3:13]
# Convert subset_elements to numeric vector
#response_vector <- as.numeric(subset_elements)
# Fit linear model
#lm_model <- lm(response_vector ~ Concentration)

# Initialize list to store models
models_list <- list()

for (i in 1:nrow(Summary_area)) {
  # Extract row data
  data_row <- Summary_area[i, 3:13]
  data_row <- as.numeric(data_row)
  
  # Create data frame for regression
  df <- data.frame(
    Concentration = Concentration,
    Response = data_row
  )
  
  # Fit linear model
  lm_fit <- lm(Response ~ Concentration, data = df)
  
  # Get the name from the first column
  model_name <- Summary_area[i, 1]
  
  # Assign to list with the name
  models_list[[model_name]] <- lm_fit
}

# Create a list to store the results
coefficients_list <- list()

for (name in names(models_list)) {
  model <- models_list[[name]]
  # Extract coefficients
  coeffs <- coef(model)
  # Store in the list with the model name
  coefficients_list[[name]] <- coeffs
}


# Transform list into dataframe
df_coeffs <- do.call(rbind, lapply(names(coefficients_list), function(name) {
  coeffs <- coefficients_list[[name]]
  # Convert to a data frame row with model name
  data.frame(
    Model = name,
    t(coeffs),
    stringsAsFactors = FALSE
  )
}))
colnames(df_coeffs)[colnames(df_coeffs) == "Concentration"] <- "Slope"
colnames(df_coeffs)[colnames(df_coeffs) == "Model"] <- "Chemicals"
#-------------------------------------------------------------------------------

Summary_area_unknown <- Summary_area[, c(1, 14:ncol(Summary_area))]
# Merge the dataframes on 'Chemicals'
df_summary_unknown <- cbind(df_coeffs, Summary_area_unknown[ , !(names(Summary_area_unknown) %in% names(df_coeffs))])
for (col in 4:ncol( df_summary_unknown)) {
  # Extract original column name
  orig_name <- names(df_summary_unknown)[col]
  
  # Remove first 12 characters and last 4 characters
  new_name <- substr(orig_name, 10, nchar(orig_name) - 4)
  
  
  #Calculate the concentrations in real samples based on the CAL curve
  new_col_name <- paste0(new_name, "_calc")
  df_summary_unknown[[new_col_name]] <- (df_summary_unknown[[col]] - df_summary_unknown[[2]]) / df_summary_unknown[[3]]
}

write.csv(df_summary_unknown, file = "Quantification.csv", row.names = FALSE)