# To do EI-MS spectra matching against NIST library in a batch mode.
# by Haotian Wang, 20250425
#-------------------------------------------------------------------------------
# Set up NIST library and NIST MS search software
# Use NIST23 in D:\Krauss\NIST23\MSSEARCH at datascience2
# Only use 2 libraries
#libraries, including:
# (1) mainlab (347100 spectra);
# (2) replib (46954 spectra);
# There are a total of 394,054 EI spectra for 347,100 compounds.
# Options->Library Search Options->Other Options->Automation, should let this on
# After several rounds of testing, the returned entries will be sorted by descending
# order of prob. when using 'Full Spectrum Search(score)' in the Spectrum Search Options.
# Details see manual (ver30Man) Page 53.
# However, the entries with highest rmf may not be ranked the top.
# Therefore, in this update, the number of hits in automation is set up at 20.
# Afterwards, the returned list is resorted by descending order of rmf, and the top 5 entries
# will be retained.
# To reduce manually manipulation of dataframe, the id numbers will also be written in text file
#-------------------------------------------------------------------------------
# Set up library
# Needs R version > 4.3
# Install 'mssearchr' from CRAN:
# https://cran.r-project.org/web/packages/mssearchr/index.html
install.packages("mssearchr")
install.packages('tidyr')
install.packages("writexl")
install.packages("openxlsx")
library(openxlsx)
library("mssearchr")
library('tidyr')
library('dplyr')
library(writexl)
#Run the function LibrarySearchUsingNistApi_change.It is in another script.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# MS-DIAL is first used to do peak picking and deconvolution.
# The peak lists are exported in a msp file format.
# To improve matching efficacy, deconvoluted spectra with S/N < 6 are removed

#Set up directory here, and put all the relevant files into this directory
setwd("D:/Haotian")

#List of MSP files to process
msp_files <- c(
  '20241202_12_CAL_1000_ppb.msp',
  '20241202_18_Procedural_Blank_1.msp',
  '20241202_19_Schmilka.msp',
  '20241202_20_Zehren.msp',
  '20241202_21_Sandau.msp',
  '20241202_22_Magdeburg.msp',
  '20241202_23_Lauenburg.msp',
  '20241202_24_Wittenberg.msp',
  '20241202_25_Torgau.msp',
  '20241202_26_Dresden.msp',
  '20241202_27_Mulde.msp',
  '20241202_28_Domitz.msp',
  '20241202_29_Procedural_Blank_2.msp'
)


#Corresponding CSV filter files for each MSP file
filter_csv_files <- c(
  '20241202_12_CAL_1000_ppb_filtered.csv',
  '20241202_18_Procedural_Blank_1_filtered.csv',
  '20241202_19_Schmilka_filtered.csv',
  '20241202_20_Zehren_filtered.csv',
  '20241202_21_Sandau_filtered.csv',
  '20241202_22_Magdeburg_filtered.csv',
  '20241202_23_Lauenburg_filtered.csv',
  '20241202_24_Wittenberg_filtered.csv',
  '20241202_25_Torgau_filtered.csv',
  '20241202_26_Dresden_filtered.csv',
  '20241202_27_Mulde_filtered.csv',
  '20241202_28_Domitz_filtered.csv',
  '20241202_29_Procedural_Blank_2_filtered.csv'
)


#Filter the msp_objs using the Unknown ID
filtered_msp <- list()
for (i in seq_along(msp_files)) {
  file_name <- msp_files[i]
  filter_csv <- filter_csv_files[i]
  
  #Read the MSP file
  msp_objs <- ReadMsp(file_name)
  
  #Extract UnknownID
  for (j in 1:length(msp_objs)) {
    name_value <- msp_objs[[j]]$name
    unknown_id <- sub("UnknownID=(\\d+)\\|.*", "\\1", name_value)
    msp_objs[[j]]$UnknownID <- as.numeric(unknown_id)
  }
  
  #Load filtered UnknownID list based on the current CSV
  filtered_ids_df <- read.csv(filter_csv)
  colnames(filtered_ids_df)[1] <- "UnknownID"
  filter_ids <- filtered_ids_df$UnknownID

  #Filter msp_objs
  filtered_msp_objs <- Filter(function(x) x$UnknownID %in% filter_ids, msp_objs)
  
  #Save or process results as needed
  filtered_msp[[file_name]] <- list(msp_objs = filtered_msp_objs)
}



# Loop through each entry in the filtered_msp list
all_hitlists <- list()
for (file_name in names(filtered_msp)) {
  
  #Access the msp_objs for the current file
  msp_objs <- filtered_msp[[file_name]]$msp_objs
  
  #Run the NIST API search
  start_time <- proc.time()
  hitlists <- LibrarySearchUsingNistApi_change(
    msp_objs,
    mssearch_dir = 'D:/Krauss/NIST23/MSSEARCH',
    temp_msp_file_dir = 'D:/Haotian',
    overwrite_spec_list = FALSE,
    comments = NULL
  )
  end_time <- proc.time()
  
  #Print the time taken
  print(paste("Processing", file_name, "took", (end_time - start_time)[3], "seconds"))
  ############################################################################################################
  #The returned RI value = X is a bitwise concatenation of RI = (X & 0x3FFF) and RI_type = (X / 0x4000)
  #assuming C/C++ language notation, all variables are integer)
  #RI_type is an ASCII code of a capital English letter:
  #A = Any Column Type (n-alkane standard),
  #N = Standard Nonpolar,
  #S =  Semi-Standard Non-Polar,
  #P = Standard Polar,
  #V = AI estimation of Semi-Standard Non-Polar
  #U = Unspecified
  ###########################################################################################################
  # Loop through each list in 'hitlists'
  for (list_idx in seq_along(hitlists)) {
    df <- hitlists[[list_idx]]
    
    # Initialize columns to store extracted RI and RI_type
    df$RI_extracted <- NA
    df$RI_type_char <- NA
    
    for (i in seq_len(nrow(df))) {
      ri_value <- as.numeric(df$ri[i])
      if (!is.na(ri_value) && length(ri_value) == 1) {
        # Extract RI (lower 14 bits)
        RI <- bitwAnd(ri_value, 0x3FFF)
        # Extract RI_type ASCII code (upper bits)
        RI_type_code <- bitwShiftR(ri_value, 14)
        # Convert ASCII code to character
        RI_type_char <- intToUtf8(RI_type_code)
        
        # Store in data frame
        df$RI_extracted[i] <- RI
        df$RI_type_char[i] <- RI_type_char
      } else {
        # NA if ri_value is missing or not scalar
        df$RI_extracted[i] <- NA
        df$RI_type_char[i] <- NA
      }
    }
    
    # Save updated data frame back to list
    hitlists[[list_idx]] <- df
  }
  all_hitlists[[file_name]] <- hitlists
}


# Save the returned list
save(all_hitlists, file = "all_hitlists.RData")
#load("all_hitlists.RData")
#------------------------------------------------------------------------------------
# Data curation for hitlists
# Loop through each list and append the $name from filtered_msp_objs to hitlists
all_reformatted_hitlists <- list()

for (file_name in names(filtered_msp)) {
  # Access the msp_objs and hitlists for the current file
  msp_objs <- filtered_msp[[file_name]]$msp_objs
  hitlists <- all_hitlists[[file_name]]
  
  #Loop through each spectrum in the current file
  for (j in seq_along(msp_objs)) {
    # Extract the spectrum's name
    name_value <- msp_objs[[j]]$name
    modelion_value <- msp_objs[[j]]$modelion
    modelionarea_value <- msp_objs[[j]]$modelionarea
    #Append 'ID' with spectrum name to the hitlist
    hitlists[[j]]$ID <- NA
    hitlists[[j]]$modelion <- NA  
    hitlists[[j]]$modelionarea <- NA
    hitlists[[j]]$ID[1] <- name_value
    hitlists[[j]]$modelion[1] <- modelion_value
    hitlists[[j]]$modelionarea[1] <- modelionarea_value
    #Reorder columns to put 'ID' first
    hitlists[[j]] <- hitlists[[j]][, c("ID", "modelion", "modelionarea", setdiff(names(hitlists[[j]]), c("ID", "modelion", "modelionarea")))]
  }
  

  #reformate each hitlist dataframe
  reformatted_list <- list()
  for (j in seq_along(hitlists)) {
    df <- hitlists[[j]]
    
    # Reformat with separate and mutate, handling the ID string
    reformatted_df <- df %>%
      separate(ID, into = c("Unknown_ID", "RT_str", "RI_str"), sep = "\\|", fill = "right") %>%
      mutate(
        RT = as.numeric(sub("RT=", "", RT_str)),
        RI = as.numeric(sub("RI=", "", RI_str)),
        `Unknown ID` = gsub("UnknownID=(\\d+)", "\\1", Unknown_ID)
      ) %>%
      select(`Unknown ID`, RT, RI, everything()) %>%
      select(-Unknown_ID, -RT_str, -RI_str)
    
    # provide a rank column
    reformatted_df$rank <- seq_len(nrow(reformatted_df))
    
    # Fill NA values with the first row's values
    if (nrow(reformatted_df) > 0) {
      first_row <- reformatted_df[1, ]
      reformatted_df <- reformatted_df %>%
        mutate(across(everything(), ~ replace_na(., first_row[[cur_column()]])))
    }
    # Save the reformatted df
    reformatted_list[[j]] <- reformatted_df
  }
  
  # Store all reformatted dataframes for this file
  all_reformatted_hitlists[[file_name]] <- reformatted_list
}


#Filter each hitlists use RI tolerance, Match factor and Reverse Match Factor
filtered_results_list <- list()
for (file_name in names(all_reformatted_hitlists)) {
  reformatted_hitlists <- all_reformatted_hitlists[[file_name]]
  
  #Combine all ranked hitlists into one data frame
  combined_df <- do.call(rbind, reformatted_hitlists)
  
  combined_df$`RI tolerance` <- abs(combined_df$RI_extracted - combined_df$RI)
  
  # Filter the features in hitlists based on multi criteria
  filtered_df <- combined_df %>%
    # Only keep rows where absolute value of RI tolerance is less then 200
    filter(!(`RI tolerance` > 200))%>%
    # From these remaining, keep rows where:
    #'mf' > 700 OR
    #'mf' < 700 AND 'rmf' > 700
    filter(mf > 700 | (mf < 700 & rmf > 700))
  
  #Save the combined hitlist (before filtering)
  filename_no_ext <- sub("\\.msp$", "", file_name)
  combined_csv <- paste0("combined_hitlist_", filename_no_ext, ".csv")
  write.csv(combined_df, combined_csv, row.names = FALSE)
  
  #Save the filtered result CSV (after filtering)
  filtered_csv <- paste0("filtered_hitlists_", filename_no_ext, ".csv")
  write.csv(filtered_df, filtered_csv, row.names = FALSE)
  
  # Save the filtered result as an Excel file
  filtered_excel <- paste0("filtered_hitlists_", filename_no_ext, ".xlsx")
  write.xlsx(filtered_df, filtered_excel)
  
  # ---- Extract and save Unknown IDs based on lib ----
  # Example: mainlib
  mainlib_ids <- filtered_df %>%
    filter(lib == "mainlib") %>%
    pull(id)
  
  writeLines(as.character(mainlib_ids), paste0("mainlib_ids_", filename_no_ext, ".txt"))
  
  # Example: replib
  replib_ids <- filtered_df %>%
    filter(lib == "replib") %>%
    pull(id)
  
  writeLines(as.character(replib_ids), paste0("replib_ids_", filename_no_ext, ".txt"))

  
  # Store the filtered df for further processing if needed
  filtered_results_list[[file_name]] <- filtered_df
}

