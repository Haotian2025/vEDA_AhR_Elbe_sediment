
# The built-in function LibrarySearchUsingNistApi has several limitations:
# (1) it can only handle 100 spectra in a msp object;
# (2) it can not handle greek characters in the name of hits. There are many names
#     containing greek characters in the mainlab, e.g., 80-75-1, Pregn-4-ene-3,20-dione, 11-hydroxy-, (11α)-
#     When the hits are these entries, the function will stop and report an error.

# Take some time to modify LibrarySearchUsingNistApi to lift the above limitations.
# Here, provide the modified function, called LibrarySearchUsingNistApi_change.
# The modified function can handle 1500 spectra each time.

LibrarySearchUsingNistApi_change <- function(msp_objs,
                                             mssearch_dir = NULL,
                                             temp_msp_file_dir = NULL,
                                             overwrite_spec_list = FALSE,
                                             comments = NULL) {
  
  #--[ Input check ]------------------------------------------------------------
  
  # System check
  if (.Platform$OS.type != "windows"){
    stop("MS Search (NIST) software is available for Windows only.")
  }
  
  # 'msp_objs'
  if (!is.list(msp_objs) || !is.list(msp_objs[[1]]) ||
      is.null(msp_objs[[1]]$mz) || is.null(msp_objs[[1]]$intst) ||
      is.null(msp_objs[[1]]$name)) {
    stop("'msp_objs' is invalid.")
  }
  if (length(msp_objs) > 2000L) {
    stop("'msp_objs' contains more than 2000 mass spectra.")     #### make a modification here from '100L' to '2000L'
  }
  
  # 'mssearch_dir'
  if (!is.null(mssearch_dir)) {
    if (!is.character(mssearch_dir) || length(mssearch_dir) != 1) {
      stop("'mssearch_dir' must be a string.")
    }
  }
  
  # 'temp_msp_file_dir'
  if (!is.null(temp_msp_file_dir)) {
    if (!is.character(temp_msp_file_dir) || length(temp_msp_file_dir) != 1) {
      stop("'temp_msp_file_dir' must be a string.")
    }
    if (!dir.exists(temp_msp_file_dir)) {
      stop("The directory ", "'", temp_msp_file_dir, "'", " does not exist.")
    }
  }
  
  # 'overwrite_spec_list'
  if (!is.logical(overwrite_spec_list) || length(overwrite_spec_list) != 1) {
    stop("'overwrite_spec_list' must be a logical value")
  }
  
  # 'comments'
  # Any R object can be passed.
  
  
  
  #--[ Getting 'mssearch_dir' ]-------------------------------------------------
  
  if (is.null(mssearch_dir)) {
    # According to the manual, 4 lines ('Path32=...', 'WorkDir32=...',
    # 'Amdis32Path=...', and  'AmdisMSPath=...') should be present in the
    # '[NISTMS]' section of the 'win.ini' file.
    
    # The method for finding the 'win.ini' file was suggested by Ivan Krylov
    win_ini_path <- file.path(Sys.getenv("WINDIR"), "win.ini")
    if (!file.exists(win_ini_path)) {
      stop("Set the 'mssearch_dir' argument manually",
           "(the 'win.ini' file was not found).")
    }
    win_ini <- readLines(win_ini_path)
    idxs <- which(win_ini == "[NISTMS]") + seq(1L, 4L)
    if (length(idxs) == 0) {
      stop("Set the 'mssearch_dir' argument manually",
           " (the '[NISTMS]' section was not found)")
    }
    temp <- grep("^Path32=", win_ini[idxs], value = TRUE)
    if (length(temp) == 0L) {
      stop("Set the 'mssearch_dir' argument manually",
           " (the 'Path32=...' line was not found")
    }
    mssearch_dir <- gsub("\\\\", "/", sub("Path32=", "", temp))
  }
  if (!dir.exists(mssearch_dir)) {
    stop("The directory ", "'", mssearch_dir, "'", " does not exist.")
  }
  
  
  
  #--[ Setting paths ]----------------------------------------------------------
  
  if (!grepl("/$", mssearch_dir)) {
    mssearch_dir <- paste0(mssearch_dir, "/")
  }
  
  nistexe_file <- paste0(mssearch_dir, "nistms$.exe")
  if (!file.exists(nistexe_file)) {
    stop("The file ", "'", basename(nistexe_file), "'", " does not exist.")
  }
  srcready_file <- paste0(mssearch_dir, "SRCREADY.TXT")
  srcreslt_file <- paste0(mssearch_dir, "SRCRESLT.TXT")
  
  # 'AUTOIMP.MSD' file (the name of this file cannot be changed)
  first_locator <- paste0(mssearch_dir, "AUTOIMP.MSD")
  
  # 'FILESPEC.FIL' file (any path and file name can be set)
  second_locator <- paste0(mssearch_dir, "FILESPEC.FIL")
  
  # 'temp.msp'
  if (is.null(temp_msp_file_dir)) {
    temp_msp_file_dir <- mssearch_dir
  }
  msp_file <- paste0(normalizePath(temp_msp_file_dir, winslash = "/"), "/",
                     "RMSSEARCH.MSP")
  
  
  
  #--[ Creating and deleting files ]--------------------------------------------
  
  # 'RMSSEARCH.MSP'
  WriteMsp(msp_objs, msp_file, fields = c("name"))
  
  # 'AUTOIMP.MSD'
  writeLines(gsub("/", "\\\\", second_locator), first_locator)
  
  # 'FILESPEC.FIL'
  if (overwrite_spec_list) {
    overwrite_par <- "OVERWRITE"
  } else {
    overwrite_par <- "APPEND"
  }
  writeLines(paste(gsub("/", "\\\\", msp_file), overwrite_par),
             second_locator)
  
  # Deleting 'SRCREADY.TXT'
  if (file.exists(srcready_file)) {
    if (!file.remove(srcready_file)) {
      stop(dirname(srcready_file), " cannot be deleted.")
    }
  }
  
 
  
  #--[ Searching ]--------------------------------------------------------------
  
  nist_par <- "/INSTRUMENT /PAR=2"
  system2(nistexe_file, nist_par)
  message("The library search is being carried out ... ", appendLF = FALSE)
  while (!file.exists(srcready_file)){
    Sys.sleep(1) # check every second
  }
  
  if (file.exists(msp_file)) {
    file.remove(msp_file)
  }
  
  temp <- readLines(srcready_file)
  n_processed_spectra <- as.integer(temp[[1]])
  if (n_processed_spectra != length(msp_objs)) {
    warning("The 'SRCREADY.TXT' contains incorrect number of spectra.")
  }
  
  res <- .ParseSrcreslt(srcreslt_file, comments)
  message("\r", "The library search has been finished.       ")
  
  #--[ Output ]-----------------------------------------------------------------
  
  return(invisible(res))

}

#==============================================================================#
#' Parse a single hitlist from a 'SRCRESLT.TXT' file
#'
#' @description
#'   Parse a single hitlist from a \emph{SRCRESLT.TXT} file.
#'
#' @param input_lines
#'   A character vector. Each element of the vector is a line from a
#'   \emph{SRCRESLT.TXT} file.
#'
#' @details
#'   It is a hidden helper, which is called from the
#'   \code{\link{.ParseSrcreslt}} function. It is supposed that some fields go
#'   in a particular order in a \emph{SRCRESLT.TXT} file: 'Name', 'Formula', any
#'   number of other fields, 'Library', and at least one more field. Moreover,
#'   the value of the 'Name', 'Formula', and 'Library' should be presented in
#'   the following format: \code{<<VALUE>>}.
#'
#' @return
#'   Return a data frame. The data frame contains the following elements (i.e.,
#'   columns): \code{name}, \code{mf}, \code{rmf}, \code{prob}, \code{lib},
#'   \code{cas}, \code{formula}, \code{mw}, \code{id}, \code{ri}. The name of an
#'   unknown compound and InLib (i.e., Compound in Library Factor) are saved as
#'   \code{unknown_name} and \code{inlib} attributes of the data frame returned.
#'
#' @noRd
#'
#==============================================================================#
.ParseSingleHitlist <- function(input_lines) {
  
  #--[ Input check ]------------------------------------------------------------
  
  # It is a hidden function. Any check is not performed.
  
  # Clean and convert input lines to a standard encoding
  input_lines <- iconv(input_lines, from = "UTF-8", to = "UTF-8", sub = "byte") #### make a modification here by adding this syntax to convert strings to a standard encoding
  
  
  #--[ Processing ]-------------------------------------------------------------
  
  field_names <- c("mf", "rmf", "prob", "cas", "mw", "id", "ri")

  
  # if (grepl("Compound in Library Factor = [0-9]+", input_lines[[1L]])) {
  #   inlib <- as.integer(sub(".*? = ([0-9]+)", "\\1", input_lines[[1L]]))
  # } else {
  #   inlib <- NA_integer_
  # }
  
  temp <- unlist(strsplit(input_lines[[1L]], "Compound in Library Factor = ",
                          fixed = TRUE))
  #> [1] "Unknown: Undecane                                                    "
  #> [2] "342"
  
  
  
  if (length(temp) == 2L) {
    unknown_name <- trimws(substring(temp[[1]], 10))
    if (temp[[2]] == "N/A") {
      inlib <- NA_integer_
    } else {
      inlib <- as.integer(temp[[2]])
    }
  } else {
    unknown_name <- NA_character_
    inlib <- NA_integer_
  }
  
  temp_res1 <- strsplit(input_lines[-1L], "<<|>>")
  #> [[1]]
  #> [1] "Hit 1  : "
  #> [2] "n-Hexane"
  #> [3] ";"
  #> [4] "C6H14"
  #> [5] "; MF: 999; RMF: 999; Prob: 89.48; CAS:110-54-3; Mw: 86; Lib: "
  #> [6] "mainlib"
  #> [7] "; Id: 27695; RI: 600."
  #> ...
  
  #print(temp_res1)
  
  temp_res2 <- lapply(temp_res1, function(x) {
    unlist(strsplit(paste0(x[[5L]], x[[7L]]), "[ ;:]+"))
    #> [1] ""         "MF"       "999"      "RMF"      "999"      "Prob"
    #> [7] "89.48"    "CAS"      "110-54-3" "Mw"       "86"       "Lib"
    #> [13] "Id"       "27695"    "RI"       "600."
  })
  
  
  
  temp_res3 <- vapply(seq_along(temp_res1), function(hit_no) {
    idxs <- match(field_names, tolower(temp_res2[[hit_no]])) + 1L
    c(temp_res1[[hit_no]][c(2L, 4L, 6L)], temp_res2[[hit_no]][idxs])
  }, character(3L + length(field_names)))
  
  
  out <- data.frame(name = temp_res3[1L, ],
                    mf = as.integer(temp_res3[4L, ]),
                    rmf = as.integer(temp_res3[5L, ]),
                    prob = as.numeric(temp_res3[6L, ]),
                    lib = temp_res3[3L, ],
                    cas = temp_res3[7L, ],
                    formula = temp_res3[2L, ],
                    mw = as.integer(temp_res3[8L, ]),
                    id = as.integer(temp_res3[9L, ]),
                    ri = as.numeric(temp_res3[10L, ]))
  

  
  attr(out, "unknown_name") <- unknown_name
  attr(out, "inlib") <- inlib
  
  
  
  #--[ Output ]-----------------------------------------------------------------
  
  return(out)
}





#==============================================================================#
#' Parse a 'SRCRESLT.TXT' file
#'
#' @description
#'   Parse a \emph{SRCRESLT.TXT} file containing the results of the library
#'   search performed using the MS Search (NIST) software.
#'
#' @param input_file
#'    A string. The full path to the the \emph{SRCRESLT.TXT} file.
#' @inheritParams LibrarySearchUsingNistApi
#'
#' @inherit LibrarySearchUsingNistApi return
#'
#' @noRd
#'
#==============================================================================#
.ParseSrcreslt <- function(input_file,
                           comments = NULL) {
  
  #--[ Input check ]------------------------------------------------------------
  
  # 'input_file'
  if (!is.character(input_file) || length(input_file) != 1) {
    stop("'input_file' must be a string.")
  }
  if (!file.exists(input_file)) {
    stop("The file ", "'", input_file, "'", " does not exist.")
  }
  
  # 'comments'
  # Any R object can be passed.
  
  
  
  #--[ Reading and preprocessing ]----------------------------------------------
  
  all_lines <- readLines(input_file)
  
  # Convert to standard encoding
  all_lines <- iconv(all_lines, from = "UTF-8", to = "UTF-8", sub = "byte")  #### make a modification here by adding this syntax to convert strings to a standard encoding
  
  first_line_idxs <- grep("^Unknown", all_lines, ignore.case = TRUE)
  if (length(first_line_idxs) == 0L) {
    stop("The 'SRCREADY.TXT' file is invalid.")
  }
  last_line_idxs <- c(first_line_idxs[-1] - 1L, length(all_lines))
  
  
  
  #--[ Parsing ]----------------------------------------------------------------
  
  hitlists <- lapply(seq_along(first_line_idxs), function(res_no) {
    idxs <- seq(first_line_idxs[[res_no]], last_line_idxs[[res_no]])
    out <- .ParseSingleHitlist(all_lines[idxs])
    return(out)
  })
  attr(hitlists, "comments") <- comments
  
  
  
  #--[ Output ]-----------------------------------------------------------------
  
  return(hitlists)
}
