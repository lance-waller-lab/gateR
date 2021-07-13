#' Prepare Flow Cytometry Standard (FCS) formatted data for gateR package
#'
#' Converts a single or a collection of objects of class 'fcs' to an object of class 'data.frame' readable by \code{\link{gating}}, \code{\link{rrs}}, and \code{\link{lotrrs}} functions.
#'
#' @param path Character string of the path to the directory containing the FCS files and a metadata file. Note: The string must end with '/'
#' @param mdata Character string of the name of the metadata file in the directory. Note: The file must be in a CSV format
#' @param bl Object of class 'list'. This is a named list with column names of FCS files to exclude from the analysis
#' 
#' @return An object of class 'data.frame'. This is a named data frame with columns (features) for each condition specified in the metadata file and an ID for each row.
#'
#' @importFrom flowCore exprs getChannelMarker identifier read.flowSet
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.csv
#' @export
#' 
#' @examples
#' if (interactive()) {
#'  print(list("x" = "THIS IS A PLACEHOLDER", "y" = "DELETE THIS LINE")) # Delete after developing more informative test(s) 
#' }
#'
fcsprocessor <- function(path, mdata, bl) {
  ## Need metadata and will need to enter a directory name here
  fcspath <- path
  
  filecsv <- suppressWarnings(utils::read.csv(paste(fcspath, mdata, sep = ""), header = TRUE))
  metadata <- data.frame(lapply(filecsv, factor))
  
  fs <- flowCore::read.flowSet(path = fcspath, pattern = ".fcs")
  blacklist <- bl
  print(paste("Blacklist", blacklist, sep = ": "))
  
  
  if (length(metadata$Filename) != length(fs)) {
    stop('Rows in metadata must correspond to number of files.')
  }
  
  #Build out a dataframe of the markers and channels so we have a reference table
  channels <- colnames(fs[[1]])
  channels <- channels[!(channels %in% blacklist)]
  
  m_dim <- length(channels)
  
  markers <- as.data.frame(matrix(nrow = 2, ncol = m_dim))
  for (ch in 1:m_dim) {
    if (!is.na(match(flowCore::getChannelMarker(fs[[1]], channels[ch])$desc, markers[1, ]))) {
      firstindex <- match(flowCore::getChannelMarker(fs[[1]], channels[ch])$desc, markers[1, ])
      suffix <- gsub("[()]", "", channels[ch])
      orig_suffix <- gsub("[()]", "", channels[firstindex])
      
      markers[1, ch] <- paste(flowCore::getChannelMarker(fs[[1]], channels[ch])$desc, suffix, sep = "-")
      markers[1, firstindex] <- paste(flowCore::getChannelMarker(fs[[1]], channels[ch])$desc,
                                      orig_suffix,
                                      sep = "-")
      markers[2, ch] <- channels[ch]
      
    }
    else {
      markers[1, ch] <- flowCore::getChannelMarker(fs[[1]], channels[ch])$desc
      markers[2, ch] <- channels[ch]
    }
    
  }
  
  
  #full_data will be the data_frame for the gating function
  #manual_colnames are the columns that we will be adding
  #If a second condition is in the file, "C2" will be appended later
  full_data <- NULL
  manual_colnames <- c("ID", "C1")
  
  ### Building out the input data frame
  
  for (frame in 1:length(fs)) {
    #We go through each of the files in the flowset and create a little dataframe with the necessary info
    #We then bind all of those dataframes together into full_data
    #fs_frame holds the data from each run of the loop and is cleared out each run
    
    #All of the file data
    #We then remove the columns that are not a part of the analysis -- listed in the blacklist
    fs_frame <- as.data.frame(flowCore::exprs(fs[[frame]]))
    fs_frame <- fs_frame[, -which(names(fs_frame) %in% blacklist)]
    
    
    #find the number of cells (rows) for the file
    fs_row <- nrow(fs_frame)
    
    #Repeat the filename for the number of rows (id_col), and bind it to the frame with the header Dataset
    #the identifier() function is from flowCore and extracts the filename from the current file
    #file_path_sans_ext is from tools and removes the .fcs extension
    id_col <- rep(flowCore::identifier(fs[[frame]]), fs_row)
    fs_frame <- cbind(fs_frame, Dataset = id_col)
    file_val <- tools::file_path_sans_ext(flowCore::identifier(fs[[frame]]))
    
    #Using the metadata file we can now extract the conditions and their names
    #There is no need to rename them, they can retain their original names
    
    #C1 is needed but C2 is optional
    #We deal with C1 first
    condition1val <- metadata$C1[metadata$Filename == file_val]
    condition1col <- rep(condition1val, fs_row)
    fs_frame <- cbind(fs_frame, C1 = condition1col)
    
    #And C2 if it is present
    if ('C2' %in% colnames(metadata)) {
      num_c <- 2
      c2_control <- metadata$C2[which(metadata$C2.Control == TRUE)]
      c2_case <- metadata$C2[which(metadata$C2.Control == FALSE)]
      condition2val <- metadata$C2[metadata$Filename == file_val]
      condition2col <- rep(condition2val, fs_row)
      fs_frame <- cbind(fs_frame, C2 = condition2col)
      if (!('C2' %in% manual_colnames)) {
        manual_colnames <- append(manual_colnames, "C2")
      }
    }
    
    full_data <- rbind(full_data, fs_frame)
    
  }
  
  #Add a column for the row numbers called "ID"
  row_numbers <- seq.int(nrow(full_data))
  full_data <- cbind(full_data, ID = row_numbers, row.names = NULL)
  full_data <- move_to_start(full_data, manual_colnames)
  
  return(full_data)
}
