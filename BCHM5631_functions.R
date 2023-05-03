import_peaks <- function(consensus_file_path = broadpeakfilepath) {
  peak_files <- list.files(consensus_file_path, full.names = T, pattern = ".broadPeak")
  dbp_name <- sapply(peak_files, function(x){
    y <-  str_extract(x, "([^\\/]+$)")
    # NOTE ALTERNATIVE Gsub
    # gsub("_peaks.broadPeak", "", y)
    paste(unlist(strsplit(y, "_"))[c(1,2)], collapse = "_") 
  })
  
  # NOW FOR LOOP TO IMPORT FILES
  # setting an empty list "peak_list" to be filled by for loop
  peak_list <- c()
  
  # the for loop !
  for(i in 1:length(peak_files)) {
    # Import peaks
    peaks <- rtracklayer::import(peak_files[i])
    # Append this GRanges object to the of the list.
    peak_list <- c(peak_list, peaks)
    # Name the list elements by their TF name.
    names(peak_list)[length(peak_list)] <- dbp_name[i]
  }
  return(peak_list)
}

#' CREATE CONSENSUS PEAKS
#' this function will take multiple replicate .broadPeak files (also narrow)
#' find peaks that overlap in all the replicates. 
#' @description 
#' input set of chipseq replicate peak files
#' this function then creates one merged file peaks in all samples
#' @param dbp
#' This will be extracted with names(GR_list) in the lapply at end of fun
#' You will need a "dbps" or some object for the lapply that has the 
#' name of each dbp in the named GRanges list
#' 
#' @param peak_list
#' Named list of GRanges for each chipseq replicate
#' peak_list can be generated using import_peaks function above
consensus_from_reduced <- function(dbp, peak_list) {
  dbp_peaks <- peak_list[grepl(as.character(dbp), names(peak_list))]
  suppressWarnings(all_peaks <- GenomicRanges::reduce(unlist(as(dbp_peaks, "GRangesList"))))
  all_peaks <- all_peaks[grepl("chr", seqnames(all_peaks))]
  
  # peak_exists <- lapply(dbp_peaks, function(x) {
  #   as.numeric(countOverlaps(all_peaks, x) > 0))
  # }) %>%
  # bind_rows() OR bind_cols()
  peak_exists <- matrix(NA, nrow = length(all_peaks), ncol = length(dbp_peaks))
  for(i in 1:length(dbp_peaks)) {
    suppressWarnings(peak_exists[,i] <- as.numeric(countOverlaps(all_peaks, dbp_peaks[[i]]) > 0))
  }
  # filter to consensus requiring peaks to be in all replicates
  dbp_consensus <- all_peaks[rowSums(peak_exists) == ncol(peak_exists)]
  # Required only two replicates == dbp_consensus <- all_peaks[rowSums(peak_exists) > 1]
  return(dbp_consensus)
}

# convert consensus peak files into ucsc format
ucsc_formating <- function(consensusFilePath = consensusFilePath, export_path = export_path) {
  
  consensus_file_list <- list.files(consensusFilePath, full.names = T, pattern = ".bed")

  dbps <- sapply(consensus_file_list, function(x) {
    y <- str_extract(x, "([^\\/]+$)")
    unlist(strsplit(y, "_"))[1]})
  
  peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))
  names(peaks) <- dbps
  print(length(peaks))
  canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")
  peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))
  
  headers <- paste0("track type=bed name=", names(peaks))
  new_filenames <- paste0(export_path, "/", names(peaks), ".bed")
  
  for(i in 1:length(peaks)) {
    # Write the header line
    writeLines(headers[[i]], new_filenames[[i]])
    # Append the broadPeak table data
    
    write.table(peaks[[i]], new_filenames[[i]],
                sep = "\t", col.names = FALSE, row.names = FALSE,
                quote = FALSE, append = TRUE)
  }
  
#  return(c("done?"))
}

#' function to summarize the number of events in features on each individual promoter. 
#' 
#' @description 
#' Take a gencode gtf to subset the biotype of promoters we want as a set of GRanges
#' 
#' @param features
#' set of genomic features as a GRanges object
#'  
#' @param peak_list
#' #list of peaks of dna binding proteins that will be intersected
#' 
#' @param type
#' Return either a matrix of counts over features or a binary occurrence matrix

count_peaks_per_feature <- function(features, peak_list, type = "counts") {
  
  if(!(type %in% c("counts", "occurrence"))) {
    stop("Type must be either occurrence or counts.")
  }
  
  peak_count <- matrix(numeric(), ncol = length(features), nrow = 0)
  
  for(j in 1:length(peak_list)) {
    suppressWarnings(ov <- countOverlaps(features, peak_list[[j]]))
    peak_count <- rbind(peak_count, ov)
    rownames(peak_count)[nrow(peak_count)] <- names(peak_list)[j]
    colnames(peak_count) <- features$gene_id
  }
  
  peak_matrix <- peak_count
  
  if(type == "occurrence") {
    peak_occurrence <- matrix(as.numeric(peak_count > 0), 
                              nrow = dim(peak_count)[1],
                              ncol = dim(peak_count)[2])
    rownames(peak_occurrence) <- rownames(peak_count)
    colnames(peak_occurrence) <- colnames(peak_count)
    peak_matrix <- peak_occurrence
  }
  
  return(peak_matrix)
  
}

#This function actually makes the request and returns the data only 
#(without the response headers) in a data.frame format.
encode_file_info <- function(experiment_accession,
                             base_url = "https://www.encodeproject.org/report.tsv?",
                             file_format = "fastq",
                             type = "File",
                             status = "released",
                             fields = c("accession", "read_count", "md5sum",
                                        "controlled_by", "paired_end",
                                        "paired_with", "replicate", "target")) {
  
  # Now we are creating a url that encode will understand
  path <- "report.tsv?"
  base_url <- modify_url("https://www.encodeproject.org/", path = path)
  url <- construct_query(experiment_accession,
                         base_url = base_url,
                         file_format,
                         type,
                         status,
                         fields)
  
  # this is now retrieving the data with GET function in httr and any error messages
  resp <- GET(url)
  if (http_error(resp)) {
    # error out message
    error_message <- content(resp, type = "text/html", encoding = "UTF-8") %>%
      xml_find_all("//p") %>%
      xml_text() %>%
      first()
    stop(
      # error out message
      sprintf(
        "ENCODE API request failed [%s]\n%s",
        status_code(resp),
        error_message
      ),
      call. = FALSE
    )
  }
  # another error out message
  if (http_type(resp) != "text/tsv") {
    stop("API did not return text/tsv", call. = FALSE)
  }
  body <- read_tsv(content(resp, "text"), skip = 1) %>%
    clean_names()
  return(body)
}
