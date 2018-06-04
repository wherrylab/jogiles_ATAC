library(tidyverse)

# For saving files
expGroup <- c("_") # Tag for saving files


### Read in file with peak locations and reads
# DESeq2 normalized counts with peak locations and IDs
# columns: chr, start, end, peakID, samples (normalized counts)

peaks_reads <- read_csv(file = "deseq_normcounts.csv", col_names = TRUE)
bed <- peaks_reads %>% select(chr, start, end, peakID)
bed <- write_tsv(bed, path = "unionPeakList_all.bed", col_names = FALSE) # To be used by psea_findOverlaps.sh

normCounts_sub <- peaks_reads %>% dplyr::select(-c(chr, start, end)) # To be used in RenamePSEA function


# Use bedtools intersect (psea_findOverlaps.sh) to find overlap between experiment peak list and peaksets of interest


# Read in the peakset lists

filenames <- list.files(path=getwd(), pattern = ".grp") 

list_peaksets <- list()

for(i in filenames){
  x <- read_tsv(file = i, col_names = FALSE)
  colnames(x) <- c("peakID")
  name <- tools::file_path_sans_ext(i)
  assign(name, x)
  list_peaksets[[i]] <- x
}


# Output files from psea_findOverlaps.sh

filenames <- list.files(path=getwd(), pattern = "_overlap.txt") 

list_expOverlap <- list()

for(i in filenames){
  x <- read_tsv(file = i, col_names = FALSE)
  colnames(x) <- c("chr", "start", "end", "exp_peakID", "ps_overlap_chr", "ps_overlap_start", "ps_overlap_end", "ps_peakID")
  name <- tools::file_path_sans_ext(i)
  assign(name, x)
  list_expOverlap[[i]] <- x
}


##########

RenamePSEA <- function(overlap, peakset, overlap_names, peakset_names){
  
  ps_rle <- rle(overlap$ps_peakID) # Create run length vector -- with the peakID and number of times is occurs --> use to generate new peakID
  
  overlap_M  <- overlap %>% mutate(ps_peakID_2 = ifelse(duplicated(ps_peakID), paste0(rep(ps_rle$values, times = ps_rle$lengths), "_", unlist(lapply(ps_rle$lengths, seq_len))), ps_peakID)) # peakset_peakIDs that appear more than once must be renamed to create unique IDs
  
  counts_renamed <- left_join(normCounts_sub, overlap_M, by = c("peakID" = "exp_peakID")) %>% dplyr::rename(exp_peakID = peakID) %>% dplyr::select(ps_peakID_2, 1:ncol(normCounts_sub)) %>% mutate(new_peakID = ifelse(!is.na(ps_peakID_2), ps_peakID_2, exp_peakID)) %>% dplyr::select(new_peakID, 1:ncol(.)) # Add these new peakset_peakIDs to the normalized count matrix -- if a peak in the sample dataset overlaps with a peakset peak --> use the peakset_peakID instead of the sample/experimental peakID
  
  
  overlap_peakIDs <- counts_renamed %>% select(new_peakID) %>% filter(grepl("ps",new_peakID))
  peakset_M <- peakset %>% dplyr::rename(new_peakID = peakID)
  peakset_M <- rbind(peakset_M, overlap_peakIDs) %>% unique() # Use the modified peakset_peakIDs in the "gene" set for gsea 
  
  # The IDs in "geneset" (.grp file) must be present in the normalized count data (first column in the gsea gct file) -- need add these into the normalized data --> give all samples 0 value since these areas were not called as peaks by by macs2
  missing_ps_peakIDs <- anti_join(peakset_M, counts_renamed)
  n <- nrow(missing_ps_peakIDs)
  
  missing_df <- matrix(data = NA, nrow = n, ncol = ncol(counts_renamed)) %>% data.frame()
  colnames(missing_df) <- colnames(counts_renamed) 
  missing_df[ ,1] <- missing_ps_peakIDs$new_peakID
  
  missing_df[is.na(missing_df[])] <- 0
  
  counts_renamed2 <- rbind(counts_renamed, missing_df)
  
  
  write_tsv(counts_renamed2, path = paste0("normCounts_", expGroup, overlap_names), col_names = TRUE) # File with all peaknames and counts -- for reference
  write_tsv(peakset_M, path = paste0("newPeakIDs_", peakset_names), col_names = FALSE)  # .grp peakset file for gsea
  
  
  
  
  
  # Make gsea-read gct file
  
  counts_renamed2_gsea <- counts_renamed2[ , -c(2:3)]
  counts_renamed2_gsea <- counts_renamed2_gsea %>% dplyr::rename(NAME = new_peakID) %>% mutate(Description = NAME) %>% select(NAME, Description, 2:ncol(counts_renamed2_gsea))
  
  gsea_addOns <- matrix(data = NA, nrow = 3, ncol = ncol(counts_renamed2_gsea)) %>% data.frame()
  gsea_addOns[1,1] <- "#1.2"
  gsea_addOns[2,1] <- nrow(counts_renamed2)
  gsea_addOns[2,2] <- ncol(counts_renamed2_gsea)-2
  gsea_addOns[3, ] <- colnames(counts_renamed2_gsea)
  colnames(gsea_addOns) <- colnames(counts_renamed2_gsea)
  
  counts_renamed2_gsea2 <- rbind(gsea_addOns, counts_renamed2_gsea)
  
  write_tsv(counts_renamed2_gsea2, path = paste0("normCounts_", expGroup, overlap_names, ".gct"), col_names = FALSE, na = "")
  
  
}
##################


mapply(RenamePSEA, list_expOverlap, list_peaksets, names(list_expOverlap), names(list_peaksets))

