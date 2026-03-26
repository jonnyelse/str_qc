#!/usr/bin/env Rscript
install.packages("/usr/bin/farver_2.0.0.tar.gz", repos = NULL, type = "source")
library(farver)
install.packages("/usr/bin/labeling_0.3.tar.gz", repos = NULL, type = "source")
library(labeling)
install.packages("/usr/bin/colorspace_1.4-0.tar.gz", repos = NULL, type = "source")
library(colorspace)
install.packages("/usr/bin/munsell_0.5.0.tar.gz", repos = NULL, type = "source")
library(munsell)
install.packages("/usr/bin/RColorBrewer_1.1-2.tar.gz", repos = NULL, type = "source")
library(RColorBrewer)
install.packages("/usr/bin/viridisLite_0.3.0.tar.gz", repos = NULL, type = "source")
library(viridisLite)
install.packages("/usr/bin/scales_1.1.0.tar.gz", repos = NULL, type = "source")
library(scales)
install.packages("/usr/bin/purrr_0.3.3.tar.gz", repos = NULL, type = "source")
library(purrr)
library(dplyr)
library(stringr)
install.packages("/usr/bin/Rcpp_1.0.5.tar.gz", repos = NULL, type = "source")
library(Rcpp)
install.packages("/usr/bin/tidyr_1.1.0.tar.gz", repos = NULL, type = "source")
library(tidyr)
install.packages("/usr/bin/gtable_0.3.0.tar.gz", repos = NULL, type = "source")
library(gtable)
install.packages("/usr/bin/ggplot2_3.3.0.tar.gz", repos = NULL, type = "source")
library(ggplot2)
install.packages("/usr/bin/ggrepel_0.9.0.tar.gz", repos = NULL, type = "source")
library(ggrepel)
library(data.table)
library(psych)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
input_vcf <- args[1]
homozygosity_threshold <-  as.numeric(args[2])
adsp_threshold <- as.numeric(args[3])
work_dir <- args[4]
adsp_filepath <- args[5]
chr <- args[6]
print(work_dir)

convertGenotypesToAlleleFreq <- function(genotypeCounts) {
  # Initialize an empty list to store allele counts
  alleleCounts <- list()
  totalAlleleCount <- 0  # To store the total number of alleles
  
  # Iterate over each genotype in the vector
  for (genotype in names(genotypeCounts)) {
    count <- genotypeCounts[genotype]
    
    # Skip missing data represented by "./."
    if (genotype == "./.") next
    
    # Split the genotype to get individual alleles
    alleles <- unlist(strsplit(genotype, split = "/"))
    
    # Add counts for each allele
    for (allele in alleles) {
      # Check if the allele already has a count and update it
      if (allele %in% names(alleleCounts)) {
        alleleCounts[[allele]] <- alleleCounts[[allele]] + count
      } else {
        alleleCounts[[allele]] <- count
      }
      totalAlleleCount <- totalAlleleCount + count  # Update total allele count
    }
  }
  
  # Calculate allele frequencies
  alleleFrequencies <- sapply(alleleCounts, function(x) x / totalAlleleCount)
  
  # Convert list to a named numeric vector and return
  return(setNames(as.numeric(unlist(alleleFrequencies)), names(alleleCounts)))
}

calculate_hwep <- function(observed_counts){
  if('./.' %in% names(observed_counts)){
    no_participants <- sum(observed_counts) - unname(observed_counts['./.'])
  }
  else{no_participants <- sum(observed_counts)}
  
  alleleFrq <- convertGenotypesToAlleleFreq(observed_counts)
  homoz_pval <- list()
  for (allele in names(alleleFrq)){
    gt <- paste0(as.character(allele), "/", as.character(allele))
    if (gt %in% names(observed_counts)) {
      observed <- as.numeric(observed_counts[gt])
    } else {
      observed <- 0
    }
    btest <- binom.test(x=observed, n=no_participants, p=unname(alleleFrq[as.character(allele)])^2)
    homoz_pval[allele] <- btest$p.value
  }
  return(min(unlist(homoz_pval)))
  
}

process_chunk_to_counts <- function(file_path, skip = 0, nrows = 100) {
  df <- fread(file_path, skip = skip, nrows = nrows, header = FALSE, fill = TRUE)

  setnames(df, c("VARID", paste0("G", 2:ncol(df))))
  
  genotype_counts_list <- lapply(1:nrow(df), function(i) {
    genotypes <- as.character(df[i, -1, with = FALSE])
    genotypes <- genotypes[genotypes != "" & !is.na(genotypes)]
    table(genotypes)
  })
  
  names(genotype_counts_list) <- df$VARID
  
  return(genotype_counts_list)
}
file_path <- paste0(work_dir, "/HWE_genotypes.tsv")
command <- paste0("bcftools query -f '%VARID[\t%GT]\n' ", input_vcf, " > ", file_path)  
print(command)
system(command)



# Setup to manage file reading
total_rows <- as.integer(system(paste("wc -l", file_path, "| awk '{print $1}'"), intern = TRUE))
chunk_size <- 100  # specify chunk size
iterations <- ceiling(total_rows / chunk_size)

lowest_p <- list()
names(lowest_p) <- NULL
for (i in seq_len(iterations)) {
  skip_rows <- (i - 1) * chunk_size  # calculate the number of rows to skip
  # Process each chunk
  genotype_data_chunk <- process_chunk_to_counts(file_path, skip = skip_rows, nrows = chunk_size)
  
  # Loop through each variant in the chunk to calculate HWE p-values
  for (variant in names(genotype_data_chunk)){
    print(variant)
    lowest_p[[variant]] <- calculate_hwep(genotype_data_chunk[[variant]])
  }
  print(i)
  print("finished")
}
# Filter to include only those elements that have non-empty names
valid_entries <- lowest_p[names(lowest_p) != ""]
df <- data.frame(variant_id = names(valid_entries), lowest_p_val = unlist(valid_entries))
failed_hmz <- df %>%
  slice_min(order_by = lowest_p_val, prop = homozygosity_threshold) %>%
  pull(variant_id)
write.table(failed_hmz, "failed_homozygosity.tsv", row.names = F, col.names = F, quote = F)

#adsp
process_row <- function(row) {
  
  # Extract genotypes and adsp_values using vectorized operations
  genotypes <- row[seq(2, length(row), by = 2)]
  adsp_values <- row[seq(3, length(row), by = 2)]
  
  # Filter out any genotypes that are "./."
  valid_indices <- genotypes != "./."
  genotypes <- genotypes[valid_indices]
  adsp_values <- adsp_values[valid_indices]
  
  # Split all genotype and adsp_values at once
  split_genotypes <- strsplit(genotypes, "/")
  split_adsp_values <- strsplit(adsp_values, "/")
  
  # Create a data frame for easier manipulation
  df <- tibble(
    genotype_1 = sapply(split_genotypes, `[`, 1),
    genotype_2 = sapply(split_genotypes, `[`, 2),
    adsp_value_1 = sapply(split_adsp_values, function(x) as.numeric(x[1])),
    adsp_value_2 = sapply(split_adsp_values, function(x) as.numeric(x[2]))
  )
  
  
  # Create forward and reverse keys
  df <- df %>%
    mutate(
      forward_key = paste(genotype_1, genotype_2, sep = "_"),
      reverse_key = paste(genotype_2, genotype_1, sep = "_")
    )
  
  # Prepare the list based on unique keys
  unique_keys <- unique(c(df$forward_key, df$reverse_key))
  adsp_list <- setNames(vector("list", length(unique_keys)), unique_keys)
  adsp_list <- lapply(adsp_list, function(x) numeric(0))
  
  # Group, summarize, and apply changes
  df_grouped <- df %>%
    group_by(forward_key) %>%
    summarise(all_values = list(adsp_value_1), .groups = 'drop')
  
  walk2(df_grouped$forward_key, df_grouped$all_values, ~{ adsp_list[[.x]] <<- c(adsp_list[[.x]], unlist(.y)) })
  
  df_grouped <- df %>%
    group_by(reverse_key) %>%
    summarise(all_values = list(adsp_value_2), .groups = 'drop')
  
  walk2(df_grouped$reverse_key, df_grouped$all_values, ~{ adsp_list[[.x]] <<- c(adsp_list[[.x]], unlist(.y)) })
  
  return(adsp_list)
}


count_lines_unix <- function(file_path) {
  output <- system(paste("wc -l", shQuote(file_path)), intern = TRUE)
  as.integer(strsplit(output, "\\s+")[[1]][1])
}


chunk_size <- 100  # Define the chunk size

nrows <- count_lines_unix(adsp_filepath) 
results_list <- list()
# Process file in chunks

for (start_row in seq(1, nrows, by = chunk_size)) {
  print(start_row)
  # Read chunk of data
  data_chunk <- fread(adsp_filepath, skip = start_row-1, nrows = chunk_size, header = FALSE, sep = "\t")
  
  # Apply process_row to each row of the chunk
  list_adsp <- lapply(split(data_chunk, row.names(data_chunk)), function(row) process_row(unlist(row)))
  
  # Further analysis (e.g., statistical tests) would be here
  # Assuming some hypothetical code for analysis:
  for (idx in seq_along(list_adsp)) {
    adsp_list <- list_adsp[[as.character(idx)]]
    locus <- data_chunk$V1[idx]
    print(locus)
    all_genotypes <- unique(unlist(strsplit(gsub("\\._\\.", "", names(adsp_list)), "_")))
    for (geno_key in names(adsp_list)) {
      split_gt <- strsplit(geno_key, "_")[[1]]
      if (split_gt[1] == split_gt[2]) {
        group1 <- adsp_list[[geno_key]] / 2
        for (allele in all_genotypes){
          if (allele == split_gt[1]){next}
          group_name <- paste(split_gt[1], allele, sep="_")
          group2 <- adsp_list[[group_name]]
          # Perform Welch's t-test if both groups have sufficient data
          if (length(group1) >= 10 && length(group2) >= 10) {
            test_result <- t.test(group1, group2, var.equal = FALSE)
            gt_name <- paste(geno_key, allele, sep="_")
            # Collect p-values
            results_list[[length(results_list) + 1]] <- tibble(locus = locus, genotype = gt_name, p_value = test_result$p.value)
            }
        }
      }
    }
  }
}
results <- bind_rows(results_list)
mins <- results %>%
  group_by(locus) %>%
  summarize(min_p = min(p_value), .groups = "drop")

cutoff <- quantile(mins$min_p, probs = adsp_threshold)

failed_adsp <- results %>%
  semi_join(mins %>% filter(min_p <= cutoff), by = "locus") %>%
  group_by(locus) %>%
  filter(p_value == min(p_value)) %>%                             # keep all ties at the min
  ungroup()

test_df <- results %>% filter(locus == "chr14_102811728_102811762")
print(test_df)


# Save the results to a CSV file
write.table(failed_adsp, "failed_adsp.tsv", row.names = FALSE, quote = F)

#make QC info
path_hi_miss <- paste0(work_dir, "/QC/", chr, "_high_miss_variants.txt")
hi_miss_data <- read.table(path_hi_miss, header=F)
df_hi_miss <- data.frame("locus"=hi_miss_data$V3, "category"="hi_miss")
df_hmz <- data.frame("locus"=failed_hmz, "category"="homozygosity")
df_adsp <- data.frame("locus"=unique(failed_adsp$locus), "category"="adsp")
all_qc <- rbind(df_hi_miss, df_hmz, df_adsp)
write.table(all_qc, "all_qc.tsv", row.names=F, quote=F)

