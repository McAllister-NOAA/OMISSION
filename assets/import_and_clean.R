#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/Users/mcallister/Desktop/03_Projects/03_Mitogenomes/gap_analysis/10_testing/" #Export Directory
#args[2]<-"/Users/mcallister/software/scripts/github/OMISSION/inputs/CalNEMO/nemesis_species_list_formatted.csv" #Path to input file
#args[3]<-"scientific_name" #Header text for query taxa names
########################################

library(tidyverse)
library(readxl)

setwd(as.character(args[1]))
system(paste0("cp ",as.character(args[2])," ./list_input.txt", sep = ""))

# Load in Data - expect tab-delim w/ headers
read_any <- function(file) {
  ext <- tools::file_ext(file)
  switch(tolower(ext),
         "csv"  = read.csv(file, header=TRUE, stringsAsFactors = FALSE),
         "tsv"  = read.delim(file, sep = "\t", header=TRUE, stringsAsFactors = FALSE),
         "txt"  = read.delim(file, sep = "\t", header=TRUE, stringsAsFactors = FALSE),
         "xlsx" = read_xlsx(file, col_names = TRUE),
         stop("Unknown/Unsupported file type: ", ext)
  )
}

raw_file <- read_any(as.character(args[2]))
names(raw_file)[names(raw_file) == make.names(args[3])] <- "Scientific.Name"

list_assembled_raw <- raw_file %>% mutate(original_Scientific.Name = Scientific.Name)

#remove EOL, uppercase, and remove whitespace
list_assembled_raw <- list_assembled_raw %>%
  mutate(across(everything(), ~ str_replace_all(.x, "[\r\n]", "")))

list_assembled_raw <- list_assembled_raw %>%
  mutate(Scientific.Name = str_to_upper(Scientific.Name))

list_assembled_raw <- list_assembled_raw %>%
  mutate(Scientific.Name = str_trim(Scientific.Name))

#Remove empty scientific name rows
list_assembled_raw <- list_assembled_raw %>%
  filter(!(is.na(Scientific.Name) | Scientific.Name == "NA"))

#Remove trailing special characters (*, ?)
list_assembled_raw <- list_assembled_raw %>%
  mutate(Scientific.Name = str_replace(Scientific.Name, "[*?]+$", ""))
list_assembled_raw <- list_assembled_raw %>%
  mutate(Scientific.Name = str_trim(Scientific.Name))

#Remove trailing sp, spp, sp., spp., ssp., ssp
list_assembled_raw <- list_assembled_raw %>%
  mutate(Scientific.Name = str_replace(Scientific.Name, "\\s*\\b(SSP|SP|SPP)(\\.|)?(\\s+([A-Za-z]|[0-9]+))?\\s*$", ""))
list_assembled_raw <- list_assembled_raw %>%
  mutate(Scientific.Name = str_trim(Scientific.Name))

#Add in Row IDs
list_assembled_raw <- list_assembled_raw %>%
  mutate(row_id = str_c("R", str_pad(row_number(), width = 7, pad = "0"))) %>%
  dplyr::select(row_id, everything())

#Export list
write.table(list_assembled_raw, file='list_cleaned.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names = TRUE)

# Identify entries with abbreviated genus names
genus_combos <- list_assembled_raw %>%
  filter(str_detect(Scientific.Name, "\\b\\w\\.\\s\\w+")) %>%
  mutate(Issue = "Abbreviated genus")

#Identify entries with " AND "
AND_combos <- list_assembled_raw %>% 
  filter(str_detect(Scientific.Name, " AND ")) %>%
  mutate(Issue = "AND in name")

#Identify entries with " SP "
sp_middle <- list_assembled_raw %>% 
  filter(str_detect(Scientific.Name, " (SP|SPP|SSP)\\.? ")) %>%
  mutate(Issue = "SP/SSP/SPP middle of name")

#Identify entries with " OR "
OR_combos <- list_assembled_raw %>% 
  filter(str_detect(Scientific.Name, " OR ")) %>%
  mutate(Issue = "OR in name")

#Identify entries with "CF"
CF_combos <- list_assembled_raw %>% 
  filter(str_detect(Scientific.Name, "\\sCF\\.?\\s")) %>%
  mutate(Issue = "CF in name")

#Identify single trailing end character
single_char_end <- list_assembled_raw %>%
  filter(str_detect(Scientific.Name, "\\b\\w$")) %>%
  mutate(Issue = "Single letter tail")

#Identify all weird non-alphanumeric entries
non_alphanumeric_combos <- list_assembled_raw %>%
  filter(str_detect(Scientific.Name, "[^a-zA-Z0-9\\s]")) %>%
  mutate(Issue = "Non-alphanumeric characters")

#Combine issues
all_weird <- bind_rows(genus_combos, CF_combos, AND_combos, OR_combos, single_char_end, non_alphanumeric_combos, sp_middle) %>%
  select(Scientific.Name, Issue) %>%
  distinct() %>%
  group_by(Scientific.Name) %>%
  summarise(Issue = str_c(unique(Issue), collapse = "; "), .groups = "drop")

# Print results to screen
if (nrow(all_weird) == 0) {
  message("No odd entries found.")
} else {
  message("Odd or special scientific names detected (please manually verify or fix and rerun):")
  print(all_weird)
}

write.table(all_weird, file='list_allWeird.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names = TRUE)

