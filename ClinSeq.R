# Install packages
# install.packages("tidyverse")
# install.packages("gridExtra")
# install.packages("patchwork")
# install.packages("caret")
# install.packages("extrafont")

# Load necessary libraries
library(tidyverse)
library(gridExtra)
library(patchwork)
library(extrafont)
library(caret)

#set wd
setwd("*SET WORKING DIRECTORY*") #i.e. C:/Users/John/ClinSeq/

# Clear the environment
rm(list = ls())

####################################################################################
# Define the variables at the start of the script

number_of_rows <- 6 # Default 5

weight_factors <- list(
  NormCladeReads = 2.0,  # Default 2
  Kmers = 1.0,           # Default 1
  KmerDuplicity = 1.0,    # Default 1
  KmerCoverage = 1.0      # Default 1
)

contaminant_threshold <- 0.90

SD_factor <- 0.6 # Default 0.6
CV_factor <- 100 # Default 100


#exclude_sources from analysis
exclude_sources <- c('A35', 1355,1356,1357,1358,1359,1360,1361,1362,1909,1912,1913,1920,1930,1931,1933,1934,1938,1943,1945,1946,
                     1947,1948,1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,3918,3928,3929,3956,4448,3400,3401,3402,3403,3404,
                     3405,3406,3407,3437,3438,3439,3440,3441,3442,3443,3444,3445,3446,3447,3523,3524,3525,3526,3527,3528,3529,3530,
                     3531,3532,3533,3534,3535,3536,3539,3540
                     )
#included in analysis but excluded from results
filter_sources <- c(1363,1364,1365,1368,1370,1910,1911,1914,1915,
                     1916,1917,1918,1919,1921,1922,1923,1924,1925,
                     1926,1927,1928,1929,4399,4400,4401,4402,4403,
                     4404,4405,4406,4407,4408,4409,4410,4411,4442,
                     4379,4446,4439,4377,4433,1366,1367,1369,1932,
                     1936,1937,1939,1940,1941,1942,4412,4413,4414,
                     4415,4416,4417,4418,4419,4420,4421,4422,4423,
                     4424,4425,4426,4427,4428,1944,4381,4382,4383,
                     4384,4385,4386,4387,4388,4434,4437,4378,4438,
                     3930,1935,4389,4390,4391,4392,4393,4394,4395,
                     4396,4397,4398,4374,4375,4376,3600,3601,3602,
                     3603,3783,3784,3791,3802,3810,3818,3828,3834,
                     3835,3841,3842,3932,4189,4190,4191,4193,4194,
                     4196,3954,3919,3921,3926,3946,3950,3951,3952,
                     3955,4198,4199,4200,4201,3935,3945,3957,3958,
                     3944,3408,3537,3596,3923,3924,3608,3609,3610,
                     3611,3782,3789,3790,3797,3801,3808,3809,3816,
                     3817,3826,3827,3833,3840,4187,4188,4192,4195,
                     4202,4203,4204,3597,3598,3788,3787,3599,3830,
                     3837,3825,3832,3604,3605,3606,3607,3785,3792,
                     3793,3803,3804,3811,3812,3819,3829,3836,3843,
                     4205,4206,3612,3613,3614,3615,3959,3960,3920,
                     3925,3927)

input_file_path <- "*INPUT FILE PATH*" #i.e. C:/Users/John/ClinSeq/clinseq.csv

####################################################################################
# target pathogens to mark i.e. pathogens known to cause meningitis

marked_pathogen_viruses <- c(
  "s_Enterovirus",
  "s_Enterovirus A",
  "s_Enterovirus B",
  "s_Enterovirus coxsackiepol",
  "s_Orthoflavivirus nilense",
  "s_Orthoflavivirus",
  "s_Orthoflavivirus louisense",
  "s_Alphavirus eastern",
  "s_Alphavirus western",
  "s_Orthoflavivirus japonicum",
  "s_Orthoflavivirus powassanense",
  "s_Varicellovirus",
  "s_Simplexvirus",
  "s_Muromegalovirus",
  "s_Simplexvirus humanalpha1",
  "s_Simplexvirus humanalpha2",
  "s_Varicellovirus humanalpha3",
  "s_Lymphocryptovirus humangamma4",
  "s_Cytomegalovirus humanbeta5",
  "s_Orthorubulavirus parotitidis",
  "s_Lentivirus humimdef1",
  "s_Lentivirus humimdef2",
  "s_Mammarenavirus choriomeningitidis",
  "s_Parechovirus ahumpari",
  "s_Parechovirus beljungani",
  "s_Parechovirus cebokele",
  "s_Parechovirus deferreti",
  "s_Parechovirus efalco",
  "s_Parechovirus feterobo",
  "s_Pasivirus agallia",
  "s_Deltaretrovirus priTlym1",
  "s_Deltaretrovirus priTlym2",
  "s_Rubivirus rubellae",
  "s_Orthoflavivirus denguei",
  "s_Alphavirus chikungunya",
  "s_Morbillivirus hominis",
  "s_Lyssavirus rabies",
  "s_Lymphocytic"
)

marked_pathogen_bacteria <- c(
  "s_Mycobacterium",
  "s_Leptospira",
  "s_Borrelia",
  "s_Treponema",
  "s_Brucella",
  "s_Mycoplasma",
  "s_Cryptococcus",
  "s_Coccidioides",
  "s_Histoplasma",
  "s_Blastomyces",
  "s_Angiostrongylus",
  "s_Gnathostoma",
  "s_Toxoplasma",
  "s_Schistosoma",
  "s_Rickettsia",
  "s_Orientia",
  "s_Nocardia",
  "s_Bartonella",
  "s_Corynebacterium",
  "s_Francisella",
  "s_Sporothrix",
  "s_Paracoccidioides",
  "s_Emergomyces",
  "s_Candida",
  "s_Baylisascaris",
  "s_Strongyloides",
  "s_Trypanosoma",
  "s_Plasmodium",
  "s_Fasciola",
  "s_Streptococcus",
  "s_Neisseria",
  "s_Haemophilus",
  "s_Listeria",
  "s_Staphylococcus",
  "s_Streptococcus",
  "s_Mycoplasma",
  "s_Kingella",
  "s_Campylobacter",
  "s_Salmonella",
  "s_Aeromonas",
  "s_Meningococcus",
  "s_Listeria",
  "s_Staphylococcus",
  "s_Streptococcus",
  "s_Mycoplasma",
  "s_Kingella",
  "s_Campylobacter",
  "s_Salmonella",
  "s_Aeromonas"
)

####################################################################################

# Extract the file name without extension for the directory name
input_file_name <- tools::file_path_sans_ext(basename(input_file_path))

# Append number_of_rows to the file name
name_wd <- paste0(input_file_name, "_", number_of_rows)

# Create the working directory path
working_dir <- file.path(getwd(), name_wd)

# Check if the directory exists, if not, create it
if (!dir.exists(working_dir)) {
  dir.create(working_dir)
  message("Created directory: ", working_dir)
} else {
  message("Directory already exists: ", working_dir)
}

# Set the working directory to the newly created folder
setwd(working_dir)
message("Working directory set to: ", getwd())


# Create a log file to store parameters
log_file <- file.path(working_dir, "parameters_log.txt")

# Write parameters to the log file
writeLines(c(
  paste("Number of Rows:", number_of_rows),
  "Weight Factors:",
  paste("  NormCladeReads =", weight_factors$NormCladeReads),
  paste("  Kmers =", weight_factors$Kmers),
  paste("  KmerDuplicity =", weight_factors$KmerDuplicity),
  paste("  KmerCoverage =", weight_factors$KmerCoverage),
  paste("SD Factor:", SD_factor),
  paste("CV Factor:", CV_factor),
  paste("Excluded Sources:", paste(exclude_sources, collapse = ", ")),
  paste("Filtered Sources:", paste(filter_sources, collapse = ", ")),
  paste("Input File Path:", input_file_path),
  paste("Working Directory:", working_dir)
), log_file)

message("Parameters log saved to: ", log_file)


# Reload the data
input_file <- read_csv(input_file_path, show_col_types = FALSE)
input_file_original <- input_file

# Exclude specified Source.Name values from the dataset
input_file <- input_file %>%
 filter(!Source.Name %in% exclude_sources)

# Check the structure of the data
str(input_file)

# Ensure the column names are correct
colnames(input_file) <- c("Source.Name", "Sample.type", "Pathogen", "Detected.pat", "Run.no", "NormCladeReads", "NormTaxonReads", "TaxRank", "TaxID", 
                          "Name", "Kmers", "KmerDuplicity", "KmerCoverage", "Depth", "TaxLineage")

# Replace NA values in 'Detected.pat' for NK, LPC, HPC with NA
input_file <- input_file %>%
  mutate(
    Detected.pat = ifelse(Sample.type %in% c("NK", "LPC", "HPC") & is.na(Detected.pat), "NA", Detected.pat)
  )

# Convert necessary columns to numeric
input_file <- input_file %>% 
  mutate(
    NormCladeReads = as.numeric(NormCladeReads),
    NormTaxonReads = as.numeric(NormTaxonReads),
    Kmers = as.numeric(Kmers),
    KmerDuplicity = as.numeric(KmerDuplicity),
    KmerCoverage = as.numeric(KmerCoverage)
  )

# Create a vector with the current values and their corresponding new names. Adjust this if original naming scheme of your file (column "Pathogen") doesn't match KrakenUniq naming scheme
pathogen_rename_map <- c(
  "Coxiella sp." = "s_Coxiella",
  "Enterobacter cloacae comp." = "s_Enterobacter cloacae",
  "Bartonella spp." = "s_Bartonella",
  "Streptococcus sp. (mitis group)" = "s_Streptococcus",
  "Staphylococcus spp." = "s_Staphylococcus",
  "Achromobacter spp." = "s_Achromobacter",
  "Babesia sp. - genom" = "s_Babesia",
  "Leptospira sp." = "s_Leptospira",
  "Human betaherpesvirus 5" = "s_Cytomegalovirus humanbeta5",
  "Anaplasma phagocytophilum" = "s_Anaplasma phagocytophilum",
  "Puumala virus" = "s_Orthohantavirus puumalaense",
  "Dengue virus" = "s_Orthoflavivirus denguei",
  "Tick-borne encephalitis virus" = "s_Orthoflavivirus encephalitidis",
  "Human gammaherpesvirus 4" = "s_Lymphocryptovirus humangamma4",
  "Parvovirus B19" = "s_Erythroparvovirus primate1",
  "Yellow fever virus" = "s_Yellow fever virus",
  "Chikungunya virus" = "s_Chikungunya virus",
  "Zika virus" = "s_Orthoflavivirus zikaense",
  "Haemophilus influenzae" = "s_Haemophilus influenzae",
  "Capnocytophaga canimorsus" = "s_Capnocytophaga canimorsus",
  "Staphylococcus aureus" = "s_Staphylococcus aureus",
  "Sphingomonas paucimobilis" = "s_Sphingomonas paucimobilis",
  "Neoehrlichia mikurensis" = "s_Neoehrlichia mikurensis",
  "Streptococcus pyogenes" = "s_Streptococcus pyogenes",
  "Escherichia coli" = "s_Escherichia coli",
  "Rothia mucilaginosa" = "s_Rothia mucilaginosa",
  "Tularemia" = "s_Francisella tularensis",
  "Dobrava virus" = "s_Orthohantavirus dobravaense",
  "Plasmodium falciparum" = "s_Plasmodium falciparum"
)

# Use tryCatch to handle potential errors and continue execution
input_file <- tryCatch({
  input_file %>%
    mutate(Pathogen = recode(Pathogen, !!!pathogen_rename_map))
}, error = function(e) {
  # If an error occurs, print a message and return the input_file without changes
  message("An error occurred while updating the 'Pathogen' column: ", e$message)
  return(input_file)
})


# Filter to keep only samples with TaxRank 'S' (optional, if needed)
input_file <- input_file %>%
  filter(TaxRank == 'S')

# Create nk_adjustments: Sum of NormCladeReads for NK samples within each Run.no and Name
nk_adjustments <- input_file %>% 
  filter(Sample.type == "NK") %>% 
  group_by(Run.no, Name) %>% 
  summarise(NK_Adjustment = sum(NormCladeReads, na.rm = TRUE), .groups = "drop")

# Define a function to calculate pathogen likelihood score for each unique organism within a run before removing NK reads
calculate_pathogen_score_prenk <- function(data, weight_factors) {
  data <- data %>% 
    group_by(Run.no, Name) %>%  # Group by Run.no and Name
    mutate(
      # Use the weight factors for NormCladeReads, Kmers, KmerDuplicity, and KmerCoverage
      
      PathogenScorePreNK = (NormCladeReads^weight_factors$NormCladeReads * 
                         log10(Kmers + 1)^weight_factors$Kmers * 
                         (1 / (KmerDuplicity^weight_factors$KmerDuplicity))  * 
                         KmerCoverage^weight_factors$KmerCoverage)
    ) %>% 
    ungroup() %>%
    mutate(
      # Round PathogenScore to 6 decimal places
      PathogenScorePreNK = round(PathogenScorePreNK, 6),
      # Add a new column with the logarithm of PathogenScore (log10 of PathogenScore)
      LogPathogenScorePreNK = log10(PathogenScorePreNK + 1)  # Adding 1 to avoid log(0)
    )
  return(data)
}

# Apply the pathogen score calculation with the weight factors this is used to get the background score according to which the low positive pathogens are filtered out
input_file <- calculate_pathogen_score_prenk(input_file, weight_factors)

# Get all unique microorganisms for each run and convert to long format
unique_microorganisms_long <- input_file %>%
  group_by(Run.no) %>%
  summarise(
    Unique_Microorganisms = list(unique(Name)),
    .groups = "drop"
  ) %>%
  unnest(cols = Unique_Microorganisms)  # Convert the list of microorganisms into long format

# Rename column to make it consistent for joining
unique_microorganisms_long <- unique_microorganisms_long %>%
  rename(Microorganism = Unique_Microorganisms)

# Cross-join microorganisms with Source.Name for each Run.no
all_combinations <- input_file %>%
  distinct(Source.Name, Run.no) %>% # Get all unique Source.Name and Run.no combinations
  inner_join(
    unique_microorganisms_long,
    by = "Run.no",
    relationship = "many-to-many" # Explicitly allow many-to-many relationships
  )

# Join with input_file to get LogPathogenScorePreNK values
complete_data <- all_combinations %>%
  left_join(
    input_file %>% select(Source.Name, Run.no, Name, LogPathogenScorePreNK), # Select relevant columns
    by = c("Source.Name", "Run.no", "Microorganism" = "Name") # Match on Source.Name, Run.no, and Microorganism
  ) %>%
  mutate(LogPathogenScorePreNK = replace_na(LogPathogenScorePreNK, 0)) # Replace missing values with 0

# Calculate generalized metrics for each microorganism on a per-run basis
generalized_metrics_per_run <- complete_data %>%
  group_by(Run.no, Microorganism) %>%
  summarise(
    Mean_PatScorePreNK = mean(LogPathogenScorePreNK, na.rm = TRUE),   # Mean abundance
    SD_PatScorePreNK = sd(LogPathogenScorePreNK, na.rm = TRUE),       # Standard deviation
    CV_PatScorePreNK = (SD_PatScorePreNK / Mean_PatScorePreNK) * 100, # Coefficient of variation
    Presence = sum(LogPathogenScorePreNK > 0),                        # Count of samples where the microorganism is detected
    Total_Samples = n(),                                              # Total number of samples in the run
    Proportion_Presence = Presence / Total_Samples                    # Proportion of samples with microorganism
  ) %>%
  mutate(
    # Z-scores for detecting outliers in mean and variability within each run
    Z_Score_Mean = (Mean_PatScorePreNK - mean(Mean_PatScorePreNK, na.rm = TRUE)) / sd(Mean_PatScorePreNK, na.rm = TRUE),
    Z_Score_SD = (SD_PatScorePreNK - mean(SD_PatScorePreNK, na.rm = TRUE)) / sd(SD_PatScorePreNK, na.rm = TRUE),
    Z_Score_CV = (CV_PatScorePreNK - mean(CV_PatScorePreNK, na.rm = TRUE)) / sd(CV_PatScorePreNK, na.rm = TRUE),
    # Additional flags for extreme behavior
    Is_High_Mean = Z_Score_Mean > 2,                   # Flag for abnormally high mean
    Is_High_SD = Z_Score_SD > 2,                       # Flag for abnormally high variability
    Is_Low_Proportion = Proportion_Presence < 0.1      # Flag for microorganisms in very few samples
  ) %>%
  arrange(Run.no, desc(Mean_PatScorePreNK))            # Sort for inspection by Run.no and mean

#copy input_file for fn and tp statistic calculation
input_file_fn <- input_file
input_file_tp <- input_file

# Calculate the total number of unique Source.names for each Run.no (without grouping by Name)
total_samples_per_run <- input_file %>%
  group_by(Run.no) %>%
  summarise(TotalSamples = n_distinct(Source.Name), .groups = "drop")

# Identify all microorganisms of dataset
all_microorganisms <- input_file %>%
  group_by(Name) %>%
  select(Name) %>%
  distinct()

# Identify the presence of each microorganism in each sample
microorganism_presence <- input_file %>%
  group_by(Run.no, Name) %>%
  summarise(
    PresenceCount = sum(LogPathogenScorePreNK > 0),  # Count detected samples
    BleedthroughThreshold = mean(LogPathogenScorePreNK, na.rm = TRUE) + 
      sd(LogPathogenScorePreNK, na.rm = TRUE),  # Calculate threshold directly
    .groups = "drop"
  )

# Join the total samples back to the microorganism presence data
microorganism_presence <- microorganism_presence %>%
  left_join(total_samples_per_run, by = "Run.no") %>%
  mutate(
    ProportionPresence = PresenceCount / TotalSamples  # Calculate the proportion of presence
  )

# Flag microorganisms as potential contaminants
potential_contaminants <- microorganism_presence %>%
  filter(ProportionPresence >= contaminant_threshold) %>%
  mutate(PotentialContaminant = "TRUE")

# Filter the data for the specified sample
potential_contaminants_filtered <- potential_contaminants %>%
  filter(ProportionPresence > 0.9)

# Join contaminant info to the input file by Run.no and Name
input_file <- input_file %>%
  left_join(potential_contaminants %>% select(Run.no, Name, PotentialContaminant, BleedthroughThreshold), 
            by = c("Run.no", "Name"))

# Create the AboveBleedthroughThreshold column to ensure it's available for manipulation
input_file <- input_file %>%
  mutate(AboveBleedthroughThreshold = NA_character_)  # Initialize with NA

# Flag contaminants as above or below the bleedthrough threshold
input_file <- input_file %>%
  mutate(AboveBleedthroughThreshold = case_when(
    PotentialContaminant == "TRUE" & LogPathogenScorePreNK < BleedthroughThreshold ~ "FALSE",
    PotentialContaminant == "TRUE" & LogPathogenScorePreNK >= BleedthroughThreshold ~ "TRUE",
    TRUE ~ AboveBleedthroughThreshold  # Keep existing value (NA) for non-contaminants
  ))

# Exact match for viruses (TRUE if match, FALSE if no match)
all_microorganisms$virus_match <- all_microorganisms$Name %in% marked_pathogen_viruses

# Extract genus (first word) for both viruses and bacteria
genus_viruses <- sub("^(s_[a-zA-Z]+).*", "\\1", all_microorganisms$Name)

# Genus-based match for viruses (check if genus is in marked_pathogen_viruses)
# Extract the genus part from the virus list as well
genus_viruses_list <- sub("^(s_[a-zA-Z]+).*", "\\1", marked_pathogen_viruses)
all_microorganisms$virus_genus_match <- genus_viruses %in% genus_viruses_list

# Check for genus-based bacteria matches
genus_bacteria <- sub("^(s_[a-zA-Z]+).*", "\\1", all_microorganisms$Name)

# Check if the genus is in the list of marked pathogens for bacteria
all_microorganisms$bacteria_match <- genus_bacteria %in% marked_pathogen_bacteria

# New column that is TRUE if any of the three match columns is TRUE
all_microorganisms$target_pathogen <- with(all_microorganisms, 
                                       virus_match | virus_genus_match | bacteria_match)

# Remove the original matching columns
all_microorganisms <- all_microorganisms[, !colnames(all_microorganisms) %in% c("virus_match", "virus_genus_match", "bacteria_match")]

# Final dataframe with matches
print(all_microorganisms)

# Assuming 'all_microorganisms' already has the 'target_pathogen' column from earlier code
# If you want to work with a filtered dataset, for example, you can filter the organisms with TRUE in 'target_pathogen'

# Flag microorganisms based on target_pathogen (TRUE/FALSE organisms)
flagged_microorganisms <- all_microorganisms %>%
  filter(target_pathogen == TRUE)  # You can change this condition to flag based on other criteria

# Join back to the input file to label microorganisms as contaminants
input_file <- input_file %>%
  left_join(flagged_microorganisms %>% select(Name, target_pathogen), 
            by = c("Name"))

# Generate separate files for each Run.no
# Create a list of unique Run.nos
run_numbers <- unique(input_file$Run.no)

# Save a separate file for each Run.no
for(run in run_numbers) {
  # Filter data for the current Run.no
  run_data <- potential_contaminants_filtered %>% filter(Run.no == run)
  
  # Create a file name based on the Run.no
  file_name <- paste0("Run_", run, "_Potential_Contaminants.csv")
  
  # Write the data to a CSV file
  write.csv(run_data, file = file_name, row.names = FALSE)
}

# View the result
head(input_file)

# Subtract NK reads from SAMPLEs within each run and ensure no negative values
input_file <- input_file %>% 
  left_join(nk_adjustments, by = c("Run.no", "Name")) %>% 
  mutate(
    NormCladeReads = ifelse(Sample.type == "SAMPLE", 
                            pmax(NormCladeReads - coalesce(NK_Adjustment, 0), 0), 
                            NormCladeReads)
  ) %>% 
  select(-NK_Adjustment)

# Filter out rows where NormCladeReads == 0
input_file <- input_file %>%
  filter(NormCladeReads != 0)

# Remove rows with Sample.type as "NK", "LPC", "HPC", or "PK"
input_file <- input_file %>%
  filter(!Sample.type %in% c("NK", "LPC", "HPC", "PK"))

# Define the prefixes to standardize
prefixes <- c(
  "s_Coxiella",
  "s_Enterobacter cloacae",
  "s_Bartonella",
  "s_Streptococcus",
  "s_Staphylococcus",
  "s_Achromobacter",
  "s_Babesia",
  "s_Leptospira"
)

# Standardize 'Name' column based on the prefixes
input_file <- input_file %>%
  mutate(
    SimplifiedName = sapply(Name, function(n) {
      match <- prefixes[sapply(prefixes, function(p) grepl(paste0("^", p), n))]
      if (length(match) > 0) match else n
    })
  )

# Matching function for specified pathogens (partial match) and others (exact match)
match_pathogen <- function(pathogen, name) {
  
  # Define the list of pathogens that require partial matching
  partial_match_pathogens <- c("s_Coxiella", "s_Enterobacter cloacae", "s_Bartonella", "s_Streptococcus", 
                               "s_Staphylococcus", "s_Achromobacter", "s_Babesia", "s_Leptospira")
  
  # Check if pathogen is in the partial match list
  if (pathogen %in% partial_match_pathogens) {
    return(grepl(pathogen, name))  # Partial match for specified pathogens
  } else {
    return(pathogen == name)  # Exact match for all other pathogens
  }
}

# Apply the matching function to check if Pathogen and Name match
input_file <- input_file %>%
  mutate(
    PathogenMatch = mapply(match_pathogen, Pathogen, Name)
  )


# Define a function to calculate pathogen likelihood score for each unique organism within a run
calculate_pathogen_score <- function(data, weight_factors) {
  data <- data %>% 
    group_by(Run.no, Name) %>%  # Group by Run.no and Name
    mutate(
      # Use the weight factors for NormCladeReads, Kmers, KmerDuplicity, and KmerCoverage
      
      PathogenScore = (NormCladeReads^weight_factors$NormCladeReads * 
                         log10(Kmers + 1)^weight_factors$Kmers * 
                         (1 / (KmerDuplicity^weight_factors$KmerDuplicity))  * 
                         KmerCoverage^weight_factors$KmerCoverage)
    ) %>% 
    ungroup() %>%
    mutate(
      # Round PathogenScore to 6 decimal places
      PathogenScore = round(PathogenScore, 6),
      # Add a new column with the logarithm of PathogenScore (log10 of PathogenScore)
      LogPathogenScore = log10(PathogenScore + 1)  # Adding 1 to avoid log(0)
    )
  return(data)
}

# Apply the pathogen score calculation with the weight factors this is used to get the background score according to which the low positive pathogens are filtered out
input_file <- calculate_pathogen_score(input_file, weight_factors)

# Calculate the background threshold per microorganism per run
generalized_metrics_per_run <- generalized_metrics_per_run %>%
  mutate(
    Threshold_Mean_SD = Mean_PatScorePreNK + SD_factor * SD_PatScorePreNK,  # Mean + 0.6×SD
    Threshold_CV = Mean_PatScorePreNK * (1 + CV_PatScorePreNK / CV_factor)  # Mean weighted by CV
  )

# Join the thresholds back to the input file
input_file <- input_file %>%
  left_join(
    generalized_metrics_per_run %>% 
      select(Run.no, Microorganism, Threshold_Mean_SD, Threshold_CV), 
    by = c("Run.no", "Name" = "Microorganism")
  )

# Categorize microorganisms as "likely" or "unlikely"
input_file <- input_file %>%
  mutate(
    Likelihood = case_when(
      LogPathogenScore > Threshold_Mean_SD ~ "LIKELY",      # Use Mean + 0.6×SD threshold
      #LogPathogenScore > Threshold_CV ~ "LIKELY",           # Optionally use CV-based threshold
      TRUE ~ "UNLIKELY"
    )
  )

# Output for top 15 pathogens based on log-transformed pathogen score
top_filtered_pathogens <- input_file %>%
  filter(Sample.type == "SAMPLE", Likelihood == "LIKELY") %>%  # Ignore "UNLIKELY" pathogens in this step
  group_by(Source.Name) %>%  # Group by sample (Source.Name)
  arrange(Source.Name, desc(LogPathogenScore)) %>%  # Sort within each sample by LogPathogenScore
  slice_head(n = number_of_rows) %>%  # Get the top 15 pathogens for each sample
  ungroup()  # Ungroup after selection

# Create a table for top 15 pathogens, showing relevant columns
top_filtered_pathogens_table <- top_filtered_pathogens %>%
  select(Run.no, Source.Name, Pathogen, Name, LogPathogenScore) %>%  # Include relevant columns
  arrange(Source.Name, desc(LogPathogenScore), Name)  # Sort by LogPathogenScore

# Create a table for top 15 pathogens, showing relevant columns
top_filtered_pathogens_table_contaminants <- top_filtered_pathogens %>%
  select(Run.no, Source.Name, Pathogen, Name, PotentialContaminant, AboveBleedthroughThreshold, LogPathogenScore) %>%  # Include relevant columns
  arrange(Source.Name, desc(LogPathogenScore), Name)  # Sort by LogPathogenScore

# Create a table for top 15 pathogens, showing relevant columns
top_filtered_pathogens_table_contaminants_target <- top_filtered_pathogens %>%
  select(Run.no, Source.Name, Pathogen, Name, PotentialContaminant, AboveBleedthroughThreshold, target_pathogen, LogPathogenScore) %>%  # Include relevant columns
  arrange(Source.Name, desc(LogPathogenScore), Name)  # Sort by LogPathogenScore

# Create a list of unique microorganism names for easier ncbi searching
unique_microorganisms_from_top <- top_filtered_pathogens %>%
  select(Name) %>%
  mutate(Name = str_remove(Name, "^s_")) %>%
  distinct() %>%
  arrange(Name)
write_csv(unique_microorganisms_from_top, "unique_microorganisms_from_top.csv")

# Save the top 15 pathogens to a CSV file
top_output_file <- "top_filtered_pathogens_by_log_score.csv"
write_csv(top_filtered_pathogens_table, top_output_file)

cat("Top 15 pathogens by log-transformed pathogen score saved to:", top_output_file, "\n")

# Save the top 15 pathogens with contaminants to a CSV file
top_output_file_contaminants <- "top_filtered_pathogens_by_log_score_contaminants.csv"
write_csv(top_filtered_pathogens_table_contaminants, top_output_file_contaminants)

cat("Top 15 pathogens by log-transformed pathogen score saved to:", top_output_file_contaminants, "\n")

# Save the top 15 pathogens with contaminants and target pathogens to a CSV file
top_output_file_contaminants_target <- "top_filtered_pathogens_by_log_score_contaminants_target.csv"
write_csv(top_filtered_pathogens_table_contaminants_target, top_output_file_contaminants_target)

cat("Top 15 pathogens by log-transformed pathogen score saved to:", top_output_file_contaminants_target, "\n")

#################################
# Graphs to help interpret data #
# and find suitable cut-off     #
#################################

# Rank samples by pathogen score
ranked_data <- top_filtered_pathogens_table %>%
  group_by(Run.no, Source.Name) %>%  # Group by 'Run.no' and 'Source.Name'
  mutate(Rank = row_number()) %>%    # Assign ranks within each group
  ungroup()  # Ungroup after ranking

# Calculate average LogPathogenScore for each rank
average_scores <- ranked_data %>%
  group_by(Rank) %>%
  summarise(AverageLogPathogenScore = mean(LogPathogenScore, na.rm = TRUE)) %>%
  ungroup()

# Create the plot and store it in a variable
ggplot(average_scores, aes(x = Rank, y = AverageLogPathogenScore)) +
  geom_line(fill = "skyblue", color = "black") +  # Bar plot
  theme_bw() +
  labs(
    title = "Average LogPathogenScore by Rank",
    x = "Rank",
    y = "Average LogPathogenScore"
  )

# Save the plot to a file (PNG format, but you can also use PDF or others)
ggsave("average_plot.png", width = 8, height = 6)

# Calculate the overall average LogPathogenScore
overall_average <- mean(ranked_data$LogPathogenScore, na.rm = TRUE)

# Calculate average LogPathogenScore for each rank and the deviation from the overall average
deviation_scores <- ranked_data %>%
  group_by(Rank) %>%
  summarise(AverageLogPathogenScore = mean(LogPathogenScore, na.rm = TRUE)) %>%
  mutate(DeviationFromAverage = pmax(AverageLogPathogenScore - overall_average)) %>%
  ungroup()

# Create the plot
ggplot(deviation_scores, aes(x = Rank, y = DeviationFromAverage)) +
  geom_line(fill = "skyblue", color = "black") +
  theme_bw() +
  labs(
    title = "Deviation from Overall Average LogPathogenScore by Rank",
    x = "Rank",
    y = "Deviation from Average LogPathogenScore"
  )

# Save the plot to a file (PNG format, but you can also use PDF or others)
ggsave("deviation_from_average_plot.png", width = 8, height = 6)

ranked_data_run <- top_filtered_pathogens_table %>%
  group_by(Run.no, Source.Name) %>%  # Group by Run.no and Source.Name
  mutate(Rank = row_number()) %>%    # Assign ranks within each group
  ungroup()

average_scores_per_run <- ranked_data_run %>%
  group_by(Run.no, Rank) %>%
  summarise(AverageLogPathogenScore = mean(LogPathogenScore, na.rm = TRUE)) %>%
  ungroup()

ggplot(average_scores_per_run, aes(x = Rank, y = AverageLogPathogenScore)) +
  geom_line(fill = "skyblue", color = "black") +
  theme_bw() +
  labs(
    title = "Average LogPathogenScore by Rank (Faceted by Run.no)",
    x = "Rank",
    y = "Average LogPathogenScore"
  ) +
  facet_wrap(~Run.no, scales = "free_x")  # Creates separate plots per Run.no


# Save the plot
ggsave("average_plot_per_run.png", width = 8, height = 6)

# Calculate the overall average LogPathogenScore for each Run.no
overall_avg_per_run <- ranked_data %>%
  group_by(Run.no) %>%
  summarise(OverallAverage = mean(LogPathogenScore, na.rm = TRUE), .groups = "drop")

# Calculate the average score per rank within each Run.no
deviation_scores_per_run <- ranked_data %>%
  group_by(Run.no, Rank) %>%
  summarise(AverageLogPathogenScore = mean(LogPathogenScore, na.rm = TRUE), .groups = "drop") %>%
  left_join(overall_avg_per_run, by = "Run.no") %>%  # Merge in overall average
  mutate(DeviationFromAverage = AverageLogPathogenScore - OverallAverage) %>%  # Compute deviation
  ungroup()

ggplot(deviation_scores_per_run, aes(x = Rank, y = DeviationFromAverage, color = as.factor(Run.no))) +
  geom_line(fill = "skyblue", color = "black") +
  theme_bw() +
  labs(
    title = "Deviation from Overall Average LogPathogenScore by Rank per Run.no",
    x = "Rank",
    y = "Deviation from Average LogPathogenScore",
    color = "Run.no"
  ) +
  facet_wrap(~Run.no)  # Creates separate plots per Run.no


# Save the plot
ggsave("deviation_from_average_plot_per_run.png", width = 8, height = 6)

# Find the first rank where deviation is less than 0
first_negative_deviation_rank <- deviation_scores %>%
  filter(DeviationFromAverage < 0) %>%
  slice(1) %>%  # Take the first row where DeviationFromAverage < 0
  pull(Rank)    # Extract the Rank column

# Print the rank
print(paste("The first rank where deviation is less than 0 is:", first_negative_deviation_rank))

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
# 
# # Specify the sample you want to plot
# specified_sample <- "3597"  # Replace with the actual sample name
# 
# # Filter the data for the specified sample
# filtered_sample_data <- ranked_data %>%
#   filter(Source.Name == specified_sample)
# 
# # # Create a bar plot with Rank on the x-axis and LogPathogenScore on the y-axis
# ggplot(filtered_sample_data, aes(x = factor(Rank), y = LogPathogenScore)) +  # factor() to ensure rank is treated as discrete
#   geom_col(fill = "skyblue", color = "black") +  # Create bars
#   geom_text(aes(label = round(LogPathogenScore, 3)),
#             vjust = 0.5, hjust = -0.5, color = "black", angle = 90) +  # Adjust text position
#   theme_minimal() +
#   labs(
#     title = paste("LogPathogenScore by Rank for", specified_sample),
#     x = "Rank",
#     y = "LogPathogenScore"
#   )
# 
# # # Create a bar plot for all samples, faceted by 'Source.Name'
# ggplot(ranked_data, aes(x = factor(Rank), y = LogPathogenScore)) +  # factor() to ensure rank is treated as discrete
#   geom_col(fill = "skyblue", color = "black") +  # Create bars
#   geom_text(aes(label = round(LogPathogenScore, 3)),
#             vjust = 0.5, hjust = -0.5, color = "black", angle = 90) +  # Adjust text position
#   theme_minimal() +
#   labs(
#     title = "LogPathogenScore by Rank for All Samples",
#     x = "Rank",
#     y = "LogPathogenScore"
#   ) +
#   facet_wrap(~ Source.Name, scales = "free_x")

###################################################################################
#  _____ _   _ ____     ___  _____   _   _ _   _ _  ___   _  _____        ___   _ #
# | ____| \ | |  _ \   / _ \|  ___| | | | | \ | | |/ / \ | |/ _ \ \      / / \ | |#
# |  _| |  \| | | | | | | | | |_    | | | |  \| | ' /|  \| | | | \ \ /\ / /|  \| |#
# | |___| |\  | |_| | | |_| |  _|   | |_| | |\  | . \| |\  | |_| |\ V  V / | |\  |#
# |_____|_| \_|____/   \___/|_|      \___/|_| \_|_|\_\_| \_|\___/  \_/\_/  |_| \_|#
###################################################################################

df_data_check <- input_file

print(class(df_data_check))

df_data_check <- as.data.frame(df_data_check)

# Stop script if 'Pathogen' column does not exist
if (!"Pathogen" %in% names(df_data_check)) {
  stop("The 'Pathogen' column does not exist in the dataset.")
}

# Stop script if 'Pathogen' column is empty or contains only NAs
if (nrow(df_data_check) == 0 || all(is.na(df_data_check$Pathogen) | df_data_check$Pathogen == "")) {
  stop("The 'Pathogen' column is empty or contains only missing values.")
}

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################

# Flag rows where Pathogen and SimplifiedName are identical
filtered_data <- input_file %>%
  mutate(PathogenMatch = ifelse(Pathogen == SimplifiedName, "Match", "No Match"))

# Refine detection logic to ensure proper detection status
filtered_data <- filtered_data %>%
  mutate(
    Detection = case_when(
      PathogenMatch == "Match" & Detected.pat == "YES" ~ "TRUE",  # Matched and detected
      PathogenMatch == "No Match" & Detected.pat == "NO" ~ "FALSE", # No match but not detected
      TRUE ~ "FALSE"  # For all other cases, mark as FALSE
    )
  )

# Separate the data into Match and No Match datasets
matched_data <- filtered_data %>%
  filter(PathogenMatch == "Match")  # Only rows with "Match"

# Identify Source.Name samples that have at least one match entry
samples_with_matches <- matched_data %>%
  select(Source.Name) %>%
  distinct()

# Sort the matched data by LogPathogenScore in descending order and keep the top match for each sample
sorted_matched_data <- matched_data %>%
  arrange(Source.Name, Run.no, desc(LogPathogenScore)) %>%  # Sort by LogPathogenScore within each sample group
  group_by(Source.Name, Run.no) %>%  # Group by Source.Name and Run.no to apply sorting per sample
  slice_head(n = 1) %>%  # Keep only the top match (highest LogPathogenScore) for each sample group
  ungroup()  # Remove the grouping to return to the regular data frame structure

# Exclude these samples from unmatched_data and sort by sample first, then pathogen score
unmatched_data <- filtered_data %>%
  filter(!(Source.Name %in% samples_with_matches$Source.Name)) %>%  # Exclude samples with any "Match" entries
  group_by(Source.Name) %>%  # Group by Source.Name
  arrange(Source.Name, desc(LogPathogenScore)) %>%  # Sort within each sample by descending LogPathogenScore
  ungroup()  # Remove grouping to return to the regular data frame structure

# For unmatched data, filter out samples with at least one match and retain only samples that had no match
unmatched_data_cleaned <- unmatched_data %>%
  filter(PathogenMatch == "No Match") %>%  # Only rows with "No Match"
  group_by(Source.Name) %>%
  filter(all(PathogenMatch == "No Match")) %>%  # Keep only those Source.Name with all "No Match"
  ungroup() %>%
  group_by(Source.Name) %>%
  slice_head(n = number_of_rows) %>%  # Keep the top 15 rows for each sample
  mutate(
    Name = NA  # Set Name to NA for No Match samples
    # LogPathogenScore = NA  # Uncomment if you also want to set LogPathogenScore to NA
  ) %>%
  ungroup()

# Combine the matched and unmatched datasets without repetition of each sample and keep the Pathogen column
final_combined_data <- bind_rows(
  sorted_matched_data,
  unmatched_data_cleaned %>%
    distinct(Source.Name, .keep_all = TRUE) %>%  # Remove duplicate rows based on Source.Name
    mutate(across(-c(Source.Name, Sample.type, Detected.pat, Run.no, Pathogen, PathogenMatch, Detection), ~NA))  # Set all columns except the specified ones to NA
) %>%
  arrange(Source.Name, Run.no)  # Sort the final combined dataset by sample and Run.no

# Select only necessary columns (optional)
final_combined_data <- final_combined_data %>%
  mutate(Name = ifelse(is.na(Name), Pathogen, Name)) %>%
  select(Source.Name, Pathogen, Name, LogPathogenScore, PathogenMatch, Detected.pat, Detection)

# Print or save the final combined data as required
# For example, saving it to a CSV file
write_csv(final_combined_data, "final_combined_match_status.csv")

cat("Final combined match status saved to: final_combined_match_status.csv\n")

# Summarize top 15 pathogens and LogPathogenScore into separate columns by Source.Name
top_grouped <- top_filtered_pathogens_table %>%
  group_by(Source.Name) %>%
  slice_head(n = number_of_rows) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  # Create columns for both LogPathogenScore and Name
  pivot_wider(
    names_from = row, 
    values_from = c(Name, LogPathogenScore), 
    names_prefix = "Pathogen_"
  )

# View the resulting top_grouped dataframe
print(top_grouped)

# Join final_combined_data with top_grouped by Source.Name
joined_data <- final_combined_data %>%
  left_join(top_grouped, by = "Source.Name")

# Clean Name and pathogen columns (remove whitespace and make lowercase)
joined_data <- joined_data %>%
  mutate(
    Name = str_trim(tolower(Name)),  # Remove extra spaces and make lowercase
    across(starts_with("Name_Pathogen_"), ~str_trim(tolower(.)))  # Apply to pathogen columns
  )

# Shorten name to insure match on a genus level
joined_data <- joined_data %>%
  mutate(
    Name = str_trim(tolower(Name)),  # Remove extra spaces and make lowercase
    Name.Short = str_split_fixed(Name, " ", 2)[, 1],  # Extract the first word before the space
    across(starts_with("Name_Pathogen_"), ~str_trim(tolower(.)))  # Apply to pathogen columns
  )

# Iterate across each Name_Pathogen column to create a new Match column for each
joined_data <- joined_data %>%
  mutate(across(
    starts_with("Name_Pathogen_"),
    ~ str_split_fixed(., " ", 2)[, 1] == Name.Short,  # Split and compare the first word
    .names = "Match_{.col}"  # Create a new column for each pathogen match
  ))

# View the resulting data to confirm new columns were added
print(joined_data)

# Write to CSV
write_csv(joined_data, "matched_pathogens_individual_matches.csv")

# Verify if the file was created successfully
cat("CSV file created: matched_pathogens_individual_matches.csv\n")

# Create a new column 'OverallMatch' based on the Match columns
joined_data <- joined_data %>%
  mutate(
    OverallMatch = case_when(
      # If any column is TRUE, set OverallMatch to TRUE
      rowSums(select(., starts_with("Match_Name_Pathogen_")) == TRUE, na.rm = TRUE) > 0 ~ TRUE,
      
      # If no TRUE values but at least one FALSE, set OverallMatch to FALSE
      rowSums(select(., starts_with("Match_Name_Pathogen_")) == FALSE, na.rm = TRUE) > 0 ~ FALSE,
      
      # If all columns are NA, set OverallMatch to FALSE
      TRUE ~ FALSE
    )
  )

# View the updated data
print(joined_data)

# Write to CSV
write_csv(joined_data, "matched_pathogens_with_overall_match.csv")

# Verify if the file was created successfully
cat("CSV file created: matched_pathogens_with_overall_match.csv\n")

# Rename 'YES' to 'Positive' and 'NO' to 'Negative'
joined_data <- joined_data %>%
  mutate(
    Detected.pat = case_when(
      Detected.pat == "YES" ~ "Positive",
      Detected.pat == "NO" ~ "Negative",
      TRUE ~ Detected.pat  # Keep other values unchanged if any
    )
  )

# Create the 'Result' column based on 'Detected.pat' and 'OverallMatch'
joined_data <- joined_data %>%
  mutate(
    Result = case_when(
      Detected.pat == "Positive" & OverallMatch == TRUE ~ "True Positive",
      Detected.pat == "Positive" & OverallMatch == FALSE ~ "False Negative",
      Detected.pat == "Negative" & OverallMatch == TRUE ~ "False Positive",
      Detected.pat == "Negative" & OverallMatch == FALSE ~ "True Negative",
      TRUE ~ NA_character_
    )
  )

# Exclude specified Source.Name values from the dataset
joined_data <- joined_data %>%
  filter(Source.Name %in% filter_sources)

# Write the updated data with 'Result' column to CSV
write_csv(joined_data, "matched_pathogens_with_result.csv")

# Calculate counts for each result type
tp_count <- sum(joined_data$Result == "True Positive", na.rm = TRUE)
fn_count <- sum(joined_data$Result == "False Negative", na.rm = TRUE)
fp_count <- sum(joined_data$Result == "False Positive", na.rm = TRUE)
tn_count <- sum(joined_data$Result == "True Negative", na.rm = TRUE)

# Calculate total count and rates, rounding to 3 decimal places
total_count <- nrow(joined_data)

tp_rate <- round(tp_count / total_count, 3)
fn_rate <- round(fn_count / total_count, 3)
fp_rate <- round(fp_count / total_count, 3)
tn_rate <- round(tn_count / total_count, 3)

# Output the rates to the console
cat("True Positive Rate (TP):", tp_rate, tp_count, "\n")
cat("False Negative Rate (FN):", fn_rate, fn_count, "\n")
cat("False Positive Rate (FP):", fp_rate, fp_count, "\n")
cat("True Negative Rate (TN):", tn_rate, tn_count, "\n")

# Create the data frame for the agreement plot
result_df <- data.frame(
  Detected.pat = c("Positive", "Negative", "Negative", "Positive"),
  OverallMatch = c("True", "True", "False", "False"),
  Freq = c(tp_count, fp_count, tn_count, fn_count),
  rate = c(tp_rate, fp_rate, tn_rate, fn_rate)
)

# View the resulting data frame
print(result_df)

# Extract counts
tp_count <- result_df$Freq[result_df$Detected.pat == "Positive" & result_df$OverallMatch == "True"]
fp_count <- result_df$Freq[result_df$Detected.pat == "Negative" & result_df$OverallMatch == "True"]
tn_count <- result_df$Freq[result_df$Detected.pat == "Negative" & result_df$OverallMatch == "False"]
fn_count <- result_df$Freq[result_df$Detected.pat == "Positive" & result_df$OverallMatch == "False"]

# Calculate total number of observations (N)
N <- sum(result_df$Freq)

# Calculate observed agreement P_o
P_o <- (tp_count + tn_count) / N

# Calculate expected agreement P_e
p_positive <- (tp_count + fp_count) / N  # Proportion of Positive cases
p_negative <- (tn_count + fn_count) / N  # Proportion of Negative cases

P_e <- (p_positive * (tp_count + fn_count) / N) + (p_negative * (fp_count + tn_count) / N)

# Calculate Cohen's Kappa
kappa <- (P_o - P_e) / (1 - P_e)

interpretation <- ""
if (kappa >= 0.81) {
  interpretation <- "Almost perfect agreement"
} else if (kappa >= 0.61) {
  interpretation <- "Substantial agreement"
} else if (kappa >= 0.41) {
  interpretation <- "Moderate agreement"
} else if (kappa >= 0.21) {
  interpretation <- "Fair agreement"
} else if (kappa >= 0.01) {
  interpretation <- "Slight agreement"
} else if (kappa == 0) {
  interpretation <- "No agreement beyond chance"
} else {
  interpretation <- "Less agreement than expected by chance"
}

# Create the full interpretation scale
interpretation_scale <- "
Interpretation Scale:
1: Perfect agreement
0.81 to 1.00: Almost perfect agreement
0.61 to 0.80: Substantial agreement
0.41 to 0.60: Moderate agreement
0.21 to 0.40: Fair agreement
0.01 to 0.20: Slight agreement
0: No agreement beyond chance
Negative values: Less agreement than expected by chance
"
# Create the text content to write to a file
output_text <- paste("Cohen's Kappa: ", kappa, "\nInterpretation: ", interpretation, "\n\n", interpretation_scale, sep = "")

# Specify the filename
filename <- "Cohens_Kappa_Interpretation.txt"

# Write the output to a text file in the working directory
write(output_text, file = filename)

# Notify the user that the file has been saved
cat("The Cohen's Kappa value, interpretation, and scale have been written to ", filename, "\n", sep = "")

# Print the Cohen's Kappa value
print(paste("Cohen's Kappa: ", kappa))

# Create the plot with counts and rates
result_df$OverallMatch <- factor(result_df$OverallMatch, levels = c("True", "False"))
agreement_plot <- ggplot(result_df, aes(x = OverallMatch, y = Detected.pat, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = paste("Number:", Freq, "\nRatio:", rate)), color = "black", size = 5, family = "Times New Roman") +
  scale_fill_gradient(low = "white", high = "lightgreen") +
  theme_bw() +
  labs(title = "Agreement plot of Detection vs Overall Match",
       x = "ClinSeq detection",
       y = "Manual review detection",
       fill = "Number") +
  theme(text = element_text(family = "Times New Roman"),
        axis.title.x = element_text(size = 14),           # Set x-axis title size
        axis.title.y = element_text(size = 14),           # Set y-axis title size
        axis.text.x = element_text(size = 12),            # Set x-axis labels size
        axis.text.y = element_text(size = 12)             # Set y-axis labels size
        )

# Save the plot with the filename based on the weight factors
file_name <- paste0("Agreement plot", 
                    weight_factors$NormCladeReads, 
                    weight_factors$Kmers, 
                    weight_factors$KmerDuplicity, 
                    weight_factors$KmerCoverage, 
                    ".png")

# Print the plot and save it
print(agreement_plot)
ggsave(file_name, plot = agreement_plot, width = 6, height = 6)

# Filter for False Negative and False Positive samples and return their source.name
false_negative <- joined_data %>%
  filter(Result %in% c("False Negative", "False Positive")) %>%
  select(Source.Name, Result)

# Display the source.name and result for False Negative and False Positive samples
print(n = 50, false_negative)

# Optionally save the False Negative and False Positive source names to a CSV file
write.csv(false_negative, "false_negative.csv", row.names = FALSE)


# Standardize 'Name' column based on the prefixes
input_file_fn <- input_file_fn %>%
  mutate(
    SimplifiedName = sapply(Name, function(n) {
      match <- prefixes[sapply(prefixes, function(p) grepl(paste0("^", p), n))]
      if (length(match) > 0) match else n
    })
  )


# Matching function for specified pathogens (partial match) and others (exact match)
match_pathogen_fn <- function(pathogen, name) {
  
  # Define the list of pathogens that require partial matching
  partial_match_pathogens <- c("s_Coxiella", "s_Enterobacter cloacae", "s_Bartonella", "s_Streptococcus", 
                               "s_Staphylococcus", "s_Achromobacter", "s_Babesia", "s_Leptospira")
  
  # Check if pathogen is in the partial match list
  if (pathogen %in% partial_match_pathogens) {
    return(grepl(pathogen, name))  # Partial match for specified pathogens
  } else {
    return(pathogen == name)  # Exact match for all other pathogens
  }
}

# Apply the matching function to check if Pathogen and Name match
input_file_fn <- input_file_fn %>%
  mutate(
    PathogenMatch = mapply(match_pathogen, Pathogen, Name)
  )

# Create file for NK
input_file_fn_nk <- input_file_fn %>%
  filter(Sample.type == 'NK')

# Filter for TRUE matches of PathogenMatch
input_file_fn <- input_file_fn %>%
  filter(PathogenMatch == 'TRUE')

# Keep only samples with Source.Name present in false_negative
filtered_data_fn <- input_file_fn %>%
  semi_join(false_negative %>% select(Source.Name), by = "Source.Name")




# Filter for False Negative and False Positive samples and return their source.name
true_positive <- joined_data %>%
  filter(Result %in% c("True Positive")) %>%
  select(Source.Name, Result)

# Display the source.name and result for False Negative and False Positive samples
print(true_positive)

# Optionally save the False Negative and False Positive source names to a CSV file
write.csv(true_positive, "true_positive.csv", row.names = FALSE)

# Standardize 'Name' column based on the prefixes
input_file_tp <- input_file_tp %>%
  mutate(
    SimplifiedName = sapply(Name, function(n) {
      match <- prefixes[sapply(prefixes, function(p) grepl(paste0("^", p), n))]
      if (length(match) > 0) match else n
    })
  )


# Filter for False Negative and False Positive samples and return their source.name
true_negative <- joined_data %>%
  filter(Result %in% c("True Negative")) %>%
  select(Source.Name, Result)

# Display the source.name and result for False Negative and False Positive samples
print(true_negative)

# Optionally save the False Negative and False Positive source names to a CSV file
write.csv(true_negative, "true_negative.csv", row.names = FALSE)



# Matching function for specified pathogens (partial match) and others (exact match)
match_pathogen_tp <- function(pathogen, name) {
  
  # Define the list of pathogens that require partial matching
  partial_match_pathogens <- c("s_Coxiella", "s_Enterobacter cloacae", "s_Bartonella", "s_Streptococcus", 
                               "s_Staphylococcus", "s_Achromobacter", "s_Babesia", "s_Leptospira")
  
  # Check if pathogen is in the partial match list
  if (pathogen %in% partial_match_pathogens) {
    return(grepl(pathogen, name))  # Partial match for specified pathogens
  } else {
    return(pathogen == name)  # Exact match for all other pathogens
  }
}

# Apply the matching function to check if Pathogen and Name match
input_file_tp <- input_file_tp %>%
  mutate(
    PathogenMatch = mapply(match_pathogen, Pathogen, Name)
  )

# Filter for TRUE matches of PathogenMatch
input_file_tp <- input_file_tp %>%
  filter(PathogenMatch == 'TRUE')

# Keep only samples with Source.Name present in false_negative
filtered_data_tp <- input_file_tp %>%
  semi_join(true_positive %>% select(Source.Name), by = "Source.Name")

# Filter out duplicates based on Source.Name and keep the row with the highest LogPathogenScorePreNK
filtered_data_tp_unique <- filtered_data_tp %>%
  group_by(Source.Name) %>%
  filter(LogPathogenScorePreNK == max(LogPathogenScorePreNK, na.rm = TRUE)) %>%
  ungroup()


# Combine the filtered data with the NK samples
final_filtered_data_fn <- bind_rows(filtered_data_fn, input_file_fn_nk)

# Filter SAMPLE rows and collect unique microorganism names per Run.no
valid_microorganisms_per_run <- filtered_data_fn %>%
  filter(Sample.type != "NK") %>%  # Exclude NK samples
  group_by(Run.no) %>%
  summarise(Valid_Names = list(unique(Name)), .groups = "drop")

# Transform valid_microorganisms_per_run to long format
valid_microorganisms_long <- valid_microorganisms_per_run %>%
  unnest(Valid_Names)  # Convert the list column to long format

# Function to filter NK samples based on valid microorganisms for each Run.no
filter_valid_nk_samples <- function(nk_samples, valid_microorganisms) {
  nk_samples %>%
    left_join(valid_microorganisms, by = "Run.no", relationship = "many-to-many") %>%
    group_by(Run.no) %>%  # Group by Run.no to compare within the same run
    filter(Name %in% Valid_Names) %>%  # Check if Name matches any valid names for the run
    ungroup() %>%  # Remove the grouping after filtering
    select(-Valid_Names)  # Drop the Valid_Names column after filtering
}

# Apply the function to filter the NK samples
filtered_nk_samples <- filter_valid_nk_samples(input_file_fn_nk, valid_microorganisms_long)

# Combine the filtered NK samples with the other filtered data
final_filtered_data_fn <- bind_rows(filtered_data_fn, filtered_nk_samples)



# Add a new column to label the samples as 'FN' (False Negative) or 'TP' (True Positive)
filtered_data_fn <- filtered_data_fn %>%
  mutate(Group = "FN")

filtered_data_tp_unique <- filtered_data_tp_unique %>%
  mutate(Group = "TP")

# Combine the two datasets into one
combined_data <- bind_rows(filtered_data_fn, filtered_data_tp_unique)


# Function to perform the statistical test and format the p-value
perform_stat_test <- function(data, metric) {
  test_result <- wilcox.test(data[[metric]] ~ data$Group, data = data)
  p_value <- test_result$p.value
  # Format p-value and return as text, with values below 0.001 shown as "< 0.001"
  p_label <- ifelse(p_value < 0.001, paste("< 0.001"), paste("p =", round(p_value, 3)))
  return(p_label)
}

# Manually calculate Wilcoxon test for each metric by group
metrics <- c("NormCladeReads", "Kmers", "KmerDuplicity", "KmerCoverage")
test_results <- lapply(metrics, function(metric) {
  # Perform Wilcoxon test for the given metric
  test <- wilcox.test(combined_data[[metric]] ~ combined_data$Group)
  # Return test result including p-value
  return(list(p_value = test$p.value, statistic = test$statistic))
})

# Print the test results for each metric
test_results


# Create boxplots for each metric with 'Group' as a factor for comparison
metrics <- c("NormCladeReads", "Kmers", "KmerDuplicity", "KmerCoverage")
p_values <- sapply(metrics, function(metric) perform_stat_test(combined_data, metric))

# Create the boxplots with legends for the first plot only
plot_norm_clade_reads_box <- ggplot(combined_data, aes(x = Group, y = NormCladeReads, fill = Group)) +
  geom_boxplot(color = "black", alpha = 1) +
  scale_fill_manual(
    values = c("FN" = "lightcoral", "TP" = "lightgreen"),
    labels = c("FN" = "CS(-) | mREV (+)", "TP" = "CS(+) | mREV (+)")
  ) +
  scale_y_log10() +
  scale_x_discrete(labels = c(
    "FN" = "CS(-)\nmREV (+)",
    "TP" = "CS(+)\nmREV (+)"
  )) +
  labs(title = "Normalized reads", 
       subtitle = bquote(italic(p) == .(p_values["NormCladeReads"])), 
       x = "Group", 
       y = expression("Count (Log"[10]*")"), 
       fill = "Group") +  # Change the legend title here
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10))


plot_kmers_box <- ggplot(combined_data, aes(x = Group, y = Kmers, fill = Group)) +
  geom_boxplot(color = "black", alpha = 1) +
  scale_fill_manual(
    values = c("FN" = "lightcoral", "TP" = "lightgreen"),
    labels = c("FN" = "CS(-) | mREV (+)", "TP" = "CS(+) | mREV (+)")
  ) +
  scale_y_log10() +
  scale_x_discrete(labels = c(
    "FN" = "CS(-)\nmREV (+)",
    "TP" = "CS(+)\nmREV (+)"
  )) +
  labs(title = expression(italic(k)*"mer count"),
       subtitle = bquote(italic(p) == .(p_values["Kmers"])), 
       x = "Group", 
       y = expression("Count (Log"[10]*")"),  # Add Log10 to the y-axis title
       fill = "Group") +
       theme_minimal() +
  theme(legend.position = "none",
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

plot_kmer_duplicity_box <- ggplot(combined_data, aes(x = Group, y = KmerDuplicity, fill = Group)) +
  geom_boxplot(color = "black", alpha = 1) +
  scale_y_log10() +
  scale_fill_manual(
    values = c("FN" = "lightcoral", "TP" = "lightgreen"),
    labels = c("FN" = "CS(-) | mREV (+)", "TP" = "CS(+) | mREV (+)")
  ) +
  scale_y_log10() +
  scale_x_discrete(labels = c(
    "FN" = "CS(-)\nmREV (+)",
    "TP" = "CS(+)\nmREV (+)"
  )) +
  labs(title = expression(italic(k)*"mer duplicity"),
       subtitle = bquote(italic(p) == .(p_values["KmerDuplicity"])), 
       x = "Group", 
       y = expression("Count (Log"[10]*")"),  # Add Log10 to the y-axis title
  fill = "Group") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

plot_kmer_coverage_box <- ggplot(combined_data, aes(x = Group, y = KmerCoverage, fill = Group)) +
  geom_boxplot(color = "black", alpha = 1) +
  scale_y_log10() +
  scale_fill_manual(
    values = c("FN" = "lightcoral", "TP" = "lightgreen"),
    labels = c("FN" = "CS(-) | mREV (+)", "TP" = "CS(+) | mREV (+)")
  ) +
  scale_y_log10() +
  scale_x_discrete(labels = c(
    "FN" = "CS(-)\nmREV (+)",
    "TP" = "CS(+)\nmREV (+)"
  )) +
  labs(title = expression(italic(k)*"mer coverage"),
       subtitle = bquote(italic(p) == .(p_values["KmerCoverage"])), 
       x = "Group", 
       y = expression("Ratio (Log"[10]*")"), # Add Log10 to the y-axis title
  fill = "Group") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

# Combine the plots into a 2x2 grid with spacing between top and bottom rows
combined_plots <- plot_norm_clade_reads_box + plot_kmers_box + 
  plot_kmer_duplicity_box + plot_kmer_coverage_box +
  plot_layout(ncol = 2, nrow = 2) +
  plot_layout(guides = 'collect') & theme(legend.position = 'right') &
  theme(plot.margin = margin(10, 10, 20, 10))  # Add space between top and bottom

# Add labels to the plots (A, B, C, D) using patchwork's annotation
statistical_significance_metrics <- combined_plots + 
  plot_annotation(
    tag_levels = 'a'  # Automatically add tags A, B, C, D to the individual plots
  )

# Display the final plot
print(statistical_significance_metrics)
ggsave("statistical_significance_metrics.png", plot = statistical_significance_metrics, width = 7, height = 6)
# 
# Calculate means for each variable grouped by Group
mean_values <- combined_data %>%
  group_by(Group) %>%
  summarise(
    Mean_NormCladeReads = mean(NormCladeReads, na.rm = TRUE),
    Mean_Kmers = mean(Kmers, na.rm = TRUE),
    Mean_KmerDuplicity = mean(KmerDuplicity, na.rm = TRUE),
    Mean_KmerCoverage = mean(KmerCoverage, na.rm = TRUE)
  )
