# ClinSeq

**ClinSeq** is a bioinformatics tool for ranking microorganisms detected in Next-Generation Sequencing (NGS) data.  
It integrates k-mer statistics from **KrakenUniq** and calculates a composite score for each detection based on:

- K-mer counts  
- Duplicity  
- Coverage  
- Normalized read counts  

ClinSeq assigns this score to each detected microorganism and produces a **ranked list** from highest to lowest score.  
This ranking helps prioritize potentially important detections in complex NGS datasets.

---

## Features

- Uses **KrakenUniq** output as input  
- Calculates a reproducible score for each detection  
- Accounts for **k-mer count, duplicity, coverage, and normalized read counts**  
- Produces an easy-to-interpret **ranked list of organisms**  
- Lightweight and designed for integration into NGS analysis workflows  

---

## Dependencies

ClinSeq was developed and tested in **R**. To ensure reproducibility, the following versions were used:

- **R version**: 4.4.0  
- **Required R packages / libraries**:
  - `tidyverse`: 2.0.0
  - `gridExtra`: 2.3  
  - `patchwork`: 1.3.0  
  - `extrafont`: 0.10 
  - `caret`: 7.0.1 

---

## Multiple Output Files

Running ClinSeq will generate **multiple output files**, including:

### Plots

- `average_plot.png` – Graph plotting average score per rank  
- `average_plot_per_run.png` – Graph plotting average score per rank per run  
- `deviation_from_average_plot.png` – Graph plotting the deviation of each score from the average  
- `deviation_from_average_plot_per_run.png` – Graph plotting the deviation per run
- `Agreement plotXXXX.png` – Plot showing `True positive`, `True negative`, `False positive`, and `False negative` detections. 
- `statistical_significance_metrics.png` – Graph plotting showing statistical significance between the `True positive` and `False negative` groups. 

### CSV Files

- `final_combined_match_status.csv` – Overall matching status of all samples  
- `matched_pathogens_individual_matches.csv` – Matches for individual ranks  
- `matched_pathogens_with_overall_match.csv` – Matches for individual ranks + overall match column  
- `matched_pathogens_with_result.csv` – Matches for individual ranks + overall match + final result  
- `Run_X_Potential_Contaminants.csv` – Potential contaminants per run (one file per run)  
- `top_filtered_pathogens_by_log_score.csv` – Top-ranked pathogens by log score  
- `top_filtered_pathogens_by_log_score_contaminants.csv` – Top-ranked pathogens + potential contaminants column  
- `top_filtered_pathogens_by_log_score_contaminants_target.csv` – Top-ranked pathogens + contaminants + target microorganisms  
- `true_positive.csv` – List of all true positive detections.  
- `true_negative.csv` – List of all true negative detections.  
- `false_positive.csv` – List of all false positive detections.  
- `false_negative.csv` – List of all false negative detections.  

### TXT Files

- `Cohens_Kappa_Interpretation.txt` – Info about Cohens Kappa when comparing ClinSeq to manual review or known targets.  
- `parameters_log.txt` – Log file with variables used for ClinSeq calculations.  


---

## Installation & Usage

Simply download the R script and adjust the input parameters at the start of the code.

Make sure to prepare the input .csv file according to the list below.

### Input .csv file column names and order

- **Source.Name** – Sample identifier (unique ID)  
- **Sample.type** – Sample type or category (e.g., SAMPLE, NK, HPC)  
- **Pathogen** – Known pathogen (if available)  
- **Detected.pat** – Whether the target microorganism is known or has been detected through manual analysis (`YES`/`NO`)  
- **Run.no** – Sequencing run number (multiple runs can be analyzed simultaneously)  
- ✦ **NormCladeReads** – Normalized reads assigned to the clade (normalized to 10 million total reads)  
- ✦ **NormTaxonReads** – Normalized reads assigned to the taxon (normalized to 10 million total reads)  
- † **TaxRank** – Taxonomic rank from KrakenUniq (species, genus, family, etc.)  
- † **TaxID** – Taxonomy ID  
- † **Name** – Scientific name of the taxon  
- † **Kmers** – Number of k-mers observed  
- † **KmerDuplicity** – Average multiplicity of k-mers  
- † **KmerCoverage** – Fraction of unique k-mers covered  
- † **Depth** – Taxonomic level  
- † **TaxLineage** – Full taxonomic lineage  

> Items marked with † are unchanged KrakenUniq outputs, ✦ are calculated from KrakenUniq outputs, and unmarked items are manually added when constructing the final “master” CSV file.

---

## Example Input

### Known Pathogens

| Source.Name | Sample.type | Pathogen      | Detected.pat | Run.no | NormCladeReads  | NormTaxonReads   | TaxRank | TaxID | Name   | Kmers     | KmerDuplicity | KmerCoverage | Depth | TaxLineage |
|------------|-------------|---------------|--------------|--------|----------------|-----------------|---------|-------|--------|-----------|---------------|--------------|-------|------------|
| 0001       | SAMPLE      | Puumala virus | NO           | 1      | 1039.574819    | 22.42737554     | -       | 1     | -_root | 38108557  | 1.91          | 0.0007985    | 0     | -_root     |
| 0002       | NK          | NK            |              | 1      | 11556.32958    | 1065.312947     | -       | 1     | -_root | 30142     | 1.84          | 6.316E-07    | 0     | -_root     |
| 0003       | HPC         | HPC           |              | 1      | 3480974.202    | 1004.772273     | -       | 1     | -_root | 1368030   | 47.9          | 0.00002866   | 0     | -_root     |

### Unknown Pathogens

| Source.Name | Sample.type | Pathogen | Detected.pat | Run.no | NormCladeReads      | NormTaxonReads        | TaxRank | TaxID | Name   | Kmers     | KmerDuplicity | KmerCoverage | Depth | TaxLineage |
|------------|-------------|----------|--------------|--------|--------------------|----------------------|---------|-------|--------|-----------|---------------|--------------|-------|------------|
| 0001        | SAMPLE      |          |              | 1      | 1350972.1497       | 2848.46096           | -       | 1     | -_root | 16705211  | 15            | 0.00035      | 0     | -_root     |
| 0002        | NK          |          |              | 1      | 2522121.4314       | 2782.56998           | -       | 1     | -_root | 259500    | 2.5           | 0.000005437  | 0     | -_root     |
| 0003       | HPC         |          |              | 1      | 3963956.9891       | 1324.26579           | -       | 1     | -_root | 1752282   | 180           | 0.00003672   | 0     | -_root     |
