# Load required packages
library(data.table)
library(dplyr)

# Function to compute pairwise Q tests for subtypes vs overall
pairwise_q <- function(df, group_name, exposures, subtypes) {
  out <- data.frame(
    Trait = character(),
    Cancer = character(),
    Subtype = character(),
    Q_statistic = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  for (exposure_trait in exposures) {
    df_sub <- df %>%
      filter(exposure == exposure_trait, outcome %in% c(group_name, subtypes))
    
    overall <- df_sub %>% filter(outcome == group_name)
    others <- df_sub %>% filter(outcome %in% subtypes)
    
    if (nrow(overall) == 1 && nrow(others) > 0) {
      others_q <- others %>%
        rowwise() %>%
        mutate(
          Q_statistic = ((b - overall$b)^2) / (se^2 + overall$se^2),
          P_value = pchisq(Q_statistic, df = 1, lower.tail = FALSE),
          Trait = gsub(" \\|\\|.*", "", exposure),  # remove ' || id:...'
          Cancer = group_name,
          Subtype = outcome
        ) %>%
        ungroup() %>%
        select(Trait, Cancer, Subtype, Q_statistic, P_value)
      out <- bind_rows(out, others_q)
    }
  }
  return(out)
}

# Define subtype groupings
cancer_groups <- list(
  "Endometrial cancer" = c("Endometrioid endometrial cancer", "Non-endometrioid endometrial cancer"),
  "Breast cancer" = c(
    "Luminal A-like breast cancer", "Luminal B-HER2-negative-like breast cancer",
    "Luminal B-like breast cancer", "HER2-enriched-like breast cancer",
    "Triple negative or basal-like breast cancer"
  ),
  "Ovarian cancer" = c(
    "High grade serous carcinoma ovarian cancer", "Low grade serous carcinoma ovarian cancer",
    "Invasive mucinous ovarian cancer", "Endometrioid ovarian cancer",
    "Clear cell ovarian cancer", "Low malignant potential ovarian cancer"
  ),
  "Colorectal cancer (overall)" = c("Colon cancer", "Proximal colon cancer", "Distal colon cancer", "Rectal cancer")
)

# Function to process one result file
process_file <- function(path) {
  df <- fread(path)
  df <- df %>% filter(method %in% c("Inverse variance weighted", "Wald ratio"))
  exposures <- unique(df$exposure)
  
  results <- data.frame(
    Trait = character(),
    Cancer = character(),
    Subtype = character(),
    Q_statistic = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cg in names(cancer_groups)) {
    subs <- cancer_groups[[cg]]
    results_sub <- pairwise_q(df, cg, exposures, subs)
    results <- bind_rows(results, results_sub)
  }
  return(results)
}

# Process both datasets
results_adiposity <- process_file("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
results_molecular <- process_file("results/LiverFatMR/Molecular_traits_cancers_results.csv")

# Combine and format final results
combined_results <- bind_rows(results_adiposity, results_molecular) %>%
  mutate(
    `Q statistic` = format(round(Q_statistic, 2), scientific = TRUE),
    `P value` = format(round(P_value, 2), scientific = TRUE)
  ) %>%
  select(Trait, Cancer, Subtype, `Q statistic`, `P value`)  # drop numeric Q_statistic/P_value

print(combined_results)

# Write output
fwrite(combined_results, "results/LiverFatMR/pairwise_Q_all_traits_cancers.csv")
