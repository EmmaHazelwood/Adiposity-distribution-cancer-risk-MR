sink("Rscript2")

library(dplyr)
library(TwoSampleMR)
library(data.table)
library(LDlinkR)
library(tidyr)
library(gwasvcf)


# Run updated proxy function ----------------------------------------------

#' @title proxy_search: search for proxies when using read_outcome_data()
#' @description
#' This function gets proxy SNPs for MR analyses when using local outcome data.
#' It takes formatted exposure and outcome data from `TwoSampleMR`, searches for
#' missing SNPs in the outcome, identifies which of these missing-SNPs are
#' available in the provided reference panel, extracts all proxy-SNPs for the
#' missing-SNPs from the reference panel, returns a data frame with the top
#' proxy-SNP for each missing-SNP. You MUST have a local reference panel. Only
#' works with rsID.
#' @return a data frame
#' @param data_exposure exposure data frame
#' @param data_outcome outcome data frame
#' @param data_outcome_path file path used for `read_outcome_data()`
#' @param data_reference reference data; bim file
#' @param data_reference_path file path for your downloaded reference panel
#' @param tag_r2 r2 for proxy SNP; from `get_ld_proxies()`; default = 0.8
#' @param tag_kb window to look for proxy SNPs; from `get_ld_proxies()`; default = 5000
#' @param tag_nsnp from `get_ld_proxies()`; default = 5000
#' @param outcome_sep separator for your outcome GWAS
#' @param outcome_phenotype phenotype column name of your GWAS
#' @param outcome_SNP SNP column name of your GWAS
#' @param outcome_BETA BETA column name of your GWAS
#' @param outcome_SE SE column name of your GWAS
#' @param outcome_P P column name of your GWAS
#' @param outcome_EA EA column name of your GWAS
#' @param outcome_OA OA column name of your GWAS
#' @param outcome_EAF EAF column name of your GWAS
#' @param outcome_N N column name of your GWAS
#' @param outcome_ID ID column name of your GWAS
#' @param outcome_CHR CHR column name of your GWAS
#' @param outcome_POS POS column name of your GWAS
proxy_search <- function(data_exposure, data_outcome, data_outcome_path, data_reference, data_reference_path,
                         tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                         outcome_sep, outcome_phenotype, outcome_SNP, outcome_BETA, outcome_SE, outcome_P,
                         outcome_EA, outcome_OA, outcome_EAF, outcome_N, outcome_ID, outcome_CHR, outcome_POS) {
  
  # Parameter Validation
  if (!file.exists(data_outcome_path) || !file.exists(paste0(data_reference))) {
    stop("Invalid file paths provided.")
  }
  
  if (!is.data.frame(data_exposure) || !is.data.frame(data_outcome)) {
    stop("data_exposure and data_outcome must be data frames.")
  }
  
  # exposure snps missing from outcome ====
  snps_missing <- setdiff(unique(as.factor(data_exposure$SNP)), unique(as.factor(data_outcome$SNP)))
  if (length(snps_missing) == 0) {
    message("No missing SNPs. Returning the original dataframe.")
    return(data_outcome)
  }
  
  # look-up missing SNPs in reference panel ====
  message(paste0("# 1. looking up ", length(unique(as.factor(snps_missing))), " missing-SNP(s) in the reference panel"))
  reference_header <- readLines(data_reference, n = 1) ## read the first row to identify the column containing "rs*"
  column_index <- grep("rs", strsplit(reference_header, "\t")[[1]]) ## get column index
  ## bash command to filter rows based on the identified column
  cmd <- sprintf("awk '$%d ~ /^(%s)$/' %s", column_index, paste(snps_missing, collapse = "|"), paste0(data_reference))
  ## use system to run the bash command and then read the result with fread
  reference <- data.table::fread(cmd = cmd, quote = "")
  if (nrow(reference) != 0) {
  snps_reference <- intersect(unique(as.factor(snps_missing)), unique(as.factor(reference$V2)))
  message(paste0("## ", length(unique(as.factor(snps_reference))), " of ", length(unique(as.factor(snps_missing))), " missing-SNP(s) available in the reference panel"))
  
  # find proxies for available SNPs ====
  message(paste0("# 2. extracting proxy-SNP(s) for the ", length(unique(as.factor(snps_reference))), " missing-SNP(s) from the reference panel"))
  gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
  proxies <- functions::get_ld_proxies(rsid = snps_missing,
                                       bfile = data_reference_path,
                                       searchspace = NULL,
                                       tag_kb = tag_kb,
                                       tag_nsnp = tag_nsnp,
                                       tag_r2 = tag_r2,
                                       threads = 1,
                                       out = tempfile())
  ## format proxy data: change column order and names, add proxy.outcome = TRUE
  proxies <- proxies %>%
    dplyr::select(target_snp.outcome = SNP_A,
                  proxy_snp.outcome = SNP_B,
                  target_a1.outcome = A1,
                  target_a2.outcome = A2,
                  proxy_a1.outcome = B1,
                  proxy_a2.outcome = B2,
                  R) %>%
    dplyr::mutate(proxy.outcome = TRUE,
                  SNP = proxy_snp.outcome) %>%
    dplyr::select(proxy.outcome, everything())
  message(paste0("## proxy-SNP(s) for ", length(unique(proxies$target_snp.outcome)), " missing-SNP(s) found; ", "proxy-SNP(s) for ", length(unique(as.factor(snps_reference))) - length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) not available (e.g., no proxy-SNP or r2 < provided)"))
  
  # extract proxies from outcome ====
  message(paste0("# 3. extracting proxy-SNP(s) from outcome"))
  proxy_snps <- proxies %>% # select unique proxy SNPs to extract
    dplyr::distinct(proxy_snp.outcome) %>%
    dplyr::pull(proxy_snp.outcome)
  proxy_snps<-c(proxy_snps,data_outcome$SNP[1])
  data_outcome_proxies <- TwoSampleMR::read_outcome_data(filename = data_outcome_path,
                                                         snps = proxy_snps,
                                                         sep = outcome_sep,
                                                         phenotype_col = outcome_phenotype,
                                                         snp_col = outcome_SNP,
                                                         beta_col = outcome_BETA,
                                                         se_col = outcome_SE,
                                                         eaf_col = outcome_EAF,
                                                         effect_allele_col = outcome_EA,
                                                         other_allele_col = outcome_OA,
                                                         pval_col = outcome_P,
                                                         samplesize_col = outcome_N,
                                                         id_col = outcome_ID,
                                                         chr_col = outcome_CHR,
                                                         pos_col = outcome_POS)
  data_outcome_proxies<-data_outcome_proxies[-nrow(data_outcome_proxies),]
  data_outcome_proxies <- dplyr::left_join(data_outcome_proxies, proxies, by = c("SNP" = "SNP"))
  message(paste0("## proxy-SNP(s) for ", length(unique(as.factor(data_outcome_proxies$target_snp.outcome))), " of ", length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) extracted"))
  
  # select proxy-SNP(s) with the highest R2 ====
  data_outcome_proxies <- data_outcome_proxies %>%
    dplyr::group_by(target_snp.outcome) %>%
    dplyr::filter(R == max(R)) %>%
    dplyr::slice(1) %>%
    dplyr::select(-R)
  
  ## Bind rows of data_outcome with data_outcome_proxies
  data_outcome <- dplyr::bind_rows(data_outcome, data_outcome_proxies)
  
  return(data_outcome)
  } else {
  return(data_outcome)
}}


results<-data.frame(matrix(nrow=0,ncol=0))


#BMI
exposure_dat <- extract_instruments(outcomes="ieu-b-40", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Esophagus adenocarcinoma
print("Esophagus adenocarcinoma")
outfilepath<-"data/Obesity_related_cancer_GWAS/Esophagus_adenocarcinoma/GCST003739_buildGRCh37.tsv"

outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele"
)

if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="variant_id", outcome_BETA="beta", outcome_SE="standard_error", outcome_P="p_value",
                            outcome_EA="effect_allele", outcome_OA="other_allele",outcome_CHR="chromosome", outcome_POS="base_pair_location",outcome_EAF="missing",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}


dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$ncase.outcome<-4112
dat$ncontrol.outcome<-17159
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Esophagus adenocarcinoma"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Colorectal
print("Colorectal")
outfilepath<-"data/gecco/annotated/overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt"

outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=" ", outcome_phenotype="exposure", outcome_SNP="SNP", outcome_BETA="Effect", outcome_SE="StdErr", outcome_P="P.value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="Chromosome", outcome_POS="Position",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-98715
dat$ncontrol.outcome<-52775
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Colorectal cancer (overall)"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Liver
print("Liver")
outfilepath<-"data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS_formatted.csv"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1")
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="MarkerName", outcome_BETA="beta", outcome_SE="se", outcome_P="P-value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="Chromosome", outcome_POS="Position",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-862
dat$ncontrol.outcome<-715717
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Liver cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



#Gallbladder
print("Gallbladder")
outfilepath<-"data/Obesity_related_cancer_GWAS/Gallbladder_cancer/Meta_analysis_Gallbladder_cancer_GWAS_formatted.csv"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="MarkerName", outcome_BETA="beta", outcome_SE="se", outcome_P="P-value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="Chromosome", outcome_POS="Position",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-279
dat$ncontrol.outcome<-715718
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Gallbladder cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



#Breast
print("Breast")
outfilepath<-"data/Obesity_related_cancer_GWAS/Breast_cancer/annotated/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_annotated.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "beta.Gwas",
  se_col = "SE.Gwas",
  effect_allele_col = "Effect.Gwas",
  other_allele_col = "Baseline.Gwas",
  eaf_col = "Freq.Gwas"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="rsID", outcome_BETA="beta.Gwas", outcome_SE="SE.Gwas", outcome_P="missing",
                            outcome_EA="Effect.Gwas", outcome_OA="Baseline.Gwas",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq.Gwas",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-133384
dat$ncontrol.outcome<-113789
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Endometrial
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006464", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-12906
dat$ncontrol.outcome<-108979
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Endometrial cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Ovarian
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1120", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-25509
dat$ncontrol.outcome<-40941
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Kidney (renal-cell)
print("Kidney")
outfilepath<-"data/Obesity_related_cancer_GWAS/Kidney_cancer/Meta_analysis_Kidney_cancer_GWAS_formatted.csv"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="MarkerName", outcome_BETA="beta", outcome_SE="se", outcome_P="P-value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-2085
dat$ncontrol.outcome<-635902
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Kidney (renal-cell) cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Meningioma
print("Meningioma")
outfilepath<-"data/Obesity_related_cancer_GWAS/Meningioma/finngen_R8_C3_MENINGIOMA_EXALLC"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="beta", outcome_SE="sebeta", outcome_P="mising",
                            outcome_EA="alt", outcome_OA="ref",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="af_alt",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1012
dat$ncontrol.outcome<-259583
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Meningioma"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Thyroid
print("Thyroid")
outfilepath<-"data/Obesity_related_cancer_GWAS/Thyroid_cancer/Meta_analysis_Thyroid_cancer_GWAS_formatted.csv"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="MarkerName", outcome_BETA="beta", outcome_SE="se", outcome_P="P-value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1904
dat$ncontrol.outcome<-715554
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Thyroid cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Multiple myeloma
print("Multiple myeloma")
outfilepath<-"data/Obesity_related_cancer_GWAS/Multiple_myeloma/Meta_analysis_Multiple_myeloma_cancer_GWAS_formatted.csv"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "MarkerName",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="MarkerName", outcome_BETA="beta", outcome_SE="se", outcome_P="P-value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1649
dat$ncontrol.outcome<-727247
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Multiple myeloma"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Pancreas cancer
print("Pancreas cancer")
outfilepath<-"data/PanScan/PanScan/PanScan_with_rsid.csv"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "V2",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="V2", outcome_BETA="BETA", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="ALT", outcome_OA="REF",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="missing",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$ncase.outcome<-4951
dat$ncontrol.outcome<-3456
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Pancreas cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Colorectal
#colon
print("colon")
outfilepath<-"data/gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "data/gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=" ", outcome_phenotype="exposure", outcome_SNP="SNP", outcome_BETA="Effect", outcome_SE="StdErr", outcome_P="P.value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="Chromosome", outcome_POS="Position",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-28736
dat$ncontrol.outcome<-43099
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Colon cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#proximal
print("proximal")
outfilepath<-"data/gecco/annotated/proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "data/gecco/annotated/proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=" ", outcome_phenotype="exposure", outcome_SNP="SNP", outcome_BETA="Effect", outcome_SE="StdErr", outcome_P="P.value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="Chromosome", outcome_POS="Position",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-14416
dat$ncontrol.outcome<-43099
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Proximal colon cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#distal
print("distal")
outfilepath<-"data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=" ", outcome_phenotype="exposure", outcome_SNP="SNP", outcome_BETA="Effect", outcome_SE="StdErr", outcome_P="P.value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="Chromosome", outcome_POS="Position",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-12879
dat$ncontrol.outcome<-43099
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Distal colon cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#rectal
print("rectal")
outfilepath<-"data/gecco/annotated/rectal_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=" ", outcome_phenotype="exposure", outcome_SNP="SNP", outcome_BETA="Effect", outcome_SE="StdErr", outcome_P="P.value",
                            outcome_EA="Allele1", outcome_OA="Allele2",outcome_CHR="Chromosome", outcome_POS="Position",outcome_EAF="Freq1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-14150
dat$ncontrol.outcome<-43099
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Rectal cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Endometrial
#Endometrioid
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006465", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-8758
dat$ncontrol.outcome<-46126
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Endometrioid endometrial cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Non-endometrioid
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006466", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1230
dat$ncontrol.outcome<-35447
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Non-endometrioid endometrial cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Breast
#Luminal A-like
print("luminal a-like")
outfilepath<-"data/Obesity_related_cancer_GWAS/Breast_cancer/annotated/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "Luminal_A_log_or_meta",
  se_col = "Luminal_A_se_meta",
  effect_allele_col = "Effect.Meta",
  other_allele_col = "Baseline.Meta",
  eaf_col = "Freq.Gwas"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="rsID", outcome_BETA="Luminal_A_log_or_meta", outcome_SE="Luminal_A_se_meta", outcome_P="missing",
                            outcome_EA="Effect.Meta", outcome_OA="Baseline.Meta",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq.Gwas",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-63767
dat$ncontrol.outcome<-91477
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Luminal A-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Luminal B/HER2-negative-like
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "Luminal_B_HER2Neg_log_or_meta",
  se_col = "Luminal_B_HER2Neg_se_meta",
  effect_allele_col = "Effect.Meta",
  other_allele_col = "Baseline.Meta",
  eaf_col = "Freq.Gwas"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="rsID", outcome_BETA="Luminal_B_HER2Neg_log_or_meta", outcome_SE="Luminal_B_HER2Neg_se_meta", outcome_P="missing",
                            outcome_EA="Effect.Meta", outcome_OA="Baseline.Meta",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq.Gwas",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-15942
dat$ncontrol.outcome<-91477
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Luminal B-HER2-negative-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Luminal B-like
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "Luminal_B_log_or_meta",
  se_col = "Luminal_B_se_meta",
  effect_allele_col = "Effect.Meta",
  other_allele_col = "Baseline.Meta",
  eaf_col = "Freq.Gwas"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="rsID", outcome_BETA="Luminal_B_log_or_meta", outcome_SE="Luminal_B_se_meta", outcome_P="missing",
                            outcome_EA="Effect.Meta", outcome_OA="Baseline.Meta",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq.Gwas",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-15942
dat$ncontrol.outcome<-91477
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Luminal B-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#HER2-enriched-like
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "HER2_Enriched_log_or_meta",
  se_col = "HER2_Enriched_se_meta",
  effect_allele_col = "Effect.Meta",
  other_allele_col = "Baseline.Meta",
  eaf_col = "Freq.Gwas"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="rsID", outcome_BETA="HER2_Enriched_log_or_meta", outcome_SE="HER2_Enriched_se_meta", outcome_P="missing",
                            outcome_EA="Effect.Meta", outcome_OA="Baseline.Meta",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq.Gwas",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-10628
dat$ncontrol.outcome<-91477
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "HER2-enriched-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Triple negative or basal-like
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "rsID",
  beta_col = "Triple_Neg_log_or_meta",
  se_col = "Triple_Neg_se_meta",
  effect_allele_col = "Effect.Meta",
  other_allele_col = "Baseline.Meta",
  eaf_col = "Freq.Gwas"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep=",", outcome_phenotype="exposure", outcome_SNP="rsID", outcome_BETA="Triple_Neg_log_or_meta", outcome_SE="Triple_Neg_se_meta", outcome_P="missing",
                            outcome_EA="Effect.Meta", outcome_OA="Baseline.Meta",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="Freq.Gwas",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$id.outcome<-outcome_dat$id.outcome[1]
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-8602
dat$ncontrol.outcome<-91477
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Triple negative or basal-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Ovarian
#High grade serous carcinoma - ieu-a-1121
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1121", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-25509
dat$ncontrol.outcome<-40941
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "High grade serous carcinoma ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Low grade serous carcinoma - ieu-a-1122
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1122", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1012
dat$ncontrol.outcome<-40941
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Low grade serous carcinoma ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Mucinous - ieu-a-1123
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1123", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1417
dat$ncontrol.outcome<-40941
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Invasive mucinous ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Endometrioid - ieu-a-1125
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1125", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-2810
dat$ncontrol.outcome<-40941
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Endometrioid ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Clear cell - ieu-a-1124
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1124", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1366
dat$ncontrol.outcome<-40941
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Clear cell ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Low malignant potential tumours - ieu-a-1233
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1233", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-3103
dat$ncontrol.outcome<-40941
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Low malignant potential ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

fwrite(results,"results/LiverFatMR/BMI_obesity_cancers_results.csv")

sink()
