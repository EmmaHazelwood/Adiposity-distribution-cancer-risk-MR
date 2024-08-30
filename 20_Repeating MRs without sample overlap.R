sink("Rscript.txt")
library(dplyr)
library(TwoSampleMR)
library(data.table)
library(LDlinkR)
library(tidyr)

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
}



results<-data.frame(matrix(nrow=0,ncol=0))


# Measures of adiposity and cancers ---------------------------------------

#Thyroid
exposure_dat <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/VAT_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$chr.exposure<-NA
exposure_dat$pos.exposure<-NA
exposure_dat$samplesize.exposure<-38965
exposure_dat$exposure<-"VAT"
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Thyroid
print("Thyroid")
outfilepath<-"data/Obesity_related_cancer_GWAS/Thyroid_cancer/193_PheCode.v1.0.fastGWA"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="SNP", outcome_BETA="BETA", outcome_SE="SE", outcome_P="P",
                            outcome_EA="A1", outcome_OA="A2",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="AF1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1525
dat$ncontrol.outcome<-259585
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Thyroid cancer (FinnGen)"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Multiple myeloma
exposure_dat <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Pancreas fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$chr.exposure<-NA
exposure_dat$pos.exposure<-NA
exposure_dat$samplesize.exposure<-25617
exposure_dat$exposure<-"Pancreas fat"
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Multiple myeloma
print("Multiple myeloma")
outfilepath<-"data/Obesity_related_cancer_GWAS/Multiple_myeloma/204.4_PheCode.v1.0.fastGWA"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="SNP", outcome_BETA="BETA", outcome_SE="SE", outcome_P="P",
                            outcome_EA="A1", outcome_OA="A2",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="AF1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1085
dat$ncontrol.outcome<-271463
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
dat$outcome <- "Multiple myeloma (FinnGen)"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Liver
#Exposures
exposure_dat1 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Liver fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1$chr.exposure<-NA
exposure_dat1$pos.exposure<-NA
exposure_dat1$samplesize.exposure<-32858
exposure_dat1$exposure<-"Liver fat"

exposure_dat2 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/ASAT_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat2$chr.exposure<-NA
exposure_dat2$pos.exposure<-NA
exposure_dat2$samplesize.exposure<-38965
exposure_dat2$exposure<-"ASAT"

exposure_dat3 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/VAT_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat3$chr.exposure<-NA
exposure_dat3$pos.exposure<-NA
exposure_dat3$samplesize.exposure<-38965
exposure_dat3$exposure<-"VAT"

exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

print("Liver cancer")
outfilepath<-"data/Obesity_related_cancer_GWAS/Liver_cancer/155_PheCode.v1.0.fastGWA"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="SNP", outcome_BETA="BETA", outcome_SE="SE", outcome_P="P",
                            outcome_EA="A1", outcome_OA="A2",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="AF1",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
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
dat$outcome <- "Liver cancer (FinnGen)"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


# Measures of adiposity to molecular traits -------------------------------

#ASAT and GFAT
exposure_dat1 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/ASAT_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1$chr.exposure<-NA
exposure_dat1$pos.exposure<-NA
exposure_dat1$samplesize.exposure<-38965
exposure_dat1$exposure<-"ASAT"

exposure_dat2 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/GFAT_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat2$chr.exposure<-NA
exposure_dat2$pos.exposure<-NA
exposure_dat2$samplesize.exposure<-38965
exposure_dat2$exposure<-"GFAT"

exposure_dat<-rbind(exposure_dat1,exposure_dat2)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200


#HDL cholesterol
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST002223", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$outcome <- "HDL cholesterol (no UK Biobank)"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#IGF-1
print("IGF-1")
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/IGF1/2952_75_IGF1.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="missing",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35560
dat <- steiger_filtering(dat)
dat$outcome <- "IGF-1 no UK Biobank"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#ASAT, GFAT, Pancreas fat
exposure_dat3 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Pancreas fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat3$chr.exposure<-NA
exposure_dat3$pos.exposure<-NA
exposure_dat3$samplesize.exposure<-25617
exposure_dat3$exposure<-"Pancreas fat"

exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Triglycerides
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST002216", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$samplesize.outcome<-94595
dat <- steiger_filtering(dat)
dat$outcome <- "Triglycerides (no UK Biobank)"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



# Molecular traits to cancers ---------------------------------------------
#HDL -> triple negative breast cancer
exposure_dat <- extract_instruments("ebi-a-GCST002223")

#Triple negative or basal-like
print("ASAT #riple negative or basal-like")
outfilepath<-"data/Obesity_related_cancer_GWAS/Breast_cancer/annotated/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt"

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
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-8602
dat$ncontrol.outcome<-91477
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
missing_af<-dat[is.na(dat$eaf.exposure),]
missing_af1<-merge(missing_af,af_ref,by.x=c("SNP","effect_allele.exposure","other_allele.exposure"),by.y=c("ID","REF","ALT"))
missing_af2<-merge(missing_af,af_ref,by.x=c("SNP","effect_allele.exposure","other_allele.exposure"),by.y=c("ID","ALT","REF"))
missing_af<-rbind(missing_af1,missing_af2)
missing_af<-dplyr::select(missing_af,SNP,effect_allele.exposure,other_allele.exposure,ALT_FREQS)
dat<-merge(dat,missing_af,by=c("SNP","effect_allele.exposure","other_allele.exposure"),all.x=TRUE)
dat$eaf.exposure[is.na(dat$eaf.exposure)]<-dat$ALT_FREQS[is.na(dat$eaf.exposure)]
still_missing<-dat[is.na(dat$eaf.exposure),]
dat<-dat[!is.na(dat$eaf.exposure),]
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
still_missing$units.outcome<-NA
still_missing$units.exposure<-NA
still_missing$rsq.exposure<-NA
still_missing$effective_n.exposure<-NA
still_missing$rsq.outcome<-NA
still_missing$steiger_dir<-TRUE
still_missing$steiger_pval<-NA
still_missing$r.outcome<-NA
dat<-rbind(dat,still_missing)  
dat$outcome <- "Triple negative or basal-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#IGF-1 -> Luminal B/Her 2 breast cancer
exposure_dat <-read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/IGF1/IGF_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$chr.exposure<-NA
exposure_dat$pos.exposure<-NA
exposure_dat$samplesize.exposure<-38965
exposure_dat$exposure<-"VAT"
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Luminal B/HER2-negative-like
print("ASAT #Luminal B/HER2-negative-like")
outfilepath<-"data/Obesity_related_cancer_GWAS/Breast_cancer/annotated/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt"

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
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-15942
dat$ncontrol.outcome<-91477
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
missing_af<-dat[is.na(dat$eaf.exposure),]
missing_af1<-merge(missing_af,af_ref,by.x=c("SNP","effect_allele.exposure","other_allele.exposure"),by.y=c("ID","REF","ALT"))
missing_af2<-merge(missing_af,af_ref,by.x=c("SNP","effect_allele.exposure","other_allele.exposure"),by.y=c("ID","ALT","REF"))
missing_af<-rbind(missing_af1,missing_af2)
missing_af<-dplyr::select(missing_af,SNP,effect_allele.exposure,other_allele.exposure,ALT_FREQS)
dat<-merge(dat,missing_af,by=c("SNP","effect_allele.exposure","other_allele.exposure"),all.x=TRUE)
dat$eaf.exposure[is.na(dat$eaf.exposure)]<-dat$ALT_FREQS[is.na(dat$eaf.exposure)]
still_missing<-dat[is.na(dat$eaf.exposure),]
dat<-dat[!is.na(dat$eaf.exposure),]
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
dat <- steiger_filtering(dat)
still_missing$units.outcome<-NA
still_missing$units.exposure<-NA
still_missing$rsq.exposure<-NA
still_missing$effective_n.exposure<-NA
still_missing$rsq.outcome<-NA
still_missing$steiger_dir<-TRUE
still_missing$steiger_pval<-NA
still_missing$r.outcome<-NA
dat<-rbind(dat,still_missing)  
dat$outcome <- "Luminal B-HER2-negative-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

results$bh_p <- p.adjust(results$pval,method = "BH")

fwrite(results,"results/LiverFatMR/Repeating_no_sample_overlap.csv")

sink()
