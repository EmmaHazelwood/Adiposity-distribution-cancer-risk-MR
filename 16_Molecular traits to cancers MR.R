sink("Rscript.txt")
library(dplyr)
library(TwoSampleMR)
library(data.table)
library(LDlinkR)
library(tidyr)
library(gwasvcf)


af_ref<-fread("data/1000GenomesReferenceFiles/EUR.afreq")

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
  data_outcome$id.outcome<-data_outcome$id.outcome[1]
  return(data_outcome)
}


results<-data.frame(matrix(nrow=0,ncol=0))


# ASAT  ------------------------------------------

exposure_dat1 <- extract_instruments(outcomes=c("ieu-b-4870","ieu-b-4871","ieu-b-4865","ieu-b-4869","ieu-b-4868","ieu-b-111","ieu-b-109"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat1$samplesize.exposure[exposure_dat1$id.exposure=="ieu-b-4870"]<-214989
exposure_dat1$samplesize.exposure[exposure_dat1$id.exposure=="ieu-b-4871"]<-185221
exposure_dat1$samplesize.exposure[exposure_dat1$id.exposure=="ieu-b-4865"]<-199569
exposure_dat1$samplesize.exposure[exposure_dat1$id.exposure=="ieu-b-111"]<-441016
exposure_dat1$samplesize.exposure[exposure_dat1$id.exposure=="ieu-b-109"]<-403943
exposure_dat1$samplesize.exposure[exposure_dat1$id.exposure=="ieu-b-4869"]<-180386
exposure_dat1$samplesize.exposure[exposure_dat1$id.exposure=="ieu-b-4868"]<-184205


exposure_dat2 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/IGFBP1/IGFBP_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat2$chr.exposure<-NA
exposure_dat2$pos.exposure<-NA
exposure_dat2$samplesize.exposure<-35560
exposure_dat2$exposure<-"IGFBP1"

exposure_dat3 <- read_exposure_data(
  filename = "data/Sinnott_Armstrong_et_al_2019/IGF_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat3$chr.exposure<-NA
exposure_dat3$pos.exposure<-NA
exposure_dat3$samplesize.exposure<-317114
exposure_dat3$exposure<-"IGF-1"

exposure_dat4 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/ADIPOQ/Adiponectin_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat4$chr.exposure<-NA
exposure_dat4$pos.exposure<-NA
exposure_dat4$samplesize.exposure<-35560
exposure_dat4$exposure<-"Adiponectin"

exposure_dat5 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/PAI1/PAI1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat5$chr.exposure<-NA
exposure_dat5$pos.exposure<-NA
exposure_dat5$samplesize.exposure<-35560
exposure_dat5$exposure<-"PAI-1"

exposure_dat6 <- read_exposure_data(
  filename = "data/FI_GWAS/FI_combined_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat6$chr.exposure<-NA
exposure_dat6$pos.exposure<-NA
exposure_dat6$samplesize.exposure<-98210
exposure_dat6$exposure<-"Fasting insulin"

exposure_dat7 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/CXCL8/CXCL8_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat7$chr.exposure<-NA
exposure_dat7$pos.exposure<-NA
exposure_dat7$samplesize.exposure<-35560
exposure_dat7$exposure<-"CXCL8"

exposure_dat8 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/RETN/Resistin_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat8$chr.exposure<-NA
exposure_dat8$pos.exposure<-NA
exposure_dat8$samplesize.exposure<-35560
exposure_dat8$exposure<-"Resistin"

exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3,exposure_dat4,exposure_dat5,exposure_dat6,exposure_dat7,exposure_dat8)
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

#Endometrial
print("endometrial")

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006464", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-12906
dat$ncontrol.outcome<-108979
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
dat$outcome <- "Endometrial cancer"
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
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-862
dat$ncontrol.outcome<-715717
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
dat$outcome <- "Liver cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Luminal B-like
print("ASAT #luminal b-like")

outfilepath<-"data/Obesity_related_cancer_GWAS/Breast_cancer/annotated/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt"

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
dat$outcome <- "Luminal B-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Ovarian
print("ASAT #ovarian")

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1120", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-25509
dat$ncontrol.outcome<-40941
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
dat$outcome <- "Ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


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
  
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$ncase.outcome<-4112
dat$ncontrol.outcome<-17159
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Esophagus adenocarcinoma"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Endometrioid - ieu-a-1125
print("ASAT endometrioid")

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1125", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-2810
dat$ncontrol.outcome<-40941
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Endometrioid ovarian cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Endometrial
#Endometrioid
print("ASAT endometrioid endometrial")

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006465", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-8758
dat$ncontrol.outcome<-46126
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Endometrioid endometrial cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)




# VAT ---------------------------------------------------------------------



# GFAT --------------------------------------------------------------------

exposure_dat1 <- extract_instruments(outcomes=c("ieu-b-109","ieu-b-4871","ieu-b-4869","ieu-b-4868"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat1$samplesize.exposure[exposure_dat1$exposure=="Bioavailable Testosterone || id:ieu-b-4868"]<-184205
exposure_dat1$samplesize.exposure[exposure_dat1$exposure=="Bioavailable Testosterone || id:ieu-b-4869"]<-180386
exposure_dat1$samplesize.exposure[exposure_dat1$exposure=="HDL cholesterol || id:ieu-b-109"]<-403943
exposure_dat1$samplesize.exposure[exposure_dat1$exposure=="Sex hormone binding globulin (SHBG) || id:ieu-b-4871"]<-185221



exposure_dat2 <- read_exposure_data(
  filename = "data/Sinnott_Armstrong_et_al_2019/IGF_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat2$chr.exposure<-NA
exposure_dat2$pos.exposure<-NA
exposure_dat2$samplesize.exposure<-317114
exposure_dat2$exposure<-"IGF-1"


exposure_dat3 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/ADIPOQ/Adiponectin_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat3$chr.exposure<-NA
exposure_dat3$pos.exposure<-NA
exposure_dat3$samplesize.exposure<-35560
exposure_dat3$exposure<-"Adiponectin"

exposure_dat4 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/IGFBP1/IGFBP_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat4$chr.exposure<-NA
exposure_dat4$pos.exposure<-NA
exposure_dat4$samplesize.exposure<-35560
exposure_dat4$exposure<-"IGFBP1"

exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3,exposure_dat4)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200


#Mucinous - ieu-a-1123
print("GFAT mucinous")

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1123", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1417
dat$ncontrol.outcome<-40941
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Invasive mucinous ovarian cancer"
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
  
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1012
dat$ncontrol.outcome<-259583
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Meningioma"
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
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-133384
dat$ncontrol.outcome<-113789
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Non-endometrioid
print("GFAT non-endometrioid")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006466", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1230
dat$ncontrol.outcome<-35447
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Non-endometrioid endometrial cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Triple negative or basal-like
print("GFAT triple neg")
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Triple negative or basal-like breast cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)






# Pancreas fat ------------------------------------------------------------

exposure_dat <- extract_instruments(outcomes=c("ieu-b-111","ieu-b-4864"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$samplesize.exposure[exposure_dat$exposure=="Total Testosterone || id:ieu-b-4864"]<-199569
exposure_dat$samplesize.exposure[exposure_dat$exposure=="triglycerides || id:ieu-b-111"]<-441016
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Endometrioid - ieu-a-1125
print("pancreas endometrioid")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1125", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-2810
dat$ncontrol.outcome<-40941
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Endometrioid ovarian cancer"
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
  
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-1649
dat$ncontrol.outcome<-727247
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Multiple myeloma"
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
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-14416
dat$ncontrol.outcome<-43099
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Proximal colon cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

# Liver fat ---------------------------------------------------------------

exposure_dat1 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/ADIPOQ/Adiponectin_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat1$chr.exposure<-NA
exposure_dat1$pos.exposure<-NA
exposure_dat1$samplesize.exposure<-35560
exposure_dat1$exposure<-"Adiponectin"

exposure_dat2 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/CXCL8/CXCL8_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat2$chr.exposure<-NA
exposure_dat2$pos.exposure<-NA
exposure_dat2$samplesize.exposure<-35560
exposure_dat2$exposure<-"CXCL8"

exposure_dat3 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/PAI1/PAI1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat3$chr.exposure<-NA
exposure_dat3$pos.exposure<-NA
exposure_dat3$samplesize.exposure<-35560
exposure_dat3$exposure<-"PAI-1"

exposure_dat4 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/FASN/FASN_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat4$chr.exposure<-NA
exposure_dat4$pos.exposure<-NA
exposure_dat4$samplesize.exposure<-35560
exposure_dat4$exposure<-"FASN"

exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3,exposure_dat4)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Liver
print("Liver liver")
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
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$ncase.outcome<-862
dat$ncontrol.outcome<-715717
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
if(nrow(still_missing)>0){
  still_missing$units.outcome<-NA
  still_missing$units.exposure<-NA
  still_missing$rsq.exposure<-NA
  still_missing$effective_n.exposure<-NA
  still_missing$rsq.outcome<-NA
  still_missing$steiger_dir<-TRUE
  still_missing$steiger_pval<-NA
  still_missing$r.outcome<-NA
  dat<-rbind(dat,still_missing)}  
dat$outcome <- "Liver cancer"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


fwrite(results,"results/LiverFatMR/Molecular_traits_cancers_results.csv")

sink()
