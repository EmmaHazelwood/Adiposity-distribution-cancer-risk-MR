sink("Rscript.txt")
library(dplyr)
library(TwoSampleMR)
library(data.table)
library(LDlinkR)
library(tidyr)
library(functions)
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


# Main script -------------------------------------------------------------

results<-data.frame(matrix(nrow=0,ncol=0))

#Exposures
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
  filename = "data/LiverFatMR/Instruments/VAT_instruments.csv",
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
exposure_dat2$exposure<-"VAT"

exposure_dat3 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/GFAT_instruments.csv",
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
exposure_dat3$exposure<-"GFAT"

exposure_dat4 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Liver fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat4$chr.exposure<-NA
exposure_dat4$pos.exposure<-NA
exposure_dat4$samplesize.exposure<-32858
exposure_dat4$exposure<-"Liver fat"

exposure_dat5 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Pancreas fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat5$chr.exposure<-NA
exposure_dat5$pos.exposure<-NA
exposure_dat5$samplesize.exposure<-25617
exposure_dat5$exposure<-"Pancreas fat"

exposure_dat6 <- extract_instruments(outcomes="ieu-b-40", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)

exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3,exposure_dat4,exposure_dat5,exposure_dat6)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Outcomes
#IL-1B
print("IL-1B")
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/IL1B/3037_62_IL1B_IL_1b.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "IL-1B"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#IL-6
print("IL6")
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/IL6/4673_13_IL6_IL_6.txt"
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "IL-6"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#TNF-a
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/TNFa/5936_53_TNF_TNF_a.txt"
print("TNFa")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "TNF-a"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#IFN-a
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/IFNA1/18389_11_IFNA1_IFNA1.txt"
print("IFNa")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "IFN-a"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#IFN-b
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/IFNB1/14127_240_IFNB1_IFN_b.txt"
print("IFNb")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "IFN-b"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Visfatin
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/Visfatin/5011_11_NAMPT_PBEF.txt"
print("visfatin")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "Visfatin"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



#Resistin
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/RETN/3046_31_RETN_resistin.txt"
print("resistin")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "Resistin"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Adiponectin
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/ADIPOQ/3554_24_ADIPOQ_Adiponectin.txt"
print("adpinectin")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "Adiponectin"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



#PAI1
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/PAI1/2925_9_SERPINE1_PAI_1.txt"
print("PAI1")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "PAI1"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#IGF2
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/IGF2/15295_81_IGF2_IGF_II.txt"
print("IGF2")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "IGF2"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



#IGFBP3
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/IGFBP3/2571_12_IGFBP3_IGFBP_3.txt"
print("IGFBP3")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "IGFBP3"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#IGFBP1
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/IGFBP1/2771_35_IGFBP1_IGFBP_1.txt"
print("IGFBP1")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "data/Ferkingstad_et_al_2021_proteins/IGFBP1/2771_35_IGFBP1_IGFBP_1.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "IGFBP1"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#FASN
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/FASN/8403_18_FASN_Fatty_acid_synthase.txt"
print("FASN")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "FASN"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#MCP1
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/MCP1/2578_67_CCL2_MCP_1.txt"
print("MCP1")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "MCP1"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



#CXCL8
outfilepath<-"data/Ferkingstad_et_al_2021_proteins/CXCL8/3447_64_CXCL8_IL_8.txt"
print("CXCL8")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsids", outcome_BETA="Beta", outcome_SE="SE", outcome_P="missing",
                            outcome_EA="effectAllele", outcome_OA="otherAllele",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="ImpMAF",outcome_N="missing",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-35559
dat <- steiger_filtering(dat)
dat$outcome <- "CXCL8"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


#Triglycerides
print("Triglycerides")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-111", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$samplesize.outcome<-441016
dat <- steiger_filtering(dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



#HDL
print("HDL")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-109", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$samplesize.outcome<-403943
dat <- steiger_filtering(dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)





#Leptin
print("lep")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST90007310", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-49909
dat <- steiger_filtering(dat)
dat$outcome <- "Leptin"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)



#Fasting insulin
outfilepath<-"data/FI_GWAS/FI_combined_1000G_density_formatted_21-03-29.txt"
print("Insulin")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = "\t",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "a2",
  other_allele_col = "a1",
  samplesize_col = "n"
)
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  outcome_dat<-proxy_search(data_exposure=exposure_dat, data_outcome=outcome_dat, data_outcome_path=outfilepath, data_reference="data/UK_Biobank_Reference/mergey.bim", data_reference_path="data/UK_Biobank_Reference/mergey",
                            tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                            outcome_sep="\t", outcome_phenotype="exposure", outcome_SNP="rsid", outcome_BETA="beta", outcome_SE="se", outcome_P="missing",
                            outcome_EA="a2", outcome_OA="a1",outcome_CHR="missing", outcome_POS="missing",outcome_EAF="missing",outcome_N="n",outcome_ID="missing")
  outcome_dat$SNP[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_snp.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$effect_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a1.outcome[!is.na(outcome_dat$target_snp.outcome)]
  outcome_dat$other_allele.outcome[!is.na(outcome_dat$target_snp.outcome)]<-outcome_dat$target_a2.outcome[!is.na(outcome_dat$target_snp.outcome)]
  
}
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-98210
dat <- steiger_filtering(dat)
dat$outcome <- "Fasting insulin"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)


# IGF1 and CRP ------------------------------------------------------------

#Get position for all SNPs in exposure
library(biomaRt)
mart<-useMart(biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp",host="grch37.ensembl.org")
ens=getBM(attributes = c('refsnp_id','chrom_start','chr_name'), filters = 'snp_filter', values = exposure_dat$SNP, mart = mart,useCache = FALSE)
exposure_dat<-merge(exposure_dat,ens,by.x="SNP",by.y="refsnp_id",all.x=TRUE)    
snps<-dplyr::select(exposure_dat,SNP,chr_name,chrom_start)
snps$chr_name<-as.numeric(snps$chr_name)
snps$chrom_start<-as.numeric(snps$chrom_start)

#CRP
out_dat<-fread("data/Sinnott_Armstrong_et_al_2019/C_reactive_protein.imp",sep = "\t")
out_dat<-merge(out_dat,snps,by.x=c("#CHROM","POS"),by.y=c("chr_name","chrom_start"))
out_dat <- data.frame(out_dat)
outcome_dat <- TwoSampleMR::format_data(
  out_dat,
  type="outcome",
  snps=exposure_dat$SNP,
  snp_col="SNP",
  beta_col="Effect",
  se_col="StdErr",
  eaf_col="MAF",
  effect_allele_col="ALT",
  other_allele_col="REF"
)

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-98210
dat <- steiger_filtering(dat)
dat$outcome <- "CRP"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#IGF1
out_dat<-fread("data/Sinnott_Armstrong_et_al_2019/IGF_1.imp",sep = "\t")
out_dat<-merge(out_dat,snps,by.x=c("#CHROM","POS"),by.y=c("chr_name","chrom_start"))
out_dat <- data.frame(out_dat)
outcome_dat <- TwoSampleMR::format_data(
  out_dat,
  type="outcome",
  snps=exposure_dat$SNP,
  snp_col="SNP",
  beta_col="Effect",
  se_col="StdErr",
  eaf_col="MAF",
  effect_allele_col="ALT",
  other_allele_col="REF"
)

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-98210
dat <- steiger_filtering(dat)
dat$outcome <- "IGF-1"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

# Sex-specific traits -----------------------------------------------------

#Female exposures

exposure_dat1 <- extract_instruments(outcomes="ieu-a-974", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)

exposure_dat2 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT_Female.uvinput",
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
exposure_dat2$samplesize.exposure<-19872
exposure_dat2$exposure<-"Female ASAT"

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

exposure_dat4 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT_Female.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat4$chr.exposure<-NA
exposure_dat4$pos.exposure<-NA
exposure_dat4$samplesize.exposure<-19872
exposure_dat4$exposure<-"Female GFAT"

exposure_dat5 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Liver fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat5$chr.exposure<-NA
exposure_dat5$pos.exposure<-NA
exposure_dat5$samplesize.exposure<-32858
exposure_dat5$exposure<-"Liver fat"

exposure_dat6 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Pancreas fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat6$chr.exposure<-NA
exposure_dat6$pos.exposure<-NA
exposure_dat6$samplesize.exposure<-25617
exposure_dat6$exposure<-"Pancreas fat"


exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3,exposure_dat4,exposure_dat5,exposure_dat6)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200


#Total testosterone female
print("total t")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-4864", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$samplesize.outcome<-199569
dat <- steiger_filtering(dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Bioavailable testosterone female
print("bio t")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-4869", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$samplesize.outcome<-180386
dat <- steiger_filtering(dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#SHBG female
print("SHBG")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-4870", rsq = 0.80)
dat$samplesize.outcome<-214989
dat <- steiger_filtering(dat)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Male exposures
exposure_dat1 <- extract_instruments(outcomes="ieu-a-785", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)

exposure_dat2 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT_Male.uvinput",
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
exposure_dat2$samplesize.exposure<-19093
exposure_dat2$exposure<-"Male ASAT"

exposure_dat3 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/VAT_Male.uvinput",
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
exposure_dat3$samplesize.exposure<-19093
exposure_dat3$exposure<-"Male VAT"

exposure_dat4 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT_Male.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat4$chr.exposure<-NA
exposure_dat4$pos.exposure<-NA
exposure_dat4$samplesize.exposure<-19093
exposure_dat4$exposure<-"Male GFAT"

exposure_dat5 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Liver fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat5$chr.exposure<-NA
exposure_dat5$pos.exposure<-NA
exposure_dat5$samplesize.exposure<-32858
exposure_dat5$exposure<-"Liver fat"

exposure_dat6 <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Pancreas fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat6$chr.exposure<-NA
exposure_dat6$pos.exposure<-NA
exposure_dat6$samplesize.exposure<-25617
exposure_dat6$exposure<-"Pancreas fat"


exposure_dat<-rbind(exposure_dat1,exposure_dat2,exposure_dat3,exposure_dat4,exposure_dat5,exposure_dat6)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200

#Total testosterone male
print("total t")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-4865", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$samplesize.outcome<-199569
dat <- steiger_filtering(dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#Bioavailable testosterone male
print("bio t")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-4868", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$samplesize.outcome<-184205
dat <- steiger_filtering(dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

#SHBG male
print("SHBG")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-4871", rsq = 0.80)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat$samplesize.outcome<-185221
dat <- steiger_filtering(dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)

fwrite(results,"results/LiverFatMR/Adiposity_measures_molecular_traits_results_new_proxies.csv")

print("fin")

sink()



