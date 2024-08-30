library(TwoSampleMR)
library(data.table)
library(dplyr)

#ASAT
out_dat<-fread("data/Measures_of_adiposity_Agrawal_2022/0321_asat_bgen_stats",sep = "\t")
out_dat<-data.frame(out_dat)
exposure_dat<-format_data(
  out_dat,
  type = "exposure",
  header = TRUE,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ"
)
exposure_dat <- subset(exposure_dat, pval.exposure<5e-8)
exposure_dat <- clump_data(exposure_dat,
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p1 = 0.00000005,
                           clump_p2 = 0.00000005,
                           pop = "EUR")
fwrite(exposure_dat,"data/Measures_of_adiposity_Agrawal_2022/ASAT_Male.uvinput")

#VAT
out_dat<-fread("data/Measures_of_adiposity_Agrawal_2022/0321_vat_bgen_stats",sep = "\t")
out_dat<-data.frame(out_dat)
exposure_dat<-format_data(
  out_dat,
  type = "exposure",
  header = TRUE,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ"
)
exposure_dat <- subset(exposure_dat, pval.exposure<5e-8)
exposure_dat <- clump_data(exposure_dat,
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p1 = 0.00000005,
                           clump_p2 = 0.00000005,
                           pop = "EUR")
fwrite(exposure_dat,"data/Measures_of_adiposity_Agrawal_2022/VAT_Male.uvinput")

#GFAT
out_dat<-fread("data/Measures_of_adiposity_Agrawal_2022/0321_gfat_bgen_stats",sep = "\t")
out_dat<-data.frame(out_dat)
exposure_dat<-format_data(
  out_dat,
  type = "exposure",
  header = TRUE,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ"
)
exposure_dat <- subset(exposure_dat, pval.exposure<5e-8)
exposure_dat <- clump_data(exposure_dat,
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p1 = 0.00000005,
                           clump_p2 = 0.00000005,
                           pop = "EUR")
fwrite(exposure_dat,"data/Measures_of_adiposity_Agrawal_2022/GFAT_Male.uvinput")
