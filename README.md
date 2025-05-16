# Adiposity Distribution and Cancer Risk MR Analyses


Code used for all analyses in the paper:

Hazelwood, E., Goudswaard, L. J., Lee, M. A., Vabistsevits, M., Pournaras, D. J., Brenner, H., Buchanan, D. D., Gruber, S. B., Gsur, A., Li, L., Vodickova, L., Grant, R. C., Samadder, N. J., Timpson, N. J., Gunter, M. J., Schuster-Böckler, B., Yarmolinsky, J., Richardson, T. G., Freisling, H., Murphy, N., & Vincent, E. E. (2025). *Adiposity distribution and risks of twelve obesity-related cancers: a Mendelian randomization analysis*. https://doi.org/10.1101/2025.01.10.25320324

---

## Table of Contents

- [Overview of Scripts](#overview-of-scripts)  
- [Step-by-Step Walkthrough](#step-by-step-walkthrough-for-adiposity-traits---cancer-risk-mr)
- [Acknowledgements](#acknowledgements)

---

## Overview of Scripts

All scripts used in the analyses described in the paper can be found in this repository. These include scripts used to:

### Format GWAS data

- [Formatting Neale lab cancer GWAS.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/2_Formatting%20Neale%20lab%20cancer%20GWAS.R)  
- [Extracting GWAS from OpenGWAS for meta-analysis with metal.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/5_Extracting%20GWAS%20from%20OpenGWAS%20for%20meta-analysis%20with%20metal.R)  
- [Getting rsID for IGF1 and CRP summary data.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/14_Getting%20rsID%20for%20IGF1%20and%20CRP%20summary%20data.R)  
- [Formatting sex-specific GWAS.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/18_Formatting%20sex-specific%20GWAS.R)
- [Formatting molecular trait GWAS.r](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/4_Formatting%20Molecular%20trait%20GWAS.R)

### Construct genetic instruments

- [Generating ASAT VAT and GFAT instruments.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/1_Generating%20ASAT%20VAT%20and%20GFAT%20instruments.R)  
- [Making male measures of adiposity instruments.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/3_Making%20male%20measures%20of%20adiposity%20instruments.R)  
- [Formatting adiposity trait genetic instruments.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/10_Formatting%20adiposity%20trait%20genetic%20instruments.R)  
- [Instruments for molecular traits.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/15_Instruments%20for%20molecular%20traits.R)  
- [F stats and r2.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/21_F%20stats%20and%20r2.R)

### Perform the meta-analysis of UK Biobank and FinnGen cancer risk GWAS

- [METAL_meta_analysis_script_1.txt](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/6_METAL_meta_analysis_script_1.txt)  
- [METAL_meta_analysis_script_2.txt](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/7_METAL_meta_analysis_script_2.txt)  
- [Filtering meta-analysis results.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/8_Filtering%20meta-analysis%20results.R)  
- [Formatting all meta-analysed SNPs.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/9_Formatting%20all%20meta-analysed%20SNPs.R)

### Conduct MR analyses

- [Adiposity measures cancer MR.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/11_Adiposity%20measures%20cancer%20MR.R)  
- [Adiposity measures to molecular traits MR.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/12_Adiposity%20measures%20to%20molecular%20traits%20MR.R)  
- [Molecular traits to cancers MR.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/16_Molecular%20traits%20to%20cancers%20MR.R)
- [MR BMI onto all cancers.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/17_MR%20BMI%20onto%20all%20cancers.R)  
- [Sex-specific MRs.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/19_Sex-specific%20MRs.R)  
- [Repeating MRs without sample overlap.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/20_Repeating%20MRs%20without%20sample%20overlap.R)

### Perform multivariable MR

- [Creating files for MVMR](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/22_Creating%20files%20for%20MVMR.R)  
- [MVMR molecular mechanisms](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/23_MVMR%20molecular%20mechanisms.R)  
- [MVMR sensitivity analyses.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/29_MVMR%20sensitivity%20analyses.R)

### Calculate genetic correlations

- [Genetic correlation.txt](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/27_Genetic%20correlation.txt)


### Make figures

- [Figures.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/24_Figures.R)  
- [Power curves.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/25_Power%20curves.R)  
- [Genetic correlations matrix](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/28_Genetic%20correlations%20matrix.R)

### Make tables

- [Tables.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/26_Tables.R)


Scripts are numbered in the order in which they were run. Most were written in R, making extensive use of the [`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/) package.

All summary genetic data generated by the meta-analyses are available via the [GWAS Catalogue](https://www.ebi.ac.uk/gwas/home) (GCST accession numbers: GCST90570370-GCST90570374).

---

## Step-by-Step Walkthrough for Adiposity Traits - Cancer Risk MR

Below is a detailed step-by-step guide for the main analysis: an MR of adiposity traits on cancer risk, as described in the paper. As noted above, see the [`TwoSampleMR` GitHub guide](https://mrcieu.github.io/TwoSampleMR/) for more detail on how to perform MR in R.

### Step 1: Generating Genetic Instruments for Adiposity Traits

Summary genetic data for ASAT, VAT, and GFAT were obtained from:

Agrawal S, Wang M, Klarqvist MDR, Smith K, Shin J, Dashti H, Diamant N, Choi SH, Jurgens SJ, Ellinor PT, Philippakis A, Claussnitzer M, Ng K, Udler MS, Batra P, Khera AV. *Inherited basis of visceral, abdominal subcutaneous and gluteofemoral fat depots*. Nat Commun. 2022 Jun 30;13(1):3771. doi: [10.1038/s41467-022-30931-2](https://doi.org/10.1038/s41467-022-30931-2)

Summary genetic data for liver fat and pancreas fat were obtained from:

Liu Y, Basty N, Whitcher B, Bell JD, Sorokin EP, van Bruggen N, Thomas EL, Cule M. *Genetic architecture of 11 organ traits derived from abdominal MRI using deep learning*. eLife. 2021 Jun 15;10:e65554. doi: [10.7554/eLife.65554](https://doi.org/10.7554/eLife.65554)

While the liver fat and pancreas fat GWAS are available in [OpenGWAS](https://gwas.mrcieu.ac.uk/) and can be read directly into R using the [`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/) package, the ASAT, VAT, and GFAT data required some reformatting.

I downloaded these data from the [Cardiovascular Disease Knowledge Portal](https://cvd.hugeamp.org/).

Then I ran script [1_Generating ASAT VAT and GFAT instruments.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/1_Generating%20ASAT%20VAT%20and%20GFAT%20instruments.R) to reformat them and selected genetic instruments based on genome-wide significance (*P* < 5 × 10<sup>−8</sup>) and LD independence (r² < 0.001), using the `format_data` and `clump_data` functions from [`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/).

For example for ASAT:
``` R
library(TwoSampleMR)
library(data.table)
library(dplyr)

#ASAT
#Read in data
out_dat<-fread("data/Measures_of_adiposity_Agrawal_2022/0321_asat_bgen_stats") 
out_dat<-data.frame(out_dat)

#Put into TwoSampleMR format
exposure_dat<-format_data(
  out_dat,
  type = "exposure",
  header = TRUE,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0"
)

#Save
fwrite(exposure_dat,"data/Measures_of_adiposity_Agrawal_2022/ASAT.allSNPs")

#Construct genetic instruments
#Limit to genome-wide significant SNPs
exposure_dat <- subset(exposure_dat, pval.exposure<5e-8)

#Clump for LD
exposure_dat <- clump_data(exposure_dat,
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p1 = 0.00000005,
                           clump_p2 = 0.00000005,
                           pop = "EUR")

#Save
fwrite(exposure_dat,"data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput")
```

I then applied further formatting and excluded SNPs with F-statistics < 10 in script [10_Formatting adiposity trait genetic instruments.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/10_Formatting%20adiposity%20trait%20genetic%20instruments.R).

For example for ASAT:
``` R
#Read in data
exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
#Set exposure name and sample size
exposure_dat$exposure<-"ASAT"
exposure_dat$samplesize.exposure<-38965

#Calculate F-statistics
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2

#Filter for any with F-stat<10
exposure_dat<-exposure_dat[exposure_dat$F>=10,]

#Save
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))
```

### Step 2: Formatting Cancer GWAS

Minor reformatting of cancer GWAS data was performed using [2_Formatting Neale lab cancer GWAS.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/2_Formatting%20Neale%20lab%20cancer%20GWAS.R) and [5_Extracting GWAS from OpenGWAS for meta-analysis with metal.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/5_Extracting%20GWAS%20from%20OpenGWAS%20for%20meta-analysis%20with%20metal.R).

### Step 3: Meta-Analysis of Cancer GWAS (UK Biobank and FinnGen)

Meta-analysis was conducted using METAL (v2011-03-25). For details, see the [METAL documentation](https://genome.sph.umich.edu/wiki/METAL_Documentation).

Two METAL scripts were used:

- [6_METAL_meta_analysis_script_1.txt](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/6_METAL_meta_analysis_script_1.txt) 

This script loads and meta-analyses the UK Biobank and FinnGen data per cancer, adjusting for population stratification.

Using liver cancer as an example:

``` bash
#Liver cancer
#Apply genomic control correction to input statistics
GENOMICCONTROL ON

#The model will use effect sizes and their standard errors
SCHEME STDERR 

#Allow for tracking allele frequencies to identify allele flips
#Average allele frequencies across studies
AVERAGEFREQ ON

#Use allele frequencies to identify strand flips
MINMAXFREQ ON

#Upload first file - UK Biobank
MARKER SNP
ALLELE A1 A2
FREQ AF1
EFFECT BETA
STDERR SE
PVAL P

#There is no sample size column in our GWAS data, so we are going to specify this for the whole GWAS - effective sample size as calculated in the metal paper
#Sample size is not provided in the data
WEIGHT DONTUSECOLUMN

#Instead use study-wide effective sample size
DEFAULTWEIGHT 855.598587

PROCESS data/Obesity_related_cancer_GWAS/Liver_cancer/155_PheCode.v1.0.fastGWA

#Upload the second file - FinnGen
MARKER rsids
ALLELE alt ref
FREQ af_alt
EFFECT beta
STDERR sebeta
PVAL pval

WEIGHT DONTUSECOLUMN
DEFAULTWEIGHT 2585.545673

PROCESS data/Obesity_related_cancer_GWAS/Liver_cancer/FinnGen/finngen_R8_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC


#Execute meta-analysis
OUTFILE data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS .csv

ANALYSE HETEROGENEITY

#Clear the data and settings to prepare for the next cancer type
CLEAR
```

Effective sample sizes were calculated as recommended in:

Willer CJ, Li Y, Abecasis GR. *METAL: fast and efficient meta-analysis of genomewide association scans*. Bioinformatics. 2010;26(17):2190–2191. doi: [10.1093/bioinformatics/btq340](https://doi.org/10.1093/bioinformatics/btq340)

with the formula:

$$
N_\text{eff} = \frac{4}{\frac{1}{N_\text{cases}} + \frac{1}{N_\text{controls}}}
$$

- [7_METAL_meta_analysis_script_2.txt](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/7_METAL_meta_analysis_script_2.txt)

This script reads in the outputs of the first script, and repeats the genomic control feature to control for population stratification.

``` R
#Liver cancer
GENOMICCONTROL ON
SCHEME STDERR 

AVERAGEFREQ ON
MINMAXFREQ ON

MARKER MarkerName
ALLELE Allele1 Allele2
FREQ Freq1
EFFECT Effect
STDERR StdErr
PVAL P-value

WEIGHT DONTUSECOLUMN
DEFAULTWEIGHT 1

PROCESS data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS1.csv

#Execute meta-analysis
OUTFILE data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS_GC .csv

ANALYSE

CLEAR
```
This two-step process is recommended in the [METAL documentation](https://genome.sph.umich.edu/wiki/METAL_Documentation).

After the meta-analysis, I filtered SNPs with *P*<sub>heterogeneity</sub> < 0.05 and reformatted the data using [8_Filtering meta-analysis results.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/8_Filtering%20meta-analysis%20results.R) and [9_Formatting all meta-analysed SNPs.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/9_Formatting%20all%20meta-analysed%20SNPs.R).

### Step 4: MR of Adiposity Traits to Cancer Risk

MR analyses were conducted in R using [`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/), and a modified version of Dr Matt Lee’s [`proxy_search`](#) function (included in [11_Adiposity measures cancer MR.R](https://github.com/EmmaHazelwood/Adiposity-distribution-cancer-risk-MR/blob/main/11_Adiposity%20measures%20cancer%20MR.R)) to identify proxy SNPs. Steiger filtering was used to remove instruments explaining more outcome than exposure variance.

For example, for an MR of ASAT to oesophageal adenocarcinoma risk:

``` R
library(dplyr)
library(TwoSampleMR)
library(data.table)
library(LDlinkR)
library(tidyr)
library(functions)
library(gwasvcf)

#Read in exposure data
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

#Set chromosome, position, sample size and exposure columns (sample size is needed for Steiger filtering, chromosome and position columns must be present to enable combining with other exposure data, but can be NA)
exposure_dat1$chr.exposure<-NA
exposure_dat1$pos.exposure<-NA
exposure_dat1$samplesize.exposure<-38965
exposure_dat1$exposure<-"ASAT"

#Read in outcome data 
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

#Check if proxies are required, and if yes use proxy_search function
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

#Harmonise exposure and outcome data
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)

#Steiger filtering
#Provide number of cases and controls for cancer data
dat$ncase.outcome<-4112
dat$ncontrol.outcome<-17159
#Calculate prevalence
dat$prevalence.outcome<-dat$ncase.outcome/(dat$ncase.outcome+dat$ncontrol.outcome)
#Use exposure allele frequency if this is not present in outcome data
dat$eaf.outcome[is.na(dat$eaf.outcome)]<-dat$eaf.exposure
#Calculate r2 for outcome
dat$r.outcome <- get_r_from_lor(lor=dat$beta.outcome,
                                af=dat$eaf.outcome,
                                ncase=dat$ncase.outcome,
                                ncontrol=dat$ncontrol.outcome,
                                prevalence=dat$prevalence.outcome)
dat$rsq.outcome <- dat$r.outcome*dat$r.outcome
#Perform the Steiger filtering
dat <- steiger_filtering(dat)

#Perform MR
dat$outcome <- "Esophagus adenocarcinoma"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
```

The same approach was used for all other MR analyses, including those involving molecular traits and sex-specific data.

---

## Acknowledgements

This GitHub walkthrough was created by Emma Hazelwood in May 2025. I would like to acknowledge my co-authors, as well as the funders and participants who made this work possible.
