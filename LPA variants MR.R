
rm(list=ls())

# #Setting working directory and packages ----------------------------------------------


setwd("C:/Users/BWAD0002/OneDrive - Region Hovedstaden/Projekt remnant cholesterol/Artikel 9 - SNPs in LPA/Working folder")

library(dplyr)
library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(mygene)
library(ggplot2)
library(eoffice)
library(cowplot)
library(forestploter)
library(grid)
library(stringr)
library(officer)
library(data.table)
library(haven)
library(rms)
library(MendelianRandomization)

#Formatting P-values function for forestploter
format_p<-function(plot,pvals,plot_cols){
  
  for (q in 1:length(pvals))
    for(i in 1:length(pvals[[q]])){
      
      p_i<-if (grepl( "e",pvals[[q]][i])){
        bquote(.(paste0(word(pvals[[q]][i], 1, sep="e")," x 10"))^.(paste0("-",sub("^0+", "", word(pvals[[q]][i], 2, sep="-")))))} else {
          str_pad(pvals[[q]][i], width=4, side="right", pad="0")}
      
      plot <- add_text(plot, p_i, part="body",just ="left", row=i,col=plot_cols[q],gp=gpar(fontsize=9))
    }
  return(plot)
}

# Variables of interest
variables<-c("sex","alder","lpanew","kiv2",
             "bmi","current_smoking",
             "nonhdl","nonhdl_mgdl",
             "systolic","dmall_pre",
             "chol_med","hypertens_med")


variables_names<-c("Women","Age, years","Lipoprotein(a), mg/dL","Sum of KIV-2 repeats, number",
                   "Body mass index, kg/m2","Current smoking",
                   "Non-HDL cholesterol, mmol/L","Non-HDL cholesterol, mg/dL",
                   "Systolic blood pressure, mmHg","Diabetes",
                   "Use of statins","Use of antihypertensives")

# Loading data ------------------------------------------------------------

# Variables of interest
variables<-c("sex","alder","lpanew","kiv2",
             "bmi","current_smoking",
             "nonhdl","nonhdl_mgdl",
             "systolic","dmall_pre",
             "chol_med","hypertens_med")


variables_names<-c("Women","Age, years","Lipoprotein(a), mg/dL","Sum of KIV-2 repeats, number",
                   "Body mass index, kg/m2","Current smoking",
                   "Non-HDL cholesterol, mmol/L","Non-HDL cholesterol, mg/dL",
                   "Systolic blood pressure, mmHg","Diabetes",
                   "Use of statins","Use of antihypertensives")


lpa_cgps<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/lpanew_obush.dta")
lpa_cchs<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/lpanew_obus3+4.dta")

kiv_cgps<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/OBUSH_kiv2.dta")
kiv_cchs<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/obus3_4_kiv2.dta")

cgps<-merge(lpa_cgps,kiv_cgps,by="obushnr")
cchs<-merge(lpa_cchs,kiv_cchs,by="lobenr")

snps_cgps_exome<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/LPA_region_CGPS_Exomearray_Imputed_20260204.dta")
snps_cchs_exome<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/LPA_region_CCHS_Exomearray_Imputed_20260204.dta")

snps_cchs_exome<-dplyr::rename(snps_cchs_exome,lobenr=qbus3nr) # changing variable id name

#Adding other necessary variables for cgps
cgps_other<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Dataset_obush.dta")

cgps<-merge(cgps,cgps_other,by="obushnr")

cgps$nonhdl<-cgps$chol-cgps$hdl
cgps$nonhdl_mgdl<-cgps$nonhdl*38.6

cgps$birthdate<-as.numeric(cgps$fodedato)/365.25+1970

cgps<-subset(cgps,duplicate_between_cgps_cchs1_4==0 & duplicate_within_cgps==0 )

cgps$dmall_pre<-with(cgps, ifelse(dm_alldato<=usdato  & dm_all==1 | glucose>11 | v27==1 |m11==1 |m12==1, 1, 0)) 
cgps$current_smoking<-ifelse(cgps$status_smoking==1,1,0)

cgps<-dplyr::rename(cgps,kiv2="KIV2",hypertens_med="m2a",chol_med="m5") #correcting variable name

cgps<-cgps[,c("obushnr",variables)] #Only keeping necessry variables

cgps<-merge(cgps,snps_cgps_exome,by="obushnr")

#Adding other necessary variables for cchs
cchs_other<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/CCHS_database_04112020.dta")
cchs<-merge(cchs_other,cchs,by="lobenr")

cchs$alder<-with(cchs,ifelse(is.na(age3),age4,age3))

cchs$nonhdl3<-cchs$chol3-cchs$hdl3
cchs$nonhdl4<-cchs$chol4-cchs$hdl4

cchs$nonhdl<-with(cchs,ifelse(is.na(nonhdl3),nonhdl4,nonhdl3))

cchs$nonhdl_mgdl<-cchs$nonhdl*38.6

cchs$us3dato[cchs$us3att==0]<-NA
cchs$us34dato<-pmin(cchs$us4dato, cchs$us3dato, na.rm = TRUE)

cchs$dmall_pre<-with(cchs, ifelse(dm_alldato<=us34dato |dm_entry3==1  | dm_entry4==1 & us3att==0 |m_11_01==1 & us3att==0 |m_10_01==1 & us3att==0 , 1, 0))
cchs$dmall_pre[is.na(cchs$dmall_pre)]<-0 

cchs<-dplyr::rename(cchs,hypertens_med="hyp_med34",chol_med="m_04_01",systolic="sbt34")
cchs$chol_med[is.na(cchs$chol_med)]<-0

cchs<-dplyr::rename(cchs,current_smoking="smoke34",bmi="bmi34") #correcting variable name

cchs$usdato<-as.Date(origin="1970-01-01",cchs$us3dato)
cchs$birthdate<-as.numeric(with(cchs,usdato - (alder * 365.25)))/365.25+1970

cchs<-cchs[,c("lobenr",variables)] #Only keeping necessry variables

cchs<-merge(cchs,snps_cchs_exome,by="lobenr")

# Changing column names for CCHS
cchs$obushnr<-cchs$lobenr+9e4 #creating unique number

df<-rbind(cgps,cchs[,colnames(cgps)])

dd<-datadist(df)
options(datadist=dd)

# Formatting genetic variants----------------------------------------------


covariates<-c("sex","alder")
# Candidate predictors: everything except the outcome
predictors <- grep("16",colnames(df),value=T)

# Convert character predictors to factor so lm() handles them
for (nm in predictors) {
  if (is.character(df[[nm]])) df[[nm]] <- as.factor(df[[nm]])
}

#Remove a all columns without genetic variation
df <- df[ , sapply(df, function(col) length(unique(col)) > 1) ]

#Changing names of variants to also contain REF and ALT alleles, and changing values to allele counts

# Helper to extract the first integer-like position from a column name
.extract_pos <- function(x) {
  m <- regexpr("[0-9]+", x)
  ifelse(m > 0, substr(x, m, m + attr(m, "match.length") - 1L), x)
}

# Recalculate to allele counts and rename columns to POS_REF_ALT
new_names <- character(length(predictors))

for (i in seq_along(predictors)) {
  col <- predictors[i]
  
  # Convert to character
  g <- as.character(df[[col]])
  
  # Split genotypes into alleles (supports "A|G" and "A/G")
  alleles <- strsplit(g, "[|/]")
  
  # Pool alleles to determine major/minor
  allele_pool <- unlist(alleles)
  allele_pool <- allele_pool[!is.na(allele_pool) & allele_pool != "" & allele_pool != "."]
  
  # If no alleles present, set whole col to NA and keep original name
  if (length(allele_pool) == 0L) {
    df[[col]] <- NA_integer_
    pos <- .extract_pos(col)
    new_names[i] <- pos  # nothing else we can infer
    next
  }
  
  # Allele frequencies
  tbl <- sort(table(allele_pool), decreasing = TRUE)
  major_allele <- names(tbl)[1]                  # REF = most common
  # ALT = most common non-major (if present), tie-break by alphabetical
  if (length(tbl) >= 2L) {
    non_major <- tbl[-1]
    # break ties lexicographically for stability
    freq_max <- max(non_major)
    alts_max <- names(non_major)[non_major == freq_max]
    minor_allele <- sort(alts_max)[1]
  } else {
    # monomorphic column: choose ALT as NA and we'll drop later if desired
    minor_allele <- NA_character_
  }
  
  # Convert each genotype to allele count (0/1/2 of non-major alleles)
  counts <- vapply(alleles, function(a) {
    if (length(a) == 0 || all(is.na(a))) return(NA_integer_)
    a <- a[!is.na(a) & a != "" & a != "."]
    if (length(a) == 0) return(NA_integer_)
    sum(a != major_allele)
  }, integer(1))
  
  df[[col]] <- as.integer(counts)
  
  # Build new name: <POSITION>_<REF>_<ALT>
  pos <- .extract_pos(col)
  if (is.na(minor_allele) || minor_allele == "") {
    # fall back if no ALT can be determined
    new_names[i] <- pos
  } else {
    new_names[i] <- paste0(pos, "_", major_allele, "_", minor_allele)
  }
}

# Apply new names back to df
names(df)[match(predictors, names(df))] <- new_names

df<-subset(df,!is.na(lpanew)& !is.na(kiv2))
fwrite(df,"L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/df_colnamescorrected.tsv")


# Selecting genetic variants P<5e-6, R2<0.2----------------------------------------------

outcomes <- c("lpanew", "kiv2")

df <- df[complete.cases(df[, c(outcomes, covariates)]), ] #cleaning data frame of missing values

pos_ref_alt_pat <- "^\\d+_[ACGT]+_[ACGT]+$"

kept<-grep(pos_ref_alt_pat, names(df),value=T)

df <- df[, c("obushnr","lpanew","kiv2",covariates,kept), drop = FALSE]

# mew candidate predictors: everything except the outcome
predictors <- setdiff(names(df), c("lpanew","kiv2","obushnr",covariates))

# Prepare results with separate columns per outcome
results <- data.frame(
  predictor = predictors,
  n_1 = NA_integer_, p_1 = NA_real_,
  n_2 = NA_integer_, p_2 = NA_real_,
  p_min = NA_real_,  outcome_min = NA_character_,
  stringsAsFactors = FALSE
)

# Helper to compute the predictor's p-value for one outcome (your original logic)
get_pred_p <- function(df, outcome, pred) {
  # Pairwise complete cases for outcome + predictor + covariates
  dat <- df[, c(outcome, pred, covariates), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  n <- nrow(dat)
  
  if (n == 0L)
    return(list(n = 0L, p = NA_real_))
  
  # Build formula safely
  f <- stats::as.formula(
    paste0("`", outcome, "` ~ `", pred, "` + ",
           paste0("`", covariates, "`", collapse = " + "))
  )
  
  # Fit model
  fit <- try(stats::lm(f, data = dat), silent = TRUE)
  if (inherits(fit, "try-error"))
    return(list(n = n, p = NA_real_))
  
  # ---- FIX 2 IMPLEMENTED HERE ----
  # Extract coefficient p-values (no drop1)
  coefs <- summary(fit)$coefficients
  # match predictor rows (handles numeric and factor terms)
  pred_rows <- grep(paste0("^`?", pred, "`?"), rownames(coefs))
  
  if (length(pred_rows) == 0) {
    p <- NA_real_
  } else {
    # numeric predictor -> 1 row
    # factor predictor -> multiple rows -> take max or min
    p <- min(coefs[pred_rows, 4], na.rm = TRUE)
  }
  # ---------------------------------
  
  list(n = n, p = p)
}

# Main loop over predictors, inner loop over two outcomes
for (i in seq_along(predictors)) {
  pred <- predictors[i]
  
  # Outcome 1
  r1 <- get_pred_p(df, outcomes[1], pred)
  results$n_1[i] <- r1$n
  results$p_1[i] <- r1$p
  
  # Outcome 2
  r2 <- get_pred_p(df, outcomes[2], pred)
  results$n_2[i] <- r2$n
  results$p_2[i] <- r2$p

}


results$p_min<-pmin(results$p_1,results$p_2)

# Decide which predictors to keep (unadjusted p)
results$keep <- !is.na(results$p_min) & results$p_min < 5e-6

# Columns to keep in the filtered data
kept <- results$predictor[results$keep]
kept_cols <- c("obushnr",outcomes, covariates, kept)

# Filtered data frame with only significant predictors
df_filtered <- df[, kept_cols, drop = FALSE]

# Now filter and select only lead variants and remove variants with R2 >0.2

# Extract p-values and ensure names match df_filtered column names
pvals <- results$p_min
names(pvals) <- results$predictor

# Keep only variants present in the genotype matrix
common_vars <- intersect(colnames(df_filtered), names(pvals))
df_geno <- df_filtered[, common_vars, drop = FALSE]
pvals <- pvals[common_vars]

# Order variants by increasing p-value
ordered_vars <- names(sort(pvals, decreasing = FALSE))

### --- COMPUTE CORRELATION MATRIX AND R2 --- ###

# Pairwise correlation of dosage-coded variants (0/1/2)
cor_mat <- cor(df_geno, use = "pairwise.complete.obs")

# Convert to R^2
r2_mat <- cor_mat^2
r2_mat[is.na(r2_mat)] <- 0   # treat NA correlations as 0


### --- LD CLUMPING: PICK LEADS AND REMOVE HIGH-LD NEIGHBOURS --- ###

r2_threshold <- 0.20

lead_variants <- c()
removed <- setNames(rep(FALSE, length(ordered_vars)), ordered_vars)

for (snp in ordered_vars) {
  if (removed[snp]) next
  
  # select this variant as a lead
  lead_variants <- c(lead_variants, snp)
  
  # find all variants in LD (R2 > 0.2) with this variant
  ld_hits <- names(which(r2_mat[snp, ] > r2_threshold))
  
  # mark them as removed
  removed[ld_hits] <- TRUE
}

lead_variants <- unique(lead_variants)


# Optional: summary table of lead variants with p-values
lead_summary <- data.frame(
  variant = lead_variants,
  p_value = pvals[lead_variants]
)[order(pvals[lead_variants]), ]


### --- OUTPUT GENOTYPE MATRIX OF LEAD VARIANTS --- ###

df_leads <- df_filtered[, c("obushnr",outcomes,covariates,lead_variants), drop = FALSE]


### --- OUTPUT OBJECTS --- ###
lead_variants     # list of selected lead SNPs
str(df_leads)          # filtered genotype matrix
lead_summary      # table with p-values



# Selecting lead variants with RSID -----------------------------


library(dplyr)
library(tidyr)
library(stringr)


variant_df <- data.frame(raw = lead_variants) %>%
  separate(raw, into = c("POS", "REF", "ALT"), sep = "_") %>%
  mutate(
    CHR = "6",
    POS = as.integer(POS)
  )

#Removing unknown variant
variant_df<-variant_df[!is.na(variant_df$POS),]

# Fetcing RSIDs
library(httr)
library(jsonlite)
library(purrr)

fetch_rsids_rest <- function(CHR, POS, REF, ALT) {
  url <- paste0(
    "https://grch37.rest.ensembl.org/overlap/region/human/",
    CHR, ":", POS, "-", POS,
    "?feature=variation"
  )
  
  res <- GET(url, add_headers("Content-Type"="application/json"))
  if (status_code(res) != 200) return(NA_character_)
  
  dat <- fromJSON(content(res, "text"))
  if (length(dat) == 0) return(NA_character_)
  
  # Filter to real rsIDs only:
  rs_only <- dat$id[grepl("^rs[0-9]+$", dat$id)]
  
  if (length(rs_only) == 0) return(NA_character_)
  
  rs_only[1]
}


# Example for vector of positions:

rsids <- pmap_chr(
  variant_df[, c("CHR", "POS", "REF", "ALT")],
  fetch_rsids_rest
)
variant_df$RSID<-rsids
variant_df$SNP <- paste(variant_df$POS,
                                variant_df$REF,
                                variant_df$ALT,
                                sep = "_")

lead_summary <- lead_summary %>%
  left_join(
    variant_df[,c("SNP", "RSID")],
    by = c("variant" = "SNP")
  )


#Selecting only the variants with an RSID
lead_summary<-lead_summary[!is.na(lead_summary$RSID),]

#Now getting the beta for the variants for their association with lpa and CNV
variants_to_use <- lead_summary$variant

#Removing variants not associated at bonferroni corrected p value <0.05 in model containing all variants
geno <- df_leads[, c("lpanew","sex","alder",variants_to_use)]
model <- ols(lpanew ~ .+sex+alder, data = geno)
model_lpa<-fastbw(model,rule="p")

#Removing variants not associated at bonferroni corrected p value <0.05 in model containing all variants
geno <- df_leads[, c("kiv2","sex","alder",variants_to_use)]
model <- ols(kiv2 ~ .+sex+alder, data = geno)
model_kiv2<-fastbw(model,rule="p")

significant_variants <- unique(c(model_lpa$names.kept,model_kiv2$names.kept))
significant_variants <- gsub("`", '', significant_variants) #Removing  ´
lead_summary <- lead_summary[ lead_summary$variant %in% significant_variants, ]

#Now getting the beta for the variants for their association with lpa and CNV


variants_to_use <- lead_summary$variant



for(s in outcomes){

results <- sapply(variants_to_use, function(g) {
  
  # g must be a column name
  stopifnot(length(g) == 1)
  
  
  f <- as.formula(paste0(s,"~`",g,"`+",paste0(covariates,collapse="+")))
  
  model <- lm(f,data=df_leads)
  sm    <- summary(model)
  
  beta  <- sm$coefficients[2, 1]   # effect estimate
  se    <- sm$coefficients[2, 2]   # standard error
  pval  <- sm$coefficients[2, 4]   # p‑value
  
  c(beta = beta, se = se, p = pval)
  
})

lead_summary[,paste0("beta_",s)] <- results["beta", ][lead_summary$variant ]
lead_summary[,paste0("se_",s)] <- results["se", ][ lead_summary$variant ]
lead_summary[,paste0("p_",s)]   <- results["p", ][lead_summary$variant ]
}

lead_summary$SNP<-lead_summary$variant
#Adding POS, REF and ALT as separate columns
lead_summary <- lead_summary %>%
  separate(variant, into = c("POS", "REF", "ALT"), sep = "_") %>%
  mutate(
    POS = as.integer(POS),
    REF = toupper(REF),
    ALT = toupper(ALT)
  )


#Saving genetic variants for genetic scores

results_isodom<-data.frame(subset(fread("Results_isodom.txt"),!is.na(beta_adj)))
results_isosum<-data.frame(subset(fread("Results_isosum.txt"),!is.na(beta_adj)))

results_isodom$SNP<- sub("^chr6_", "", results_isodom$SNP) #Changing names to match data frame columns
results_isosum$SNP<- sub("^chr6_", "", results_isosum$SNP)


vars<-c("SNP","RSID","POS","REF","ALT")
#Creating list containing genetic variants for the scores
Genetic_score_variants<-list("LPA"=subset(lead_summary,p_lpanew<5e-6)[,vars],
                       "KIV2"=subset(lead_summary,p_kiv2<5e-6)[,vars],
                       "isodom"=results_isodom[,vars],
                       "isosum"=results_isosum[,vars])


#Replacing genetic variants that are in LD with R quared >99%
Genetic_score_variants$LPA[Genetic_score_variants$LPA$RSID %in% c("rs118039278","rs535797769"),]<-subset(Genetic_score_variants$isodom,RSID %in% c("rs10455872","rs189149041"))
Genetic_score_variants$KIV2[Genetic_score_variants$KIV2$RSID %in% c("rs118039278","rs535797769"),]<-subset(Genetic_score_variants$isodom,RSID %in% c("rs10455872","rs189149041"))

saveRDS(Genetic_score_variants,"Genetic_score_variants.R")


#Saving data frames separately for DST
df_selected_cchs<-df_selected[df_selected$obushnr%in% (cchs$lobenr+9e4) ,]
df_selected_cgps<-df_selected[!df_selected$obushnr%in% (cchs$lobenr+9e4) ,]
df_selected_cchs$lobenr<-df_selected_cchs$obushnr-9e4
df_selected_cchs<-data.frame(df_selected_cchs)[,!colnames(df_selected_cchs)%in%"obushnr"]
df_selected_cchs<-df_selected_cchs[,c("lobenr",colnames(df_selected_cchs)[1:length(colnames(df_selected_cchs))-1])]

fwrite(df_selected_cchs,"L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/Til Danmark statistik/CCHS34_LPA_SNPs.csv")
fwrite(df_selected_cgps,"L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/Til Danmark statistik/CGPS_LPA_SNPs.csv")



# SNP associations for all sNPS in genetic scores-----------------------------


Genetic_score_variants<-readRDS("Genetic_score_variants.R")

all_variants <- do.call(rbind, Genetic_score_variants)
all_variants<-subset(all_variants,!duplicated(all_variants))

##Aligning the SNPs to the data frame (all were already aligned)
df<-fread("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/df_colnamescorrected.tsv")

colnames(df) <- sub("^X", "", colnames(df))
               
cn <- grep("1",colnames(df),value=T)

parts <- do.call(rbind, strsplit(cn, "_"))

geno_map <- data.frame(
  POS = as.integer(parts[, 1]),
  REF = parts[, 2],
  ALT = parts[, 3],
  SNP = cn,
  stringsAsFactors = FALSE
)

setDT(geno_map)
setDT(all_variants)

ead_aligned <- all_variants[
  geno_map,
  on = .(POS, REF, ALT),
  nomatch = 0
]

#Converting to data frame and removing redunant column
ead_aligned<-data.frame(ead_aligned)
ead_aligned<-ead_aligned[, !colnames(ead_aligned) %in% "i.SNP"]
all_variants<-ead_aligned


#Estimating SNP-lpa associations for the isoform SNPs


df$lpanew_log<-log(df$lpanew+1) #Log transformation due to right-skweness, adding +1 since some values were 0
df$kiv2_log<-log(df$kiv2) #same for KIV2 repeats for comparability
outcomes<-c("lpanew","lpanew_log","kiv2","kiv2_log")

for(s in outcomes){
  
  results <- sapply(all_variants$SNP, function(g) {
    
    # g must be a column name
    stopifnot(length(g) == 1)
    
    f <- as.formula(paste0(s,"~`",g,"`+",paste0(covariates,collapse="+")))
    
    model <- lm(f,data=df)
    sm    <- summary(model)
    
    beta  <- sm$coefficients[2, 1]   # effect estimate
    se    <- sm$coefficients[2, 2]   # standard error
    pval  <- sm$coefficients[2, 4]   # p‑value
    
    c(beta = beta, se = se, p = pval)
    
  })
  
  all_variants[,paste0("beta_",s)] <- results["beta", ][all_variants$SNP ]
  all_variants[,paste0("se_",s)] <- results["se", ][ all_variants$SNP ]
  all_variants[,paste0("p_",s)]   <- results["p", ][all_variants$SNP ]
}

#File for one-sample MR
fwrite(all_variants,"Summary_lpa.tsv")

##Calculating correlation between cariants

# Identify SNP columns as those that look like genotype counts
variants_to_use<-all_variants$SNP[all_variants$SNP %in% colnames(df)]
G <- as.matrix(as.data.frame(df)[, variants_to_use])
storage.mode(G) <- "numeric"

# Center (scaling not needed)
G_centered <- scale(G, center = TRUE, scale = FALSE)

# LD correlation matrix
R <- cor(G_centered, use = "pairwise.complete.obs")

#Changing names to RSID
snps<-subset(all_variants,SNP %in%colnames(R))

colnames(R)<-snps$RSID
rownames(R)<-snps$RSID

#Aquared
save(R,file="correlation.Rdata")



# Creating Table 1 ----------------------------------------------------------



df<-fread("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/df_colnamescorrected.tsv")




#Preparation
library(tangram)

# Customizing formatting of tables
EHJ<-nejm
EHJ[['Cell']][["n"]]<-function (n, class = NULL, hdr = FALSE, possible = NULL, ...) 
{
  n <- formatC(n, format="d", big.mark=",")
  cell(n, class = c("cell_n", class), ...)
}
EHJ[['Cell']][['fraction']] <- function(numerator, denominator, format=3, ...)
{ paste0(formatC(numerator, format="d", big.mark=",")," (",render_f(round(100*numerator/denominator, 0), 0),'%)') }
EHJ[['Cell']][['iqr']] <-function (x, format = NA, na.rm = TRUE, names = FALSE, type = 8, 
                                   msd = FALSE, quant = c(0.25, 0.5, 0.75), ...) {
  if (length(quant)%%2 == 0) 
    stop("hmisc_iqr quant argument must be an odd length")
  m <- median(1:length(quant))
  y <- quantile(x, quant, na.rm=TRUE, names, type)
  
  if (median(x, na.rm=T)<10) {
    y<-round(y,1)
  } else {
    y<-round(y,0)
  }
  
  if (is.na(format)) 
    format <- format_guess(y)
  ql <- "-"
  if (sum(!is.na(y)) > 0) {
    ql <- sapply(y, function(x) as.character(x))
    ql <- paste0(ql[m], " (", paste0(ql[1:(m - 1)], 
                                     collapse = "-"), "-", paste0(ql[(m + 
                                                                        1):length(quant)], collapse = "-"), ")")
    if (msd) 
      ql <- paste0(ql, "  ", render_f(mean(x, na.rm = TRUE), 
                                      format), "?", render_f(sd(x, na.rm = TRUE), 
                                                             format))
  }
  cell(ql, ...)
}

EHJ[['Numerical']][['Categorical']] <-function (table, row, column, cell_style, pformat = NULL, msd = FALSE, 
                                                quant = c(0.25, 0.5, 0.75), overall = NULL, test = FALSE, 
                                                useNA = "no", ...) {
  overall_label <- if (is.character(overall)) 
    overall
  else "Overall"
  overall <- column$value != "1" && (isTRUE(overall) || 
                                       is.character(overall))
  datar <- row$data
  datac <- as.categorical(column$data)
  categories <- if (overall) 
    c(levels(datac), overall_label)
  else levels(datac)
  categories <- if (length(categories) == 1) 
    overall_label
  else categories
  format <- ifelse(is.na(row$format), format_guess(datar), 
                   row$format)
  useNA <- useNA == "always" || (sum(is.na(datar)) > 
                                   0 && useNA == "ifany")
  subN <- lapply(levels(datac), FUN = function(cat) {
    cell_style[["n"]](length(datac[datac == cat & !is.na(datac)]), 
                      subcol = cat, hdr = TRUE, possible = length(datac), 
                      ...)
  })
  if (overall) 
    subN[[length(subN) + 1]] <- cell_style[["n"]](sum(!is.na(column$data)), 
                                                  hdr = TRUE, subcol = "Overall", possible = length(column$data), 
                                                  ...)
  stat <- if (length(categories) == 1) {
    ""
  }
  else {
    tst <- suppressWarnings(spearman2(datac, datar, na.action = na.retain))
    cell_style[["fstat"]](f = render_f(tst["F"], 
                                       "%.2f"), df1 = tst["df1"], df2 = tst["df2"], 
                          p = cell_style[["p"]](tst["P"], pformat))
  }
  tbl <- table %>% row_header(derive_label(row))
  if (msd) 
    tbl <- tbl %>% row_header("   Mean?SD")
  if (useNA) 
    tbl <- tbl %>% row_header("   Missing (%)")
  tbl <- if (test) {
    col_header(tbl,  categories, "Test Statistic") %>% 
      col_header("", subN, "")
  }
  else {
    col_header(tbl,  categories) %>% col_header("", 
                                                subN)
  }
  
  tbl <- table_apply(tbl, categories, function(tbl, category) {
    x <- if (category == overall_label) 
      datar
    else datar[datac == category]
    tbl <- tbl %>% add_row(cell_style[["iqr"]](x, 
                                               row$format, subcol = category, na.rm = TRUE, msd = FALSE, 
                                               quant = quant))
    if (msd) {
      tbl <- add_row(tbl, cell(paste0(render_f(mean(x, 
                                                    na.rm = TRUE), row$format), "?", render_f(sd(x, 
                                                                                                 na.rm = TRUE), row$format))))
    }
    if (useNA) {
      tbl <- add_row(tbl, cell_style[["fraction"]](sum(is.na(x)), 
                                                   length(x), format = row$format, subcol = category))
    }
    tbl %>% new_col()
  })
  tbl <- home(tbl) %>% cursor_right(length(categories) + 1)
  if (test) 
    tbl <- add_col(tbl, stat)
  tbl
}

EHJ[['Categorical']][['Categorical']]<-function (table, row, column, cell_style, pformat = NULL, collapse_single = TRUE, 
                                                 overall = NULL, test = FALSE, row_percents = FALSE, useNA = "no", ...) 
{
  grid <- table(as.categorical(row$data), as.categorical(column$data), 
                useNA = useNA)
  if (is.na(colnames(grid)[ncol(grid)])) 
    grid <- grid[, 1:(ncol(grid) - 1)]
  validcol <- which(!apply(grid, 2, FUN = function(x) {
    all(x == 0)
  }))
  validrow <- which(!apply(grid, 1, FUN = function(x) {
    all(x == 0)
  }))
  stat <- if (length(validrow) < 2 || length(validcol) < 2) 
    NA
  else suppressWarnings(chisq.test(grid[validrow, validcol], 
                                   correct = FALSE))
  ncol <- dim(grid)[2]
  nrow <- dim(grid)[1]
  denominators <- if (row_percents) 
    matrix(rep(rowSums(grid), ncol), ncol = ncol, byrow = FALSE)
  else matrix(rep(colSums(grid), nrow), ncol = ncol, byrow = TRUE)
  rowlabels <- rownames(grid)
  subN <- lapply(colnames(grid), FUN = function(cat) cell_style[["n"]](sum(column$data == 
                                                                             cat, na.rm = TRUE), subcol = cat, possible = length(column$data), 
                                                                       hdr = TRUE, ...))
  if (!is.null(overall) && (!is.logical(overall) || overall)) {
    denominators <- cbind(denominators, rep(sum(grid), nrow))
    grid <- cbind(grid, rowSums(grid))
    colnames(grid)[ncol + 1] <- if (is.character(overall)) 
      overall
    else "Overall"
    subN[[ncol + 1]] <- cell_style[["n"]](sum(!is.na(column$data)), 
                                          possible = length(column$data), subcol = "Overall", 
                                          hdr = TRUE, ...)
    ncol <- ncol + 1
  }
  if (collapse_single && dim(grid)[1] <= 2) {
    name <- row$name()
    try({
      l2 <- attr(row$data, "label")
      if (!is.null(l2)) {
        name <- l2
      }
    })
    pos <- dim(grid)[1]
    x <- matrix(grid[pos, ], nrow = 1)
    colnames(x) <- colnames(grid)
    rownames(x) <- paste(name, ":", rownames(grid)[pos])
    grid <- x
    denominators <- matrix(denominators[pos, ], nrow = 1)
    nrow <- 1
  }
  else {
    if (is.na(rownames(grid)[nrow(grid)])) {
      rownames(grid)[nrow(grid)] <- "Missing (%)"
    }
    rownames(grid) <- lapply(rownames(grid), FUN = function(x) paste("  ", 
                                                                     x))
  }
  if (inherits(test, "function")) {
    test_result <- test(row, column, cell_style, ...)
    test <- TRUE
  }
  else if (test) {
    test_result <- if (any(is.na(stat))) 
      cell("NA")
    else cell_style[["chi2"]](render_f(stat$statistic, 
                                       2), stat$parameter, cell_style[["p"]](stat$p.value, 
                                                                             pformat))
  }
  if (test) {
    table <- col_header(table,  colnames(grid), 
                        "Test Statistic")
    table <- col_header(table, subN, "")
  }
  else {
    table <- col_header(table,  colnames(grid))
    table <- col_header(table,  subN)
  }
  if (nrow > 1) 
    table <- row_header(table, derive_label(row, ...))
  for (nm in rownames(grid)) table <- row_header(table, nm)
  for (j in 1:ncol) {
    if (nrow > 1) 
      table <- add_row(table, "")
    format <- if (is.na(row$format) || is.null(row$format)) 
      format_guess(as.vector(grid/denominators))
    else row$format
    for (i in 1:nrow) {
      table <- if (denominators[i, j] == 0) 
        add_row(table, "")
      else add_row(table, cell_style[["fraction"]](grid[i, 
                                                        j], denominators[i, j], format = format, subcol = colnames(grid)[i], 
                                                   subrow = rownames(grid)[j]))
    }
    table <- new_col(table)
  }
  if (test) {
    table <- add_row(table, test_result)
    if (nrow > 1) 
      table <- add_row(table, rep("", nrow))
  }
  table
}

df<-data.frame(df)

#Formatting factors
for ( i in c("sex","chol_med","hypertens_med","current_smoking","dmall_pre")){
  df[,i]<-factor(df[,i])
}


# PRS remnant
all <- tangram(paste0("1~",paste0(variables,collapse="+")),
               data = df, test=F, id="override", transforms = EHJ, msd=FALSE, digits = 1)

#Setting titles and number of individuals
all[1][[1]][[2]]<-cell_header("Copenhagen studies")

#Changing to variable names from row names
all[2][[1]][[1]]<-cell_header("No. of individuals")

for(i in 2+1:length(variables_names)){
  all[i][[1]][[1]]<-cell_header(variables_names[i-2])
}

csv(all, file="Temp.csv")


# Association of PRS with lp(a) and KIV-2 repeats -----------------------------------------------

df<-fread("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/df_colnamescorrected.tsv")
colnames(df) <- sub("^X", "", colnames(df))
df<-as.data.frame(df)

Genetic_score_variants<-readRDS("Genetic_score_variants.R")

#Claculating PRS


#Calculating lpa PRS
variants_to_use <- Genetic_score_variants$LPA$SNP
geno <- df[, variants_to_use]
model_lpa <- lm(df$lpanew ~ ., data = geno)
df$lpa_PRS<-predict(model_lpa,newdata=df)

#Calculating kiv2 PRS
variants_to_use <- Genetic_score_variants$KIV2$SNP
geno <- df[, variants_to_use]
model_kiv2 <- lm(df$kiv2 ~ ., data = geno)
df$kiv2_PRS<-predict(model_kiv2,newdata=df)



#For Claudia:

coefs_lpa <- model_lpa$coefficients[-1]
coefs_kiv2 <- model_kiv2$coefficients[-1]

coefs_lpa <- data.frame(SNP = names(coefs_lpa),beta_adj = unname(coefs_lpa), row.names = NULL)
coefs_kiv2 <- data.frame(SNP = names(coefs_kiv2), beta_adj = unname(coefs_kiv2),row.names = NULL)

coefs_kiv2$SNP <- gsub("[`']", "", coefs_kiv2$SNP)
coefs_lpa$SNP  <- gsub("[`']", "", coefs_lpa$SNP)


Genetic_scores_adj<-list()
Genetic_scores_adj[["KIV2"]] <- merge(
  Genetic_score_variants$KIV2,
  coefs_kiv2[, c("SNP", "beta_adj")],
  by = "SNP"
)
Genetic_scores_adj[["LPA"]] <- merge(
  Genetic_score_variants$LPA,
  coefs_lpa[, c("SNP", "beta_adj")],
  by = "SNP"
)
Genetic_scores_adj

for (nm in names(Genetic_scores_adj)) {
  write.table(
    Genetic_scores_adj[[nm]],
    file = file.path( paste0(nm, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

#Genetic scores for isoforms

list_iso<-list()
for(i in c("isodom","isosum")){
  list_iso[[i]]<-data.frame(subset(fread(paste0("Results_",i,".txt")),!is.na(beta_adj)))
  
  list_iso[[i]]$SNP<- sub("^chr6_", "", list_iso[[i]]$SNP) #Changing names to match data frame columns
  
  snps_in_df <- list_iso[[i]]$SNP[list_iso[[i]]$SNP %in% colnames(df)]
  
  weights <- list_iso[[i]]$beta_adj[match(snps_in_df, list_iso[[i]]$SNP)]
  
  G <- as.matrix(df[, snps_in_df])
  
  df[,paste0("PRS",i)] <- G %*% weights
}



cor(df$lpa_PRS,df$kiv2_PRS)
cor(df$lpa_PRS,df$PRSisosum)
cor(df$lpa_PRS,df$PRSisodom)

#Loading data for German CKD study
dat_gckd<-fread("variables_for_figure_GCKD.txt")
dat_gckd<-dplyr::rename(dat_gckd,lpa_PRS="PRS.LPA",kiv2_PRS="PRS.KIV")

dat_gckd$isodom<-dat_gckd$isodom-9
dat_gckd$isosum<-dat_gckd$isosum-18

fig<-list()
for  (q in c("lpanew","kiv2","isodom","isosum"))
  for (i in c("lpa_PRS","kiv2_PRS")){
    
    col<-if (i=="lpa_PRS"){"darkcyan"} else if (i=="kiv2_PRS") {"darkorchid4"} else {"darkorchid4"}
    
    dat<-data.frame(if(q %in% c("lpanew","kiv2")){df} else {dat_gckd})
    dat[,i]<-scale(dat[,i])
    
    #Calclating meeans in groups of 5
    dat_binned <- dat %>%
      arrange(.data[[i]]) %>%
      mutate(group = ceiling(row_number() / 5)) %>%
      group_by(group) %>%
      summarise(
        x_mean = mean(.data[[i]]),
        y_mean = mean(.data[[q]])
      )
    
    model_linear<-ols(as.formula(paste0(q,"~",i)),data=dat)
    model_nonlinear<-ols(as.formula(paste0(q,"~rcs(",i,",4)")),data=dat)
    
    y<-if(q=="lpanew"){200} else {80} 
    
    gg<-ggplot(data=dat, aes_string(x=i,y=q))+theme_nothing()+
      theme(axis.line = element_line(), axis.title.x = element_text(),axis.title.y = element_text(angle = 90),axis.ticks = element_line(),axis.text = element_text(),axis.ticks.length = unit(.2, "cm"))+
      #geom_point(data=dat_binned,aes(x=x_mean,y=y_mean),color="grey30",alpha=0.2)+
      geom_smooth(color=col,fill=col,method=lm,formula = y ~ rcs(x,4))+
      geom_smooth(color="black",fill="black",linetype="dashed",alpha=0.5,method=lm,formula = y ~ x)+
      xlab(i)+ylab(q)+
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) +
      coord_cartesian(xlim=c(-3,3),ylim=c(0,y))+
      annotate("text",x = -2, y = y * 0.8,
        label =  paste0("Linear R2: ",round(model_linear$stats["R2"] * 100, 0), "%",
            "\nNon-linear R2: ",round(model_nonlinear$stats["R2"] * 100, 0),  "%") , size = 5  )
    
    gg_dens<-ggplot(data=dat) + 
      geom_density(aes_string(x=i),color=col, fill=col, alpha=0.4)+
      theme_nothing()+scale_x_continuous(limits=c(-3,3))
    
    fig[paste(i,q)]<- list(plot_grid(gg_dens,gg,ncol=1,nrow=2,rel_heights= c(0.3,1)))
    
  }


#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 2
formals(plot_grid) <- fargs

figure<-do.call(plot_grid,fig)


gg_dens_PRS_lpanew<-ggplot(data=df) + 
  geom_density(aes_string(x="lpanew"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ coord_flip()+scale_x_continuous(limits=c(0,200))

gg_dens_PRS_kiv2<-ggplot(data=df) + 
  geom_density(aes_string(x="kiv2"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ coord_flip()+scale_x_continuous(limits=c(0,80))

gg_dens_PRS_isodom<-ggplot(data=dat_gckd) + 
  geom_density(aes_string(x="isodom"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ coord_flip()+scale_x_continuous(limits=c(0,80))

gg_dens_PRS_isosum<-ggplot(data=dat_gckd) + 
  geom_density(aes_string(x="isosum"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ coord_flip()+scale_x_continuous(limits=c(0,80))

figure<-plot_grid(figure,plot_grid(gg_dens_PRS_lpanew,gg_dens_PRS_kiv2,gg_dens_PRS_isodom,gg_dens_PRS_isosum,ncol=1),ncol=2, rel_widths = c(1,0.1))

rm(gg_dens_PRS_ldl,gg_dens_PRS_remncholcalc)
figure

# Exporting plots to powerpoint
doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =figure), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



#Now for isoform scores (supplemental figures)

fig<-list()
for  (q in c("lpanew","kiv2","isodom","isosum"))
  for (i in c("PRSisodom","PRSisosum")){
    
    col<-if (i=="PRSisodom"){"darkorchid4"} else if (i=="PRSisosum") {"darkorchid4"} else {"darkorchid4"}
    
    dat<-data.frame(if(q %in% c("lpanew","kiv2")){df} else {dat_gckd})
    dat[,i]<-scale(dat[,i])
    
    #Calculating meeans in groups of 5
    dat_binned <- dat %>%
      arrange(.data[[i]]) %>%
      mutate(group = ceiling(row_number() / 5)) %>%
      group_by(group) %>%
      summarise(
        x_mean = mean(.data[[i]]),
        y_mean = mean(.data[[q]])
      )
    
    model_linear<-ols(as.formula(paste0(q,"~",i)),data=dat)
    model_nonlinear<-ols(as.formula(paste0(q,"~rcs(",i,",4)")),data=dat)
    
    y<-if(q=="lpanew"){200} else {80}
    
    gg<-ggplot(data=dat, aes_string(x=i,y=q))+theme_nothing()+
      theme(axis.line = element_line(), axis.title.x = element_text(),axis.title.y = element_text(angle = 90),axis.ticks = element_line(),axis.text = element_text(),axis.ticks.length = unit(.2, "cm"))+
      #geom_point(data=dat_binned,aes(x=x_mean,y=y_mean),color="grey30",alpha=0.2)+
      geom_smooth(color=col,fill=col,method=lm,formula = y ~ rcs(x,4))+
      geom_smooth(color="black",fill="black",linetype="dashed",alpha=0.5,method=lm,formula = y ~ x)+
      xlab(i)+ylab(q)+
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) +
      coord_cartesian(xlim=c(-3,3),ylim=c(0,y))+
      annotate("text",x = -2, y = y * 0.8,
               label =  paste0("Linear R2: ",round(model_linear$stats["R2"] * 100, 0), "%",
                               "\nNon-linear R2: ",round(model_nonlinear$stats["R2"] * 100, 0),  "%") , size = 5  )
    
    gg_dens<-ggplot(data=dat) + 
      geom_density(aes_string(x=i),color=col, fill=col, alpha=0.4)+
      theme_nothing()+scale_x_continuous(limits=c(-3,3))
    
    fig[paste(i,q)]<- list(plot_grid(gg_dens,gg,ncol=1,nrow=2,rel_heights= c(0.3,1)))
    
  }


#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 2
formals(plot_grid) <- fargs

figure<-do.call(plot_grid,fig)


gg_dens_PRS_lpanew<-ggplot(data=df) + 
  geom_density(aes_string(x="lpanew"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ coord_flip()+scale_x_continuous(limits=c(0,200))

gg_dens_PRS_kiv2<-ggplot(data=df) + 
  geom_density(aes_string(x="kiv2"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ coord_flip()+scale_x_continuous(limits=c(0,80))

gg_dens_PRS_isodom<-ggplot(data=dat_gckd) + 
  geom_density(aes_string(x="isodom"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ coord_flip()+scale_x_continuous(limits=c(0,80))

gg_dens_PRS_isosum<-ggplot(data=dat_gckd) + 
  geom_density(aes_string(x="isosum"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ coord_flip()+scale_x_continuous(limits=c(0,80))

figure<-plot_grid(figure,plot_grid(gg_dens_PRS_lpanew,gg_dens_PRS_kiv2,gg_dens_PRS_isodom,gg_dens_PRS_isosum,ncol=1),ncol=2, rel_widths = c(1,0.1))

rm(gg_dens_PRS_ldl,gg_dens_PRS_remncholcalc)




# Exporting plots to powerpoint
doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =figure), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")




# Observational analysis cardiovasular outcomes -----------------------------------------------


lpa_cgps<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/lpanew_obush.dta")
lpa_cchs<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/lpanew_obus3+4.dta")

kiv_cgps<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/OBUSH_kiv2.dta")
kiv_cchs<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/obus3_4_kiv2.dta")

cgps<-merge(lpa_cgps,kiv_cgps,by="obushnr")
cchs<-merge(lpa_cchs,kiv_cchs,by="lobenr")

snps_cgps_exome<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/LPA_region_CGPS_Exomearray_Imputed_20260204.dta")
snps_cchs_exome<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/LPA_region_CCHS_Exomearray_Imputed_20260204.dta")

snps_cchs_exome<-dplyr::rename(snps_cchs_exome,lobenr=qbus3nr) # changing variable id name

#Adding other necessary variables for cgps
cgps_other<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Dataset_obush.dta")

cgps<-merge(cgps,cgps_other,by="obushnr")
cgps<-merge(cgps,snps_cgps_exome[,1],by="obushnr")

cgps$nonhdl<-cgps$chol-cgps$hdl
cgps$nonhdl_mgdl<-cgps$nonhdl*38.6

cgps$birthdate<-as.numeric(cgps$fodedato)/365.25+1970

cgps<-subset(cgps,duplicate_between_cgps_cchs1_4==0 & duplicate_within_cgps==0 )

cgps$dmall_pre<-with(cgps, ifelse(dm_alldato<=usdato  & dm_all==1 | glucose>11 | v27==1 |m11==1 |m12==1, 1, 0)) 
cgps$current_smoking<-ifelse(cgps$status_smoking==1,1,0)

cgps<-dplyr::rename(cgps,kiv2="KIV2",hypertens_med="m2a",chol_med="m5") #correcting variable name

cgps<-cgps[,c("obushnr",variables)] #Only keeping necessry variables


cchs_other<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/CCHS_database_04112020.dta")
cchs<-merge(cchs_other,cchs,by="lobenr")
cchs<-merge(cchs,snps_cchs_exome[,1],by="lobenr")

cchs$alder<-with(cchs,ifelse(is.na(age3),age4,age3))

cchs$nonhdl3<-cchs$chol3-cchs$hdl3
cchs$nonhdl4<-cchs$chol4-cchs$hdl4

cchs$nonhdl<-with(cchs,ifelse(is.na(nonhdl3),nonhdl4,nonhdl3))

cchs$nonhdl_mgdl<-cchs$nonhdl*38.6

cchs$us3dato[cchs$us3att==0]<-NA
cchs$us34dato<-pmin(cchs$us4dato, cchs$us3dato, na.rm = TRUE)

cchs$dmall_pre<-with(cchs, ifelse(dm_alldato<=us34dato |dm_entry3==1  | dm_entry4==1 & us3att==0 |m_11_01==1 & us3att==0 |m_10_01==1 & us3att==0 , 1, 0))
cchs$dmall_pre[is.na(cchs$dmall_pre)]<-0 

cchs<-dplyr::rename(cchs,hypertens_med="hyp_med34",chol_med="m_04_01",systolic="sbt34")
cchs$chol_med[is.na(cchs$chol_med)]<-0

cchs<-dplyr::rename(cchs,current_smoking="smoke34",bmi="bmi34") #correcting variable name

cchs$usdato<-as.Date(origin="1970-01-01",cchs$us3dato)
cchs$birthdate<-as.numeric(with(cchs,usdato - (alder * 365.25)))/365.25+1970

cchs<-cchs[,c("lobenr",variables)] #Only keeping necessry variables

# Changing column names for CCHS
cchs$obushnr<-cchs$lobenr+9e4 #creating unique number

df<-rbind(cgps,cchs[,colnames(cgps)])


df_cchs_outcomes<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Dataset_obus_outcomes.dta")
df_cgps_outcomes<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Dataset_obush.dta")

df_cchs_outcomes$obushnr<-df_cchs_outcomes$lobenr+9e4
df_cgps_outcomes<-subset(df_cgps_outcomes,duplicate_between_cgps_cchs1_4==0 &duplicate_within_cgps==0 & obushnr<9e4)

df_outcomes<-rbind(df_cchs_outcomes[,c("obushnr","ihd","ami")],df_cgps_outcomes[,c("obushnr","ihd","ami")])

df<-merge(df,df_outcomes[,c("obushnr","ihd","ami")],by="obushnr",all.x=T,all.y=F)
df<-subset(df,!is.na(kiv2) & !is.na(lpanew) & !is.na(ihd))

distrib<- datadist(df)
options(datadist = "distrib")


# Doing Cox regression 
observational_analysis<-data.frame()
effects<-data.frame()
for (i in c("ihd")){
  
  
  model_lpanew<-lrm(as.formula(paste0(i,"~lpanew+alder*sex")),data=df)
  model_kiv2<-lrm(as.formula(paste0(i,"~kiv2+alder*sex")),data=df)
  
  #Number of individuals and events
  observational_analysis[paste(c("LPA","KIV2"),i),c("\nNo. of\ncases")]<-formatC(model_lpanew$freq["1"],format="d", big.mark=",")
  observational_analysis[paste(c("LPA","KIV2"),i),c("\nNo. of\nnon-cases")]<-formatC(model_lpanew$freq["0"],format="d", big.mark=",")

  #Exposure label
  observational_analysis[paste(c("LPA","KIV2"),i),"\n \nExposure"]<-c("Lipoprotein(a)","Sum of KIV-2           ")
  
  observational_analysis[paste(c("LPA","KIV2"),i),"\n \nUnit change"]<-c("Per 50 mg/dL higher levels","Per 15 fewer repeats")
  
  
  #Results on hazard ratio scale
  effect_lpa<-exp(summary(model_lpanew,lpanew=c(0,50),est.all=F)[c("lpanew"),c("Effect","Lower 0.95","Upper 0.95")])
  effect_kiv2<-exp(summary(model_kiv2,kiv2=c(0,-15),est.all=F)[c("kiv2"),c("Effect","Lower 0.95","Upper 0.95")])
  
  #Hazard ratios for forest plot
  effects[paste(c("LPA"),i),c("Effect_uni","Lower 0.95_uni","Upper 0.95_uni")]<-effect_lpa
  
  #Hazard ratios for forest plot
  effects[paste(c("KIV2"),i),c("Effect_uni","Lower 0.95_uni","Upper 0.95_uni")]<-effect_kiv2

  #Hazard ratios for text
  observational_analysis[paste(c("LPA"),i),"Age and sex adjusted odds\nratio (95% CI) for risk of\ncoronary artery disease"]<-  sprintf("%.2f (%.2f-%.2f)",  effect_lpa[1], effect_lpa[2], effect_lpa[3])
  observational_analysis[paste(c("KIV2"),i),"Age and sex adjusted odds\nratio (95% CI) for risk of\ncoronary artery disease"]<-  sprintf("%.2f (%.2f-%.2f)",  effect_kiv2[1], effect_kiv2[2], effect_kiv2[3])
  
  #Creating column with space for forest plot
  observational_analysis[paste(c("LPA","KIV2"),i),"  "]<-paste(rep(" ", 20), collapse = " ")
  
  #P-values
  effects[paste(c("LPA"),i),"\n\nP "]<-as.numeric(pnorm(abs(model_lpanew$coef[2])/sqrt(model_lpanew$var[2,2]),lower.tail=F)*2)
  effects[paste(c("KIV2"),i),"\n\nP "]<-as.numeric(pnorm(abs(model_kiv2$coef[2])/sqrt(model_kiv2$var[2,2]),lower.tail=F)*2)
                                                     

  observational_analysis[paste(c("LPA","KIV2"),i),"\n\nP "]<-"            "
  
  #Mutually adjusted model
  model<-lrm(as.formula(paste0(i,"~lpanew+kiv2+alder*sex")),data=df)

  #Results on hazard ratio scale
  effect<-exp(summary(model,lpanew=c(0,50),kiv2=c(0,-15),est.all=F)[c("lpanew","kiv2"),c("Effect","Lower 0.95","Upper 0.95")])
  
  #Hazard ratios for forest plot
  effects[paste(c("LPA","KIV2"),i),c("Effect","Lower 0.95","Upper 0.95")]<-effect
  
  #Hazard ratios for text
  observational_analysis[paste(c("LPA","KIV2"),i),"Mutually adjusted odds\nratio (95% CI) for risk of\ncoronary artery disease"]<-c(
    sprintf("%.2f (%.2f-%.2f)",  effect["lpanew",1], 
            effect["lpanew",2], effect["lpanew",3]),
    sprintf("%.2f (%.2f-%.2f)",  effect["kiv2",1], 
            effect["kiv2",2], effect["kiv2",3]))
  
  #Creating column with space for forest plot
  observational_analysis[paste(c("LPA","KIV2"),i)," "]<-paste(rep(" ", 20), collapse = " ")
 
  #P-values
  effects[paste(c("LPA","KIV2"),i),"\n\nP"]<-c(pnorm(abs(model$coef[2])/sqrt(model$var[2,2]),lower.tail=F)*2,
                                                    pnorm(abs(model$coef[3])/sqrt(model$var[3,3]),lower.tail=F)*2)
  observational_analysis[paste(c("LPA","KIV2"),i),"\n\nP"]<-"            "
  
  
}

#Setting theme for plot
tm <- forest_theme(base_size = 9,base_family="arial",
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c("white","white","gray90","gray90"))))


p<-forestploter::forest(observational_analysis,
          est = list(effects$Effect_uni,effects$Effect),
          lower =list(effects$`Lower 0.95_uni`,effects$`Lower 0.95`), 
          upper = list(effects$`Upper 0.95_uni`,effects$`Upper 0.95`),
          sizes = 1,
          ref_line = 1,
          ci_column = c(6,9),
          xlab = c("\n     Age and sex adjusted odds ratio (95% CI)\n for risk of coronary artery disease",
                   "\n     Mutually adjusted odds ratio (95% CI)\nfor risk of coronary artery disease"),
          xlim = c(0.9,1.5),
          ticks_at = seq(0.9,1.5,0.1),
          title=paste0("Observational analysis results"),
          theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#Bold font for first column
p <- edit_plot(p, col = c(3), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(1:2), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))


# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df))[c(F,T)],col = c(1:2), which="text" ,label="")

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(3),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(3),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(6,9), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(6,9), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

#Adding P-values with superscript
effects[,ncol(effects)]<-sprintf(c("%3.1g", "%3.2g")[cut(effects[,ncol(effects)], c(-Inf, 0.10, Inf))], effects[,ncol(effects)])
effects[,4]<-sprintf(c("%3.1g", "%3.2g")[cut(effects[,4], c(-Inf, 0.10, Inf))], effects[,4])


p<-format_p(plot=p,pvals=list(effects[,4],effects[,ncol(effects)]),plot_cols = c(7,10))



# Exporting plots to powerpoint
doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")


# Extracting outcome data for CAD-----------------------------------------------------------

all_variants <- read_exposure_data("Summary_lpa.tsv",sep=",",snp_col = "RSID", beta_col = "beta_lpanew",se_col="se_lpanew",pval_col = "p_lpanew",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")

IHD<-list()
IHD[["ukbb"]]<-extract_outcome_data(all_variants$SNP,"ebi-a-GCST90013864",proxies=F)
IHD[["finngen"]]<-extract_outcome_data(all_variants$SNP,"finn-b-I9_IHD",proxies=F)
IHD[["c4d"]]<-extract_outcome_data(all_variants$SNP,"ieu-a-7",proxies=F)
IHD[["mvp"]] <- read_outcome_data(
  snps = all_variants$SNP,
  filename = "C:/Users/BWAD0002/OneDrive - Region Hovedstaden/Projekt remnant cholesterol/GWAS summary/dbGaP/GWAS results/MVP_CAD_Meta_analysis_EUR.withrsID.txt",
  sep = "\t",
  snp_col = "rsmid",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "pvalue",
  units_col = "",
  samplesize_col = "TotalSampleSize",
  chr_col="CHR",
  pos_col="POS_b37")

IHD[["mvp"]] <-dplyr::rename(IHD[["mvp"]],chr=chr.outcome,pos=pos.outcome) #Correcting names


#meta-analyzing data frames

#First harmonising data
for( i in names(IHD)){
  
  IHD[[i]]<-harmonise_data(all_variants,IHD[[i]])[,colnames(IHD[["c4d"]])[!colnames(IHD[["c4d"]]) %in% c("originalname.outcome","outcome.deprecated")]]
}


#Adding number of cases and controls
numbers<-c("cases","controls")
for (i in 1:length(numbers)){
  
  IHD[["ukbb"]][,numbers[i]]<-c(34541, 261984)[i]
  IHD[["finngen"]][,numbers[i]]<-c(30952,187840)[i]
  IHD[["c4d"]][,numbers[i]]<-c(60801, 123504)[i]
  IHD[["mvp"]][,numbers[i]]<-c(95151,197287)[i]

}



library(dplyr)
library(purrr)

### 1. Combine list into long dataset
long <- IHD %>%
  imap(~ .x %>%
         mutate(
           study = .y,
           chr = as.character(chr),          # FIX TYPE CONFLICT
           pos = as.numeric(pos),            # also good practice
           beta.outcome = as.numeric(beta.outcome),
           se.outcome   = as.numeric(se.outcome)
         )) %>%
  bind_rows()


### 3. Fixed‑effect function
fe_meta <- function(b, se){
  w <- 1 / (se^2)
  beta_fe <- sum(w * b) / sum(w)
  se_fe <- sqrt(1 / sum(w))
  z <- beta_fe / se_fe
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  tibble(beta.outcome = beta_fe, se.outcome = se_fe, pval.outcome = p)
}

### 4. Apply meta‑analysis per SNP
meta <- long %>%
  group_by(SNP, effect_allele.outcome, other_allele.outcome) %>%
  summarise(
    cases=sum(cases),
    controls=sum(controls),
    k = n(),
    studies = paste(study, collapse = ", "),
    fe_meta(beta.outcome, se.outcome),
    .groups = "drop"
  )
meta<-data.frame(meta)
meta$outcome<-"Coronary artery disease"
meta$id.outcome<-"CAD"
meta$eaf.outcome<-NA

IHD[["meta"]]<-meta

IHD[["meta"]]<-subset(IHD[["meta"]],!duplicated(SNP))

saveRDS(IHD,"IHD.R")


# Loading prepared data ---------------------------------------------------
IHD<-readRDS("IHD.R")
IHD[["meta"]]<-subset(IHD[["meta"]],!duplicated(SNP))

#Reading exposure data
df_lpa<-read_exposure_data("Summary_lpa.tsv",sep=",",snp_col = "RSID", beta_col = "beta_lpanew",se_col="se_lpanew",pval_col = "p_lpanew",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")
df_lpa_log<-read_exposure_data("Summary_lpa.tsv",sep=",",snp_col = "RSID", beta_col = "beta_lpanew_log",se_col="se_lpanew_log",pval_col = "p_lpanew_log",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")
df_kiv2<-read_exposure_data("Summary_lpa.tsv",sep=",",snp_col = "RSID", beta_col = "beta_kiv2",se_col="se_kiv2",pval_col = "p_kiv2",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")
df_isodom<-read_exposure_data("results_isodom.txt",sep="\t",snp_col = "RSID", beta_col = "beta",se_col="se",pval_col = "p",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")
df_isosum<-read_exposure_data("results_isosum.txt",sep="\t",snp_col = "RSID", beta_col = "beta",se_col="se",pval_col = "p",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")

#Data for the lpa score
df_isodom_lpa<-read_exposure_data("results_isoform_lpa.txt",sep="\t",snp_col = "RSID", beta_col = "beta_isodom",se_col="se_isodom",pval_col = "p_isodom",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")
df_isosum_lpa<-read_exposure_data("results_isoform_lpa.txt",sep="\t",snp_col = "RSID", beta_col = "beta_isosum",se_col="se_isosum",pval_col = "p_isosum",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")

df_isodom<-subset(rbind(df_isodom,df_isodom_lpa),!duplicated(SNP))
df_isosum<-subset(rbind(df_isosum,df_isosum_lpa),!duplicated(SNP))

load(file="correlation.Rdata") #Loading correlation as R
df_lpa$exposure<-"Lipoprotein(a)"
df_lpa_log$exposure<-"Lipoprotein(a)"
df_kiv2$exposure<-"KIV2 repeats"
df_isodom$exposure<-"Dominant apo(a) isoform repeats"
df_isosum$exposure<-"Sum of apo(a) isoform repeats"
df_lpa$id.exposure<-"a"
df_lpa_log$id.exposure<-"a"
df_kiv2$id.exposure<-"b"
df_isodom$id.exposure<-"b"
df_isosum$id.exposure<-"b"

list_comb<-list()
list_comb[["KIV2"]]<-rbind(df_lpa,df_kiv2)
list_comb[["isodom"]]<-rbind(df_lpa,df_isodom)
list_comb[["isosum"]]<-rbind(df_lpa,df_isosum)
list_comb[["KIV2_log"]]<-rbind(df_lpa_log,df_kiv2)
list_comb[["isodom_log"]]<-rbind(df_lpa_log,df_isodom)
list_comb[["isosum_log"]]<-rbind(df_lpa_log,df_isosum)

Genetic_score_variants<-readRDS("Genetic_score_variants.R")

#Univariable MR data frames
list_uni<-list()
list_uni[["LPA"]]<-subset(df_lpa,SNP %in% Genetic_score_variants$LPA$RSID)
list_uni[["LPA_log"]]<-subset(df_lpa_log,SNP %in% Genetic_score_variants$LPA$RSID)
list_uni[["KIV2"]]<-subset(df_kiv2,SNP %in% Genetic_score_variants$KIV2$RSID)
list_uni[["isodom"]]<-subset(df_isodom,SNP %in% Genetic_score_variants$isodom$RSID)
list_uni[["isosum"]]<-subset(df_isosum,SNP %in% Genetic_score_variants$isosum$RSID)

#Filtering SNPs in multivariable score data frames
for( i in names(list_comb)){
  list_comb[[i]]<-list_comb[[i]][list_comb[[i]]$SNP %in% unique(c(list_uni[["LPA"]]$SNP,list_uni[[i]]$SNP)),] #Only keeping variants from scores
}


#genetic variants in each model:
library(dplyr)
library(purrr)
library(tibble)
variants_long_unique <- imap_dfr(
  Genetic_score_variants,
  ~ mutate(.x, exposure = .y)
) %>%
  distinct(SNP, .keep_all = TRUE)



SNP_per_model<-list()
for(s in names(list_comb)[1:3]) 
  for (i in names(IHD)){
mvdat<-mv_harmonise_data(list_comb[[s]],IHD[[i]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data

snps_kept <- rownames(mvdat$exposure_beta)
SNP_per_model[[paste0(s,i)]]<-subset(variants_long_unique,RSID %in%snps_kept)

}


#Calculating conditional F-statistics
saveRDS(SNP_per_model,"SNP_per_model.R")

df<-fread("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/df_colnamescorrected.tsv")
colnames(df) <- sub("^X", "", colnames(df))
df<-as.data.frame(df)

#Conditional F-stastistics in CGPS
library(OneSampleMR)
library(ivreg)

# snp_vars: column names of the instrument dosages used in the MR model
covs <- c("sex", "alder")

fit_condF <- function(d, x1, x2, snp_vars) {
  d <- na.omit(as.data.frame(d)[, c(x1, x2, snp_vars, covs)])
  snp_vars <- snp_vars[sapply(d[snp_vars], var) > 0]   # drop monomorphic SNPs
  
  d$y_dummy <- rnorm(nrow(d)) # placeholder; does not affect the conditional F
  
  z <- paste0("`", snp_vars, "`")                      # backticks protect non-standard names
  f <- as.formula(paste("y_dummy ~", paste(c(x1, x2, covs), collapse = " + "),
                        "|", paste(c(z, covs), collapse = " + ")))
  fsw(ivreg(f, data = d))
}

F_conditional<-list()
for (i in grep("KIV2",names(SNP_per_model),value=T)){
res<-fit_condF(df, "lpanew", "kiv2", SNP_per_model[[i]]$SNP )
F_conditional[[i]]<-res$fswres[,1]
}

#And data from the GCKD study:
F_conditional <- c(F_conditional, list(
  isodomukbb    = c(lpa = 45.72162, isodom = 43.01584),
  isodomfinngen = c(lpa = 37.65903, isodom = 35.65099),
  isodomc4d     = c(lpa = 52.91560, isodom = 49.26431),
  isodommvp     = c(lpa = 42.25804, isodom = 40.04541),
  isodommeta    = c(lpa = 37.43773, isodom = 35.41028),
  
  isosumukbb    = c(lpa = 49.60657, isosum = 61.73595),
  isosumfinngen = c(lpa = 38.53562, isosum = 50.11764),
  isosumc4d     = c(lpa = 55.07026, isosum = 67.23680),
  isosummvp     = c(lpa = 40.45150, isosum = 53.12037),
  isosummeta    = c(lpa = 35.70876, isosum = 46.90351)
))



# Doing multivariable MR by cohort --------------------------------------------------

#Accounting for correlated variants using MendelianRandomization package


df_results<-data.frame()
df_forest<-data.frame()
q<- c("Lipoprotein(a)","Sum of KIV-2")
for(i in names(IHD)){
  
data<-dat_to_MRInput(harmonise_data(list_uni[["LPA"]],IHD[[i]]))
data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
res_lpa<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs,F=data@Fstat)

data<-dat_to_MRInput(harmonise_data(list_uni[["KIV2"]],IHD[[i]]))
data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
res_kiv2<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs,F=data@Fstat)

#Filtering correlation matrix with kept variants
mvdat<-mv_harmonise_data(list_comb[["KIV2"]],IHD[[i]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data

snps_kept <- rownames(mvdat$exposure_beta)
R_filtered<-R[snps_kept, snps_kept, drop = FALSE]

mv_input<-mr_mvinput(
  bx = mvdat$exposure_beta,     # SNP–exposure associations matrix
  by  = mvdat$outcome_beta,     # SNP–outcome associations vector
  bxse = mvdat$exposure_se,
  byse  = mvdat$outcome_se,
  snps   = rownames(mvdat$exposure_beta),
  exposure=mvdat$expname$exposure,
  correlation=R_filtered)

res_mv<-mr_mvivw(mv_input)

df_results[paste(i,q)," "]<-if(grepl("ukbb",i)){"UK Biobank"} else if(
  grepl("fin",i)){"FinnGen"}else if(
    grepl("mvp",i)){"Million Veteran Program"}else if(
      grepl("c4d",i)){"CARDIoGRAMplusC4D"}else if(
        grepl("bbj",i)){"Biobank Japan"}else if(
      grepl("meta",i)){"Meta-analysis"}

df_results[paste(i,q),"\nNo. of\ncases"]<-formatC(IHD[[i]]$cases[1],format="d", big.mark=",")
df_results[paste(i,q),"\nNo. of\nnon-cases"]<-formatC(IHD[[i]]$controls[1],format="d", big.mark=",")

df_results[paste(i,q),"\nGenetic\nscore"]<-q

df_results[paste(i,q),"\n\nUnit change"]<-c("Per 50 mg/dL higher levels","Per 15 fewer repeats")

df_results[paste(i,q),"\nNo. of\nSNPs"]<-formatC(c(res_lpa$nsnp,res_kiv2$nsnp),format="d", big.mark=",")

df_results[paste(i,q),"\nF-\nstatistic"]<-round(c(res_lpa$F,res_kiv2$F))

res_lpa_exp<-exp(c(res_lpa$b,res_lpa$b-(1.959*res_lpa$se),res_lpa$b+(1.959*res_lpa$se))*50) #per 50 mg/dl increment
res_kiv2_exp<-exp(c(res_kiv2$b,res_kiv2$b+(1.959*res_kiv2$se),res_kiv2$b-(1.959*res_kiv2$se))*-15) #per 20 repeats increment


df_results[paste(i,q),
   "Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  c(res_lpa_exp[1],res_kiv2_exp[1]), c(res_lpa_exp[2],res_kiv2_exp[2]), c(res_lpa_exp[3],res_kiv2_exp[3]))
df_forest[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-c(rbind(res_lpa_exp,res_kiv2_exp))

#Creating column with space for forest plot
df_results[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")

df_results[paste(i,q),"\n\nP"]<-"                   "
df_forest[paste(i,q),"\nP"]<-c(res_lpa$pval,res_kiv2$pval)

res_mv_exp<-exp(c(res_mv@Estimate,res_mv@CILower,res_mv@CIUpper)*c(50,-15))

res_mv_exp[c(length(res_mv_exp) - 2, length(res_mv_exp))] <-res_mv_exp[c(length(res_mv_exp), length(res_mv_exp) - 2)] #oorrecting order of confidence intervals


df_results[paste(i,q),"\nConditional\nF-statistic"]<-round(F_conditional[[paste0("KIV2",i)]])

df_results[paste(i,q),
   "Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv_exp[1:2],res_mv_exp[3:4], res_mv_exp[5:6])
df_forest[paste(i,q),c("Effect_mult","Lower 0.95_mult","Upper 0.95_mult")]<-res_mv_exp

#Creating column with space for forest plot
df_results[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")

df_results[paste(i,q),"\n\nP "]<-"                   "
df_forest[paste(i,q),"\nP "]<-res_mv@Pvalue



}

#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c(rep("white",1),rep("gray90",1)))))

xlab_1<-"\n
Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

xlab_2<-"\n
Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

p<-forestploter::forest(df_results,
                        est = list(df_forest$Effect,df_forest$Effect_mult),
                        lower =list(df_forest$`Lower 0.95`,df_forest$`Lower 0.95_mult`), 
                        upper = list(df_forest$`Upper 0.95`,df_forest$`Upper 0.95_mult`), 
                        sizes = 1,
                        ref_line = 1,
                        ci_column = c(9,13),
                        xlab = c(xlab_1,xlab_2),
                        xlim = c(0.9,1.8),ticks_at=c(0.9,seq(1,1.8,0.2)),
                        theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line before last cohort
p <- add_border(p, part = "header", row = 9, where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1,4), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2:3), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))

# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df_results))[c(FALSE,TRUE)],col = c(1:3), which="text" ,label="")

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(4),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(4),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(9,13), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(9,13), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))


#Adding P-values with superscript
df_forest[,4]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,4], c(-Inf, 0.10, Inf))], df_forest[,4])
df_forest[,8]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,8], c(-Inf, 0.10, Inf))], df_forest[,8])

p<-format_p(plot=p,pvals=list(df_forest[,4],df_forest[,8]),plot_cols=c(10,14))

p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[9]]<-unit(4,"cm")
p[["widths"]][[13]]<-unit(4,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# Doing multivariable MR by kiv-2 and isoform --------------------------------------------------

df_results<-data.frame()
df_forest<-data.frame()

for(i in c("KIV2","isodom","isosum")){
  
  data<-dat_to_MRInput(harmonise_data(list_uni[["LPA"]],IHD[["meta"]]))
  data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
  data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
  res_lpa <- data.frame(b = data@Estimate, se = data@StdError, pval = data@Pvalue,
                        nsnp = data@SNPs, F = data@Fstat)
  
  data<-dat_to_MRInput(harmonise_data(list_uni[[i]], IHD[["meta"]]))
  data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
  data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
  res_kiv2 <- data.frame(b = data@Estimate, se = data@StdError, pval = data@Pvalue,
                         nsnp = data@SNPs, F = data@Fstat)
  
  #Filtering correlation matrix with kept variants
  mvdat<-mv_harmonise_data(list_comb[[i]],IHD[["meta"]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data
  
  snps_kept <- rownames(mvdat$exposure_beta)
  
  R_filtered<-R[snps_kept, snps_kept, drop = FALSE]
  mv_input<-mr_mvinput(
    bx = mvdat$exposure_beta,     # SNP–exposure associations matrix
    by  = mvdat$outcome_beta,     # SNP–outcome associations vector
    bxse = mvdat$exposure_se,
    byse  = mvdat$outcome_se,
    snps   = rownames(mvdat$exposure_beta),
    exposure=mvdat$expname$exposure,
    correlation=R_filtered)
  
  res_mv<-mr_mvivw(mv_input)
  
  q<- c("Lipoprotein(a)", if(grepl("KIV",i)){"Sum of KIV-2"} else if(
    grepl("isodom",i)){"Dominant apo(a) isoform"}else if(
      grepl("isosum",i)){"Sum of apo(a) isoforms"})
  
  df_results[paste(i,q),"\nGenetic\nscore"]<-q
  
  df_results[paste(i,q),"\n\nUnit change"]<-c("Per 50 mg/dL higher levels","Per 15 fewer repeats")
  
  df_results[paste(i,q),"\nNo. of\nSNPs"]<-formatC(c(res_lpa$nsnp,res_kiv2$nsnp),format="d", big.mark=",")
  
  df_results[paste(i,q),"\nF-\nstatistic"]<-round(c(res_lpa$F,res_kiv2$F))
  
  res_lpa_exp<-exp(c(res_lpa$b,res_lpa$b-(1.959*res_lpa$se),res_lpa$b+(1.959*res_lpa$se))*50) #per 50 mg/dl increment
  res_kiv2_exp<-exp(c(res_kiv2$b,res_kiv2$b+(1.959*res_kiv2$se),res_kiv2$b-(1.959*res_kiv2$se))*-15) #per 15 repeats increment
  
  df_results[paste(i,q),
             "Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  c(res_lpa_exp[1],res_kiv2_exp[1]), c(res_lpa_exp[2],res_kiv2_exp[2]), c(res_lpa_exp[3],res_kiv2_exp[3]))
  df_forest[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-c(rbind(res_lpa_exp,res_kiv2_exp))
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
  
  df_results[paste(i,q),"\n\nP"]<-"                   "
  df_forest[paste(i,q),"\nP"]<-c(res_lpa$pval,res_kiv2$pval)
  
  res_mv_exp<-exp(c(res_mv@Estimate,res_mv@CILower,res_mv@CIUpper)*c(50,-15))
  
  res_mv_exp[c(length(res_mv_exp) - 2, length(res_mv_exp))] <-res_mv_exp[c(length(res_mv_exp), length(res_mv_exp) - 2)] #oorrecting order of confidence intervals
  
  df_results[paste(i,q),"\nConditional\nF-statistic"]<-round(F_conditional[[paste0(i,"meta")]])
  
  df_results[paste(i,q),
             "Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv_exp[1:2],res_mv_exp[3:4], res_mv_exp[5:6])
  df_forest[paste(i,q),c("Effect_mult","Lower 0.95_mult","Upper 0.95_mult")]<-res_mv_exp
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
  
  df_results[paste(i,q),"\n\nP "]<-"                   "
  df_forest[paste(i,q),"\nP "]<-res_mv@Pvalue
  
}

#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c(rep("white",1),rep("gray90",1)))))

xlab_1<-"\n
Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

xlab_2<-"\n
Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

p<-forestploter::forest(df_results,
                        est = list(df_forest$Effect,df_forest$Effect_mult),
                        lower =list(df_forest$`Lower 0.95`,df_forest$`Lower 0.95_mult`), 
                        upper = list(df_forest$`Upper 0.95`,df_forest$`Upper 0.95_mult`), 
                        sizes = 1,
                        ref_line = 1,
                        ci_column = c(6,10),
                        xlab = c(xlab_1,xlab_2),
                        xlim = c(0.9,1.6),ticks_at=c(0.9,seq(1,1.6,0.2)),
                        title=paste0("Mendelian randomization for coronary artery disease"),
                        theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line between results
p <- add_border(p, part = "header", row = c(3,5), where = "bottom",gp = grid::gpar(lty=1,lwd=1))

#Bold font for first column
p <- edit_plot(p, col = c(1), gp =  grid::gpar(fontface = "bold"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(1),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(1),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(6,10), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(6,10), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

#Adding P-values with superscript
df_forest[,4]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,4], c(-Inf, 0.10, Inf))], df_forest[,4])
df_forest[,8]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,8], c(-Inf, 0.10, Inf))], df_forest[,8])

p<-format_p(plot=p,pvals=list(df_forest[,4],df_forest[,8]),plot_cols=c(7,11))

p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[6]]<-unit(4,"cm")
p[["widths"]][[10]]<-unit(4,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")

# Fraction mediated by kiv-2 and isoform --------------------------------------------------



prop_mediated_ci <- function(beta_total, se_total,
                             beta_direct, se_direct,
                             alpha = 0.05) {
  
  # Proportion NOT mediated (direct effect)
  Q <- beta_direct / beta_total
  
  # Delta-method variance
  var_Q <- ( (beta_direct / beta_total^2)^2 * se_total^2 ) +
    ( (1 / beta_total)^2 * se_direct^2 )
  
  se_Q <- sqrt(var_Q)
  
  # Wald Z
  z_val <- Q / se_Q
  
  # Numerically stable two-sided p-value
  log_p <- log(2) + pnorm(abs(z_val),
                          lower.tail = FALSE,
                          log.p = TRUE)
  p_val <- exp(log_p)
  
  # Confidence interval
  z_alpha <- qnorm(1 - alpha / 2)
  CI <- c(
    lower = Q - z_alpha * se_Q,
    upper = Q + z_alpha * se_Q
  )
  
  # Convert to percent and cap at 100
  CI_pct <- pmin(CI * 100, 100)
  
  list(
    P    = Q * 100,
    SE   = se_Q,
    CI   = CI_pct,
    Z    = z_val,
    p    = p_val,
    logp = log_p
  )
}



#Accounting for correlated variants using MendelianRandomization package

df_results<-data.frame()
df_forest<-data.frame()

for(i in c("KIV2","isodom","isosum")){
  IHD[["meta"]]<-subset(IHD[["meta"]],!duplicated(SNP))
  
  data<-dat_to_MRInput(harmonise_data(list_uni[["LPA"]],IHD[["meta"]]))
  data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
  data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
  res_lpa<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs)
  
  data<-dat_to_MRInput(harmonise_data(list_uni[[i]], IHD[["meta"]]))
  data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
  data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
  res_kiv2<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs)
  
  #Filtering correlation matrix with kept variants
  mvdat<-mv_harmonise_data(list_comb[[i]],IHD[["meta"]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data
  
  snps_kept <- rownames(mvdat$exposure_beta)
  R_filtered<-R[snps_kept, snps_kept, drop = FALSE]
  mv_input<-mr_mvinput(
    bx = mvdat$exposure_beta,     # SNP–exposure associations matrix
    by  = mvdat$outcome_beta,     # SNP–outcome associations vector
    bxse = mvdat$exposure_se,
    byse  = mvdat$outcome_se,
    snps   = rownames(mvdat$exposure_beta),
    exposure=mvdat$expname$exposure,
    correlation=R_filtered)
  
  res_mv<-mr_mvivw(mv_input)
  
  
  # Univariable estimates (log-OR scale)
  beta_lpa_uv <- res_lpa$b
  se_lpa_uv   <- res_lpa$se
  
  beta_kiv_uv <- res_kiv2$b
  se_kiv_uv   <- res_kiv2$se
  
  # Multivariable estimates
  exp_names  <- mv_input@exposure
  beta_mv    <- res_mv@Estimate
  se_mv      <- res_mv@StdError
  
  beta_lpa_mv <- beta_mv[exp_names == "Lipoprotein(a)"]
  se_lpa_mv   <- se_mv[exp_names == "Lipoprotein(a)"]
  
  beta_kiv_mv <- beta_mv[exp_names != "Lipoprotein(a)"]
  se_kiv_mv   <- se_mv[exp_names != "Lipoprotein(a)"]
  
  med_lpa <- prop_mediated_ci(
    beta_total  = beta_lpa_uv,
    se_total    = se_lpa_uv,
    beta_direct = beta_lpa_mv,
    se_direct   = se_lpa_mv
  )
  
  P_lpa  <- med_lpa$P
  CI_lpa <- med_lpa$CI
  p_lpa  <- med_lpa$p
  
  # Fraction mediated by isoform length
  P_med  <- 100 - P_lpa
  cI_med<-100-CI_lpa
  
  q<- c("Lipoprotein(a)")
  
  df_results[paste(i,q),"\n\n"]<-" "
  
  df_results[paste(i,q),"\n\nAnalysis"]<- paste0("Lipoprotein(a) + ",if(grepl("KIV",i)){"Sum of KIV-2"} else if(
    grepl("isodom",i)){"Dominant apo(a) isoform"}else if(
      grepl("isosum",i)){"Sum of apo(a) isoforms"})
  
  df_results[paste(i,q),"\nNo. of\nSNPs"]<-paste0(res_lpa$nsnp," + ",formatC(res_kiv2$nsnp,format="d", big.mark=","))
  
  F_c <- round(F_conditional[[paste0(i,"meta")]])
  df_results[paste(i,q),"\nConditional\nF-statistics"] <- paste0(F_c[1], " + ", F_c[2])
  
  # Fraction NOT mediated (keep exactly as before)
  df_results[paste(i,q),
             "Fraction of coronary artery\ndisease risk explained by\nhigh lipoprotein(a) (95% CI)"] <-
    sprintf("%.f%% (%.f–%.f%%)",
            P_lpa, CI_lpa["lower"], CI_lpa["upper"])
  
  # Fraction mediated (point estimate only)
  df_results[paste(i,q),
             "Fraction of coronary artery\ndisease risk\nexplained by short apo(a) (95% cI)"] <-
    sprintf("%.f%% (%.f–%.f%%)", P_med,cI_med["upper"],cI_med["lower"])
  
  res_lpa<-c(P_lpa,CI_lpa)
  
  df_forest[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-rbind(res_lpa)
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
  

  
  
}


library(gridExtra)
library(gtable)
library(grid)




plot_df <- data.frame(
  label    = rownames(df_results),
  est_not  = df_forest$Effect,          # NOT mediated
  est_med  = 100 - df_forest$Effect,    # mediated
  row      = seq_len(nrow(df_forest))
)


# Optional: grouping for alternating colours
plot_df$group <- ifelse(plot_df$row %% 2 == 1, "odd", "even")

#Function for bar plot
bar_grob <- function(value_not,
                     value_med,
                     xmin = 0,
                     xmax = 100,
                     col_not = "darkcyan",
                     col_med = "darkorchid4") {
  
  # Clamp
  value_not <- max(min(value_not, xmax), xmin)
  value_med <- max(min(value_med, xmax - value_not), xmin)
  
  w_not <- value_not / (xmax - xmin)
  w_med <- value_med / (xmax - xmin)
  
  gTree(children = gList(
    
    ## NOT mediated
    rectGrob(
      x = w_not / 2,
      y = 0.5,
      width  = w_not,
      height = 0.6,
      gp = gpar(fill = col_not, col = NA)
    ),
    
    ## Mediated (stacked to the right)
    rectGrob(
      x = w_not + w_med / 2,
      y = 0.5,
      width  = w_med,
      height = 0.6,
      gp = gpar(fill = col_med, col = NA)
    )
    
  ))
}


#Creating table


tbl <- tableGrob(
  df_results,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 9,
    core = list(
      fg_params = list(
        hjust = 0,        # left-align text
        x = 0.02          # small left padding
      )
    ),
    colhead = list(
      fg_params = list(
        hjust = 0,
        x = 0.02,
        fontface = "bold"
      )
    )
  )
)

##Adding barcols to table
bar_col <- 7 
tbl$widths[[bar_col]] <- unit(4, "cm")

n_bar <- nrow(plot_df)

for (i in seq_len(n_bar)) {
  
  g <- bar_grob(
    value_not = plot_df$est_not[i],
    value_med = plot_df$est_med[i]
  )
  
  idx <- which(
    tbl$layout$t == (i + 1) &
      tbl$layout$l == bar_col &
      grepl("-fg$", tbl$layout$name)
  )
  
  stopifnot(length(idx) == 1)
  
  tbl$grobs[[idx]] <- g
}


# Create x-axis grob (0–100)
x_axis <- xaxisGrob(
  at = seq(0, 1, by = 0.2),
  label = seq(0, 100, by = 20),
  gp = gpar(fontsize = 9)
)



# Add a new row at the bottom
tbl <- gtable_add_rows(
  tbl,
  heights =  unit(0, "mm")
)

# Index of the new row
axis_row <- length(tbl$heights)

# Add x-axis ONLY under the bar column
tbl <- gtable_add_grob(
  tbl,
  grobs = x_axis,
  t = axis_row,
  b = axis_row,   # ✅ THIS is the missing piece
  l = bar_col,
  r = bar_col,
  clip = "off"
)


#Adding veritcal line at 0 

first_data_row <- 2
last_data_row  <- first_data_row + nrow(plot_df) - 1


vline_0 <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(0, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(1, "npc"),
  gp = gpar(
    col = "black",
    lwd = 1,
    lty = c(3, 1)
  )
)



tbl <- gtable_add_grob(
  tbl,
  grobs = vline_0,
  t = first_data_row,
  b = last_data_row,
  l = bar_col,
  r = bar_col,
  clip = "off"
)


#Applying line under header

header_bottom_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(1, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(0, "npc"),
  gp = gpar(col = "black", lwd = 2)
)

tbl <- gtable_add_grob(
  tbl,
  grobs = header_bottom_line,
  t = 1,
  b = 1,
  l = 1,
  r = ncol(tbl),
  clip = "off"
)


#Adding other lines
header_bottom_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(1, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(0, "npc"),
  gp = gpar(col = "black", lwd = 1)
)


for (i in c(4)){
  tbl <- gtable_add_grob(
    tbl,
    grobs = header_bottom_line,
    t = i,
    b = 1,
    l = 1,
    r = ncol(tbl),
    clip = "off"
  )
}




format_p_expr <- function(p, digits = 2) {
  if (is.na(p)) return(expression(NA))
  
  if (p >0.001) {
    return(signif(p,2))
  }
  
  expo <- floor(log10(p))
  base <- signif(p / 10^expo,1)
  base <- round(base, digits)
  
  # Return plotmath expression: base %*% 10^expo
  bquote(.(base)~ "×" ~ 10^.(expo))
}

p_exprs <- lapply(df_forest[,"\nP"], format_p_expr)


p_col <- 8

stopifnot(p_col <= ncol(tbl))

n_p <- length(p_exprs)

for (i in seq_len(n_p)) {
  
  p_grob <- textGrob(
    label = p_exprs[[i]],
    gp = gpar(fontsize = 9),
    hjust = 0,                 # left-aligned
    x = unit(0.02, "npc")         # 🔴 critical for superscripts
  )
  
  idx <- which(
    tbl$layout$t == (i + 1) &  # +1 because row 1 is header
      tbl$layout$l == p_col &
      grepl("-fg$", tbl$layout$name)
  )
  
  stopifnot(length(idx) >= 1)
  
  tbl$grobs[[idx[1]]] <- p_grob
}


xlab_grob <- textGrob(
  "\n\nFraction of coronary artery disease risk\nCyan = directly due to elevated lipoprotein(a)\nPurple = indirectly due to short apo(a)",
  gp = gpar(fontsize = 10),
  just = "center"
)


tbl <- gtable_add_rows(
  tbl,
  heights = unit(22, "mm")   # dedicated label row
)

label_row <- length(tbl$heights)


tbl <- gtable_add_grob(
  tbl,
  grobs = xlab_grob,
  t = label_row,
  b = label_row,
  l = bar_col,
  r = bar_col,
  clip = "off"
)

# Shade background of the middle data row only

middle_row <- 2                    # data row index
table_row  <- middle_row + 1        # +1 because row 1 is header

bg_grob <- rectGrob(
  gp = gpar(fill = "grey90", col = NA)
)

tbl <- gtable_add_grob(
  tbl,
  grobs = bg_grob,
  t = table_row,
  b = table_row,
  l = 1,
  r = ncol(tbl),
  z = -Inf,        # 🔑 ensures background is BEHIND text and bars
  clip = "off"
)

tbl$widths[[1]] <- unit(0, "cm")

grid.newpage()
grid.draw(tbl)

doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(tbl)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# Explaining sorting of SNPs ----------------------------------------------


kiv2_outcome<-read_exposure_data("Summary_lpa.tsv",sep=",",snp_col = "RSID", beta_col = "beta_kiv2",se_col="se_kiv2",pval_col = "p_kiv2",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")

#Filtering correlation matrix with kept variants
kiv2_outcome <- format_data(
 kiv2_outcome,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize")

df_lpa<-read_exposure_data("Summary_lpa.tsv",sep=",",snp_col = "RSID", beta_col = "beta_lpanew",se_col="se_lpanew",pval_col = "p_lpanew",effect_allele_col = "ALT",other_allele_col = "REF", pos_col = "POS")


# # Removing SNPS that are too highly correlated
r2_mat <- R^2
diag(r2_mat) <- 0

threshold <- 0.2
snps <- colnames(r2_mat)

keep <- rep(TRUE, length(snps))
names(keep) <- snps

while (TRUE) {
  sub_r2 <- r2_mat[keep, keep, drop = FALSE]
  
  # Find any pairs violating threshold
  bad_pairs <- which(sub_r2 >= threshold, arr.ind = TRUE)
  
  if (nrow(bad_pairs) == 0) break
  
  # Take the first problematic pair
  i <- bad_pairs[1, 1]
  j <- bad_pairs[1, 2]
  
  snp_i <- rownames(sub_r2)[i]
  snp_j <- colnames(sub_r2)[j]
  
  # Remove one SNP (arbitrarily remove j)
  keep[snp_j] <- FALSE
}

pruned_snps <- names(keep)[keep]

df_lpa<-subset(df_lpa,SNP %in%pruned_snps)

ggplot_data<-harmonise_data(df_lpa,kiv2_outcome,action=1) 
ggplot_data<-subset(ggplot_data,SNP %in%IHD[["meta"]]$SNP) #Removing SNPS not in outcome data.

ggplot_data[ggplot_data$beta.exposure<0,"beta.outcome"]<-(-ggplot_data[ggplot_data$beta.exposure<0,"beta.outcome"])
ggplot_data[ggplot_data$beta.exposure<0,"beta.exposure"]<-(-ggplot_data[ggplot_data$beta.exposure<0,"beta.exposure"]) #Flipping alleles if negative values

ggplot_data<-with(ggplot_data,data.frame(ggplot_data,
           x_min=beta.exposure-se.exposure*1.96,x_max=beta.exposure+se.exposure*1.96,
           y_min=beta.outcome-se.outcome*1.96, y_max=beta.outcome+se.outcome*1.96))

ggplot_data<-data.frame(ggplot_data,ratio=(ggplot_data$beta.outcome)/(ggplot_data$beta.exposure))

i<-3
ggplot_data$Cluster<-factor(Hmisc::cut2(ggplot_data$ratio,g=i),labels=c("Apolipoprotein(a) length dominant","Intermediate","Lipoprotein(a) level dominant"))

p<-ggplot(ggplot_data,aes(x=beta.exposure,y=beta.outcome,xmin=x_min,xmax=x_max,ymin=y_min,ymax=y_max,color=Cluster))+
  geom_point()+geom_errorbar()+geom_errorbarh()+geom_smooth(method = "lm", formula = y ~ x + 0,se=F)+
  scale_color_manual(values=c("darkorchid4", "grey","darkcyan"))+
  theme_classic()+labs(x="Difference in lipoprotein(a) (95% CI), mg/dL",y="Difference in sum of KIV-2 repeats (95% CI), number")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location(newlabel = "") )
nrow(ggplot_data)
print(doc,target="temp.pptx")

# Two-sample MR: lipoprotein(a) by change in KIV2-repeats -----------------


df_results<-data.frame()
df_forest<-data.frame()
for (i in rev(levels(ggplot_data$Cluster))){
  

df<-subset(df_lpa, SNP %in% subset(ggplot_data,Cluster==i)$SNP )
data<-dat_to_MRInput(harmonise_data(df,kiv2_outcome))
data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
res_kiv2<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs)

data<-dat_to_MRInput(harmonise_data(df,IHD[["meta"]]))
data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
res_IHD<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs)

df_results[paste(i),"\n\nGenetic variant cluster"]<-i

effect_kiv2<-c(res_kiv2$b,res_kiv2$b-(1.959*res_kiv2$se),res_kiv2$b+(1.959*res_kiv2$se))*50

df_results[paste(i),"\nNo. of\nSNPs"]<-formatC(res_IHD$nsnp,format="d", big.mark=",")

df_results[paste(i),"\nΔ in sum of KIV-2 repeats (95% CI) per\n50 mg/dL higher lipoprotein(a) "]<-sprintf("%.0f (%.0f to %.0f)",  effect_kiv2[1], effect_kiv2[2], effect_kiv2[3])

res_IHD_exp<-exp(c(res_IHD$b,res_IHD$b-(1.959*res_IHD$se),res_IHD$b+(1.959*res_IHD$se))*50) #per 50 mg/dl increment

df_results[paste(i),"Mendelian randomization odds ratio (95% CI)\nfor coronary artery disease per\n50 mg/dL higher lipoprotein(a)"]<-sprintf("%.2f (%.2f-%.2f)",  res_IHD_exp[1], res_IHD_exp[2], res_IHD_exp[3])

df_forest[paste(i),c("Effect","Lower 0.95","Upper 0.95")]<-res_IHD_exp

#Creating column with space for forest plot
df_results[paste(i),"    "]<-paste(rep(" ", 20), collapse = " ")

df_results[paste(i),"\n\nP"]<-"                  "
df_forest[paste(i),"P"]<-res_IHD$pval

#For calculating P-for trend
df_forest[paste(i),"se"]<-res_IHD$se
df_forest[paste(i),"effect_kiv2"]<-effect_kiv2[1]
df_forest[paste(i),"se_kiv2"]<-res_kiv2$se

}


#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = "white")))

xlab_1<-"\n\nMendelian randomization odds ratio (95% CI)\nfor coronary artery disease per\n50 mg/dL higher lipoprotein(a)"

p<-forestploter::forest(df_results,
                        est = list(df_forest$Effect),
                        lower =list(df_forest$`Lower 0.95`), 
                        upper = list(df_forest$`Upper 0.95`), 
                        sizes = 1,
                        ref_line = 1,
                        ci_column = c(5),
                        xlab = c(xlab_1),
                        xlim = c(1,1.8),ticks_at=c(seq(1,1.8,0.2)),
                        title=paste0("Mendelian randomization for coronary artery disease"),
                        theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))


#Bold font for first column
p <- edit_plot(p, col = c(1), gp =  grid::gpar(fontface = "bold"))
p <- edit_plot(p, col = c(2), which="text" ,hjust = unit(1, "npc"),
               x = unit(0.9, "npc"))
# Deleging repeated text
p <- edit_plot(p, row=c(1:nrow(df_results)-1),col = 7, which="text" ,label="")


colors<-c("darkcyan","grey","darkorchid4")

for (i in 1:length(colors)){
  #For text
p <- edit_plot(p, col = c(1) ,row=seq(i,nrow(df_results),3) , which="text",gp = gpar(col = colors[i]))
#For confidence interval
p<-edit_plot(p, col = c(5), row = seq(i,nrow(df_results),3), which = "ci", gp = gpar(col = colors[i], fill = colors[i]))

}


#Adding P-values with superscript
df_forest[,4]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,4], c(-Inf, 0.10, Inf))], df_forest[,4])

p<-format_p(plot=p,pvals=list(df_forest[,4]),plot_cols=c(6))

p[["widths"]][[2]]<-unit(6,"cm")
p[["widths"]][[5]]<-unit(4,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")


# Calculating instrument strength for two-sample MR -----------------------


#Calculating F-statistics

mvdat<-mv_harmonise_data(clump_data(list_comb[["isosum"]],clump_r2=0.01),IHD[["meta"]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data

#Calculating conditional F-statistics
library(MVMR)


mv_input <- format_mvmr(
  BXGs   = mvdat$exposure_beta,
  BYG    = mvdat$outcome_beta,
  seBXGs = mvdat$exposure_se,
  seBYG  = mvdat$outcome_se,
  RSID   =  rownames(mvdat$exposure_beta)
)


betaX <- as.numeric(mvdat$exposure_beta[, 2])
seX   <- as.numeric(mvdat$exposure_se[, 2])

LD<-R[rownames(mvdat$exposure_beta),rownames(mvdat$exposure_beta)]


D  <- diag(seX)
Sigma <- D %*% LD %*% D

K <- length(betaX)

F_ld <- as.numeric(
  t(betaX) %*% solve(Sigma) %*% betaX / K
)


F_naive <- mean( (betaX / seX)^2 )

c(F_naive = F_naive,
  F_LD_adjusted = F_ld)

mv_input


# Calculate conditional F-statistics
strength_results <- strength_mvmr(mv_input,gencov=phenocov_mvmr(Pcov=matrix(c(1,0.47,0.47,1),nrow=2,ncol=2),seBXGs=mv_input[,c("sebetaX1","sebetaX2")]))



#Instrument strength for low LD sensitivity analysis 

q<-"KIV2"

mvdat<-mv_harmonise_data(list_comb[[q]],IHD[["meta"]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data

LD<-R[rownames(mvdat$exposure_beta),rownames(mvdat$exposure_beta)]



# # Removing SNPS that are too highly correlated
r2_mat <- LD^2
diag(r2_mat) <- 0

threshold <- 0.01
snps <- colnames(r2_mat)

keep <- rep(TRUE, length(snps))
names(keep) <- snps

while (TRUE) {
  sub_r2 <- r2_mat[keep, keep, drop = FALSE]
  
  # Find any pairs violating threshold
  bad_pairs <- which(sub_r2 >= threshold, arr.ind = TRUE)
  
  if (nrow(bad_pairs) == 0) break
  
  # Take the first problematic pair
  i <- bad_pairs[1, 1]
  j <- bad_pairs[1, 2]
  
  snp_i <- rownames(sub_r2)[i]
  snp_j <- colnames(sub_r2)[j]
  
  # Remove one SNP (arbitrarily remove j)
  keep[snp_j] <- FALSE
}

pruned_snps <- names(keep)[keep]


df_lpa<-subset(list_comb[[q]],SNP %in%pruned_snps)


mvdat<-mv_harmonise_data(df_lpa,IHD[["meta"]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data



mv_input <- format_mvmr(
  BXGs   = mvdat$exposure_beta,
  BYG    = mvdat$outcome_beta,
  seBXGs = mvdat$exposure_se,
  seBYG  = mvdat$outcome_se,
  RSID   =  rownames(mvdat$exposure_beta)
)

strength_results <- strength_mvmr(mv_input)



# Sensitivity analysis: restricted to variants in low LD--------------------------------------------------


r2_prune <- 0.01

ld_prune_snps <- function(snps, R, r2_threshold = 0.01) {
  
  R_sub <- R[snps, snps, drop = FALSE]
  r2 <- R_sub^2
  
  keep <- rep(TRUE, length(snps))
  names(keep) <- snps
  
  for (i in seq_along(snps)) {
    if (!keep[i]) next
    high_ld <- which(r2[i, ] > r2_threshold)
    high_ld <- setdiff(high_ld, i)
    keep[high_ld] <- FALSE
  }
  
  return(names(keep)[keep])
}

# Helper: prune a harmonised univariable data set, strongest variants first
prune_uni <- function(h, R, r2_threshold = r2_prune) {
  h <- h[h$mr_keep, ]
  h <- h[!duplicated(h$SNP), ]
  snps_ordered <- h$SNP[order(h$pval.exposure)]
  keep <- ld_prune_snps(snps_ordered, R, r2_threshold)
  h[h$SNP %in% keep, ]
}

# Helper: prune a harmonised multivariable data set, ordered by smallest p across exposures
prune_mv <- function(mvdat, R, r2_threshold = r2_prune) {
  snps  <- rownames(mvdat$exposure_beta)
  min_p <- apply(mvdat$exposure_pval, 1, min)
  keep  <- ld_prune_snps(snps[order(min_p)], R, r2_threshold)
  idx   <- match(keep, snps)                     # positions of kept SNPs
  
  mvdat$exposure_beta <- mvdat$exposure_beta[idx, , drop = FALSE]
  mvdat$exposure_se   <- mvdat$exposure_se[idx, , drop = FALSE]
  mvdat$exposure_pval <- mvdat$exposure_pval[idx, , drop = FALSE]
  mvdat$outcome_beta  <- mvdat$outcome_beta[idx]
  mvdat$outcome_se    <- mvdat$outcome_se[idx]
  mvdat$outcome_pval  <- mvdat$outcome_pval[idx]
  mvdat
}

df_results<-data.frame()
df_forest<-data.frame()

for(i in c("KIV2","isodom","isosum")){
  
  
  dat_lpa <- harmonise_data(list_uni[["LPA"]], IHD[["meta"]])
  
  snps_lpa_pruned <- prune_uni(
    h = dat_lpa,
    R = R)$SNP
  
  dat_lpa <- dat_lpa[dat_lpa$SNP %in% snps_lpa_pruned, ]
  
  data <- dat_to_MRInput(dat_lpa)
  data[[1]]@correlation <- R[snps_lpa_pruned, snps_lpa_pruned]
  
  data <- MendelianRandomization::mr_ivw(data[[1]], correl = TRUE)
  
  res_lpa <- data.frame(
    b = data@Estimate,
    se = data@StdError,
    pval = data@Pvalue,
    nsnp = length(snps_lpa_pruned)
  )
  
  dat_i <- harmonise_data(list_uni[[i]], IHD[["meta"]])
  
  snps_i_pruned <- prune_uni(
   dat_i,
    R = R
  )$SNP
  
  dat_i <- dat_i[dat_i$SNP %in% snps_i_pruned, ]
  
  data <- dat_to_MRInput(dat_i)
  data[[1]]@correlation <- R[snps_i_pruned, snps_i_pruned]
  
  data <- MendelianRandomization::mr_ivw(data[[1]], correl = TRUE)
  
  res_kiv2 <- data.frame(
    b = data@Estimate,
    se = data@StdError,
    pval = data@Pvalue,
    nsnp = length(snps_i_pruned)
  )
  #Filtering correlation matrix with kept variants
  mvdat<-mv_harmonise_data(list_comb[[i]],IHD[["meta"]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data
  
  
  mvdat <- prune_mv( mvdat,
    R = R
  )
  
  snps_mv_pruned <- rownames(mvdat$exposure_beta)
  
  R_filtered <- R[snps_mv_pruned, snps_mv_pruned, drop = FALSE]
  
  mv_input <- mr_mvinput(
    bx        = mvdat$exposure_beta,
    by        = mvdat$outcome_beta,
    bxse      = mvdat$exposure_se,
    byse      = mvdat$outcome_se,
    snps      = snps_mv_pruned,
    exposure  = mvdat$expname$exposure,
    correlation = R_filtered
  )
  
  res_mv <- mr_mvivw(mv_input)
  
  q<- c("Lipoprotein(a)", if(grepl("KIV",i)){"Sum of KIV-2"} else if(
    grepl("isodom",i)){"Dominant apo(a) isoform"}else if(
      grepl("isosum",i)){"Sum of apo(a) isoforms"})
  
  df_results[paste(i,q),"\nGenetic\nscore"]<-q
  
  df_results[paste(i,q),"\n\nUnit change"]<-c("Per 50 mg/dL higher levels","Per 15 fewer repeats")
  
  df_results[paste(i,q),"\nNo. of\nSNPs"]<-formatC(c(res_lpa$nsnp,res_kiv2$nsnp),format="d", big.mark=",")
  
  res_lpa_exp<-exp(c(res_lpa$b,res_lpa$b-(1.959*res_lpa$se),res_lpa$b+(1.959*res_lpa$se))*50) #per 50 mg/dl increment
  res_kiv2_exp<-exp(c(res_kiv2$b,res_kiv2$b+(1.959*res_kiv2$se),res_kiv2$b-(1.959*res_kiv2$se))*-15) #per 20 repeats increment
  
  df_results[paste(i,q),
             "Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  c(res_lpa_exp[1],res_kiv2_exp[1]), c(res_lpa_exp[2],res_kiv2_exp[2]), c(res_lpa_exp[3],res_kiv2_exp[3]))
  df_forest[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-c(rbind(res_lpa_exp,res_kiv2_exp))
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
  
  df_results[paste(i,q),"\n\nP"]<-"                   "
  df_forest[paste(i,q),"\nP"]<-c(res_lpa$pval,res_kiv2$pval)
  
  res_mv_exp<-exp(c(res_mv@Estimate,res_mv@CILower,res_mv@CIUpper)*c(50,-15))
  
  res_mv_exp[c(length(res_mv_exp) - 2, length(res_mv_exp))] <-res_mv_exp[c(length(res_mv_exp), length(res_mv_exp) - 2)] #oorrecting order of confidence intervals
  
  
  df_results[paste(i,q),
             "Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv_exp[1:2],res_mv_exp[3:4], res_mv_exp[5:6])
  df_forest[paste(i,q),c("Effect_mult","Lower 0.95_mult","Upper 0.95_mult")]<-res_mv_exp
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
  
  df_results[paste(i,q),"\n\nP "]<-"                   "
  df_forest[paste(i,q),"\nP "]<-res_mv@Pvalue
  
  
  
}

#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c(rep("white",1),rep("gray90",1)))))

xlab_1<-"\n
Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

xlab_2<-"\n
Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

p<-forestploter::forest(df_results,
                        est = list(df_forest$Effect,df_forest$Effect_mult),
                        lower =list(df_forest$`Lower 0.95`,df_forest$`Lower 0.95_mult`), 
                        upper = list(df_forest$`Upper 0.95`,df_forest$`Upper 0.95_mult`), 
                        sizes = 1,
                        ref_line = 1,
                        ci_column = c(5,8),
                        xlab = c(xlab_1,xlab_2),
                        xlim = c(0.9,1.8),ticks_at=c(0.9,seq(1,1.8,0.2)),
                        title=paste0("Mendelian randomization for coronary artery disease"),
                        theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line between results
p <- add_border(p, part = "header", row = c(3,5), where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1), gp =  grid::gpar(fontface = "bold"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(1),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(1),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(5,8), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(5,8), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))


#Adding P-values with superscript
df_forest[,4]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,4], c(-Inf, 0.10, Inf))], df_forest[,4])
df_forest[,8]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,8], c(-Inf, 0.10, Inf))], df_forest[,8])

p<-format_p(plot=p,pvals=list(df_forest[,4],df_forest[,8]),plot_cols=c(6,9))

p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[5]]<-unit(4,"cm")
p[["widths"]][[8]]<-unit(4,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")




# Sensitivity analysis: two-sample MR with log-transformed values ---------


df_results<-data.frame()
df_forest<-data.frame()

for(i in c("KIV2","isodom","isosum")){
  
  
  data<-dat_to_MRInput(harmonise_data(list_uni[["LPA_log"]],IHD[["meta"]]))
  data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
  data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
  res_lpa<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs)
  
  data<-dat_to_MRInput(harmonise_data(list_uni[[i]], IHD[["meta"]]))
  data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
  data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
  res_kiv2<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs)

  #Filtering correlation matrix with kept variants
  mvdat<-mv_harmonise_data(list_comb[[paste0(i,"_log")]],IHD[["meta"]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data
  
  snps_kept <- rownames(mvdat$exposure_beta)
  R_filtered<-R[snps_kept, snps_kept, drop = FALSE]
  mv_input<-mr_mvinput(
    bx = mvdat$exposure_beta,     # SNP–exposure associations matrix
    by  = mvdat$outcome_beta,     # SNP–outcome associations vector
    bxse = mvdat$exposure_se,
    byse  = mvdat$outcome_se,
    snps   = rownames(mvdat$exposure_beta),
    exposure=mvdat$expname$exposure,
    correlation=R_filtered)
  
  res_mv<-mr_mvivw(mv_input)
  
  q<- c("Lipoprotein(a)", if(grepl("KIV",i)){"Sum of KIV-2"} else if(
    grepl("isodom",i)){"Dominant apo(a) isoform"}else if(
      grepl("isosum",i)){"Sum of apo(a) isoforms"})
  
  df_results[paste(i,q),"\nGenetic\nscore"]<-q
  
  df_results[paste(i,q),"\n\nUnit change"]<-c("Per doubling in levels","Per 15 fewer repeats")
  
  df_results[paste(i,q),"\nNo. of\nSNPs"]<-formatC(c(res_lpa$nsnp,res_kiv2$nsnp),format="d", big.mark=",")
  
  res_lpa_exp<-exp(c(res_lpa$b,res_lpa$b-(1.959*res_lpa$se),res_lpa$b+(1.959*res_lpa$se))*log(2)) #per 50 mg/dl increment
  res_kiv2_exp<-exp(c(res_kiv2$b,res_kiv2$b+(1.959*res_kiv2$se),res_kiv2$b-(1.959*res_kiv2$se))*-15) #per 15 repeats increment
  
  df_results[paste(i,q),
             "Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  c(res_lpa_exp[1],res_kiv2_exp[1]), c(res_lpa_exp[2],res_kiv2_exp[2]), c(res_lpa_exp[3],res_kiv2_exp[3]))
  df_forest[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-c(rbind(res_lpa_exp,res_kiv2_exp))
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
  
  df_results[paste(i,q),"\n\nP"]<-"                   "
  df_forest[paste(i,q),"\nP"]<-c(res_lpa$pval,res_kiv2$pval)
  
  res_mv_exp<-exp(c(res_mv@Estimate,res_mv@CILower,res_mv@CIUpper)*c(log(2),-15))
  res_mv_exp[c(length(res_mv_exp) - 2, length(res_mv_exp))] <-res_mv_exp[c(length(res_mv_exp), length(res_mv_exp) - 2)] #oorrecting order of confidence intervals
  
  
  df_results[paste(i,q),
             "Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv_exp[1:2],res_mv_exp[3:4], res_mv_exp[5:6])
  df_forest[paste(i,q),c("Effect_mult","Lower 0.95_mult","Upper 0.95_mult")]<-res_mv_exp
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
  
  df_results[paste(i,q),"\n\nP "]<-"                   "
  df_forest[paste(i,q),"\nP "]<-res_mv@Pvalue
  
  
  
}

#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c(rep("white",1),rep("gray90",1)))))

xlab_1<-"\n
Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

xlab_2<-"\n
Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

p<-forestploter::forest(df_results,
                        est = list(df_forest$Effect,df_forest$Effect_mult),
                        lower =list(df_forest$`Lower 0.95`,df_forest$`Lower 0.95_mult`), 
                        upper = list(df_forest$`Upper 0.95`,df_forest$`Upper 0.95_mult`), 
                        sizes = 1,
                        ref_line = 1,
                        ci_column = c(5,8),
                        xlab = c(xlab_1,xlab_2),
                        xlim = c(0.9,1.8),ticks_at=c(0.9,seq(1,1.8,0.2)),
                        title=paste0("Mendelian randomization for coronary artery disease"),
                        theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line between results
p <- add_border(p, part = "header", row = c(3,5), where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1), gp =  grid::gpar(fontface = "bold"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(1),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(1),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(5,8), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(5,8), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))


#Adding P-values with superscript
df_forest[,4]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,4], c(-Inf, 0.10, Inf))], df_forest[,4])
df_forest[,8]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,8], c(-Inf, 0.10, Inf))], df_forest[,8])

p<-format_p(plot=p,pvals=list(df_forest[,4],df_forest[,8]),plot_cols=c(6,9))

p[["widths"]][[2]]<-unit(4.5,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")






# Association of genetic scores with variables ----------------------------------------------


#Removing all columns with missing values
valid_cols <- vapply(df, function(x) {
  x <- x[!is.na(x)]
  length(x) >= 2 && length(unique(x)) >= 2
}, logical(1))

df <- df[, valid_cols, drop = FALSE]

distrib <- datadist(df)

# Plotting effects of genes on metabolomics
df_summary<-data.frame()
df_forest<-list()

for (k in c("lpa_PRS","kiv2_PRS")){
  
  df_forest[[k]]<-data.frame()
  
  # Getting the change in PRS needed to change remnant cholesterol or LDL cholesterol by 1 mmol/L
  model<-ols(as.formula(paste0(if (grepl("lpa",k)) {"lpanew"} else {"kiv2"},
                               "~",k,"+alder*sex")),data=df)
  
  distrib<- datadist(df)
  
  # for difference scaled to 50 mg/dL higher lipoprotein(a): 
  distrib$limits[k][c("Low:effect","High:effect"),]<-c(0,(if(k=="lpa_PRS"){50} else if (k=="kiv2_PRS"){-15}/model$coefficients[2]))
  options(datadist = "distrib")
  
  # Getting the number of sd increments needed for a change corresponding to 1 mmol/L remnant cholesterol
  
  
  for (i in variables[ -which(variables %in% c("sex","alder","nonhdl_mgdl"))]){
    
    #Setting variable name
    df_summary[i,"\n "]<-gsub("\\,.*","",variables_names[ which(variables %in% i)])
    
    regress_formula<-as.formula(paste0(i,"~",k,"+sex*alder"))
    
    
    unit<-if(
      i %in% c("remncholcalc_direct","ldl_direct","trig","nonhdl")){"mmol/L"} else if(
        i %in% c("systolic","diastolic")){"mmHg"} else if(
          i %in% c("bmi")){"kg/m2"} else if(
            i %in% c("startage")){"years"} else if(
              i %in% c("cum_smoking")){"pack-years"} else if(
                i %in% c("creatinine")){"µmol/L"} else if(
                  i %in% c("hba1c")){"mmol/mol"} else if(
                    i %in% c("apob","lpanew")){"mg/dL"} else if(
                      i %in% c("hscrp")){"mg/L"} else if(
                        i %in% c("physical")){"MET hours/week"} else if(
                          i %in% c("kiv2")){"no."} 
    
    
    
    if (length(unique(df[!is.na(df[,i]),i]))>2) { # For continuous variables
      

      # log-model for left skewed data
      #if(mean(ukbb_50[,i],na.rm=T)-median(ukbb_50[,i],na.rm=T)>0.2*sd(ukbb_50[,i],na.rm=T)){
      #dat[,i]<-log2(ukbb_50[,i]+0.001)
      #unit<-"doublings"}
      
      
      model<-ols(regress_formula,data=df) #regression
      effect<-summary(model)[k,c("Effect","Lower 0.95","Upper 0.95")]/sd(df[,i],na.rm=T)#effect on SD increments
      
      df_forest[[k]][i,c("Effect","Lower 0.95","Upper 0.95")]<-effect  
      #Creating column with space for forest plot
      df_summary[i," "]<-paste(rep(" ", 25), collapse = " ")
      
      #effect on standard units
      effect<-effect*sd(df[,i],na.rm=T)
      
      if (unit!="doublings"){
        unit<-if(i %in% c("nonhdl")){
          paste0("mmol/L, ",sprintf("%.0f (%.0f to %.0f)",  effect[1]*38.67, effect[2]*38.67, effect[3]*38.67)," mg/dL")} else if (i %in% c("trig")) {
            paste0("mmol/L, ",sprintf("%.0f (%.0f to %.0f)",  effect[1]*88.57, effect[2]*88.57, effect[3]*88.57)," mg/dL")}else if (i %in% c("creatinine")) {
              paste0("µmol/L, ",sprintf("%.0f (%.0f to %.0f)",  effect[1]*88.42, effect[2]*88.42, effect[3]*88.42)," mg/dL")}else {unit} # in mg per dL for cholesterol and triglycerides
      }
      
      df_summary[i, paste0("Change per ",if (grepl("lpa",k)){"50 mg/dL\nhigher lipoprotein(a)"} else {"15 fewer\nKIV-2 repeats"}," (95% CI)")]<-paste0(sprintf("%.1f (%.1f to %.1f)",  effect[1], effect[2], effect[3]) ," ",unit )
      #Creating column with space for forest plot
      df_summary[i,"  "]<-paste(rep(" ", 25), collapse = " ")
      
    } else { # For binary variables
      
      
      model<-lrm(regress_formula,data=df)
      effect<-summary(model)[k,c("Effect","Lower 0.95","Upper 0.95")]
      
      df_forest[[k]][i,c("Effect","Lower 0.95","Upper 0.95")]<-effect #effect log-odds ratio
      
      #Creating column with space for forest plot and for extra space after
      df_summary[i," "]<-paste(rep(" ", 25), collapse = " ")
      
      effect<-exp(effect)
      df_summary[i, paste0("Change per ",if (grepl("lpa",k)){"50 mg/dL\nhigher lipoprotein(a)"} else {"15 fewer\nKIV-2 repeats"}," (95% CI)")]<-paste0("OR: ",sprintf("%.1f (%.1f to %.1f)",  effect[1], effect[2], effect[3])  ) #effect odds ratio
      
      #Creating column with space for forest plot
      df_summary[i,"  "]<-paste(rep(" ", 25), collapse = " ")
      
    }
    
 
  }
  
}



#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = rep(c(rep("white",1),rep("gray95",1)),nrow(df_summary)/2)),
                     fg_params = list( fontsize = 9, fontfamily = "arial")))



#Formatting shape of points for odds ratios
OR_rows<-c(4,7,8,9)

size_c<-rep(1,nrow(df_summary))
size_c[OR_rows]<-1.001

pointtype<-function (est, lower, upper, pch, sizes = 1, gp = gpar(), t_height = NULL, 
                     xlim = c(0, 1), nudge_y = 0, name = NULL) 
{
  gTree(est = est, lower = lower, upper = upper, pch = if(sizes==1.001){18} else{20}, 
        size = sizes, gp = gp, t_height = t_height, xlim = xlim, 
        nudge_y = nudge_y, name = name, cl = "makeci")
}


p<-forestploter::forest(df_summary,
          est = list(df_forest[[1]]$Effect,df_forest[[2]]$Effect),
          lower =list(df_forest[[1]]$`Lower 0.95`,df_forest[[2]]$`Lower 0.95`), 
          upper = list(df_forest[[1]]$`Upper 0.95`,df_forest[[2]]$`Upper 0.95`),
          sizes = list(size_c,size_c),
          fn_ci =pointtype,
          ref_line = 0,
          ci_column = c(2,4),  
          xlab = list("\n\nChange per 50 mg/dL higher lipoprotein(a),\n   standard deviations   log-OR (95% CI)",
                      "\n\nChange per 15 fewer KIV-2 repeats,\n   standard deviations   log-OR (95% CI)"),
          xlim = list(c(-2,2),c(-2,2)),
          ticks_at = list(seq(-2,2,1),seq(-2,2,1)),
          theme=tm)


#adding line under headerhttp://127.0.0.1:35275/graphics/plot_zoom_png?width=2165&height=1131
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))


#blue font for odds ratios

col<-"gray40"

p<-edit_plot(p, col = c(2,4), row = OR_rows, which = "ci", gp = gpar( fill =col))

#changing width of column
p[["widths"]][[2]]<-unit(5,"cm")
p[["widths"]][[3]]<-unit(4.5,"cm")
p[["widths"]][[4]]<-unit(7,"cm")
p[["widths"]][[5]]<-unit(4.5,"cm")
p[["widths"]][[6]]<-unit(7,"cm")


for (i in 4:length(p[["heights"]])){ #Changing heights of rows
  p[["heights"]][[i]] <- unit(0.5,"cm")
}


# Exporting plots to powerpoint
doc <- read_pptx("widescreen.pptx")

doc<-add_slide(doc)
doc <- ph_with(x = doc, value=rvg::dml(ggobj=ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# Supplementary tables for summary data -----------------------------------


library(openxlsx)
library(dplyr)

# Prepare CAD outcome table (only needed columns)
cad_df <- IHD[["meta"]] %>%
  dplyr::select(
    SNP,
    beta.outcome,
    se.outcome,
    pval.outcome
  ) %>%
  dplyr::rename(
    "Beta for CAD" = beta.outcome,
    "Standard error for CAD" = se.outcome,
    "P-value for CAD" = pval.outcome
  )

# And cluster columns
ratio_cluster_df <- ggplot_data %>%
  dplyr::select(SNP, ratio, Cluster) %>%
  dplyr::rename(
    "Ratio KIV-2/lipoprotein(a)" = ratio,
    "Cluster" = Cluster)


# Create Excel workbook
wb <- createWorkbook()

# Helper function: clean exposure table + merge CAD info
format_df_with_cad <- function(df, cad_df, ratio_cluster_df) {
  
  exposure_df <- df %>%
    dplyr::select(
      SNP,
      pos.exposure,
      effect_allele.exposure,
      other_allele.exposure,
      beta.exposure,
      se.exposure,
      pval.exposure
    ) %>%
    dplyr::rename(
      "Genetic variant"  = SNP,
      "Position"         = pos.exposure,
      "Effect allele"    = effect_allele.exposure,
      "Reference allele" = other_allele.exposure,
      "Beta"             = beta.exposure,
      "Standard error"   = se.exposure,
      "P-value"          = pval.exposure
    )
  
  full_df <- exposure_df %>%
    left_join(cad_df, by = c("Genetic variant" = "SNP")) %>%
    left_join(ratio_cluster_df, by = c("Genetic variant" = "SNP"))
  
  return(full_df)
}

# Loop over dataframes in list_uni
for (name in names(list_uni)) {
  
  addWorksheet(wb, sheetName = name)
  
  df_out <- format_df_with_cad(
    list_uni[[name]],
    cad_df,
    ratio_cluster_df)
  
  
  writeDataTable(
    wb,
    sheet = name,
    x = df_out,
    tableStyle = "TableStyleMedium2"
  )
  
  setColWidths(
    wb,
    sheet = name,
    cols = 1:ncol(df_out),
    widths = "auto"
  )
}



# Save file
saveWorkbook(
  wb,
  file = "list_uni_genetic_variants_with_CAD.xlsx",
  overwrite = TRUE
)


# MR using lipoprotein(a) data from the UK Biobank ------------------------
# Replace lipoprotein(a) estimates with UK Biobank (ukb-d-30790_raw), nmol/L --

# All Lp(a) SNPs across data frames
lpa_snps <- unique(c(
  list_uni[["LPA"]]$SNP[list_uni[["LPA"]]$exposure == "Lipoprotein(a)"],
  unlist(lapply(list_comb, function(d) d$SNP[d$exposure == "Lipoprotein(a)"]))
))

# One query; no proxies (rare LPA variants should not be proxied)
ukb_lpa <- ieugwasr::associations(lpa_snps, "ukb-d-30790_raw", proxies = 0)
ukb_lpa <- ukb_lpa[!duplicated(ukb_lpa$rsid), ]

# SNPs not available in UK Biobank
missing_ukb <- setdiff(lpa_snps, ukb_lpa$rsid)
length(missing_ukb); missing_ukb

replace_lpa <- function(d) {
  is_lpa <- d$exposure == "Lipoprotein(a)"
  u      <- ukb_lpa[match(d$SNP[is_lpa], ukb_lpa$rsid), ]
  
  ea <- d$effect_allele.exposure[is_lpa]
  oa <- d$other_allele.exposure[is_lpa]
  
  same <- !is.na(u$ea) & u$ea == ea & u$nea == oa
  flip <- !is.na(u$ea) & u$ea == oa & u$nea == ea
  sgn  <- ifelse(same, 1, ifelse(flip, -1, NA))
  
  d$beta.exposure[is_lpa]  <- sgn * u$beta
  d$se.exposure[is_lpa]    <- ifelse(is.na(sgn), NA, u$se)
  d$pval.exposure[is_lpa]  <- ifelse(is.na(sgn), NA, u$p)
  d$eaf.exposure[is_lpa]   <- ifelse(same, u$eaf, ifelse(flip, 1 - u$eaf, NA))
  d$samplesize.exposure[is_lpa] <- u$n
  d$data_source.exposure[is_lpa] <- "ukb-d-30790_raw"
  
  # Drop Lp(a) rows not found in UKB or with mismatching alleles
  dropped <- d$SNP[is_lpa][is.na(sgn)]
  if (length(dropped) > 0) message("Dropped Lp(a) rows: ", paste(dropped, collapse = ", "))
  d[!(is_lpa & is.na(d$beta.exposure)), ]
}

list_uni_ukbb<-list_uni
list_comb_ukbb<-list_comb

list_uni_ukbb[["LPA"]] <- replace_lpa(list_uni_ukbb[["LPA"]])
list_comb_ukbb         <- lapply(list_comb_ukbb, replace_lpa)

#remove all snps with missing data for lpa from other exposure too
list_comb_ukbb <- lapply(list_comb_ukbb, function(d) {
  keep_snps <- d$SNP[d$exposure == "Lipoprotein(a)" & !is.na(d$beta.exposure)]
  d[d$SNP %in% keep_snps, ]
})



df_results<-data.frame()
df_forest<-data.frame()

for(i in c("KIV2","isodom","isosum")){
  
  
  data<-dat_to_MRInput(harmonise_data(list_uni_ukbb[["LPA"]],IHD[["meta"]]))
  data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
  data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
  res_lpa<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs)
  
  data<-dat_to_MRInput(harmonise_data(list_uni_ukbb[[i]], IHD[["meta"]]))
  data[[1]]@correlation<-R[data[[1]]@snps,data[[1]]@snps]
  data<-MendelianRandomization::mr_ivw(data[[1]],correl=T) 
  res_kiv2<-data.frame(b=data@Estimate,se=data@StdError,pval=data@Pvalue,nsnp= data@SNPs)
  
  #Filtering correlation matrix with kept variants
  mvdat<-mv_harmonise_data(list_comb_ukbb[[i]],IHD[["meta"]],harmonise_strictness =1) #Harmonizing data, keeping only variants present in outcome data
  
  snps_kept <- rownames(mvdat$exposure_beta)
  R_filtered<-R[snps_kept, snps_kept, drop = FALSE]
  mv_input<-mr_mvinput(
    bx = mvdat$exposure_beta,     # SNP–exposure associations matrix
    by  = mvdat$outcome_beta,     # SNP–outcome associations vector
    bxse = mvdat$exposure_se,
    byse  = mvdat$outcome_se,
    snps   = rownames(mvdat$exposure_beta),
    exposure=mvdat$expname$exposure,
    correlation=R_filtered)
  
  res_mv<-mr_mvivw(mv_input)
  
  q<- c("Lipoprotein(a)", if(grepl("KIV",i)){"Sum of KIV-2"} else if(
    grepl("isodom",i)){"Dominant apo(a) isoform"}else if(
      grepl("isosum",i)){"Sum of apo(a) isoforms"})
  
  df_results[paste(i,q),"\nGenetic\nscore"]<-q
  
  df_results[paste(i,q),"\n\nUnit change"]<-c("Per 100 nmol/L higher levels","Per 15 fewer repeats")
  
  df_results[paste(i,q),"\nNo. of\nSNPs"]<-formatC(c(res_lpa$nsnp,res_kiv2$nsnp),format="d", big.mark=",")
  
  res_lpa_exp<-exp(c(res_lpa$b,res_lpa$b-(1.959*res_lpa$se),res_lpa$b+(1.959*res_lpa$se))*100) #per 100 nmol/l increment
  res_kiv2_exp<-exp(c(res_kiv2$b,res_kiv2$b+(1.959*res_kiv2$se),res_kiv2$b-(1.959*res_kiv2$se))*-15) #per 20 repeats increment
  
  df_results[paste(i,q),
             "Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  c(res_lpa_exp[1],res_kiv2_exp[1]), c(res_lpa_exp[2],res_kiv2_exp[2]), c(res_lpa_exp[3],res_kiv2_exp[3]))
  df_forest[paste(i,q),c("Effect","Lower 0.95","Upper 0.95")]<-c(rbind(res_lpa_exp,res_kiv2_exp))
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"    "]<-paste(rep(" ", 20), collapse = " ")
  
  df_results[paste(i,q),"\n\nP"]<-"                   "
  df_forest[paste(i,q),"\nP"]<-c(res_lpa$pval,res_kiv2$pval)
  
  res_mv_exp<-exp(c(res_mv@Estimate,res_mv@CILower,res_mv@CIUpper)*c(100,-15))
  
  res_mv_exp[c(length(res_mv_exp) - 2, length(res_mv_exp))] <-res_mv_exp[c(length(res_mv_exp), length(res_mv_exp) - 2)] #oorrecting order of confidence intervals
  
  
  df_results[paste(i,q),
             "Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"]<-sprintf("%.2f (%.2f-%.2f)",  res_mv_exp[1:2],res_mv_exp[3:4], res_mv_exp[5:6])
  df_forest[paste(i,q),c("Effect_mult","Lower 0.95_mult","Upper 0.95_mult")]<-res_mv_exp
  
  #Creating column with space for forest plot
  df_results[paste(i,q),"     "]<-paste(rep(" ", 20), collapse = " ")
  
  df_results[paste(i,q),"\n\nP "]<-"                   "
  df_forest[paste(i,q),"\nP "]<-res_mv@Pvalue
  
  
  
}

#Setting theme for plot
tm <- forest_theme(base_size = 9,
                   ci_pch = 20,ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd=2,
                   title_col = "gray40",title_fontface = "plain",
                   core=list(
                     bg_params=list(fill = c(rep("white",1),rep("gray90",1)))))

xlab_1<-"\n
Univariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

xlab_2<-"\n
Multivariable Mendelian\nrandomization odds ratio (95% CI)\nfor risk of coronary artery disease"

p<-forestploter::forest(df_results,
                        est = list(df_forest$Effect,df_forest$Effect_mult),
                        lower =list(df_forest$`Lower 0.95`,df_forest$`Lower 0.95_mult`), 
                        upper = list(df_forest$`Upper 0.95`,df_forest$`Upper 0.95_mult`), 
                        sizes = 1,
                        ref_line = 1,
                        ci_column = c(5,8),
                        xlab = c(xlab_1,xlab_2),
                        xlim = c(0.9,1.6),ticks_at=c(0.9,seq(1,1.6,0.2)),
                        title=paste0("Mendelian randomization for coronary artery disease"),
                        theme=tm)

#adding line under header
p <- add_border(p, part = "header", row = 1, where = "bottom",gp = grid::gpar(lty=1,lwd=2))

#adding line between results
p <- add_border(p, part = "header", row = c(3,5), where = "bottom",gp = grid::gpar(lty=1,lwd=1))


#Bold font for first column
p <- edit_plot(p, col = c(1), gp =  grid::gpar(fontface = "bold"))

# Changing colors
uneven<-subset(1:200,1:200 %% 2 == 1)
even<-subset(1:200, 1:200 %% 2 == 0)

colors<-c("darkcyan","darkorchid4")

#For text
p <- edit_plot(p, col = c(1),row =uneven , which="text",gp = gpar(col = colors[1]))
p <- edit_plot(p, col = c(1),row =even , which="text",gp = gpar(col = colors[2]))

#For confidence interval
p<-edit_plot(p, col = c(5,8), row = uneven, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p<-edit_plot(p, col = c(5,8), row = even, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))


#Adding P-values with superscript
df_forest[,4]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,4], c(-Inf, 0.10, Inf))], df_forest[,4])
df_forest[,8]<-sprintf(c("%3.1g", "%3.2g")[cut(df_forest[,8], c(-Inf, 0.10, Inf))], df_forest[,8])

p<-format_p(plot=p,pvals=list(df_forest[,4],df_forest[,8]),plot_cols=c(6,9))

p[["widths"]][[2]]<-unit(4.5,"cm")
p[["widths"]][[5]]<-unit(4,"cm")
p[["widths"]][[8]]<-unit(4,"cm")


doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(p)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# Supplementary figure: rs2457564 and lp(a) levels in Copenhagen and UKBB ------------------------



df<-fread("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/LPA_SNPs/df_colnamescorrected.tsv")

df_cchs_outcomes<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Dataset_obus_outcomes.dta")
df_cgps_outcomes<-read_dta("L:/Auditdata/Forskning KBA/Benjamin Nilsson Wadstrom/Dataset_obush.dta")

df_cchs_outcomes$obushnr<-df_cchs_outcomes$lobenr+9e4
df_cgps_outcomes<-subset(df_cgps_outcomes,duplicate_between_cgps_cchs1_4==0 &duplicate_within_cgps==0 & obushnr<9e4)

df_outcomes<-rbind(df_cchs_outcomes[,c("obushnr","ihd","ami")],df_cgps_outcomes[,c("obushnr","ihd","ami")])

df<-merge(df,df_outcomes[,c("obushnr","ihd","ami")],by="obushnr",all.x=T,all.y=F)
df<-subset(df,!is.na(kiv2) & !is.na(lpanew) & !is.na(ihd))

distrib<- datadist(df)
options(datadist = "distrib")


df$rs2457564<-df$`160877640_G_A`
df$rs2504926<-df$`160876290_T_C`
df$rs10455872<-df$`161010118_A_G`

r <- cor(df$rs2457564, df$rs2504926, use = "complete.obs")
r^2

# rs2457564, lipoprotein(a) and coronary artery disease (ieu-a-7) 
# Effect allele: G. All lipoprotein(a) results in mg/dL.

effect_allele <- "G"
nmol_per_mgdl <- 2.15   # mg/dL = nmol/L / 2.15
per           <- 50     # OR per 50 mg/dL higher Lp(a)

# Helpers
align_to <- function(beta, ea, nea, allele = effect_allele) {
  if (ea == allele) return(beta)
  if (nea == allele) return(-beta)
  stop("Effect allele not found in ea/nea")
}
get_p <- function(beta, se) 2 * pnorm(-abs(beta / se))
format_pval <- function(p) {
  ifelse(p == 0, "<1e-300",
         sprintf(c("%3.1g", "%3.2g")[cut(p, c(-Inf, 0.10, Inf))], p))
}

# Lipoprotein(a): Copenhagen studies
model <- ols(lpanew ~ rs2457564 + sex + alder, data = df)

beta_cph <- -coef(model)["rs2457564"]                    # flipped to G allele
se_cph   <- sqrt(diag(vcov(model)))["rs2457564"]
n_cph    <- model$stats["n"]
eaf_cph  <- 1 - mean(df$rs2457564[!is.na(df$lpanew)], na.rm = TRUE) / 2  # assumes dosage of A allele

# Lipoprotein(a): UK Biobank 
ukb <- ieugwasr::associations("rs2457564", "ukb-d-30790_raw")

beta_ukb <- align_to(ukb$beta, ukb$ea, ukb$nea) / nmol_per_mgdl
se_ukb   <- ukb$se / nmol_per_mgdl
n_ukb    <- ukb$n
eaf_ukb  <- ifelse(ukb$ea == effect_allele, ukb$eaf, 1 - ukb$eaf)

# CAD: ieu-a-7 (same outcome for both rows) 
cad      <- ieugwasr::associations("rs2457564", "ieu-a-7")
cad_info <- ieugwasr::gwasinfo("ieu-a-7")

beta_cad <- align_to(cad$beta, cad$ea, cad$nea)          # log odds per G allele
se_cad   <- cad$se

# Wald ratios: log OR for CAD per 1 mg/dL higher Lp(a)
b_wr  <- beta_cad / c(beta_cph, beta_ukb)
se_wr <- se_cad   / abs(c(beta_cph, beta_ukb))

# Effects for forest plot 
effects_lpa <- data.frame(
  Effect = c(beta_cph, beta_ukb),
  Lower  = c(beta_cph - 1.96 * se_cph, beta_ukb - 1.96 * se_ukb),
  Upper  = c(beta_cph + 1.96 * se_cph, beta_ukb + 1.96 * se_ukb)
)
effects_lpa$P <- format_pval(get_p(c(beta_cph, beta_ukb), c(se_cph, se_ukb)))

effects_cad <- data.frame(
  Effect = exp(per * b_wr),
  Lower  = exp(per * (b_wr - 1.96 * se_wr)),
  Upper  = exp(per * (b_wr + 1.96 * se_wr))
)
effects_cad$P <- format_pval(get_p(b_wr, se_wr))

# Text table
snp_table <- data.frame(
  "\nLipoprotein(a)\ndata"         = c("Copenhagen studies", "UK Biobank"),
  "\nNo. of\nindividuals"          = formatC(c(n_cph, n_ukb), format = "d", big.mark = ","),
  "\nEffect allele\nfrequency"     = sprintf("%.2f", c(eaf_cph, eaf_ukb)),
  "\nChange (95% CI) in\nlipoprotein(a) per allele, mg/dL" =
    sprintf("%.2f (%.2f-%.2f)", effects_lpa$Effect, effects_lpa$Lower, effects_lpa$Upper),
  "  "                             = paste(rep(" ", 20), collapse = " "),
  "\n\nP"                          = "            ",
  "Odds ratio (95% CI)\nfor risk of coronary artery disease\nper 50 mg/dL higher lipoprotein(a)" =
    sprintf("%.2f (%.2f-%.2f)", effects_cad$Effect, effects_cad$Lower, effects_cad$Upper),
  " "                              = paste(rep(" ", 20), collapse = " "),
  "\n\nP "                         = "            ",
  check.names = FALSE
)

# Plot 
tm <- forest_theme(base_size = 9, base_family = "arial",
                   ci_pch = 20, ci_lwd = 2, ci_Theight = 0.4,
                   xaxis_lwd = 2,
                   title_col = "gray40", title_fontface = "plain",
                   core = list(bg_params = list(fill = c("white", "white"))))

p <- forestploter::forest(snp_table,
                          est   = list(effects_lpa$Effect, effects_cad$Effect),
                          lower = list(effects_lpa$Lower,  effects_cad$Lower),
                          upper = list(effects_lpa$Upper,  effects_cad$Upper),
                          sizes = 1,
                          ref_line = c(0, 1),
                          ci_column = c(5, 8),
                          xlab = c("\nChange (95% CI) in\nlipoprotein(a) per allele, mg/dL",
                                   "\n\nOdds ratio (95% CI)\nfor risk of coronary artery disease\nper 50 mg/dL higher lipoprotein(a)"),
                          xlim = list(c(0, 8), c(0.9, 2)),
                          ticks_at = list(seq(0, 8, 2), seq(1, 2, 0.2)),
                          title = "rs2457564, lipoprotein(a) and coronary artery disease",
                          theme = tm)

# Line under header
p <- add_border(p, part = "header", row = 1, where = "bottom", gp = grid::gpar(lty = 1, lwd = 2))

# Bold first column
p <- edit_plot(p, col = 1, gp = grid::gpar(fontface = "bold"))

# Colours
colors <- c("darkcyan", "darkorchid4")
p <- edit_plot(p, col = 1, row = 1, which = "text", gp = gpar(col = colors[1]))
p <- edit_plot(p, col = 1, row = 2, which = "text", gp = gpar(col = colors[2]))
p <- edit_plot(p, col = c(5, 8), row = 1, which = "ci", gp = gpar(col = colors[1], fill = colors[1]))
p <- edit_plot(p, col = c(5, 8), row = 2, which = "ci", gp = gpar(col = colors[2], fill = colors[2]))

# P-values with superscript (your function)
p <- format_p(plot = p, pvals = list(effects_lpa$P, effects_cad$P), plot_cols = c(6, 9))


# Export to PowerPoint 
doc <- read_pptx("widescreen.pptx")
doc <- add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj = ggplotify::as.ggplot(p)),
               location = ph_location_fullsize(newlabel = ""))
print(doc, target = "temp_snp.pptx")


## LD in South asian and European


library(LDlinkR)

token <- "080ca00c0a14"   # or paste the string directly

# r2 and D' for the pair across every 1000 Genomes population
data<-LDpop(var1  = "rs2457564",
          var2  = "rs10455872",
          pop   = c("PJL","EUR"),
          r2d   = "r2",
          token = token)

LDpair("rs2457564", "rs10455872", pop = "EUR", token = token, output = "text")

# Fraction mediated by kiv-2 and isoform using product method --------------------------------------------------


# Product method with Monte Carlo confidence interval
prop_mediated_product <- function(b_xm, se_xm,        # X -> M  (Lp(a) -> apo(a) size)
                                  b_my, se_my,        # M -> Y | X (multivariable apo(a) estimate)
                                  b_total, se_total,  # X -> Y total (univariable Lp(a))
                                  n_sim = 1e6, alpha = 0.05, seed = 2026) {
  set.seed(seed)
  sim_xm    <- rnorm(n_sim, b_xm,    se_xm)
  sim_my    <- rnorm(n_sim, b_my,    se_my)
  sim_total <- rnorm(n_sim, b_total, se_total)
  
  prop_med  <- (sim_xm * sim_my) / sim_total      # proportion mediated by apo(a) size
  prop_lpa  <- 1 - prop_med                       # proportion due to Lp(a)
  
  est_med <- (b_xm * b_my) / b_total
  ci_med  <- quantile(prop_med, c(alpha / 2, 1 - alpha / 2))
  
  # P-value for fraction due to Lp(a) (H0: 0), from simulated distribution
  se_lpa <- sd(prop_lpa)
  p_lpa  <- 2 * pnorm(-abs((1 - est_med) / se_lpa))
  
  list(
    P_med  = est_med * 100,
    CI_med = pmin(setNames(ci_med * 100, c("lower", "upper")), 100),
    P_lpa  = (1 - est_med) * 100,
    CI_lpa = pmin(setNames(rev(1 - ci_med) * 100, c("lower", "upper")), 100),
    p_lpa  = p_lpa
  )
}


# Analysis 
df_results <- data.frame()
df_forest  <- data.frame()

for (i in c("KIV2", "isodom", "isosum")) {
  IHD[["meta"]] <- subset(IHD[["meta"]], !duplicated(SNP))
  
  # Univariable Lp(a) -> CAD (total effect)
  data <- dat_to_MRInput(harmonise_data(list_uni[["LPA"]], IHD[["meta"]]))
  data[[1]]@correlation <- R[data[[1]]@snps, data[[1]]@snps]
  data <- MendelianRandomization::mr_ivw(data[[1]], correl = TRUE)
  res_lpa <- data.frame(b = data@Estimate, se = data@StdError, pval = data@Pvalue, nsnp = data@SNPs)
  
  # Univariable apo(a) measure (only used for No. of SNPs column)
  data <- dat_to_MRInput(harmonise_data(list_uni[[i]], IHD[["meta"]]))
  data[[1]]@correlation <- R[data[[1]]@snps, data[[1]]@snps]
  data <- MendelianRandomization::mr_ivw(data[[1]], correl = TRUE)
  res_kiv2 <- data.frame(b = data@Estimate, se = data@StdError, pval = data@Pvalue, nsnp = data@SNPs)
  
  # Multivariable: apo(a) measure -> CAD | Lp(a)
  mvdat <- mv_harmonise_data(list_comb[[i]], IHD[["meta"]], harmonise_strictness = 1)
  snps_kept  <- rownames(mvdat$exposure_beta)
  R_filtered <- R[snps_kept, snps_kept, drop = FALSE]
  
  mv_input <- mr_mvinput(
    bx = mvdat$exposure_beta, by = mvdat$outcome_beta,
    bxse = mvdat$exposure_se, byse = mvdat$outcome_se,
    snps = snps_kept, exposure = mvdat$expname$exposure,
    correlation = R_filtered)
  res_mv <- mr_mvivw(mv_input)
  
  exp_names <- mv_input@exposure
  is_lpa    <- exp_names == "Lipoprotein(a)"
  
  b_my  <- res_mv@Estimate[!is_lpa]
  se_my <- res_mv@StdError[!is_lpa]
  
  # Lp(a) -> apo(a) measure, using the same variants and correlation matrix
  xm_input <- mr_input(
    bx   = mvdat$exposure_beta[,  is_lpa],
    bxse = mvdat$exposure_se[,    is_lpa],
    by   = mvdat$exposure_beta[, !is_lpa],
    byse = mvdat$exposure_se[,   !is_lpa],
    snps = snps_kept,
    correlation = R_filtered)
  res_xm <- mr_ivw(xm_input, correl = TRUE)
  
  # Product method
  med <- prop_mediated_product(
    b_xm    = res_xm@Estimate, se_xm    = res_xm@StdError,
    b_my    = b_my,            se_my    = se_my,
    b_total = res_lpa$b,       se_total = res_lpa$se)
  
  q <- "Lipoprotein(a)"
  
  df_results[paste(i, q), "\n\n"] <- " "
  df_results[paste(i, q), "\n\nAnalysis"] <- paste0("Lipoprotein(a) + ",
                                                    if (grepl("KIV", i)) "Sum of KIV-2" else if (grepl("isodom", i)) "Dominant apo(a) isoform"
                                                    else "Sum of apo(a) isoforms")
  
  df_results[paste(i, q), "\nNo. of\nSNPs"] <- paste0(res_lpa$nsnp, " + ",
                                                      formatC(res_kiv2$nsnp, format = "d", big.mark = ","))
  
  df_results[paste(i, q),
             "Fraction of coronary artery\ndisease risk due to\nhigh lipoprotein(a) (95% CI)"] <-
    sprintf("%.f%% (%.f–%.f%%)", med$P_lpa, med$CI_lpa["lower"], med$CI_lpa["upper"])
  
  df_results[paste(i, q),
             "Fraction of coronary artery\ndisease risk\ndue to short apo(a) (95% CI)"] <-
    sprintf("%.f%% (%.f–%.f%%)", med$P_med, med$CI_med["lower"], med$CI_med["upper"])
  
  df_forest[paste(i, q), c("Effect", "Lower 0.95", "Upper 0.95")] <-
    c(med$P_lpa, med$CI_lpa["lower"], med$CI_lpa["upper"])
  
  df_results[paste(i, q), "    "] <- paste(rep(" ", 20), collapse = " ")
  df_results[paste(i, q), "P for fraction\ndue to\nlipoprotein(a)"] <- "                   "
  df_forest[paste(i, q), "\nP"] <- med$p_lpa
}



library(gridExtra)
library(gtable)
library(grid)




plot_df <- data.frame(
  label    = rownames(df_results),
  est_not  = df_forest$Effect,          # NOT mediated
  est_med  = 100 - df_forest$Effect,    # mediated
  row      = seq_len(nrow(df_forest))
)


# Optional: grouping for alternating colours
plot_df$group <- ifelse(plot_df$row %% 2 == 1, "odd", "even")

#Function for bar plot
bar_grob <- function(value_not,
                     value_med,
                     xmin = 0,
                     xmax = 100,
                     col_not = "darkcyan",
                     col_med = "darkorchid4") {
  
  # Clamp
  value_not <- max(min(value_not, xmax), xmin)
  value_med <- max(min(value_med, xmax - value_not), xmin)
  
  w_not <- value_not / (xmax - xmin)
  w_med <- value_med / (xmax - xmin)
  
  gTree(children = gList(
    
    ## NOT mediated
    rectGrob(
      x = w_not / 2,
      y = 0.5,
      width  = w_not,
      height = 0.6,
      gp = gpar(fill = col_not, col = NA)
    ),
    
    ## Mediated (stacked to the right)
    rectGrob(
      x = w_not + w_med / 2,
      y = 0.5,
      width  = w_med,
      height = 0.6,
      gp = gpar(fill = col_med, col = NA)
    )
    
  ))
}


#Creating table


tbl <- tableGrob(
  df_results,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 9,
    core = list(
      fg_params = list(
        hjust = 0,        # left-align text
        x = 0.02          # small left padding
      )
    ),
    colhead = list(
      fg_params = list(
        hjust = 0,
        x = 0.02,
        fontface = "bold"
      )
    )
  )
)

##Adding barcols to table
bar_col <- 6 
tbl$widths[[bar_col]] <- unit(4, "cm")

n_bar <- nrow(plot_df)

for (i in seq_len(n_bar)) {
  
  g <- bar_grob(
    value_not = plot_df$est_not[i],
    value_med = plot_df$est_med[i]
  )
  
  idx <- which(
    tbl$layout$t == (i + 1) &
      tbl$layout$l == bar_col &
      grepl("-fg$", tbl$layout$name)
  )
  
  stopifnot(length(idx) == 1)
  
  tbl$grobs[[idx]] <- g
}


# Create x-axis grob (0–100)
x_axis <- xaxisGrob(
  at = seq(0, 1, by = 0.2),
  label = seq(0, 100, by = 20),
  gp = gpar(fontsize = 9)
)



# Add a new row at the bottom
tbl <- gtable_add_rows(
  tbl,
  heights =  unit(0, "mm")
)

# Index of the new row
axis_row <- length(tbl$heights)

# Add x-axis ONLY under the bar column
tbl <- gtable_add_grob(
  tbl,
  grobs = x_axis,
  t = axis_row,
  b = axis_row,   # ✅ THIS is the missing piece
  l = bar_col,
  r = bar_col,
  clip = "off"
)


#Adding veritcal line at 0 

first_data_row <- 2
last_data_row  <- first_data_row + nrow(plot_df) - 1


vline_0 <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(0, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(1, "npc"),
  gp = gpar(
    col = "black",
    lwd = 1,
    lty = c(3, 1)
  )
)



tbl <- gtable_add_grob(
  tbl,
  grobs = vline_0,
  t = first_data_row,
  b = last_data_row,
  l = bar_col,
  r = bar_col,
  clip = "off"
)


#Applying line under header

header_bottom_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(1, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(0, "npc"),
  gp = gpar(col = "black", lwd = 2)
)

tbl <- gtable_add_grob(
  tbl,
  grobs = header_bottom_line,
  t = 1,
  b = 1,
  l = 1,
  r = ncol(tbl),
  clip = "off"
)


#Adding other lines
header_bottom_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(1, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(0, "npc"),
  gp = gpar(col = "black", lwd = 1)
)


for (i in c(4)){
  tbl <- gtable_add_grob(
    tbl,
    grobs = header_bottom_line,
    t = i,
    b = 1,
    l = 1,
    r = ncol(tbl),
    clip = "off"
  )
}




format_p_expr <- function(p, digits = 2) {
  if (is.na(p)) return(expression(NA))
  
  if (p >0.001) {
    return(signif(p,2))
  }
  
  expo <- floor(log10(p))
  base <- signif(p / 10^expo,1)
  base <- round(base, digits)
  
  # Return plotmath expression: base %*% 10^expo
  bquote(.(base)~ "×" ~ 10^.(expo))
}

p_exprs <- lapply(df_forest[,"\nP"], format_p_expr)


p_col <- 7

stopifnot(p_col <= ncol(tbl))

n_p <- length(p_exprs)

for (i in seq_len(n_p)) {
  
  p_grob <- textGrob(
    label = p_exprs[[i]],
    gp = gpar(fontsize = 9),
    hjust = 0,                 # left-aligned
    x = unit(0.02, "npc")         # 🔴 critical for superscripts
  )
  
  idx <- which(
    tbl$layout$t == (i + 1) &  # +1 because row 1 is header
      tbl$layout$l == p_col &
      grepl("-fg$", tbl$layout$name)
  )
  
  stopifnot(length(idx) >= 1)
  
  tbl$grobs[[idx[1]]] <- p_grob
}


xlab_grob <- textGrob(
  "\n\nFraction of coronary artery disease risk\nCyan = directly due to elevated lipoprotein(a)\nPurple = indirectly due to short apo(a)",
  gp = gpar(fontsize = 10),
  just = "center"
)


tbl <- gtable_add_rows(
  tbl,
  heights = unit(22, "mm")   # dedicated label row
)

label_row <- length(tbl$heights)


tbl <- gtable_add_grob(
  tbl,
  grobs = xlab_grob,
  t = label_row,
  b = label_row,
  l = bar_col,
  r = bar_col,
  clip = "off"
)

# Shade background of the middle data row only

middle_row <- 2                    # data row index
table_row  <- middle_row + 1        # +1 because row 1 is header

bg_grob <- rectGrob(
  gp = gpar(fill = "grey90", col = NA)
)

tbl <- gtable_add_grob(
  tbl,
  grobs = bg_grob,
  t = table_row,
  b = table_row,
  l = 1,
  r = ncol(tbl),
  z = -Inf,        # 🔑 ensures background is BEHIND text and bars
  clip = "off"
)

tbl$widths[[1]] <- unit(0, "cm")

grid.newpage()
grid.draw(tbl)

doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(tbl)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# Fraction mediated by kiv-2 and isoform using variants in low LD --------------------------------------------------



prop_mediated_ci <- function(beta_total, se_total,
                             beta_direct, se_direct,
                             alpha = 0.05) {
  
  # Proportion NOT mediated (direct effect)
  Q <- beta_direct / beta_total
  
  # Delta-method variance
  var_Q <- ( (beta_direct / beta_total^2)^2 * se_total^2 ) +
    ( (1 / beta_total)^2 * se_direct^2 )
  
  se_Q <- sqrt(var_Q)
  
  # Wald Z
  z_val <- Q / se_Q
  
  # Numerically stable two-sided p-value
  log_p <- log(2) + pnorm(abs(z_val),
                          lower.tail = FALSE,
                          log.p = TRUE)
  p_val <- exp(log_p)
  
  # Confidence interval
  z_alpha <- qnorm(1 - alpha / 2)
  CI <- c(
    lower = Q - z_alpha * se_Q,
    upper = Q + z_alpha * se_Q
  )
  
  # Convert to percent and cap at 100
  CI_pct <- pmin(CI * 100, 100)
  
  list(
    P    = Q * 100,
    SE   = se_Q,
    CI   = CI_pct,
    Z    = z_val,
    p    = p_val,
    logp = log_p
  )
}




# Sensitivity analysis: LD-pruned instruments (r2 < 0.01), delta method 

r2_prune <- 0.01

ld_prune_snps <- function(snps, R, r2_threshold = 0.01) {
  
  R_sub <- R[snps, snps, drop = FALSE]
  r2 <- R_sub^2
  
  keep <- rep(TRUE, length(snps))
  names(keep) <- snps
  
  for (i in seq_along(snps)) {
    if (!keep[i]) next
    high_ld <- which(r2[i, ] > r2_threshold)
    high_ld <- setdiff(high_ld, i)
    keep[high_ld] <- FALSE
  }
  
  return(names(keep)[keep])
}

# Helper: prune a harmonised univariable data set, strongest variants first
prune_uni <- function(h, R, r2_threshold = r2_prune) {
  h <- h[h$mr_keep, ]
  h <- h[!duplicated(h$SNP), ]
  snps_ordered <- h$SNP[order(h$pval.exposure)]
  keep <- ld_prune_snps(snps_ordered, R, r2_threshold)
  h[h$SNP %in% keep, ]
}

# Helper: prune a harmonised multivariable data set, ordered by smallest p across exposures
prune_mv <- function(mvdat, R, r2_threshold = r2_prune) {
  snps  <- rownames(mvdat$exposure_beta)
  min_p <- apply(mvdat$exposure_pval, 1, min)
  keep  <- ld_prune_snps(snps[order(min_p)], R, r2_threshold)
  idx   <- match(keep, snps)                     # positions of kept SNPs
  
  mvdat$exposure_beta <- mvdat$exposure_beta[idx, , drop = FALSE]
  mvdat$exposure_se   <- mvdat$exposure_se[idx, , drop = FALSE]
  mvdat$exposure_pval <- mvdat$exposure_pval[idx, , drop = FALSE]
  mvdat$outcome_beta  <- mvdat$outcome_beta[idx]
  mvdat$outcome_se    <- mvdat$outcome_se[idx]
  mvdat$outcome_pval  <- mvdat$outcome_pval[idx]
  mvdat
}

df_results <- data.frame()
df_forest  <- data.frame()

IHD[["meta"]] <- subset(IHD[["meta"]], !duplicated(SNP))

for (i in c("KIV2", "isodom", "isosum")) {
  
  # Univariable Lp(a) (total effect)
  h_lpa <- prune_uni(harmonise_data(list_uni[["LPA"]], IHD[["meta"]]), R)
  data  <- dat_to_MRInput(h_lpa)
  data[[1]]@correlation <- R[data[[1]]@snps, data[[1]]@snps]
  data  <- MendelianRandomization::mr_ivw(data[[1]], correl = TRUE)
  res_lpa <- data.frame(b = data@Estimate, se = data@StdError, pval = data@Pvalue, nsnp = data@SNPs)
  
  # Univariable apo(a) measure
  h_apoa <- prune_uni(harmonise_data(list_uni[[i]], IHD[["meta"]]), R)
  data   <- dat_to_MRInput(h_apoa)
  data[[1]]@correlation <- R[data[[1]]@snps, data[[1]]@snps]
  data   <- MendelianRandomization::mr_ivw(data[[1]], correl = TRUE)
  res_kiv2 <- data.frame(b = data@Estimate, se = data@StdError, pval = data@Pvalue, nsnp = data@SNPs)
  
  # Multivariable
  mvdat <- mv_harmonise_data(list_comb[[i]], IHD[["meta"]], harmonise_strictness = 1)
  mvdat <- prune_mv(mvdat, R)
  
  snps_kept  <- rownames(mvdat$exposure_beta)
  R_filtered <- R[snps_kept, snps_kept, drop = FALSE]
  
  mv_input <- mr_mvinput(
    bx = mvdat$exposure_beta, by = mvdat$outcome_beta,
    bxse = mvdat$exposure_se, byse = mvdat$outcome_se,
    snps = snps_kept, exposure = mvdat$expname$exposure,
    correlation = R_filtered)
  res_mv <- mr_mvivw(mv_input)
  
  # Estimates
  exp_names   <- mv_input@exposure
  beta_lpa_mv <- res_mv@Estimate[exp_names == "Lipoprotein(a)"]
  se_lpa_mv   <- res_mv@StdError[exp_names == "Lipoprotein(a)"]
  
  med_lpa <- prop_mediated_ci(
    beta_total  = res_lpa$b,
    se_total    = res_lpa$se,
    beta_direct = beta_lpa_mv,
    se_direct   = se_lpa_mv
  )
  
  P_lpa  <- med_lpa$P
  CI_lpa <- med_lpa$CI
  p_lpa  <- med_lpa$p
  P_med  <- 100 - P_lpa
  cI_med <- 100 - CI_lpa
  
  q <- "Lipoprotein(a)"
  
  df_results[paste(i, q), "\n\n"] <- " "
  df_results[paste(i, q), "\n\nAnalysis"] <- paste0("Lipoprotein(a) + ",
                                                           if (grepl("KIV", i)) "Sum of KIV-2" else if (grepl("isodom", i)) "Dominant apo(a) isoform"
                                                           else "Sum of apo(a) isoforms")
  
  df_results[paste(i, q), "\nNo. of\nSNPs"] <-
    paste0(res_lpa$nsnp, " + ", formatC(res_kiv2$nsnp, format = "d", big.mark = ","))
  
  df_results[paste(i, q),
                    "Fraction of coronary artery\ndisease risk due to\nhigh lipoprotein(a) (95% CI)"] <-
    sprintf("%.f%% (%.f–%.f%%)", P_lpa, CI_lpa["lower"], CI_lpa["upper"])
  
  df_results[paste(i, q),
                    "Fraction of coronary artery\ndisease risk\ndue to short apo(a) (95% CI)"] <-
    sprintf("%.f%% (%.f–%.f%%)", P_med, cI_med["upper"], cI_med["lower"])
  
  df_forest[paste(i, q), c("Effect", "Lower 0.95", "Upper 0.95")] <-
    c(P_lpa, CI_lpa["lower"], CI_lpa["upper"])
  
  df_results[paste(i, q), "    "] <- paste(rep(" ", 20), collapse = " ")
  df_results[paste(i, q), "P for fraction\ndue to\nlipoprotein(a)"] <- "                   "
  df_forest[paste(i, q), "\nP"] <- p_lpa
}



library(gridExtra)
library(gtable)
library(grid)




plot_df <- data.frame(
  label    = rownames(df_results),
  est_not  = df_forest$Effect,          # NOT mediated
  est_med  = 100 - df_forest$Effect,    # mediated
  row      = seq_len(nrow(df_forest))
)


# Optional: grouping for alternating colours
plot_df$group <- ifelse(plot_df$row %% 2 == 1, "odd", "even")

#Function for bar plot
bar_grob <- function(value_not,
                     value_med,
                     xmin = 0,
                     xmax = 100,
                     col_not = "darkcyan",
                     col_med = "darkorchid4") {
  
  # Clamp
  value_not <- max(min(value_not, xmax), xmin)
  value_med <- max(min(value_med, xmax - value_not), xmin)
  
  w_not <- value_not / (xmax - xmin)
  w_med <- value_med / (xmax - xmin)
  
  gTree(children = gList(
    
    ## NOT mediated
    rectGrob(
      x = w_not / 2,
      y = 0.5,
      width  = w_not,
      height = 0.6,
      gp = gpar(fill = col_not, col = NA)
    ),
    
    ## Mediated (stacked to the right)
    rectGrob(
      x = w_not + w_med / 2,
      y = 0.5,
      width  = w_med,
      height = 0.6,
      gp = gpar(fill = col_med, col = NA)
    )
    
  ))
}


#Creating table


tbl <- tableGrob(
  df_results,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 9,
    core = list(
      fg_params = list(
        hjust = 0,        # left-align text
        x = 0.02          # small left padding
      )
    ),
    colhead = list(
      fg_params = list(
        hjust = 0,
        x = 0.02,
        fontface = "bold"
      )
    )
  )
)

##Adding barcols to table
bar_col <- 6 
tbl$widths[[bar_col]] <- unit(4, "cm")

n_bar <- nrow(plot_df)

for (i in seq_len(n_bar)) {
  
  g <- bar_grob(
    value_not = plot_df$est_not[i],
    value_med = plot_df$est_med[i]
  )
  
  idx <- which(
    tbl$layout$t == (i + 1) &
      tbl$layout$l == bar_col &
      grepl("-fg$", tbl$layout$name)
  )
  
  stopifnot(length(idx) == 1)
  
  tbl$grobs[[idx]] <- g
}


# Create x-axis grob (0–100)
x_axis <- xaxisGrob(
  at = seq(0, 1, by = 0.2),
  label = seq(0, 100, by = 20),
  gp = gpar(fontsize = 9)
)



# Add a new row at the bottom
tbl <- gtable_add_rows(
  tbl,
  heights =  unit(0, "mm")
)

# Index of the new row
axis_row <- length(tbl$heights)

# Add x-axis ONLY under the bar column
tbl <- gtable_add_grob(
  tbl,
  grobs = x_axis,
  t = axis_row,
  b = axis_row,   # ✅ THIS is the missing piece
  l = bar_col,
  r = bar_col,
  clip = "off"
)


#Adding veritcal line at 0 

first_data_row <- 2
last_data_row  <- first_data_row + nrow(plot_df) - 1


vline_0 <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(0, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(1, "npc"),
  gp = gpar(
    col = "black",
    lwd = 1,
    lty = c(3, 1)
  )
)



tbl <- gtable_add_grob(
  tbl,
  grobs = vline_0,
  t = first_data_row,
  b = last_data_row,
  l = bar_col,
  r = bar_col,
  clip = "off"
)


#Applying line under header

header_bottom_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(1, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(0, "npc"),
  gp = gpar(col = "black", lwd = 2)
)

tbl <- gtable_add_grob(
  tbl,
  grobs = header_bottom_line,
  t = 1,
  b = 1,
  l = 1,
  r = ncol(tbl),
  clip = "off"
)


#Adding other lines
header_bottom_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(1, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(0, "npc"),
  gp = gpar(col = "black", lwd = 1)
)


for (i in c(4)){
  tbl <- gtable_add_grob(
    tbl,
    grobs = header_bottom_line,
    t = i,
    b = 1,
    l = 1,
    r = ncol(tbl),
    clip = "off"
  )
}




format_p_expr <- function(p, digits = 2) {
  if (is.na(p)) return(expression(NA))
  
  if (p >0.001) {
    return(signif(p,2))
  }
  
  expo <- floor(log10(p))
  base <- signif(p / 10^expo,1)
  base <- round(base, digits)
  
  # Return plotmath expression: base %*% 10^expo
  bquote(.(base)~ "×" ~ 10^.(expo))
}

p_exprs <- lapply(df_forest[,"\nP"], format_p_expr)


p_col <- 7

stopifnot(p_col <= ncol(tbl))

n_p <- length(p_exprs)

for (i in seq_len(n_p)) {
  
  p_grob <- textGrob(
    label = p_exprs[[i]],
    gp = gpar(fontsize = 9),
    hjust = 0,                 # left-aligned
    x = unit(0.02, "npc")         # 🔴 critical for superscripts
  )
  
  idx <- which(
    tbl$layout$t == (i + 1) &  # +1 because row 1 is header
      tbl$layout$l == p_col &
      grepl("-fg$", tbl$layout$name)
  )
  
  stopifnot(length(idx) >= 1)
  
  tbl$grobs[[idx[1]]] <- p_grob
}


xlab_grob <- textGrob(
  "\n\nFraction of coronary artery disease risk\nCyan = directly due to elevated lipoprotein(a)\nPurple = indirectly due to short apo(a)",
  gp = gpar(fontsize = 10),
  just = "center"
)


tbl <- gtable_add_rows(
  tbl,
  heights = unit(22, "mm")   # dedicated label row
)

label_row <- length(tbl$heights)


tbl <- gtable_add_grob(
  tbl,
  grobs = xlab_grob,
  t = label_row,
  b = label_row,
  l = bar_col,
  r = bar_col,
  clip = "off"
)

# Shade background of the middle data row only

middle_row <- 2                    # data row index
table_row  <- middle_row + 1        # +1 because row 1 is header

bg_grob <- rectGrob(
  gp = gpar(fill = "grey90", col = NA)
)

tbl <- gtable_add_grob(
  tbl,
  grobs = bg_grob,
  t = table_row,
  b = table_row,
  l = 1,
  r = ncol(tbl),
  z = -Inf,        # 🔑 ensures background is BEHIND text and bars
  clip = "off"
)

tbl$widths[[1]] <- unit(0, "cm")

grid.newpage()
grid.draw(tbl)

doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(tbl)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")



# Fraction mediated by kiv-2 and isoform using variants in low LD and with Q-statistic minimisation --------------------------------------------------


#Corrlation of lipoprotein(a) with other exposures
rho_val<-c(KIV2=-0.4691151, isodom=-0.5102852, isosum=-0.4576141)


prop_mediated_ci <- function(beta_total, se_total,
                             beta_direct, se_direct,
                             alpha = 0.05) {
  
  # Proportion NOT mediated (direct effect)
  Q <- beta_direct / beta_total
  
  # Delta-method variance
  var_Q <- ( (beta_direct / beta_total^2)^2 * se_total^2 ) +
    ( (1 / beta_total)^2 * se_direct^2 )
  
  se_Q <- sqrt(var_Q)
  
  # Wald Z
  z_val <- Q / se_Q
  
  # Numerically stable two-sided p-value
  log_p <- log(2) + pnorm(abs(z_val),
                          lower.tail = FALSE,
                          log.p = TRUE)
  p_val <- exp(log_p)
  
  # Confidence interval
  z_alpha <- qnorm(1 - alpha / 2)
  CI <- c(
    lower = Q - z_alpha * se_Q,
    upper = Q + z_alpha * se_Q
  )
  
  # Convert to percent and cap at 100
  CI_pct <- pmin(CI * 100, 100)
  
  list(
    P    = Q * 100,
    SE   = se_Q,
    CI   = CI_pct,
    Z    = z_val,
    p    = p_val,
    logp = log_p
  )
}


# Sensitivity analysis: LD-pruned instruments (r2 < 0.01), delta method 

r2_prune <- 0.01

ld_prune_snps <- function(snps, R, r2_threshold = r2_prune) {
  
  R_sub <- R[snps, snps, drop = FALSE]
  r2 <- R_sub^2
  
  keep <- rep(TRUE, length(snps))
  names(keep) <- snps
  
  for (i in seq_along(snps)) {
    if (!keep[i]) next
    high_ld <- which(r2[i, ] > r2_threshold)
    high_ld <- setdiff(high_ld, i)
    keep[high_ld] <- FALSE
  }
  
  return(names(keep)[keep])
}

# Helper: prune a harmonised univariable data set, strongest variants first
prune_uni <- function(h, R, r2_threshold = r2_prune) {
  h <- h[h$mr_keep, ]
  h <- h[!duplicated(h$SNP), ]
  snps_ordered <- h$SNP[order(h$pval.exposure)]
  keep <- ld_prune_snps(snps_ordered, R, r2_threshold)
  h[h$SNP %in% keep, ]
}

# Helper: prune a harmonised multivariable data set, ordered by smallest p across exposures
prune_mv <- function(mvdat, R, r2_threshold = r2_prune) {
  snps  <- rownames(mvdat$exposure_beta)
  min_p <- apply(mvdat$exposure_pval, 1, min)
  keep  <- ld_prune_snps(snps[order(min_p)], R, r2_threshold)
  idx   <- match(keep, snps)                     # positions of kept SNPs
  
  mvdat$exposure_beta <- mvdat$exposure_beta[idx, , drop = FALSE]
  mvdat$exposure_se   <- mvdat$exposure_se[idx, , drop = FALSE]
  mvdat$exposure_pval <- mvdat$exposure_pval[idx, , drop = FALSE]
  mvdat$outcome_beta  <- mvdat$outcome_beta[idx]
  mvdat$outcome_se    <- mvdat$outcome_se[idx]
  mvdat$outcome_pval  <- mvdat$outcome_pval[idx]
  mvdat
}

library(MVMR)
df_results <- data.frame()
df_forest  <- data.frame()

IHD[["meta"]] <- subset(IHD[["meta"]], !duplicated(SNP))

for (i in c("KIV2", "isodom", "isosum")) {
  
  # Univariable Lp(a) (total effect)
  h_lpa <- prune_uni(harmonise_data(list_uni[["LPA"]], IHD[["meta"]]), R)
  data  <- dat_to_MRInput(h_lpa)
  data[[1]]@correlation <- R[data[[1]]@snps, data[[1]]@snps]
  data  <- MendelianRandomization::mr_ivw(data[[1]], correl = TRUE)
  res_lpa <- data.frame(b = data@Estimate, se = data@StdError, pval = data@Pvalue, nsnp = data@SNPs)
  
  # Univariable apo(a) measure
  h_apoa <- prune_uni(harmonise_data(list_uni[[i]], IHD[["meta"]]), R)
  data   <- dat_to_MRInput(h_apoa)
  data[[1]]@correlation <- R[data[[1]]@snps, data[[1]]@snps]
  data   <- MendelianRandomization::mr_ivw(data[[1]], correl = TRUE)
  res_kiv2 <- data.frame(b = data@Estimate, se = data@StdError, pval = data@Pvalue, nsnp = data@SNPs)
  
  # Multivariable
  mvdat <- mv_harmonise_data(list_comb[[i]], IHD[["meta"]], harmonise_strictness = 1)
  mvdat <- prune_mv(mvdat, R)
  
  mvmr_in <- format_mvmr(
    BXGs   = mvdat$exposure_beta,
    BYG    = mvdat$outcome_beta,
    seBXGs = mvdat$exposure_se,
    seBYG  = mvdat$outcome_se,
    RSID   = rownames(mvdat$exposure_beta)
  )
  
  pcor <- matrix(c(1, rho_val[i], rho_val[i], 1), 2, 2)
  
  res_q <- qhet_mvmr(mvmr_in, pcor = pcor, CI = TRUE, iterations = 1000)
  
  ci_num <- t(sapply(strsplit(res_q$`95% CI`, "(?<=[0-9])-", perl = TRUE), as.numeric))
  colnames(ci_num) <- c("lower", "upper")
  
  res_q$se <- (ci_num[, "upper"] - ci_num[, "lower"]) / (2 * qnorm(0.975))
  
  
  # Estimates
  exp_names   <- mv_input@exposure
  beta_lpa_mv <-  res_q$`Effect Estimates`[1]
  se_lpa_mv   <- res_q$se[1]
  
  med_lpa <- prop_mediated_ci(
    beta_total  = res_lpa$b,
    se_total    = res_lpa$se,
    beta_direct = beta_lpa_mv,
    se_direct   = se_lpa_mv
  )
  
  P_lpa  <- med_lpa$P
  CI_lpa <- med_lpa$CI
  p_lpa  <- med_lpa$p
  P_med  <- 100 - P_lpa
  cI_med <- 100 - CI_lpa
  
  q <- "Lipoprotein(a)"
  
  df_results[paste(i, q), "\n\n"] <- " "
  df_results[paste(i, q), "\n\nAnalysis"] <- paste0("Lipoprotein(a) + ",
                                                    if (grepl("KIV", i)) "Sum of KIV-2" else if (grepl("isodom", i)) "Dominant apo(a) isoform"
                                                    else "Sum of apo(a) isoforms")
  
  df_results[paste(i, q), "\nNo. of\nSNPs"] <-
    paste0(res_lpa$nsnp, " + ", formatC(res_kiv2$nsnp, format = "d", big.mark = ","),
           " (", length(snps_kept), ")")
  
  df_results[paste(i, q),
             "Fraction of coronary artery\ndisease risk due to\nhigh lipoprotein(a) (95% CI)"] <-
    sprintf("%.f%% (%.f–%.f%%)", P_lpa, CI_lpa["lower"], CI_lpa["upper"])
  
  df_results[paste(i, q),
             "Fraction of coronary artery\ndisease risk\ndue to short apo(a) (95% CI)"] <-
    sprintf("%.f%% (%.f–%.f%%)", P_med, cI_med["upper"], cI_med["lower"])
  
  df_forest[paste(i, q), c("Effect", "Lower 0.95", "Upper 0.95")] <-
    c(P_lpa, CI_lpa["lower"], CI_lpa["upper"])
  
  df_results[paste(i, q), "    "] <- paste(rep(" ", 20), collapse = " ")
  df_results[paste(i, q), "P for fraction\ndue to\nlipoprotein(a)"] <- "                   "
  df_forest[paste(i, q), "\nP"] <- p_lpa
}



library(gridExtra)
library(gtable)
library(grid)




plot_df <- data.frame(
  label    = rownames(df_results),
  est_not  = df_forest$Effect,          # NOT mediated
  est_med  = 100 - df_forest$Effect,    # mediated
  row      = seq_len(nrow(df_forest))
)


# Optional: grouping for alternating colours
plot_df$group <- ifelse(plot_df$row %% 2 == 1, "odd", "even")

#Function for bar plot
bar_grob <- function(value_not,
                     value_med,
                     xmin = 0,
                     xmax = 100,
                     col_not = "darkcyan",
                     col_med = "darkorchid4") {
  
  # Clamp
  value_not <- max(min(value_not, xmax), xmin)
  value_med <- max(min(value_med, xmax - value_not), xmin)
  
  w_not <- value_not / (xmax - xmin)
  w_med <- value_med / (xmax - xmin)
  
  gTree(children = gList(
    
    ## NOT mediated
    rectGrob(
      x = w_not / 2,
      y = 0.5,
      width  = w_not,
      height = 0.6,
      gp = gpar(fill = col_not, col = NA)
    ),
    
    ## Mediated (stacked to the right)
    rectGrob(
      x = w_not + w_med / 2,
      y = 0.5,
      width  = w_med,
      height = 0.6,
      gp = gpar(fill = col_med, col = NA)
    )
    
  ))
}


#Creating table


tbl <- tableGrob(
  df_results,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 9,
    core = list(
      fg_params = list(
        hjust = 0,        # left-align text
        x = 0.02          # small left padding
      )
    ),
    colhead = list(
      fg_params = list(
        hjust = 0,
        x = 0.02,
        fontface = "bold"
      )
    )
  )
)

##Adding barcols to table
bar_col <- 6 
tbl$widths[[bar_col]] <- unit(4, "cm")

n_bar <- nrow(plot_df)

for (i in seq_len(n_bar)) {
  
  g <- bar_grob(
    value_not = plot_df$est_not[i],
    value_med = plot_df$est_med[i]
  )
  
  idx <- which(
    tbl$layout$t == (i + 1) &
      tbl$layout$l == bar_col &
      grepl("-fg$", tbl$layout$name)
  )
  
  stopifnot(length(idx) == 1)
  
  tbl$grobs[[idx]] <- g
}


# Create x-axis grob (0–100)
x_axis <- xaxisGrob(
  at = seq(0, 1, by = 0.2),
  label = seq(0, 100, by = 20),
  gp = gpar(fontsize = 9)
)



# Add a new row at the bottom
tbl <- gtable_add_rows(
  tbl,
  heights =  unit(0, "mm")
)

# Index of the new row
axis_row <- length(tbl$heights)

# Add x-axis ONLY under the bar column
tbl <- gtable_add_grob(
  tbl,
  grobs = x_axis,
  t = axis_row,
  b = axis_row,   # ✅ THIS is the missing piece
  l = bar_col,
  r = bar_col,
  clip = "off"
)


#Adding veritcal line at 0 

first_data_row <- 2
last_data_row  <- first_data_row + nrow(plot_df) - 1


vline_0 <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(0, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(1, "npc"),
  gp = gpar(
    col = "black",
    lwd = 1,
    lty = c(3, 1)
  )
)



tbl <- gtable_add_grob(
  tbl,
  grobs = vline_0,
  t = first_data_row,
  b = last_data_row,
  l = bar_col,
  r = bar_col,
  clip = "off"
)


#Applying line under header

header_bottom_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(1, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(0, "npc"),
  gp = gpar(col = "black", lwd = 2)
)

tbl <- gtable_add_grob(
  tbl,
  grobs = header_bottom_line,
  t = 1,
  b = 1,
  l = 1,
  r = ncol(tbl),
  clip = "off"
)


#Adding other lines
header_bottom_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(1, "npc"),
  y0 = unit(0, "npc"),
  y1 = unit(0, "npc"),
  gp = gpar(col = "black", lwd = 1)
)


for (i in c(4)){
  tbl <- gtable_add_grob(
    tbl,
    grobs = header_bottom_line,
    t = i,
    b = 1,
    l = 1,
    r = ncol(tbl),
    clip = "off"
  )
}




format_p_expr <- function(p, digits = 2) {
  if (is.na(p)) return(expression(NA))
  
  if (p >0.001) {
    return(signif(p,2))
  }
  
  expo <- floor(log10(p))
  base <- signif(p / 10^expo,1)
  base <- round(base, digits)
  
  # Return plotmath expression: base %*% 10^expo
  bquote(.(base)~ "×" ~ 10^.(expo))
}

p_exprs <- lapply(df_forest[,"\nP"], format_p_expr)


p_col <- 7

stopifnot(p_col <= ncol(tbl))

n_p <- length(p_exprs)

for (i in seq_len(n_p)) {
  
  p_grob <- textGrob(
    label = p_exprs[[i]],
    gp = gpar(fontsize = 9),
    hjust = 0,                 # left-aligned
    x = unit(0.02, "npc")         # 🔴 critical for superscripts
  )
  
  idx <- which(
    tbl$layout$t == (i + 1) &  # +1 because row 1 is header
      tbl$layout$l == p_col &
      grepl("-fg$", tbl$layout$name)
  )
  
  stopifnot(length(idx) >= 1)
  
  tbl$grobs[[idx[1]]] <- p_grob
}


xlab_grob <- textGrob(
  "\n\nFraction of coronary artery disease risk\nCyan = directly due to elevated lipoprotein(a)\nPurple = indirectly due to short apo(a)",
  gp = gpar(fontsize = 10),
  just = "center"
)


tbl <- gtable_add_rows(
  tbl,
  heights = unit(22, "mm")   # dedicated label row
)

label_row <- length(tbl$heights)


tbl <- gtable_add_grob(
  tbl,
  grobs = xlab_grob,
  t = label_row,
  b = label_row,
  l = bar_col,
  r = bar_col,
  clip = "off"
)

# Shade background of the middle data row only

middle_row <- 2                    # data row index
table_row  <- middle_row + 1        # +1 because row 1 is header

bg_grob <- rectGrob(
  gp = gpar(fill = "grey90", col = NA)
)

tbl <- gtable_add_grob(
  tbl,
  grobs = bg_grob,
  t = table_row,
  b = table_row,
  l = 1,
  r = ncol(tbl),
  z = -Inf,        # 🔑 ensures background is BEHIND text and bars
  clip = "off"
)

tbl$widths[[1]] <- unit(0, "cm")

grid.newpage()
grid.draw(tbl)

doc <- read_pptx("widescreen.pptx")
doc<-add_slide(doc)
doc <- ph_with(x = doc, value = rvg::dml(ggobj =ggplotify::as.ggplot(tbl)), 
               location = ph_location_fullsize(newlabel = "") )

print(doc,target="temp.pptx")


