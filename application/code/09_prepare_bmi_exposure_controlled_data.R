library(MVMR)
library(ieugwasr)
library(TwoSampleMR)

exp_data <- list.files(Sys.getenv("ERACA_CANDIDATE_DIR", unset = file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "cov_data")),full.names = TRUE)
load(Sys.getenv("ERACA_VARIANTS_RDATA", unset = file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "variants1.Rdata")))
BMI <- fread(Sys.getenv("ERACA_BMI_GWAS", unset = file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "exp_data", "21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz")))
BMI <- cbind(BMI,variants1[,-1])

# clump -------------------------------------------------------------------------------
BMI$pval <- as.numeric(BMI$pval)
BMI_P <- subset(BMI,pval < 5e-8)
BMI_P <- as.data.frame(BMI_P)

exp_iv <- BMI_P[, c("rsid", "pval")]
colnames(exp_iv) <- c("rsid", "pval")

if (nrow(exp_iv) > 1) {
  a_clump <- try(ld_clump_local(
    dat = exp_iv, clump_kb = 10000, clump_r2 = 0.001, clump_p=1,
    bfile = Sys.getenv("ERACA_LD_BFILE", unset = file.path("controlled_data", "ld_reference", "reference")),
    plink_bin = Sys.getenv("ERACA_PLINK", unset = unname(Sys.which("plink")))
  ), silent = TRUE)
  
  if (!inherits(a_clump, "try-error") && !is.null(a_clump) && nrow(a_clump) > 0) {
    BMI_P_clump <- BMI_P[BMI_P$rsid %in% a_clump$rsid, ]
    
  } else {
    BMI_P_clump <- BMI_P
  }
}

BMI_P_clump <- format_data(BMI_P_clump,
                           type = 'exposure',
                           snp_col = "rsid",
                           beta_col = "beta",
                           se_col = "se",
                           effect_allele_col = "alt",
                           other_allele_col = "ref",
                           pval_col = "pval",
                           samplesize_col = "n_complete_samples",
                           chr_col = "chr",
                           pos_col = "pos",
                           eaf_col = "minor_AF")
head(BMI_P_clump)
BMI_P_clump$exposure <- c('BMI')
head(BMI_P_clump)

save(BMI_P_clump, file = file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "BMI_10000_0.0001.RData"))
