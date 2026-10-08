#########################################################################################
# 先整理数据。mvmr做1个协变量和n个协变量，用SNP的并集。都是UKB数据
# BMI的5e-8.Rdata做clump,某个协变量的5e-8.Rdata做clump,取并集之后再做clump
# 随后再去原gwas里找到数据,整理成mvmr数据格式
# 命名规则：例如BMI和HDL做MVMR，就是先读取BMI处理好的Rdata，之后处理HDL，处理好后命名为BMI_1COV_HDL.Rdata
# 第一步：处理好每一个cov的5e-8及clump的Rata
# 第二步：处理好把bmi的snp和处理好cov的snp提取出来取并集再做clump
# 第三步：把做完clump的并集snp去bmi原始数据和cov原始数据中查找snp，输出命名为BMI_1COV_HDL.Rdata
# 第四步：之后再做mvmr
# step1:load进来都是covariate_clumped
# step3:load进来都是mvmr_data_combined，但这个data是竖表
# 日期: 2026-01-20
#########################################################################################

#########################################################################################
# Step 1: 批量处理协变量 (Filter P < 5e-8 & Clump)
# 对应任务: 处理好每一个cov的5e-8及clump的Rata
# load进来都是covariate_clumped
#########################################################################################

library(data.table)
library(ieugwasr)
library(TwoSampleMR)
library(stringr)

# ================= 1. 路径设置 (请根据实际情况修改) =================
# 原始 GWAS 数据文件夹
cov_dir <- Sys.getenv("ERACA_CANDIDATE_DIR", unset = file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "cov_data"))

# variants1.Rdata 路径 (必须先加载，用于cbind)
variants_path <- Sys.getenv("ERACA_VARIANTS_RDATA", unset = file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "variants1.Rdata"))

# PLINK 设置
plink_bin <- Sys.getenv("ERACA_PLINK", unset = unname(Sys.which("plink")))
bfile_path <- Sys.getenv("ERACA_LD_BFILE", unset = file.path("controlled_data", "ld_reference", "reference"))

# 结果保存路径 (建议新建一个文件夹存放第一步的结果)
save_dir <- file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "Step1_cov_clumped_data")
if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# ================= 2. 建立 ID 与 变量名 的映射关系 =================
# Neale Lab UK Biobank Round 2 phenotype IDs used for the candidate-trait files.
trait_map <- list(
  "21001_irnt" = "BMI",
  "48_irnt"    = "Waist_Circumference",
  "49_irnt"    = "Hip_Circumference",
  "30690_irnt" = "Total_Cholesterol",
  "30760_irnt" = "HDL",
  "30780_irnt" = "LDL",
  "30870_irnt" = "Triglycerides",
  "30740_irnt" = "Fasting_Glucose",
  "30750_irnt" = "HbA1c",
  "4080_irnt"  = "Systolic_BP",
  "4079_irnt"  = "Diastolic_BP",  
  "20116_2"    = "Smoking",       
  "20117_2"    = "Alcohol"
)

# ================= 3. 加载辅助数据 =================
cat(">>> Loading variants info...\n")
load(variants_path) 
# 确保 variants1 存在于环境中
if(!exists("variants1")) stop("variants1.Rdata 加载失败或变量名不对，请检查！")

# 获取所有 .bgz 文件
files <- list.files(cov_dir, full.names = TRUE, pattern = "\\.bgz$")
cat(">>> Found", length(files), "files to process.\n")

# ================= 4. 循环处理 =================

for (f in files) {
  
  # --- 4.1 识别变量名 ---
  fname <- basename(f)
  # 提取 ID 部分 (去掉 .gwas.imputed...)
  file_id <- str_split(fname, "\\.gwas")[[1]][1]
  
  # 从映射表中获取名称，如果没有匹配到，就用文件名 ID
  if (file_id %in% names(trait_map)) {
    trait_name <- trait_map[[file_id]]
  } else {
    trait_name <- file_id
    warning(paste("未在映射表中找到", file_id, "使用原始ID作为名称"))
  }
  
  cat("\n--------------------------------------------------------------\n")
  cat("Processing:", trait_name, "(File:", fname, ")\n")
  
  # --- 4.2 读取数据并合并 Variant ---
  # 使用 fread 读取
  dat <- fread(f)
  
  # Append variant annotations; input rows must be aligned.
  # 注意：variant1[,-1] 意味着去掉第一列，请确保逻辑与你之前做BMI时一致
  if(nrow(dat) != nrow(variants1)){
    warning(paste("Skipping", trait_name, "- Row count mismatch between data and variants!"))
    next
  }
  dat <- cbind(dat, variants1[,-1])
  
  # --- 4.3 筛选 P < 5e-8 ---
  dat$pval <- as.numeric(dat$pval)
  dat_P <- subset(dat, pval < 5e-8)
  dat_P <- as.data.frame(dat_P)
  
  cat("  Significant SNPs (P<5e-8):", nrow(dat_P), "\n")
  
  if (nrow(dat_P) == 0) {
    cat("  No significant SNPs found. Skipping clumping.\n")
    next
  }
  
  # --- 4.4 Clumping ---
  exp_iv <- dat_P[, c("rsid", "pval")]
  colnames(exp_iv) <- c("rsid", "pval")
  
  # 定义默认结果为原始筛选数据 (防错)
  dat_clumped_res <- dat_P
  
  if (nrow(exp_iv) > 1) {
    cat("  Clumping...\n")
    a_clump <- try(ld_clump_local(
      dat = exp_iv, 
      clump_kb = 10000, 
      clump_r2 = 0.001, 
      clump_p = 1,
      bfile = bfile_path,
      plink_bin = plink_bin
    ), silent = TRUE)
    
    if (!inherits(a_clump, "try-error") && !is.null(a_clump) && nrow(a_clump) > 0) {
      dat_clumped_res <- dat_P[dat_P$rsid %in% a_clump$rsid, ]
      cat("  Clumping finished. SNPs left:", nrow(dat_clumped_res), "\n")
    } else {
      cat("  Clumping failed or returned no result. Keeping original significant SNPs.\n")
    }
  }
  
  # --- 4.5 Format Data (TwoSampleMR 格式) ---
  final_dat <- format_data(dat_clumped_res,
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
                           eaf_col = "minor_AF") # Allele-frequency column in the source GWAS.
  
  final_dat$exposure <- trait_name
  
  # --- 4.6 保存 RData ---
  # 命名格式: 变量名_clumped.RData (例如 HDL_clumped.RData)
  out_file <- paste0(save_dir, trait_name, "_clumped.RData")
  
  # 为了后续方便，我们将对象名统一改为 covariate_clumped 再保存，
  # 或者直接保存 final_dat 对象
  covariate_clumped <- final_dat
  save(covariate_clumped, file = out_file)
  
  cat("  Saved to:", out_file, "\n")
  
  # 清理内存
  rm(dat, dat_P, dat_clumped_res, final_dat, covariate_clumped)
  gc()
}

cat("\n>>> Step 1 Complete! All covariates processed.\n")


#########################################################################################
# MVMR Step 2 & 3: Union, Re-clump, and Extract Original Data
# 日期: 2026-01-20
#########################################################################################

library(data.table)
library(ieugwasr)
library(TwoSampleMR)
library(dplyr)
library(stringr)

# ================= 1. 路径设置 =================
# 1.1 基础路径
base_dir <- file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"))
cov_raw_dir <- paste0(base_dir, "cov_data/")
plink_bin <- Sys.getenv("ERACA_PLINK", unset = unname(Sys.which("plink")))
bfile_path <- Sys.getenv("ERACA_LD_BFILE", unset = file.path("controlled_data", "ld_reference", "reference"))

# 1.2 输入文件
variants_path <- paste0(base_dir, "variants1.Rdata")
bmi_clump_path <- paste0(base_dir, "BMI_10000_0.0001.RData") # BMI已处理好的Clump数据
bmi_raw_file <- paste0(base_dir, "cov_data/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz")

# 1.3 第一步生成的协变量 Clump 文件夹
step1_dir <- paste0(base_dir, "Step1_cov_clumped_data/")

# 1.4 输出路径 (Step 3 结果)
output_dir <- paste0(base_dir, "Step3_mvmr_ready_data/")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ================= 2. 准备工作：建立映射与加载 BMI =================

# 2.1 变量名 -> 原始文件名 ID 的映射表 (用于Step 3反向查找原始数据)
# 注意：这里需要根据第一步的命名反推。
name_to_id_map <- list(
  "Waist_Circumference" = "48_irnt",
  "Hip_Circumference" = "49_irnt",
  "Total_Cholesterol" = "30690_irnt",
  "HDL" = "30760_irnt",
  "LDL" = "30780_irnt",
  "Triglycerides" = "30870_irnt",
  "Fasting_Glucose" = "30740_irnt",
  "HbA1c" = "30750_irnt",
  "Systolic_BP" = "4080_irnt",
  "Diastolic_BP" = "4079_irnt",
  "Smoking" = "20116_2",
  "Smoking_Prev" = "20116_0",
  "Alcohol" = "20117_2"
)

# 2.2 加载 Variants 信息 (只需加载一次)
cat(">>> Loading variants info...\n")
load(variants_path) # 加载 variants1
if(!exists("variants1")) stop("variants1 加载失败")

# 2.3 加载 BMI 原始数据 (只加载一次到内存，提速关键)
cat(">>> Loading Raw BMI Data into memory (Please wait)...\n")
bmi_raw_dt <- fread(bmi_raw_file)
# 合并 Variants 得到 rsid (按照你之前的逻辑)
if(nrow(bmi_raw_dt) == nrow(variants1)){
  bmi_raw_dt <- cbind(bmi_raw_dt, variants1[,-1])
} else {
  stop("BMI Raw Data 行数与 variants1 不匹配！")
}
cat("    BMI Raw Data Loaded. Rows:", nrow(bmi_raw_dt), "\n")

# 2.4 加载 BMI Clumped 数据
load(bmi_clump_path) # 加载 BMI_P_clump
# 统一列名方便合并 (保证有 rsid 和 pval)
bmi_ivs <- BMI_P_clump[, c("SNP", "pval.exposure")]
colnames(bmi_ivs) <- c("rsid", "pval")


# ================= 3. 循环处理每个协变量 =================

# 获取第一步生成的 RData 文件
cov_files <- list.files(step1_dir, pattern = "_clumped\\.RData$", full.names = TRUE)

for (f in cov_files) {
  
  # --- 3.1 准备环境 ---
  fname <- basename(f)
  trait_name <- gsub("_clumped\\.RData", "", fname) # 提取变量名，例如 HDL
  
  cat("\n=======================================================\n")
  cat("Processing:", trait_name, "\n")
  
  # 加载协变量的 Clump 数据
  # 注意：第一步保存的对象名可能是 covariate_clumped 或 final_dat，这里动态获取
  env_temp <- new.env()
  load(f, envir = env_temp)
  cov_obj_name <- ls(env_temp)[1]
  cov_clumped <- env_temp[[cov_obj_name]]
  
  
  # --- Step 2: 取并集 & Re-clump ---
  cat(">>> Step 2: Union and Re-clumping...\n")
  
  # 提取协变量 IVs
  cov_ivs <- cov_clumped[, c("SNP", "pval.exposure")]
  colnames(cov_ivs) <- c("rsid", "pval")
  
  # 合并 BMI 和 协变量 的 IVs
  combined_ivs <- rbind(bmi_ivs, cov_ivs)
  
  # *关键逻辑*: 
  # MVMR Re-clump 需要保留该 SNP 在所有暴露中最小的 P 值
  # 先按 P 值排序，然后去重 (keep first = keep smallest P)
  combined_ivs <- combined_ivs[order(combined_ivs$pval), ]
  combined_ivs_unique <- combined_ivs[!duplicated(combined_ivs$rsid), ]
  
  cat("    Total unique SNPs before re-clump:", nrow(combined_ivs_unique), "\n")
  
  # 执行 Re-clump
  final_ivs <- try(ld_clump_local(
    dat = combined_ivs_unique, 
    clump_kb = 10000, 
    clump_r2 = 0.001, 
    clump_p = 1, # 设为1，因为已经筛选过了，只做去连锁
    bfile = bfile_path,
    plink_bin = plink_bin
  ), silent = TRUE)
  
  if (inherits(final_ivs, "try-error") || is.null(final_ivs) || nrow(final_ivs) == 0) {
    warning("Re-clumping failed for ", trait_name)
    next
  }
  
  target_snps <- final_ivs$rsid
  cat("    Final Target SNPs after Re-clump:", length(target_snps), "\n")
  
  
  # --- Step 3: 回溯原始数据提取信息 ---
  cat(">>> Step 3: Extracting from Raw Data...\n")
  
  # A. 从 BMI 原始数据提取 (利用内存中的 bmi_raw_dt)
  # ------------------------------------------------
  cat("    Extracting from BMI Raw...\n")
  bmi_extracted <- bmi_raw_dt[bmi_raw_dt$rsid %in% target_snps, ]
  bmi_extracted <- as.data.frame(bmi_extracted)
  
  # 格式化 BMI
  bmi_formatted <- format_data(bmi_extracted, type = "exposure",
                               snp_col = "rsid", beta_col = "beta", se_col = "se",
                               effect_allele_col = "alt", other_allele_col = "ref",
                               pval_col = "pval", eaf_col = "minor_AF",
                               samplesize_col = "n_complete_samples",
                               chr_col = "chr", pos_col = "pos")
  bmi_formatted$exposure <- "BMI"
  
  
  # B. 从 协变量 原始数据提取 (需要读取对应文件)
  # ------------------------------------------------
  # 1. 找到原始文件名
  if (trait_name %in% names(name_to_id_map)) {
    raw_id <- name_to_id_map[[trait_name]]
    # 在原始文件夹中模糊匹配包含该ID的文件
    raw_file_pattern <- list.files(cov_raw_dir, pattern = raw_id, full.names = TRUE)
    # 确保只匹配到一个 .bgz 文件
    target_raw_file <- raw_file_pattern[grep("\\.bgz$", raw_file_pattern)][1]
  } else {
    warning("无法在映射表中找到 ", trait_name, " 的原始文件ID，跳过 extraction")
    next
  }
  
  if (is.na(target_raw_file) || !file.exists(target_raw_file)) {
    warning("原始文件不存在: ", raw_id)
    next
  }
  
  cat("    Reading Covariate Raw File:", basename(target_raw_file), "...\n")
  
  # 2. 读取并合并 (这里只能读一遍，内存吃紧的话可以在循环末尾 gc)
  cov_raw_dt <- fread(target_raw_file)
  if(nrow(cov_raw_dt) == nrow(variants1)){
    cov_raw_dt <- cbind(cov_raw_dt, variants1[,-1])
  } else {
    warning("行数不匹配，跳过 ", trait_name)
    next
  }
  
  # 3. 提取
  cov_extracted <- cov_raw_dt[cov_raw_dt$rsid %in% target_snps, ]
  cov_extracted <- as.data.frame(cov_extracted)
  
  # 4. 格式化 Covariate
  cov_formatted <- format_data(cov_extracted, type = "exposure",
                               snp_col = "rsid", beta_col = "beta", se_col = "se",
                               effect_allele_col = "alt", other_allele_col = "ref",
                               pval_col = "pval", eaf_col = "minor_AF",
                               samplesize_col = "n_complete_samples",
                               chr_col = "chr", pos_col = "pos")
  cov_formatted$exposure <- trait_name
  
  
  # --- 保存结果 ---
  # 将两个处理好的 exposure 数据合并保存 (MVMR 分析时可以直接提取)
  mvmr_data_combined <- rbind(bmi_formatted, cov_formatted)
  
  out_name <- paste0(output_dir, "BMI_1COV_", trait_name, ".RData")
  save(mvmr_data_combined, file = out_name)
  
  cat("    Saved:", out_name, "\n")
  
  # 清理循环内的临时大对象 (保留 bmi_raw_dt 和 variants1)
  rm(cov_raw_dt, cov_extracted, cov_formatted, bmi_formatted, combined_ivs_unique)
  gc() 
}

cat("\n>>> All steps finished! Files are ready in:", output_dir, "\n")
