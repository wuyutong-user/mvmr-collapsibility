#########################################################################################
# 根据整理好的bmi_cov数据，来做mvmr。
# 读取的out数据还是my_out_dir <- "controlled_data/outcome_gwas"地址下的，还是需要逐条chr处理
# 第四步：之后再做mvmr
# load进来都是mvmr_data_combined，但这个data是竖表，转换完横表再做mvmr。
# 日期: 2026-01-20
#########################################################################################

#########################################################################################
# MVMR Analysis Step 4 (Optimized): Run MVMR (BMI + 1 Covariate)
# 优化逻辑: Out层在最外 -> 预先提取所有Exp的SNP并集 -> 读一次Out做多次MVMR
# 日期: 2026-01-20
#########################################################################################

library(TwoSampleMR)
library(MVMR)
library(data.table)
library(stringr)

# =======================================================================================
# 1. 路径设置
# =======================================================================================

# 暴露数据路径 (BMI + Covariate 的 RData)
exp_data_dir <- file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "Step3_mvmr_ready_data")

# Outcome 原始数据路径 (ZIP)
outcome_dir_base <- Sys.getenv("ERACA_OUTCOME_DIR", unset = file.path("controlled_data", "outcome_gwas"))

# 结果保存路径
save_dir <- file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "MVMR_1COV_RESULT")
if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# 获取 Outcome ZIP 文件列表
outcome_files <- list.files(outcome_dir_base, pattern = "\\.zip$", full.names = FALSE)
# outcome_files <- outcome_files[1:2] # 测试用

# =======================================================================================
# 2. 优化后的核心分析函数
# =======================================================================================

run_mvmr_optimized <- function(exp_dir, out_dir, out_files, output_path) {
  
  # --- 1. 获取所有暴露文件 ---
  exp_files <- list.files(exp_dir, pattern = "\\.RData$", full.names = TRUE)
  if(length(exp_files) == 0) stop("未找到 Exposure RData 文件！")
  
  
  # --- 2. 预处理：构建全局 SNP 索引 (Union of all SNPs) ---
  cat(">>> Step 0: Pre-scanning Exposure files to build Global SNP List...\n")
  
  all_needed_snps <- c()
  
  # 快速遍历一遍 Exp 文件，只为了收集 SNP ID
  for (f in exp_files) {
    env_temp <- new.env()
    load(f, envir = env_temp)
    # 假设变量名是 mvmr_data_combined
    if("mvmr_data_combined" %in% ls(env_temp)){
      tmp_dat <- env_temp[["mvmr_data_combined"]]
      all_needed_snps <- c(all_needed_snps, tmp_dat$SNP)
    }
  }
  
  all_needed_snps <- unique(all_needed_snps)
  cat(paste0("    Total unique SNPs required across all covariates: ", length(all_needed_snps), "\n"))
  
  
  # --- 3. 外层循环: 遍历 Outcome (最耗时的IO操作只做一次) ---
  for (j in 1:length(out_files)) {
    
    current_zip <- out_files[j]
    full_zip_path <- file.path(out_dir, current_zip)
    outcome_name <- gsub("\\.zip$", "", current_zip)
    
    cat(paste0("\n===================================================================\n"))
    cat(paste0(">>> Processing Outcome ", j, "/", length(out_files), ": ", outcome_name, "\n"))
    cat(paste0("===================================================================\n"))
    
    if (!file.exists(full_zip_path)) { warning("File missing:", full_zip_path); next }
    
    # --- 3.1 读取并立刻筛选 Outcome ---
    # 这里的策略是：读取染色体 -> 用全局SNP表筛选 -> 存入内存
    # 这样内存里只有几千行数据，而不是几千万行
    
    zip_contents <- unzip(full_zip_path, list = TRUE)
    csv_files <- grep("\\.csv$", zip_contents$Name, value = TRUE, ignore.case = TRUE)
    
    outcome_list <- list()
    cat("    Reading & Filtering Chromosomes: ")
    
    for (csv_f in csv_files) {
      unzip(full_zip_path, files = csv_f, exdir = ".")
      
      # fread 读取
      temp_dt <- fread(csv_f)
      
      if ("ID" %in% names(temp_dt)) {
        # 【核心优化】直接用全局 SNP 列表筛选
        # 只保留那些在任意一个协变量分析中会用到的 SNP
        temp_subset <- temp_dt[ID %in% all_needed_snps, ]
        if(nrow(temp_subset) > 0){
          outcome_list[[csv_f]] <- temp_subset
        }
      }
      file.remove(csv_f) # 删掉解压的临时文件
      cat(".") 
    }
    cat("\n")
    
    outcome_raw_global <- rbindlist(outcome_list, fill = TRUE)
    
    if (nrow(outcome_raw_global) == 0) {
      cat("    No overlapping SNPs found for ANY exposure. Skipping this Outcome.\n")
      next
    }
    
    # --- 3.2 统一清洗 Outcome 数据 ---
    # 数据量现在很小，处理起来极快
    cat("    Cleaning Outcome Data (Global)...\n")
    
    outcome_raw_global <- as.data.frame(outcome_raw_global)
    outcome_raw_global$OR <- as.numeric(as.character(outcome_raw_global$OR))
    outcome_raw_global$`LOG(OR)_SE` <- as.numeric(as.character(outcome_raw_global$`LOG(OR)_SE`))
    outcome_raw_global$P <- as.numeric(as.character(outcome_raw_global$P))
    
    # 移除无效数据
    outcome_raw_global <- outcome_raw_global[!is.na(outcome_raw_global$OR) & outcome_raw_global$OR > 0, ]
    
    if (nrow(outcome_raw_global) == 0) next
    
    # 计算 Beta 和 Other Allele
    outcome_raw_global$BETA_CALC <- log(outcome_raw_global$OR)
    outcome_raw_global$OTHER_ALLELE <- ifelse(outcome_raw_global$A1 == outcome_raw_global$ALT, 
                                              outcome_raw_global$REF, outcome_raw_global$ALT)
    
    # 此时 outcome_raw_global 是一个干净的、包含了所有可能用到 SNP 的小表
    
    
    # --- 4. 内层循环: 遍历 Exposure (内存操作，极快) ---
    cat("    Starting Exposure Loop (Internal Memory Operations)...\n")
    
    for (exp_f in exp_files) {
      
      # 加载暴露数据
      env_temp <- new.env()
      load(exp_f, envir = env_temp)
      mvmr_long <- env_temp[["mvmr_data_combined"]]
      
      # 提取协变量名
      fname <- basename(exp_f)
      cov_name <- gsub("BMI_1COV_|\\.RData", "", fname)
      
      # 获取当前这个 Exposure 需要的 SNP
      current_target_snps <- unique(mvmr_long$SNP)
      
      # --- 4.1 从全局 Outcome 表中提取当前需要的子集 ---
      # 这一步是在内存中做 subset，毫秒级
      outcome_subset <- outcome_raw_global[outcome_raw_global$ID %in% current_target_snps, ]
      
      if (nrow(outcome_subset) == 0) {
        # cat("      No SNPs for", cov_name, "- skipping.\n") 
        next
      }
      
      # --- 4.2 Format Data ---
      out_dat <- format_data(
        outcome_subset,
        type = "outcome",
        snps = current_target_snps,
        snp_col = "ID",              
        beta_col = "BETA_CALC",      
        se_col = "LOG(OR)_SE",       
        effect_allele_col = "A1",    
        other_allele_col = "OTHER_ALLELE", 
        pval_col = "P",              
        samplesize_col = "OBS_CT",   
        chr_col = "#CHROM",
        pos_col = "POS"
      )
      out_dat$outcome <- outcome_name
      
      
      # --- 4.3 MVMR Analysis (复用之前的逻辑) ---
      # cat(paste0("      Running MVMR: BMI + ", cov_name, "...\n"))
      
      # Harmonise
      mv_dat <- mv_harmonise_data(exposure_dat = mvmr_long, outcome_dat = out_dat)
      
      if (is.null(mv_dat$exposure_beta) || nrow(mv_dat$exposure_beta) < 3) next
      
      # Run MVMR
      res_mvmr <- mv_multiple(mv_dat)
      res_final <- res_mvmr$result
      res_final$outcome <- outcome_name
      res_final$covariate_adjustment <- cov_name 
      
      # F-stat
      tryCatch({
        mvmr_input <- format_mvmr(
          BXGs = mv_dat$exposure_beta, BYG = mv_dat$outcome_beta,
          seBXGs = mv_dat$exposure_se, seBYG = mv_dat$outcome_se,
          RSID = rownames(mv_dat$exposure_beta)
        )
        f_stat_res <- strength_mvmr(mvmr_input, gencov = 0)
        
        res_final$conditional_f_stat <- NA
        if("BMI" %in% rownames(f_stat_res)) res_final$conditional_f_stat[res_final$exposure == "BMI"] <- f_stat_res["BMI", "F_stat"]
        if(cov_name %in% rownames(f_stat_res)) res_final$conditional_f_stat[res_final$exposure == cov_name] <- f_stat_res[cov_name, "F_stat"]
        if(all(is.na(res_final$conditional_f_stat))) {
          res_final$conditional_f_stat[1] <- f_stat_res[1,1]; res_final$conditional_f_stat[2] <- f_stat_res[2,1]
        }
      }, error = function(e) {})
      
      # Save Result
      res_final <- generate_odds_ratios(res_final)
      out_csv_name <- paste0(output_path, "MVMR_BMI_", cov_name, "_vs_", outcome_name, ".csv")
      write.csv(res_final, out_csv_name, row.names = FALSE)
      
    } # End Inner Loop (Exposure)
    
    cat(paste0("    Done with Outcome: ", outcome_name, "\n"))
    
    # 强制 GC 释放内存，为下一个大文件做准备
    rm(outcome_raw_global); gc()
    
  } # End Outer Loop (Outcome)
  
  cat("\nAll optimized MVMR analyses finished!\n")
}


# =======================================================================================
# 3. 开始运行
# =======================================================================================

run_mvmr_optimized(
  exp_dir = exp_data_dir,
  out_dir = outcome_dir_base,
  out_files = outcome_files,
  output_path = save_dir
)
