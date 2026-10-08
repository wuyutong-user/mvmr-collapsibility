#########################################################################################
# 通用 MR 分析代码 (适配 ZIP 多染色体 CSV 格式 + BMI RData 暴露)
# 日期: 2026-01-20
#########################################################################################

# 1. 定义核心分析函数
run_mr_analysis <- function(
    exp_rdata_path,      # 暴露数据的完整路径 (RData)
    out_dir,             # Outcome zip 文件所在的文件夹路径
    out_files,           # Outcome zip 文件名列表 (例如 c("AAA outcome1.zip"))
    save_path,           # 结果保存路径
    p_val_thresh = 5e-8, # 暴露筛选 P 值 (如果 RData 已筛选则此参数仅作记录)
    # meths_use,
    meths_use = c("mr_wald_ratio","mr_ivw", "mr_egger_regression"), # 默认方法
    MRPRESSO = FALSE     # 是否进行 MR-PRESSO (耗时，根据需要开启)
) {
  
  # --- 环境检查 ---
  packages <- c("TwoSampleMR", "data.table", "stringr", "MRPRESSO", "ggplot2", "readr")
  for(p in packages){
    if(!require(p, character.only = T)){
      install.packages(p)
      library(p, character.only = T)
    }
  }
  
  # 创建保存目录
  if(!dir.exists(save_path)) dir.create(save_path, recursive = T)
  
  # --- 第一步：加载并格式化暴露数据 (只做一次) ---
  cat(">>> Step 1: Loading Exposure Data (BMI)...\n")
  
  # 加载 RData (变量名通常是 BMI_P_clump，这里动态获取变量名以防万一)
  env_temp <- new.env()
  load(exp_rdata_path, envir = env_temp)
  exp_var_name <- ls(env_temp)[1] # 假设 RData 里只有一个主要对象
  exp_raw <- env_temp[[exp_var_name]]
  
  # 确保暴露数据是 data.frame 并且有 SNP 列
  exp_dat <- as.data.frame(exp_raw)
  
  # The input RData contains the clumped exposure data in TwoSampleMR format.
  # Reuse the existing formatted exposure columns.
  # 如果是原始 GWAS 格式，需要 format_data。这里假设是处理好的，只需要保证有 SNP 列。
  
  # *重要*：为了后续筛选 Outcome，我们需要 SNP 列表
  target_snps <- exp_dat$SNP
  cat("    Exposure loaded. Total SNPs: ", length(target_snps), "\n")
  
  
  # --- 第二步：循环处理每一个 Outcome (ZIP 包) ---
  for (j in 1:length(out_files)) {
    
    current_zip <- out_files[j]
    full_zip_path <- file.path(out_dir, current_zip)
    outcome_name <- gsub("\\.zip$", "", current_zip) # 去掉 .zip 后缀作为名字
    
    cat(paste0("\n>>> Step 2: Processing Outcome ", j, "/", length(out_files), ": ", outcome_name, "\n"))
    
    # --- 2.1 智能读取 ZIP 中的 CSV ---
    if (!file.exists(full_zip_path)) {
      warning(paste("File not found:", full_zip_path))
      next
    }
    
    # 列出 zip 内容
    zip_contents <- unzip(full_zip_path, list = TRUE)
    # 找到所有的 csv 文件 (排除文件夹等)
    csv_files <- grep("\\.csv$", zip_contents$Name, value = TRUE, ignore.case = TRUE)
    
    cat("    Found", length(csv_files), "chromosome files in zip. Reading and filtering...\n")
    
    outcome_list <- list()
    
    # 循环读取每条染色体 (为了省内存，读一条，筛一条，存一条)
    for (csv_f in csv_files) {
      # 解压单个文件到临时目录
      unzip(full_zip_path, files = csv_f, exdir = ".")
      
      # 读取 (使用 fread 加速)
      temp_dt <- fread(csv_f)
      
      # ***关键步骤***：立刻筛选！只保留暴露中存在的 SNP
      # The source ID column supplies SNP identifiers.
      if ("ID" %in% names(temp_dt)) {
        temp_subset <- temp_dt[ID %in% target_snps, ]
        outcome_list[[csv_f]] <- temp_subset
      } else {
        warning(paste("Column 'ID' not found in", csv_f))
      }
      
      # 删除解压出来的临时文件，保持目录干净
      file.remove(csv_f)
      
      # 进度条效果
      cat(".") 
    }
    cat("\n")
    
    # 合并所有染色体的筛选结果
    outcome_raw <- rbindlist(outcome_list, fill = TRUE)
    
    if (nrow(outcome_raw) == 0) {
      cat("    No overlapping SNPs found for this outcome. Skipping.\n")
      next
    }
    
    # --- 2.2 格式化 Outcome Data ---
    cat("    Formatting Outcome Data...\n")
    
    # 数据预处理：强制转换为 numeric，防止读入为字符型导致报错
    outcome_raw <- as.data.frame(outcome_raw)
    
    # 【新增】强制转换关键列为数字，非数字会被转为 NA
    outcome_raw$OR <- as.numeric(as.character(outcome_raw$OR))
    outcome_raw$`LOG(OR)_SE` <- as.numeric(as.character(outcome_raw$`LOG(OR)_SE`)) # 注意列名里的特殊字符
    outcome_raw$P <- as.numeric(as.character(outcome_raw$P))
    
    # 【新增】移除转换后 OR 为 NA 或 OR <= 0 的行（log无法处理负数或0）
    outcome_raw <- outcome_raw[!is.na(outcome_raw$OR) & outcome_raw$OR > 0, ]
    
    if (nrow(outcome_raw) == 0) {
      cat("    No valid numeric OR data left after cleaning. Skipping.\n")
      next
    }
    
    # 计算 Beta = log(OR)
    outcome_raw$BETA_CALC <- log(outcome_raw$OR)
    
    # 推导 Other Allele
    # 假设是双等位基因。如果 A1 == ALT, 则 Other = REF; 否则 Other = ALT
    outcome_raw$OTHER_ALLELE <- ifelse(outcome_raw$A1 == outcome_raw$ALT, 
                                       outcome_raw$REF, 
                                       outcome_raw$ALT)
    
    out_dat <- format_data(
      outcome_raw,
      type = "outcome",
      snps = target_snps,
      snp_col = "ID",              
      beta_col = "BETA_CALC",      
      se_col = "LOG(OR)_SE",       # 注意这里必须和原数据列名完全一致
      effect_allele_col = "A1",    
      other_allele_col = "OTHER_ALLELE", 
      pval_col = "P",              
      # eaf_col = "A1_FREQ",       # 如果没有频率列，保持注释
      samplesize_col = "OBS_CT",   
      chr_col = "#CHROM",
      pos_col = "POS"
    )
    
    out_dat$outcome <- outcome_name
    
    # --- 2.3 Harmonise (协同) ---
    cat("    Harmonising data...\n")
    dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2)
    
    # 剔除无法协同的 SNP
    dat <- dat[dat$mr_keep == TRUE, ]
    
    if (nrow(dat) == 0) {
      cat("    No SNPs left after harmonisation. Skipping.\n")
      next
    }
    
    # --- 2.4 Perform MR ---
    cat("    Performing MR Analysis...\n")
    
    # 计算 F 统计量和 R2 (如果你需要手动算)
    dat$r2 <- (get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure))^2
    dat$F_stat <- dat$r2 * (dat$samplesize.exposure - 2) / (1 - dat$r2)
    
    # 运行 MR
    res <- mr(dat, method_list = meths_use)
    
    # 生成 OR 值 (Exponetiate results)
    res_or <- generate_odds_ratios(res)
    
    # --- 2.5 敏感性分析 (Heterogeneity & Pleiotropy) ---
    res_heter <- mr_heterogeneity(dat)
    res_pleio <- mr_pleiotropy_test(dat)
    
    # --- 2.6 MR-PRESSO (可选) ---
    if (MRPRESSO & nrow(dat) > 3) {
      cat("    Running MR-PRESSO (this may take time)...\n")
      tryCatch({
        presso_out <- run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
        # 将 PRESSO 结果简单合并或者保存 (此处简化处理，只保存主要结果)
        # Save MR-PRESSO results separately from the main MR summary.
        write.csv(presso_out[[1]]$`Main MR results`, 
                  paste0(save_path, outcome_name, "_MR_PRESSO.csv"), row.names = F)
      }, error = function(e) { cat("    MR-PRESSO failed: ", e$message, "\n") })
    }
    
    # --- 2.7 保存结果 ---
    cat("    Saving results...\n")
    
    # 构造文件前缀
    prefix <- paste0(save_path, "BMI_vs_", outcome_name)
    
    # 1. MR 主结果
    write.csv(res_or, paste0(prefix, "_res.csv"), row.names = F)
    # 2. Harmonised Data (原始数据)
    write.csv(dat, paste0(prefix, "_dat.csv"), row.names = F)
    # 3. 异质性
    if(nrow(res_heter) > 0) write.csv(res_heter, paste0(prefix, "_heter.csv"), row.names = F)
    # 4. 多效性
    if(nrow(res_pleio) > 0) write.csv(res_pleio, paste0(prefix, "_pleio.csv"), row.names = F)
    # 5. Leave-one-out (如果 SNP 够多)
    if(nrow(dat) > 2) {
      res_loo <- mr_leaveoneout(dat)
      write.csv(res_loo, paste0(prefix, "_loo.csv"), row.names = F)
    }
    
    cat("    Done with ", outcome_name, "\n")
  }
  
  cat("\nAll outcomes processed!\n")
}


# =======================================================================================
# 运行配置区 (请修改这里)
# =======================================================================================

# 1. 暴露数据路径 (RData)
my_exp_path <- file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "BMI_10000_0.0001.RData")

# 2. Outcome 文件夹路径 (存放那些 AAA outcome1.zip 的地方)
my_out_dir <- Sys.getenv("ERACA_OUTCOME_DIR", unset = file.path("controlled_data", "outcome_gwas"))

# 3. 获取该文件夹下所有的 ZIP 文件列表
# 这样你就不用手写 c("AAA outcome1.zip", ...) 了，自动读取所有 zip
my_out_files <- list.files(my_out_dir, pattern = "\\.zip$", full.names = FALSE)

# 如果你只想测其中几个，可以手动指定，例如：
# my_out_files <- c("AAA outcome6.zip", "AAA outcome7.zip", "AAA outcome8.zip", "AAA outcome9.zip")

# 4. 结果保存路径
my_save_path <- file.path(Sys.getenv("ERACA_CONTROLLED_ROOT", unset = "controlled_data"), "UVMR_RESULT")

# =======================================================================================
# 开始运行
# =======================================================================================

run_mr_analysis(
  exp_rdata_path = my_exp_path,
  out_dir = my_out_dir,
  out_files = my_out_files,
  save_path = my_save_path,
  MRPRESSO = FALSE  # 如果跑得太慢，保持 FALSE
)
