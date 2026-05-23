# ==============================================================================
# miRNA/sRNA Target Prediction & Statistical Evaluation Pipeline
# ==============================================================================
# ==============================================================================
# Script Name: miRNA_Target_Poisson_Predictor.R
# Author: 旭旭旭 (GitHub: pengxu0907-cpu)
# Year: 2026
# 
# Description: 
# A custom pipeline for sRNA/miRNA target prediction utilizing a 1st-order 
# Markov background model with Laplace smoothing and Poisson distribution testing.
#
# Acknowledgments:
# The foundational structural logic was inspired by community bioinformatics 
# discussions, subsequently heavily optimized and refactored for strict 
# statistical symmetry (Poisson independence), floating-point underflow 
# prevention, and full-universe FDR correction.
#
# License: MIT License (Free to use, modify, and distribute with attribution)
# ==============================================================================

library(Biostrings)

# ==========================================
# 1. Background Modeling
# ==========================================
build_markov_model <- function(fasta_file) {
  seqs <- readDNAStringSet(fasta_file)
  alpha <- 1 
  mono_counts <- oligonucleotideFrequency(seqs, width = 1, as.prob = FALSE)
  mono_sums <- colSums(mono_counts)
  p_mono <- (mono_sums + alpha) / (sum(mono_sums) + 4 * alpha)
  
  di_counts <- oligonucleotideFrequency(seqs, width = 2, as.prob = FALSE)
  di_sums <- colSums(di_counts)
  
  p_transition <- list()
  bases <- c("A", "C", "G", "T")
  for (n1 in bases) {
    n1_starts <- paste0(n1, bases)
    sum_n1_starts <- sum(di_sums[n1_starts])
    for (n2 in bases) {
      dinuc <- paste0(n1, n2)
      p_transition[[dinuc]] <- (di_sums[dinuc] + alpha) / (sum_n1_starts + 4 * alpha)
    }
  }
  return(list(p_mono = p_mono, p_transition = p_transition))
}

# ==========================================
# 2. Expected Probability Estimation
# ==========================================
calculate_expected_prob <- function(seed_seq, model) {
  seed_seq <- toupper(trimws(seed_seq))
  seed_chars <- strsplit(seed_seq, "")[[1]]
  if (!all(seed_chars %in% c("A", "C", "G", "T"))) {
    stop("CRITICAL ERROR: Seed sequence contains invalid characters. Only A, C, G, T are permitted.")
  }
  prob <- model$p_mono[seed_chars[1]]
  for (i in 1:(length(seed_chars) - 1)) {
    dinuc <- paste0(seed_chars[i], seed_chars[i+1])
    prob <- prob * model$p_transition[[dinuc]]
  }
  return(prob)
}

# ==========================================
# 3. Motif Scanning & Poisson Framework
# ==========================================
predict_targets <- function(fasta_file, seed_seq, expected_prob) {
  seed_seq <- toupper(trimws(seed_seq))
  seqs <- readDNAStringSet(fasta_file)
  seed_len <- nchar(seed_seq) 
  
  results_list <- lapply(seq_along(seqs), function(i) {
    gene_id <- names(seqs)[i]
    seq_obj <- seqs[[i]] 
    seq_len <- length(seq_obj) 
    windows_n <- seq_len - seed_len + 1
    if (windows_n <= 0) return(NULL) 
    
    raw_lambda <- windows_n * expected_prob
    expected_lambda <- ifelse(raw_lambda == 0, .Machine$double.xmin, raw_lambda)
    observed_count <- countPattern(seed_seq, seq_obj)
    
    if (observed_count == 0) {
      p_val <- 1.0
    } else {
      p_val <- ppois(q = observed_count - 1, lambda = expected_lambda, lower.tail = FALSE)
    }
    
    return(data.frame(Gene_ID = gene_id, Sequence_Length = seq_len,
                      Expected_Lambda = expected_lambda, Observed_Count = observed_count,
                      P_value = p_val, stringsAsFactors = FALSE))
  })
  
  results_df <- do.call(rbind, results_list)
  if (is.null(results_df)) {
    results_df <- data.frame(Gene_ID=character(), Sequence_Length=integer(), 
                             Expected_Lambda=numeric(), Observed_Count=integer(), P_value=numeric())
  }
  return(results_df)
}

# ==========================================
# Main Execution Entry (主运行区)
# ==========================================
fasta_path <- "test_target_regions.fasta"  

target_seed <- "GGTCCG"            

cat("1. 正在评估细菌基因组马尔可夫背景...\n")
model <- build_markov_model(fasta_path)

cat(sprintf("2. 正在计算序列的背景概率: %s...\n", target_seed))
exp_prob <- calculate_expected_prob(target_seed, model)

cat("3. 正在全基因组扫描并计算泊松分布P值...\n")
results_df <- predict_targets(fasta_path, target_seed, exp_prob)

if (nrow(results_df) == 0) {
  cat("\n运行结束，未找到有效序列。\n")
} else {
  cat("4. 正在进行FDR多重检验校正...\n")
  results_df$FDR_Q_value <- p.adjust(results_df$P_value, method = "BH")
  results_df$Significant <- results_df$FDR_Q_value < 0.05
  results_df <- results_df[order(results_df$FDR_Q_value), ]
  
  output_file <- sprintf("6mer_Target_Prediction_%s.csv", target_seed)
  write.csv(results_df, output_file, row.names = FALSE)
  
  sig_num <- sum(results_df$Significant, na.rm = TRUE)
  cat("\n==============================================================\n")
  cat(sprintf("【分析报告】:\n"))
  cat(sprintf("  - 具有统计学显著性的潜在靶标基因数 (FDR < 0.05): %d\n", sig_num))
  cat("==============================================================\n")
  cat(sprintf("完整结果表格已保存到你的文件夹，文件名为: %s\n\n", output_file))
}