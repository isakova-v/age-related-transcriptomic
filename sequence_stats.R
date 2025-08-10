library(readr)
library(tidyverse)
library(dplyr)

setwd('~/Studies/Schools/Public.Health.Hackathon/')

# ----- DATA FROM THE ARTICLE -----

library(readxl)
invar_article <- read_excel("aging-v17i1-206192-supplementary-material-SD2.xls", 
                                                            sheet = "9. Final Genelits", skip = 1) %>% as.data.frame()

rownames(invar_article) <- invar_article[[1]]
invar_article <- invar_article[ , -1]

only_invar_article <- lapply(invar_article, function(col) {
  rownames(invar_article)[which(col)]
})

only_var_article <- lapply(invar_article, function(col) {
  rownames(invar_article)[which(!col)]
})

all_tissues_article <- rownames(invar_article)[rowSums(invar_article) == ncol(invar_article)]



# ----- OUR DATA -----

# dataframe of genes and their most representative transcripts
repr_tr <- read_csv('representative_transcript.csv')

# dataframe of all genes (in ENSEMBL and SYMBOL) and corresponding transcripts
tr_gene_info <- read_tsv('gencode_vM22.transcript2gene.tsv')
tr_gene_info$gene_ens_id <- sub("\\..*$", "", tr_gene_info$gene_id)

# number of transcripts for all genes
tr_counts <- tr_gene_info %>%
  group_by(gene_ens_id, gene_name) %>%
  summarise(amount_of_transcripts = n_distinct(transcript_id), .groups = "drop")

# tmm <- tmm %>%
#   dplyr::left_join(
#     tr_gene_info %>%
#       dplyr::select(gene_ens_id, gene_name = gene_name) %>%
#       dplyr::distinct(),
#     by = c("feature" = "gene_ens_id")
#   )
# 
# all_genes_our <- unique(c(tmm$gene_name, tpm$feature))


# # dataframe: gene_id, gene_name, longest transcript exon length
# gene_data <- read_tsv('gencode_vM22.gene_data.tsv')


# obtain GC content for every transcript
library(Biostrings)
tr_fasta <- readDNAStringSet("gencode.vM22.transcripts.fa")
names(tr_fasta) <- sub("\\|.*$", "", names(tr_fasta))

gc_counts <- letterFrequency(tr_fasta, letters = c("G", "C"), as.prob = FALSE)
tr_len <- width(tr_fasta)

gc_table <- data.frame(
  transcript_id = names(tr_fasta),
  length_nt = tr_len,
  GC_percent = rowSums(gc_counts) / tr_len * 100
)

# collect all info:
# gene_ens_id, repr_transcript, gene name, amount_of_transcripts, tr_len, GC content
repr_tr_full <- repr_tr %>%
  left_join(tr_counts, by = c("gene_id" = "gene_ens_id")) %>%
  left_join(gc_table, by = c("representative_transcript" = "transcript_id"))


# # all filtered genes in our data
# all_genes_our_with_id <- data.frame(gene_name = all_genes_our) %>%
#   dplyr::left_join(
#     tr_gene_info %>%
#       dplyr::select(gene_id, gene_name) %>%
#       dplyr::distinct(),
#     by = "gene_name"
#   )

# Our invariative genes
invar_our <- read_csv(
  "final_genes_all_tissues.csv",
  quote = '"'
)

only_invar_our <- setNames(
  lapply(seq_len(nrow(invar_our)), function(i) {
    tibble(
      gene_id = strsplit(invar_our$Ensembl_Genes[i], ",")[[1]],
      gene_name = strsplit(invar_our$Gene_Symbols[i], ",")[[1]]
    ) %>%
      left_join(repr_tr_full, by = "gene_id")
  }),
  invar_our$Tissue
)

# FIX SPLEEN TABLE
spleen <- read_csv("tpm_invariant_genes_Spleen.csv")
colnames(spleen)[1] <- "gene_id"

new_spleen <- spleen['gene_id'] %>%
  # rename(gene_id = feature) %>%  # переименовываем колонку для join
  left_join(repr_tr_full, by = "gene_id") 

new_spleen$gene_name.x = new_spleen$gene_name
new_spleen$gene_name.y = new_spleen$gene_name
new_spleen$gene_name <- NULL

only_invar_our$Spleen <- new_spleen


non_invar_our <- lapply(only_invar_our, function(tbl) {
  repr_tr_full %>%
    anti_join(tbl %>% select(gene_id), by = "gene_id") %>%
    select(gene_id, gene_name) %>%
    left_join(repr_tr_full, by = "gene_id")
})


# Compare genes

# Summary
get_summary <- function(df, type_str) {
  df %>%
    summarise(
      mean_GC = mean(GC_percent, na.rm = TRUE),
      median_GC = median(GC_percent, na.rm = TRUE),
      mean_length = mean(length_nt, na.rm = TRUE),
      median_length = median(length_nt, na.rm = TRUE),
      mean_transcripts = mean(amount_of_transcripts, na.rm = TRUE),
      median_transcripts = median(amount_of_transcripts, na.rm = TRUE),
      count = n()
    ) %>%
    mutate(type = type_str)
}

summary_invar_noninvar <- function(tissue_name) {
  summary_invar <- get_summary(only_invar_our[[tissue_name]], "invariant")
  summary_noninvar <- get_summary(non_invar_our[[tissue_name]], "noninvariant")
  
  bind_rows(summary_invar, summary_noninvar) %>%
    mutate(tissue = tissue_name) %>%
    select(tissue, everything())
}

all_summaries <- lapply(names(only_invar_our), summary_invar_noninvar) %>%
  bind_rows()


# Check normal distribution
library(tidyr)

ks_norm_test <- function(x) {
  if(length(x) < 3) return(NA)
  mu <- mean(x)
  sigma <- sd(x)
  if(sigma == 0) return(NA) # если нет вариации
  
  # Удаляем повторы, чтобы избежать warning
  x_unique <- unique(x)
  if(length(x_unique) < 3) return(NA)  # если после удаления дублей мало данных
  
  ks.test(x_unique, "pnorm", mean = mu, sd = sigma)$p.value
}

distr_results <- list()
vars <- c("GC_percent", "length_nt", "amount_of_transcripts")
tissues <- intersect(names(only_invar_our), names(non_invar_our))

for(tissue in tissues) {
  invar_df <- only_invar_our[[tissue]]
  noninvar_df <- non_invar_our[[tissue]]
  
  for(var in vars) {
    invar_vec <- na.omit(invar_df[[var]])
    noninvar_vec <- na.omit(noninvar_df[[var]])
    
    p_invar <- ks_norm_test(invar_vec)
    p_noninvar <- ks_norm_test(noninvar_vec)
    
    distr_results[[paste(tissue, var, sep = "_")]] <- tibble(
      tissue = tissue,
      variable = var,
      group = c("invariant", "noninvariant"),
      p_value = c(p_invar, p_noninvar)
    )
  }
}

distr_results_df <- bind_rows(distr_results)


# Compare
library(broom)

compare_groups <- function(tissue_name, var_name) {
  invar_df <- only_invar_our[[tissue_name]]
  noninvar_df <- non_invar_our[[tissue_name]]
  
  invar_vec <- na.omit(invar_df[[var_name]])
  noninvar_vec <- na.omit(noninvar_df[[var_name]])
  
  p_norm_invar <- ks_norm_test(invar_vec)
  p_norm_noninvar <- ks_norm_test(noninvar_vec)
  
  if (!is.na(p_norm_invar) && !is.na(p_norm_noninvar) && p_norm_invar > 0.05 && p_norm_noninvar > 0.05) {
    test_res <- t.test(invar_vec, noninvar_vec)
    test_used <- "t-test"
  } else {
    test_res <- wilcox.test(invar_vec, noninvar_vec)
    test_used <- "wilcox-test"
  }
  
  tidy_res <- broom::tidy(test_res) %>%
    mutate(
      tissue = tissue_name,
      variable = var_name,
      test = test_used,
      p_value = p.value
    )
  
  if("estimate" %in% colnames(tidy_res)) {
    tidy_res <- tidy_res %>% select(tissue, variable, test, p_value, estimate, statistic, method)
  } else {
    tidy_res <- tidy_res %>% select(tissue, variable, test, p_value, statistic, method)
  }
  
  return(tidy_res)
}


vars <- c("GC_percent", "length_nt", "amount_of_transcripts")
tissues <- intersect(names(only_invar_our), names(non_invar_our))

all_tests <- purrr::map_dfr(tissues, function(tiss) {
  purrr::map_dfr(vars, function(v) {
    compare_groups(tiss, v)
  })
})


# ----- PLOTS -----

library(ggplot2)
library(dplyr)
library(tidyr)

# Функция для вычисления p-value между группами (используем t-test или wilcox.test)
get_p_value <- function(df_long) {
  invar_vec <- na.omit(df_long %>% filter(group == "invariant") %>% pull(value))
  noninvar_vec <- na.omit(df_long %>% filter(group == "noninvariant") %>% pull(value))
  
  p_norm_invar <- ks_norm_test(invar_vec)
  p_norm_noninvar <- ks_norm_test(noninvar_vec)
  
  if (!is.na(p_norm_invar) && !is.na(p_norm_noninvar) && p_norm_invar > 0.05 && p_norm_noninvar > 0.05) {
    test_res <- t.test(invar_vec, noninvar_vec)
  } else {
    test_res <- wilcox.test(invar_vec, noninvar_vec)
  }
  test_res$p.value
}

plot_list <- list()
vars <- c("GC_percent", "length_nt", "amount_of_transcripts")

for(tiss in intersect(names(only_invar_our), names(non_invar_our))) {
  plot_data <- bind_rows(
    only_invar_our[[tiss]] %>% mutate(group = "invariant"),
    non_invar_our[[tiss]] %>% mutate(group = "noninvariant")
  )
  
  plot_data_long <- plot_data %>%
    select(group, all_of(vars)) %>%
    pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value")
  
  # Для каждой переменной берём только её подтаблицу и считаем p-value
  pval_df <- lapply(vars, function(v) {
    df_sub <- plot_data_long %>% filter(variable == v)
    data.frame(variable = v, p_value = get_p_value(df_sub))
  }) %>% bind_rows()
  
  pval_df <- pval_df %>%
    mutate(label = paste0("p = ", signif(p_value, 3)))
  
  p <- ggplot(plot_data_long, aes(x = group, y = value, fill = group)) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.5, aes(color = group)) +  # Добавлен jitter с цветом по группе
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0.8) +
    facet_wrap(~ variable, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(title = paste("Tissue:", tiss), x = "Group", y = "Value") +
    geom_text(
      data = pval_df,
      aes(x = 1.5, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.5,
      size = 4
    )
  
  plot_list[[tiss]] <- p
}

# Пример вывода графика по одной ткани:
print(plot_list[["Spleen"]])
