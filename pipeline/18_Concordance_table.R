library(tidyverse)
library(ggplot2)
library(boot)
library(ggtext)

# STEP 1:
setwd("[Local_computer_pathway]/data/Concordance")

# STEP 2:
bootstrap_ci_counts <- function(counts) {
  data <- rep(1:length(counts), round(counts))
  count_func <- function(data, indices) {
    sample_data <- data[indices]
    sapply(1:length(counts), function(x) sum(sample_data == x))
  }
  results <- boot(data, statistic = count_func, R = 1000)
  tibble(
    lower_ci_count = apply(results$t, 2, quantile, probs = 0.025),
    upper_ci_count = apply(results$t, 2, quantile, probs = 0.975)
  )
}

create_heatmap <- function(branch_id, conc_vectors) {
  clade_data <- conc_vectors %>%
    filter(ID == branch_id) %>%
    select(gene_psi1, gene_psi2, gene_psi3, gene_psi4, 
           site_psi1, site_psi2, site_psi3, site_psi4, 
           quartet_psi1, quartet_psi2, quartet_psi3, quartet_psi4,
           gene_psi1_N, gene_psi2_N, gene_psi3_N, gene_psi4_N, 
           site_psi1_N, site_psi2_N, site_psi3_N, site_psi4_N, 
           quartet_psi1_N, quartet_psi2_N, quartet_psi3_N, quartet_psi4_N)
  
  long_data <- clade_data %>%
    select(-ends_with("_N")) %>%
    pivot_longer(cols = everything(), names_to = "type_psi", values_to = "value") %>%
    separate(type_psi, into = c("type", "psi"), sep = "_psi") %>%
    mutate(type = factor(type, levels = c("quartet", "site", "gene")))
  
  long_data_N <- clade_data %>%
    select(ends_with("_N")) %>%
    pivot_longer(cols = everything(), names_to = "type_psi_N", values_to = "count") %>%
    separate(type_psi_N, into = c("type", "psi_N"), sep = "_psi") %>%
    mutate(type = factor(type, levels = c("quartet", "site", "gene")))
  
  long_data_N <- long_data_N %>%
    group_by(type) %>%
    group_modify(~ {
      counts <- .x$count
      total_counts <- sum(counts)
      ci <- bootstrap_ci_counts(counts) %>%
        mutate(lower_CI = (lower_ci_count / total_counts) * 100,
               upper_CI = (upper_ci_count / total_counts) * 100,
               total_counts = total_counts)
      bind_cols(.x, ci)
    }) %>%
    ungroup()
  
  long_data <- long_data %>%
    mutate(psi = paste0(psi, "_N")) %>%
    left_join(long_data_N %>% select(type, psi_N, lower_CI, upper_CI, total_counts), 
              by = c("type", "psi" = "psi_N")) %>%
    mutate(psi = str_replace(psi, "_N", ""),
           main_label = scales::number(value, accuracy = 0.001),
           ci_label = paste0("(", scales::number(lower_CI, accuracy = 0.001), ", ", 
                             scales::number(upper_CI, accuracy = 0.001), ")"))
  
  y_labels <- long_data %>%
    select(type, total_counts) %>%
    distinct() %>%
    mutate(y_label = paste0(type, "<br><span style='font-size:8pt;'>(n=", total_counts, ")</span>")) %>%
    arrange(factor(type, levels = c("quartet", "site", "gene"))) %>%
    pull(y_label)
  
  plot <- ggplot(long_data, aes(x = psi, y = type, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = main_label), color = "black", size = 6, vjust = -0.5) +
    geom_text(aes(label = ci_label), color = "black", size = 3, vjust = 1.5) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 100), name = "proportion") +
    scale_x_discrete(position = 'top', labels = c(bquote(Psi[1]), bquote(Psi[2]), bquote(Psi[3]), bquote(Psi[4]))) +
    scale_y_discrete(labels = y_labels) +
    theme_minimal() +
    ggtitle(paste("Concordance factors for branch ID", branch_id),
            subtitle = "Values are percentages, numbers in brackets are bootstrap 95% CIs") +
    theme(
      plot.title = element_text(size = 24),
      axis.title = element_blank(),
      axis.text.y = element_markdown(size = 16),
      axis.text.x = element_text(size = 16),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "bottom"
    ) +
    guides(fill = guide_colourbar(
      title = "Concordance or Discordance factor (%)",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm"),
      label.position = "bottom",
      label.theme = element_text(size = 14),
      ticks.colour = "black"
    ))
  
  output_data <- long_data %>%
    select(type, psi, value, lower_CI, upper_CI)
  
  return(list(plot = plot, data = output_data))
}

# STEP 3:
# Load data and initialize
concordance_vectors <- read_csv("concordance_vectors.csv")
all_branch_ids <- unique(concordance_vectors$ID)

# Collect all data and plots into one file each
all_output_data <- list()

pdf("all_concordance_plots.pdf", width = 8, height = 5)

for (branch_id in all_branch_ids) {
  cat("Generating plot for branch ID:", branch_id, "\n")
  
  output <- create_heatmap(branch_id, concordance_vectors)
  
  print(output$plot)
  
  # Store data
  output_data <- output$data %>% mutate(ID = branch_id)
  all_output_data[[as.character(branch_id)]] <- output_data
}

dev.off()

# Combine and save all output data tables
combined_data <- bind_rows(all_output_data)
write_csv(combined_data, "all_concordance_data.csv")

cat("Done: all_concordance_plots.pdf and all_concordance_data.csv saved.\n")

