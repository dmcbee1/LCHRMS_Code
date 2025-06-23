# Load necessary library
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify) 
library(readxl)
library(ggrepel)

#Setwd To save files created
setwd('/Users/dmcbee/Desktop/Current MS Data/XCMS_Script/XCMS_Results')

# Import the CSV file
lcms_data <- read.csv('/Users/dmcbee/Desktop/Current MS Data/XCMS_Script/XCMS_Results/metaboanalyst_table.csv')

group1 <- "EB"
group2 <- "WT"
group3 <- "yumC"
group4 <- "ispC"

extract_lcms_features <- function(lcms_data, pellet_mass_file = NULL, internal_standard_mz = NULL, internal_standard_RT = NULL) {
  
  # Extract grouping information (first row, excluding first column)
  grouping <- lcms_data[1, -1] %>% t() %>% as.data.frame()
  colnames(grouping) <- c("Group")
  grouping$Sample <- colnames(lcms_data)[-1]
  
  # Remove the first row from the data to get feature data
  feature_data <- lcms_data[-1, ]
  colnames(feature_data)[1] <- "mz_RT"  # Rename first column
  
  # Separate 'mz_RT' into 'mz' and 'RT'
  feature_data <- feature_data %>%
    separate(mz_RT, into = c("mz", "RT"), sep = "/") %>%
    mutate(mz = as.numeric(mz),
           RT = as.numeric(RT))
  
  # Add a feature ID column
  feature_data <- feature_data %>%
    mutate(ID = row_number())
  
  # Convert all feature intensity values to numeric
  feature_data <- feature_data %>% mutate(across(-c(mz, RT, ID), as.numeric))
  
  # Reshape the data into long format and set the value column name to "Value"
  long_data <- feature_data %>%
    pivot_longer(cols = -c(mz, RT, ID), names_to = "Sample", values_to = "Value")
  
  # Merge with grouping information
  cleaned_lcms_data <- long_data %>%
    left_join(grouping, by = "Sample")
  
  # Replace NA values with 0
  stats_lcms_data <- cleaned_lcms_data %>% mutate(Value = replace_na(Value, 0))
  
  # Reorder columns
  stats_lcms_data <- stats_lcms_data %>% select(ID, mz, RT, Sample, Value, Group)
  
  return(as.data.frame(stats_lcms_data))
}
normalize_data <- function(data, value_col = "Value", norm_method = c("median", "mean", "internal"), internal_peak = NULL) {
  library(dplyr)
  
  norm_method <- match.arg(norm_method)
  
  if (norm_method %in% c("median", "mean")) {
    norm_data <- data %>%
      group_by(Group) %>%
      mutate(
        norm_factor = if (norm_method == "median") {
          median(.data[[value_col]], na.rm = TRUE)
        } else {
          mean(.data[[value_col]], na.rm = TRUE)
        },
        Normalized_Value = .data[[value_col]] / norm_factor
      ) %>%
      ungroup()
    
  } else if (norm_method == "internal") {
    if (is.null(internal_peak)) {
      stop("For internal normalization, please provide an internal_peak identifier (a value in the 'ID' column).")
    }
    # Use tidy evaluation to select the proper column and rename it to internal_value
    internal_factors <- data %>%
      filter(ID == internal_peak) %>%
      dplyr::select(Sample, internal_value = !!rlang::sym(value_col))
    
    norm_data <- data %>%
      left_join(internal_factors, by = "Sample") %>%
      mutate(
        Normalized_Value = if_else(internal_value == 0, NA_real_, .data[[value_col]] / internal_value)
      )
  }
  return(norm_data)
}
scale_data <- function(data,
                       value_col   = "Normalized_Value",
                       scale_method = c("none", "mean_center", "auto", "pareto", "range", "log")) {
  # Validate method
  scale_method <- match.arg(scale_method)
  
  # Internal scaler
  apply_scaling <- function(x) {
    x <- as.numeric(x)
    
    # 1) No scaling
    if (scale_method == "none") {
      return(x)
    }
    # 2) Log scaling
    if (scale_method == "log") {
      return(log(x + 1))
    }
    
    # 3) All other methods need summary stats
    x_mean <- mean(x, na.rm = TRUE)
    x_sd   <- sd(x,   na.rm = TRUE)
    x_min  <- min(x,  na.rm = TRUE)
    x_max  <- max(x,  na.rm = TRUE)
    
    # Guard against NA or zero
    if (is.na(x_sd) || x_sd == 0)       x_sd  <- 1
    if (is.na(x_min) || x_max == x_min) x_max <- x_min + 1e-16
    
    # Switch on method
    switch(scale_method,
           
           mean_center = x - x_mean,
           
           auto        = (x - x_mean) / x_sd,        # unit variance
           
           pareto      = (x - x_mean) / sqrt(x_sd),  # Pareto scaling
           
           range       = (x - x_min)  / (x_max - x_min)
    )
  }
  
  # Apply per feature (ID)
  data %>%
    group_by(ID) %>%
    mutate(Scaled_Value = apply_scaling(.data[[value_col]])) %>%
    ungroup()
}

remove_extraction_blank <- function(df, blank_threshold = 1, value_col = "Normalized_Value") {
  library(dplyr)
  
  # Calculate the average blank intensity per feature using the specified column
  blank_summary <- df %>%
    filter(Group == "EB") %>%
    group_by(ID) %>%
    summarise(ave_blank = mean(.data[[value_col]], na.rm = TRUE)) %>%
    ungroup()
  
  # Calculate the average intensity for each non-blank group per feature
  group_summary <- df %>%
    filter(Group != "EB") %>%
    group_by(ID, Group) %>%
    summarise(ave_group = mean(.data[[value_col]], na.rm = TRUE)) %>%
    ungroup()
  
  # Left join so that features missing a blank (NA) are kept.
  # Calculate ratio (if no blank or blank is zero, assign Inf)
  group_summary <- group_summary %>%
    left_join(blank_summary, by = "ID") %>%
    mutate(ratio = if_else(is.na(ave_blank) | ave_blank == 0,
                           Inf, 
                           ave_group / ave_blank))
  
  # For each feature, check if any non-blank group meets the threshold
  features_to_keep <- group_summary %>%
    group_by(ID) %>%
    summarise(keep_feature = any(ratio >= blank_threshold)) %>%
    filter(keep_feature) %>%
    pull(ID)
  
  # Filter the original dataframe to retain only these features
  df_filtered <- df %>%
    filter(ID %in% features_to_keep)
  
  return(df_filtered)
}
merge_lcms_data <- function(stat_input) {
  library(dplyr)
  
  # 1 & 2. Summarize all mz and RT Normalizeds in list columns
  summarized_lcms <- stats_lcms_data %>%
    group_by(ID) %>%
    summarise(
      mz = median(mz),
      rt = median(RT),
      intensity = mean(Normalized_Value),
      .groups = "drop"
    )
  
  # 3. Join to fold_change so that each ID is one row,
  #    but includes all mz and RT values in list columns
  result_df <- stat_input %>%
    left_join(summarized_lcms, by = "ID")
  
  return(result_df)
}
pca <- function(data, Group = "Group", Value = "Value") {
  # Prepare the data for PCA by pivoting it to a wide format
  pca_data <- data %>%
    select(Sample, Group, ID, Value) %>%
    pivot_wider(names_from = ID, values_from = Value, values_fn = mean, values_fill = 0)
  
  # Set the Sample column as row names and remove it from the data frame
  rownames(pca_data) <- pca_data$Sample
  pca_data <- pca_data[, -1]  # Remove the Sample column
  
  # Extract the Group information
  groups <- pca_data$Group
  pca_data <- pca_data[, -1]  # Remove Group column before PCA
  
  # Perform PCA
  pca_result <- prcomp(pca_data, scale. = TRUE)
  
  # Create a data frame for plotting
  pca_scores <- as.data.frame(pca_result$x)
  pca_scores$Group <- groups
  pca_scores$Sample <- rownames(pca_scores)
  
  # Interactive plot with plotly
  p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +
    labs(title = "PCA Plot",
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal()
  
  return(p)
}
t_test <- function(data, sig = 0.05, group1, group2, Value = "Value") {
  # Filter data for the two groups of interest
  filtered_data <- data %>%
    filter(Group %in% c(group1, group2))
  
  # Perform T-tests for each ID and store results in a dataframe
  t_test_results <- filtered_data %>%
    group_by(ID) %>%
    summarize(
      p_value = t.test(Value ~ Group)$p.value,
      mz = mean(mz),   # Take the first mz value (should be consistent within ID)
      RT = mean(RT),   # Take the first RT value (should be consistent within ID)
      !!group1 := round(mean(Value[Group == group1], na.rm = TRUE)),  # Mean intensity for group 1, rounded
      !!group2 := round(mean(Value[Group == group2], na.rm = TRUE))   # Mean intensity for group 2, rounded
    )
  
  # Filter for significant p-values
  t_test_results <- t_test_results %>%
    mutate(
      log_p_value = -log10(p_value),
      significance = if_else(p_value < sig, "Significant", "Not Significant")
    )
  
  # Visualization of significant T-test results with log-transformed p-values
  tplot <- ggplot(t_test_results, aes(x = RT, y = log_p_value, color = significance)) +
    geom_point() +
    labs(
      title = paste("Significantly Different Intensities Between", group1, "and", group2),
      x = "RT",
      y = "Negative Log10 P-value",
      color = "Significance"
    ) +
    theme_classic()+
    geom_text_repel(
      data = subset(t_test_results, p_value < sig),
      aes(label = ID),
      size = 3
    )
  
  
  # Display the plot
  print(tplot)
  
  #Save files dynamically
  # 1) CSV of significant results
  csv_name <- paste0('TTestResults_',group1, "_vs_", group2, ".csv")
  write.csv(t_test_results, file = csv_name, row.names = FALSE)
  
  # 2) Static PNG of the plot
  #    Note: ggsave() works with the ggplot object (tplot), not the plotly object.
  plot_name <- paste0("TTestPlot_",group1, "_vs_", group2, ".pdf")
  ggsave(filename = plot_name, plot = tplot, width = 6, height = 4, dpi = 300)
  
  # Return the significant results for further use
  return(t_test_results)
}
volcano_plot <- function(data, group1, group2, n = 2, value_col = "Value") {
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  
  # Filter data for the two groups
  filtered_data <- data %>%
    filter(Group %in% c(group1, group2))
  
  volcano <- filtered_data %>%
    group_by(ID) %>%
    summarize(
      mz = first(mz),
      RT = first(RT),
      count_group1 = sum(.data[[value_col]][Group == group1] > 0, na.rm = TRUE),
      count_group2 = sum(.data[[value_col]][Group == group2] > 0, na.rm = TRUE),
      mean_group1  = mean(.data[[value_col]][Group == group1], na.rm = TRUE),
      mean_group2  = mean(.data[[value_col]][Group == group2], na.rm = TRUE),
      p_value = {
        if (count_group1 >= n | count_group2 >= n) {
          tryCatch(
            t.test(.data[[value_col]][Group == group1],
                   .data[[value_col]][Group == group2])$p.value,
            error = function(e) NA_real_
          )
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      fold_change = if_else(
        count_group1 >= n & count_group2 >= n,
        mean_group2 / mean_group1,
        NA_real_
      ),
      log2FC = if_else(
        count_group1 >= n & count_group2 >= n,
        log2(fold_change),
        NA_real_
      ),
      log_p_value = -log10(p_value)
    )
  
  # Rename columns to incorporate the provided group names
  volcano <- volcano %>%
    rename(
      !!paste0("count_", group1) := count_group1,
      !!paste0("count_", group2) := count_group2,
      !!paste0("mean_", group1)  := mean_group1,
      !!paste0("mean_", group2)  := mean_group2
    )
  
  # Create volcano plot, automatically dropping NA values in aesthetics
  volcano_plot <- ggplot(volcano, aes(x = log2FC, y = log_p_value)) +
    geom_point(aes(color = p_value < 0.05), alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = c("black", "red")) +
    labs(
      title = paste("Volcano Plot:", group1, "vs.", group2),
      x = "Log2 Fold Change",
      y = "-Log10 P-value",
      color = "p < 0.05?"
    ) +
    theme_classic() +
    geom_text_repel(
      data = subset(volcano, p_value < 0.05),
      aes(label = mz),
      size = 3
    )
  
  # Print the plot in the console / RStudio viewer
  print(volcano_plot)
  
  # Save the results data frame as a CSV
  csv_name <- paste0("Volcano_", group1, "_vs_", group2, ".csv")
  write.csv(volcano, file = csv_name, row.names = FALSE)
  
  # Save the plot as a PDF
  plot_name <- paste0("VolcanoPlot_", group1, "_vs_", group2, ".pdf")
  ggsave(filename = plot_name, plot = volcano_plot, width = 6, height = 4, dpi = 300)
  
  return(volcano)
}
preform_anova <- function(data) {
  library(ggplot2)
  library(dplyr)
  library(plotly)
  
  # Perform ANOVA for each Feature_ID
  anova_results <- data %>%
    group_by(ID) %>%
    summarize(p_value = {
      model <- lm(Normalized_Value ~ Group, data = .)  
      anova(model)$`Pr(>F)`[1]  
    }, .groups = "drop")  # Drop grouping after summarization
  
  # Filter for significant results (p-value < 0.05)
  significant_results <- anova_results %>%
    filter(p_value < 0.05) %>%
    mutate(log_p_value = -log10(p_value))  # Log-transform p-values
  
  # Plot the results
  aplot <- ggplot(significant_results, aes(x = ID, y = log_p_value)) +
    geom_point() +
    labs(title = "Significantly Different Intensities Between Groups (ANOVA)",
         x = "Feature ID",
         y = "Negative Log10 P-value") +
    theme_minimal() +
    # ---- Add dashed lines ----
    # Horizontal line at p = 0.05 => -log10(0.05) ~ 1.3
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    # Vertical lines at log2FC = ±1 (example cutoffs)
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")
  
  
  # Display the plot
  print(aplot)
  
  # Return the anova_results dataframe
  return(anova_results)
} #Need to work on this function- do not use yet
revert_lcms_features <- function(stats_lcms_data) {
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  # Create the combined 'mz_RT' string from mz and RT values
  stats_lcms_data <- stats_lcms_data %>%
    mutate(mz_RT = paste(mz, RT, sep = "/"))
  
  # Pivot the long-format data back into wide format.
  # Use ID and mz_RT as identifiers, and Sample as the column names.
  wide_data <- stats_lcms_data %>%
    select(ID, mz_RT, Sample, Value) %>%
    pivot_wider(names_from = Sample, values_from = Value) %>%
    arrange(ID)
  
  # Remove the ID column since the original file did not include it
  wide_data <- wide_data %>% select(-ID)
  
  # Extract unique grouping information for each sample.
  grouping <- stats_lcms_data %>%
    distinct(Sample, Group)
  
  # Determine the sample column order (all columns except the first, mz_RT)
  sample_columns <- colnames(wide_data)[-1]
  
  # Create a grouping row.
  # The first column is left blank, and subsequent cells contain the group label
  # corresponding to each sample (in the order of sample_columns).
  grouping_row <- c("" , grouping$Group[match(sample_columns, grouping$Sample)])
  names(grouping_row) <- colnames(wide_data)
  
  # Convert all columns to character so that the grouping row can be combined.
  wide_data <- wide_data %>%
    mutate(across(everything(), as.character))
  
  # Convert grouping_row to a one-row tibble using as.list
  grouping_df <- as_tibble(as.list(grouping_row))
  
  # Combine the grouping row with the wide_data
  output_data <- bind_rows(grouping_df, wide_data)
  
  return(output_data)
}


# Rework the data to preform stats; normalize if pelletmasses or standards are provided
stats_lcms_data <- extract_lcms_features(lcms_data, pellet_mass_file = pellet_mass_excel)
stats_lcms_data <- remove_extraction_blank(stats_lcms_data, blank_threshold = 10, value_col = "Value") %>%
  filter(Group != "EB")%>%
  filter(Value > 5000)
stats_lcms_data <- scale_data(stats_lcms_data, value_col = "Value", scale_method = "auto")
stats_lcms_data<-  normalize_data(stats_lcms_data,norm_method = "mean", value_col = "Scaled_Value")

length(unique(stats_lcms_data$ID))

#PCA
pca(stats_lcms_data)

ggsave('PCA_plot_EB_Removed.pdf')

# T_Test
t_test_results <- t_test(stats_lcms_data, sig = 0.05, group2, group3)

#fold cange with volcano plot
volcano_output <- volcano_plot(stats_lcms_data, group2, group3, n=1, value_col = "Value")

sig_peaks <- volcano_output %>%
  filter(
    p_value < 0.1,
    fold_change < 0.5
  )

#ANOVA
preform_anova <- anova(stats_lcms_data) #Need to work on this function- do not use yet

metaboanaylist_export <- revert_lcms_features(blank_removed) 
write.csv(metaboanaylist_export, file.path( "metaboanalyst_post_workup.csv"), row.names = FALSE)



library(xcms)
library(MSnbase)
library(dplyr)
library(RColorBrewer)

# — 0) pick which groups you actually want to plot —
keep_groups <- c("WT", "yumC","ispC")   # ← edit this as you like

# — 1) get your files & derive group labels —
parent_dir   <- '/Users/dmcbee/Desktop/Current MS Data/XCMS_Script'
setwd(parent_dir)
files        <- list.files(pattern="\\.mzML$", full.names=TRUE, recursive=TRUE)
group_labels <- basename(dirname(files))   # one label per file

# — 2) filter to only your chosen groups —
sel          <- group_labels %in% keep_groups
files_sel    <- files[sel]
group_labels <- group_labels[sel]         # now only WT, yumC, ispC, etc.

# — 3) build your palette & map colours —
uniq_groups  <- keep_groups               # in the order you want them in the legend
palette_cols <- brewer.pal(n=length(uniq_groups), name="Set1")
names(palette_cols) <- uniq_groups
group_colors <- palette_cols[group_labels]  # one colour per file

# — 4) load *only* the selected files —
raw_data     <- readMSData(files_sel, msLevel=1, mode="onDisk")

# — 5) make sure your RT is in seconds —
sig_peaks    <- sig_peaks %>% rename(RT_sec = RT)

# — 6) plotting loop (unchanged except for using *_sel variables) —
mz_tol_ppm    <- 10
rt_window_sec <- 60
ppm_to_da     <- function(mz, ppm) mz * ppm / 1e6

for (i in seq_len(nrow(sig_peaks))) {
  mz_val <- sig_peaks$mz[i]
  rt_val <- sig_peaks$RT[i]
  tol_da  <- ppm_to_da(mz_val, mz_tol_ppm)
  
  eic <- MSnbase::chromatogram(raw_data,
                               mz = c(mz_val - tol_da, mz_val + tol_da),
                               rt = c(rt_val - rt_window_sec,
                                      rt_val + rt_window_sec))
  
  plot(eic,
       type = "l",
       col  = group_colors,
       lwd  = 2,
       xlab = "Retention Time (s)",
       ylab = "Intensity",
       main = sprintf("Feature ID %d — m/z %.4f ±%d ppm",
                      sig_peaks$ID[i], mz_val, mz_tol_ppm))
  
  legend("topright",
         legend = uniq_groups,
         col    = palette_cols,
         lty    = 1,
         lwd    = 2,
         bty    = "n")
}


