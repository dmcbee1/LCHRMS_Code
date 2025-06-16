# ========== Libraries ==========
if (!requireNamespace("xcms", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("xcms")
}
if (!requireNamespace("MSnbase", quietly = TRUE)) {
  BiocManager::install("MSnbase")
}
if (!requireNamespace("mzR", quietly = TRUE)) {
  BiocManager::install("mzR")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("MsExperiment", quietly = TRUE)) {
  install.packages("MsExperiment")
}

library(xcms)
library(MSnbase)
library(mzR)
library(MsExperiment)

# ========== Parameters ==========
parent_dir <- '/Users/dmcbee/Desktop/Current MS Data/XCMS_Script'
output_dir <- file.path(parent_dir, "XCMS_Results")  # Directory for results

# CentWave parameters for peak detection
ppm <- 10
peakwidth <- c(2, 40)
snthresh <- 2
noise <- 250  # ensure noise is defined

# Retention time correction
rt_correction <- TRUE
rt_correction_method <- "Obiwarp"

# Peak grouping parameters
group_bandwidth <- 5

# ========== Script ==========
# Set working directory
setwd(parent_dir)

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
cat("All outputs will be saved in: ", output_dir, "\n")

# List all mzML files in the directory recursively (so that files in subfolders are included)
files <- list.files(path = parent_dir, pattern = "\\.mzML$", full.names = TRUE, recursive = TRUE)

# Create a vector of group labels from the parent folder name for each file
group_labels <- basename(dirname(files))

# Update metadata to include unique groups (if desired)
metadata <- list(
  ParentDirectory = parent_dir,
  OutputDirectory = output_dir,
  CentWaveParam = list(ppm = ppm, peakwidth = peakwidth, snthresh = snthresh, noise = noise),
  RetentionTimeCorrection = if (rt_correction) rt_correction_method else "None",
  GroupingParam = list(sampleGroups = unique(group_labels), bandwidth = group_bandwidth),
  FillMissingPeaks = TRUE,
  InputFiles = paste(basename(files), collapse = "; ")
)

# Check if files are detected
if (length(files) == 0) {
  stop("No mzML files found in the specified directory.")
}

# Create XCMSnExp object
cat("Reading data files...\n")
raw_data <- readMSData(files, mode = "onDisk")

# Peak detection
cat("Performing peak detection...\n")
cwp <- CentWaveParam(ppm = ppm, peakwidth = peakwidth, snthresh = snthresh, noise = noise, fitgauss = TRUE)
serialParam <- SerialParam()

# Force Single Core Processing
register(SerialParam())
xdata <- findChromPeaks(raw_data, param = cwp, BPPARAM = serialParam)

# Retention time correction
register(SerialParam())  # Force serial processing on macOS
if (rt_correction && length(files) > 1) {
  cat("Correcting retention time...\n")
  xdata <- adjustRtime(xdata, param = ObiwarpParam())
} else if (length(files) > 1) {
  cat("Skipping retention time correction as it is disabled.\n")
} else {
  cat("Only one file detected. Skipping retention time correction.\n")
}

# Group peaks across samples using group labels from folder names
cat("Grouping peaks...\n")
grouped_peaks <- groupChromPeaks(
  xdata,
  param = PeakDensityParam(sampleGroups = group_labels, bw = group_bandwidth)
)

# Fill in missing peaks
cat("Filling in missing peaks...\n")
filled_peaks <- fillChromPeaks(grouped_peaks)

# Extract grouped feature metadata
cat("Extracting grouped feature metadata...\n")
feature_definitions <- featureDefinitions(filled_peaks)

# Extract intensity matrix for grouped features
cat("Extracting feature intensities...\n")
intensity_matrix <- featureValues(filled_peaks, method = "maxint")

# Generate MetaboSeek-like output
metaboseek_output <- data.frame(
  mz = feature_definitions$mzmed,
  mzmin = feature_definitions$mzmin,
  mzmax = feature_definitions$mzmax,
  RT = feature_definitions$rtmed,
  rtmin = feature_definitions$rtmin,
  rtmax = feature_definitions$rtmax,
  Intensity = rowSums(intensity_matrix, na.rm = TRUE),
  Max_Intensity = apply(intensity_matrix, 1, max, na.rm = TRUE),
  SN = feature_definitions$npeaks
)

# Save the output as a CSV file
write.csv(metaboseek_output, file.path(output_dir, "feature-table.csv"), row.names = FALSE)
cat("MetaboSeek-like output saved as 'feature-table.csv'.\n")

# Save metadata as a CSV file
metadata_df <- as.data.frame(t(unlist(metadata)), stringsAsFactors = FALSE)
write.csv(metadata_df, file.path(output_dir, "metadata.csv"), row.names = FALSE)
cat("Metadata saved as 'metadata.csv'.\n")

# Generate MetaboAnalyst-ready output
cat("Generating MetaboAnalyst-ready output...\n")
metaboanalyst_data <- data.frame(
  Sample = paste0(feature_definitions$mzmed, "/", feature_definitions$rtmed),
  intensity_matrix
)
colnames(metaboanalyst_data) <- c("Sample", basename(files))

# Add sample group row (if needed, you might modify this to include group info)
sample_groups <- c("Label", group_labels)
metaboanalyst_data <- rbind(sample_groups, metaboanalyst_data)

# Save MetaboAnalyst-formatted data
write.csv(metaboanalyst_data, file.path(output_dir, "metaboanalyst_table.csv"), row.names = FALSE)
cat("MetaboAnalyst-ready output saved as 'metaboanalyst_table.csv'.\n")

# Optional: Plot data for quality check
cat("Generating diagnostic plots...\n")
# Plot TIC for each sample
tic <- chromatogram(raw_data, aggregationFun = "sum")
png(file.path(output_dir, "tic_plot.png"))
plot(tic, main = "Total Ion Chromatogram (TIC)")
dev.off()

# Plot retention time correction (if applicable)
if (rt_correction && length(files) > 1) {
  png(file.path(output_dir, "rt_correction.png"))
  plotAdjustedRtime(filled_peaks)
  dev.off()
}

cat("XCMS processing completed. Results saved in 'XCMS_Results' folder.\n")

#### Data Checks ######

# Define your filtering ranges (adjust these values as needed)
rtr <- c(305, 340)         # Retention time range (in seconds)
mzr <- c(652.3900, 652.4100) # m/z range

# Extract sample names (full file paths) from raw_data and derive group labels
# Each group is assumed to be the parent folder's name.
sampleGroups <- basename(dirname(fileNames(raw_data)))

# Create a color palette for the groups
unique_groups <- unique(sampleGroups)
group_colors <- setNames(rainbow(length(unique_groups)), unique_groups)

#### Data Checks ######

# Define your filtering ranges (adjust these values as needed)
rtr <- c(40, 45)         # Retention time range (in seconds)
mzr <- c(774.5700, 774.5900) # m/z range

# Derive group labels from the original 'files' vector (each file’s parent folder is its group)
group_labels <- basename(dirname(files))
print(unique(group_labels))  # Check the unique groups

# Create a color palette for the groups based on these labels
unique_groups <- unique(group_labels)
group_colors <- setNames(rainbow(length(unique_groups)), unique_groups)

# ----- Plot Raw Chromatograms -----
# Extract the chromatogram for the specified m/z and RT ranges
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
# Plot raw chromatograms using group colors: each sample is colored by its group.
plot(chr_raw, col = group_colors[group_labels], main = "Raw Chromatograms")
grid()

# ----- Plot Extracted Ion Chromatogram (XIC) -----
# Filter the raw data and plot the XIC with colors mapped by group.
raw_data |>
  filterRt(rt = rtr) |>
  filterMz(mz = mzr) |>
  plot(type = "XIC", col = group_colors[group_labels], main = "XIC Plot")

# ----- Compare Raw and Adjusted Chromatograms -----
# Set up a 2-row plotting area for comparison
par(mfrow = c(2, 1))
# Plot the raw chromatograms again with group-specific colors
plot(chr_raw, col = group_colors[group_labels], main = "Raw Chromatograms")
grid()

# Extract the chromatogram from the adjusted object (xdata) and plot it
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
plot(chr_adj, peakType = "none", col = group_colors[group_labels],
     main = "Adjusted Chromatograms")
grid()

# ----- Peak Density Plot -----
# Define PeakDensityParam using the group_labels vector (each sample is now assigned its group)
pdp <- PeakDensityParam(
  sampleGroups = group_labels, bw = group_bandwidth
)

# Extract the chromatogram for the specified m/z range from the grouped peaks object
chr_mzr <- chromatogram(grouped_peaks, mz = mzr)

# Plot the peak density using the group colors
plotChromPeakDensity(
  chr_mzr,
  col = group_colors[group_labels],
  param = pdp,
  peakPch = 16,
  main = "Peak Density by Group"
)


## Calculate m/z width of features for your data
mzw <- feature_definitions$mzmax - feature_definitions$mzmin

## Plot the median m/z vs. m/z width for your features
plot(feature_definitions$mzmed, mzw,
     xlab = "m/z", ylab = "m/z width", pch = 21,
     col = "#ff000020", bg = "#ff000010",
     main = "Feature m/z Width")


intensity_matrix <- featureValues(filled_peaks, method = "maxint")
## Extract the features and log2 transform them
## Extract the feature intensities and log₂ transform them
ft_ints <- log2(intensity_matrix)
# Remove features with any missing values
ft_ints <- ft_ints[complete.cases(ft_ints), ]

## Perform PCA (transposing so that samples become observations)
pc <- prcomp(t(ft_ints), center = TRUE)
pcSummary <- summary(pc)

## Prepare the sample colors from your pre-defined group_colors and group_labels.
# (Make sure that the length of group_labels equals the number of samples)
sample_colors <- group_colors[group_labels]

## Plot the PCA (PC1 vs. PC2)
plot(pc$x[, 1], pc$x[, 2],
     pch = 21, 
     main = "PCA of Feature Intensities",
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100, digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100, digits = 3), " % variance"),
     col = "black",           # Border color for points
     bg = sample_colors,        # Fill colors according to groups
     cex = 2)
grid()

## Add sample labels (using column names from the intensity matrix)
text(pc$x[, 1], pc$x[, 2], labels = colnames(intensity_matrix),
     col = "darkgrey", pos = 3, cex = 0.8)








