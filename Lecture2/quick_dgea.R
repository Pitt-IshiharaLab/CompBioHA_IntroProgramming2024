# Install packages if you don't have them already.
BiocManager::install("edgeR")
install.packages('cpm')

# Load necessary libraries.
load("edgeR")
load("cpm")

setwd("Downloads")  # Move to whichever directory the raw counts data is in.
dir()               # Check that you're in the right directory.

# Load in data
data_raw <- read.csv("airway_scaledcounts.csv", header = TRUE)
head(data_raw)

# Clean up unnecessary columns
data_clean <- data_raw[,2:(ncol(data_raw))]
head(data_clean)

# Load in metadata
metadata <- read.csv("airway_metadata.csv", header = TRUE)
metadata

# Use the control/treated column as the group
y <- edgeR::DGEList(counts = data_clean, group = metadata[2][[1]])
# Differential Gene Expression Analysis 
y <- edgeR::calcNormFactors(y)
y$samples
y <- edgeR::estimateDisp(y)
et <- edgeR::exactTest(y)
results_edgeR <- edgeR::topTags(et, n = nrow(data_clean), sort.by = "none")

head(results_edgeR$table)                       # Check DGEA results.
write.csv(results_edgeR,"results_edgeR.csv")    # Export DGEA results.

# The resulting .csv may not have ENSGENE ids in the index column. You can fix
# this by opening it in Excel and copy-pasting the ids from the raw counts csv
# into the new one.