# Load necessary library
library(dplyr)

# Read the input file
data <- read.table("leafcutter_ds_cluster_significance.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 1. Extract rows where p.adjust is NA
padj_na <- data %>%
  filter(is.na(p.adjust)) %>%
  select(cluster, genes)

# 2. Extract rows where loglr != 0 and p.adjust < 0.05 (significant)
significant <- data %>%
  filter(!is.na(loglr) & loglr != 0 & !is.na(p.adjust) & p.adjust < 0.05) %>%
  select(cluster, genes)

# 3. Extract rows where loglr != 0 and p.adjust >= 0.05 (non-significant)
non_significant <- data %>%
  filter(!is.na(loglr) & loglr != 0 & !is.na(p.adjust) & p.adjust >= 0.05) %>%
  select(cluster, genes)

# Write the extracted data to separate files
write.table(padj_na, "padj_na.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(significant, "significant.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(non_significant, "non_significant.txt", sep = "\t", row.names = FALSE, quote = FALSE)

