# Load necessary libraries if needed
# library(dplyr)

# Load dxr1 from the saved RDS file
output_path <- "/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/alb_x_hyp_Testis/"
dxr1 <- readRDS(file = paste0(output_path, "dxr1.rds"))

# Extract significant exons
significant_exons <- dxr1[!is.na(dxr1$padj) & dxr1$padj < 0.05 & dxr1$log2fold_Ficedula_hypoleuca_Ficedula_albicollis != 0, c("groupID", "featureID")]

# Extract NA exons
na_exons <- dxr1[is.na(dxr1$padj), c("groupID", "featureID")]

# Extract non-significant exons
non_significant_exons <- dxr1[!is.na(dxr1$padj) & dxr1$padj >= 0.05, c("groupID", "featureID")]

# Identify other exons
# This assumes that other exons are those not present in significant, NA, or non-significant exons
significant_and_na_ids <- unique(c(significant_exons$featureID, na_exons$featureID))
non_significant_ids <- unique(non_significant_exons$featureID)
other_exons <- dxr1[!(dxr1$featureID %in% c(significant_and_na_ids, non_significant_ids)), c("groupID", "featureID")]

# Diagnostic checks
cat("Number of significant exons:", nrow(significant_exons), "\n")
cat("Number of NA exons:", nrow(na_exons), "\n")
cat("Number of non-significant exons:", nrow(non_significant_exons), "\n")
cat("Number of other exons:", nrow(other_exons), "\n")

# Write to files
write.table(significant_exons, file = paste0(output_path, "significant_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(na_exons, file = paste0(output_path, "na_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(non_significant_exons, file = paste0(output_path, "non_significant_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(other_exons, file = paste0(output_path, "other_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

-------------------------------------------------------------------------------


# Load necessary libraries if needed
# library(dplyr)

# Load dxr1 from the saved RDS file
output_path <- "/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/alb_x_f1_Testis/"
dxr1 <- readRDS(file = paste0(output_path, "dxr1.rds"))

# Extract significant exons
significant_exons <- dxr1[!is.na(dxr1$padj) & dxr1$padj < 0.05 & dxr1$log2fold_Ficedula_albicollis_F1_hybrid != 0, c("groupID", "featureID")]

# Extract NA exons
na_exons <- dxr1[is.na(dxr1$padj), c("groupID", "featureID")]

# Extract non-significant exons
non_significant_exons <- dxr1[!is.na(dxr1$padj) & dxr1$padj >= 0.05, c("groupID", "featureID")]

# Identify other exons
# This assumes that other exons are those not present in significant, NA, or non-significant exons
significant_and_na_ids <- unique(c(significant_exons$featureID, na_exons$featureID))
non_significant_ids <- unique(non_significant_exons$featureID)
other_exons <- dxr1[!(dxr1$featureID %in% c(significant_and_na_ids, non_significant_ids)), c("groupID", "featureID")]

# Diagnostic checks
cat("Number of significant exons:", nrow(significant_exons), "\n")
cat("Number of NA exons:", nrow(na_exons), "\n")
cat("Number of non-significant exons:", nrow(non_significant_exons), "\n")
cat("Number of other exons:", nrow(other_exons), "\n")

# Write to files
write.table(significant_exons, file = paste0(output_path, "significant_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(na_exons, file = paste0(output_path, "na_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(non_significant_exons, file = paste0(output_path, "non_significant_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(other_exons, file = paste0(output_path, "other_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

--------------------------------------------------------------------------------------------


# Load necessary libraries if needed
# library(dplyr)

# Load dxr1 from the saved RDS file
output_path <- "/home/sfadi/Bureau/DSA_all_species_Testis/DEXSeq_output/hyp_x_f1_Testis/"
dxr1 <- readRDS(file = paste0(output_path, "dxr1.rds"))

# Extract significant exons
significant_exons <- dxr1[!is.na(dxr1$padj) & dxr1$padj < 0.05 & dxr1$log2fold_Ficedula_hypoleuca_F1_hybrid != 0, c("groupID", "featureID")]

# Extract NA exons
na_exons <- dxr1[is.na(dxr1$padj), c("groupID", "featureID")]

# Extract non-significant exons
non_significant_exons <- dxr1[!is.na(dxr1$padj) & dxr1$padj >= 0.05, c("groupID", "featureID")]

# Identify other exons
# This assumes that other exons are those not present in significant, NA, or non-significant exons
significant_and_na_ids <- unique(c(significant_exons$featureID, na_exons$featureID))
non_significant_ids <- unique(non_significant_exons$featureID)
other_exons <- dxr1[!(dxr1$featureID %in% c(significant_and_na_ids, non_significant_ids)), c("groupID", "featureID")]

# Diagnostic checks
cat("Number of significant exons:", nrow(significant_exons), "\n")
cat("Number of NA exons:", nrow(na_exons), "\n")
cat("Number of non-significant exons:", nrow(non_significant_exons), "\n")
cat("Number of other exons:", nrow(other_exons), "\n")

# Write to files
write.table(significant_exons, file = paste0(output_path, "significant_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(na_exons, file = paste0(output_path, "na_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(non_significant_exons, file = paste0(output_path, "non_significant_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(other_exons, file = paste0(output_path, "other_exons.txt"), sep = "\t", row.names = FALSE, quote = FALSE)



