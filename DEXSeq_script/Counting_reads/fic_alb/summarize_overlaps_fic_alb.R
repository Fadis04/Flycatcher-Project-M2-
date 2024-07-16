library(Rsamtools)
library(GenomicAlignments)
library(SummarizedExperiment)

# Load flattened annotation
flattenedAnnotation <- readRDS("/beegfs/home/fslimi/scripts/DEXSeq_prep/flattenedAnnotation.rds")

# Set the path to your BAM files
bam_path <- "/beegfs/data/fslimi/fic_bam/fic_alb_bam/"

# Get a list of all BAM files in the specified path
bam_files <- list.files(path = bam_path, pattern = "\\.bam$", full.names = TRUE)

# Check if BAM files are found
if (length(bam_files) == 0) {
  stop("No BAM files found in the specified path.")
}

# Create a BamFileList object
bamFiles <- BamFileList(bam_files)

# Summarize overlaps
se <- summarizeOverlaps(
  features = flattenedAnnotation,
  reads = bamFiles,
  mode = "Union",
  singleEnd = FALSE,
  fragments = TRUE,
  ignore.strand = TRUE
)

# Extract the read counts matrix
read_counts_fic_alb <- assay(se)

# Save read counts to a file in the desired path
write.table(read_counts_fic_alb , file = "/beegfs/data/fslimi/read_counts_fic_alb.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

