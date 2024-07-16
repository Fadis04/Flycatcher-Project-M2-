# List of input read count files
input_files <- list(
  f1_hyb = "/beegfs/data/fslimi/read_counts_f1_hyb.txt",
  fic_alb = "/beegfs/data/fslimi/read_counts_fic_alb.txt",
  fic_hyp = "/beegfs/data/fslimi/read_counts_fic_hyp.txt"
)

# Initialize a list to store read counts
all_read_counts <- list()

# Loop through each file, read the data and store in the list
for (sample in names(input_files)) {
  file_path <- input_files[[sample]]
  
  # Read the read count file
  read_counts <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t")
  
  # Store the read counts in the list
  all_read_counts[[sample]] <- read_counts
}

# Combine all read counts into a single data frame
final_read_counts <- do.call(cbind, all_read_counts)

# Save the final combined read counts to a file
write.table(final_read_counts, file = "/beegfs/data/fslimi/final_read_counts_all_samples.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

