# Lab 1: DNA Barcoding Sashimi
setwd("C:/Users/nelin/Desktop/GWU Information/GWU Fall 2026/Princip_Bioinformatics/LAB1/scripts")

# Read FASTA
fasta_lines <- readLines("C:/Users/nelin/Desktop/GWU Information/GWU Fall 2026/Princip_Bioinformatics/LAB1/data/query_sashimi.fasta")

# Drop the header line, join the 7 sequence lines into one string, convert to uppercase
query_seq <- toupper(paste(fasta_lines[-1], collapse = ""))

# Count the nucleotides, should be 481
nchar(query_seq) #481

# Reverse complement function
revcomp <- function(s) {
  s <- toupper(s) # Standardize case
  comp <- chartr("ACGT", "TGCA", s) # Swap each base for its partner: A<->T, C<->G
  paste(rev(strsplit(comp, "")[[1]]), collapse = "") # Split into single letters, reverse the order, then rejoin
}

# Primer Verification
# Primers picked by Primer3Plus (product size forced to 481-481)
fwd_primer <- "CGGGCAGAACTAAGCCAGC" # Left primer, 19 bp, Start 1
rev_primer <- "ATAGAAGGAGTAGTACTGCGGT" # Right primer, 22 bp, Start 481

# Store each primer's length so the code adapts if the primers change
fwd_len <- nchar(fwd_primer)
rev_len <- nchar(rev_primer)
seq_len <- nchar(query_seq)

# Forward primer check:
# Identical to the start of the seq, take bases 1-19
first_bases <- substr(query_seq, 1, fwd_len)
first_bases
identical(first_bases, fwd_primer) # Should print TRUE

# Reverse primer check:
# The reverse complement of the end of the sequence, take last 22 bases (460-481) and reverse-complement them
last_bases <-substr(query_seq, seq_len - rev_len + 1, seq_len)
last_bases
revcomp(last_bases)
identical(revcomp(last_bases), rev_primer)



# ----- Build combined FASTA of query + reference sequences -----
data_dir <- "C:/Users/nelin/Desktop/GWU Information/GWU Fall 2026/Princip_Bioinformatics/LAB1/data"

# Helper function: read a single-sequence FASTA and return the sequence as one uppercase string
read_one_seq <-function(path) {
  lines <- readLines(path) # each line of the file becomes one vector element
  s <- toupper(paste(lines[-1], collapse = "")) # drop the ">" header, join the sequence lines
  stopifnot(length(s) == 1) # safety check
  s # <- return the sequence 
}

# --- Extract COX1 from the Hamachi mitochondrial genome ---
hamachi_mito <- read_one_seq(file.path(data_dir, "s_quinqueradiata.fasta"))
nchar(hamachi_mito) #16537

cox1_start <- 5508 # Start coordinate of COX1 from the GenBank FEATURES section
cox1_end <- 7058 # End coordinate of COX1

hamachi_cox1 <- substr(hamachi_mito, cox1_start, cox1_end) # Keep only the COX1 gene
nchar(hamachi_cox1) # 1550

# --- Collect all sequences, each with short label ---
# Short labels for ape - uses the header as the row/column name in the distance matrix
# Full GenBank header would make the matrix unreadable
seqs <- c(
  Query_sashimi = query_seq,
  Tilapia_Oniloticus_MN756471 = read_one_seq(file.path(data_dir, "o_niloticus.fasta")),
  Hirame_Polivaceus_MH032483  = read_one_seq(file.path(data_dir, "p_olivaceus.fasta")),
  Sake_Ssalar_MZ407797        = read_one_seq(file.path(data_dir, "s_salar.fasta")),
  Hamachi_Squinq_NC016868     = hamachi_cox1,
  Maguro_Talbacares_PQ812466  = read_one_seq(file.path(data_dir, "t_albacares.fasta"))
)
nchar(seqs) # Length of each sequence
# Result:

# --- Write into one FASTA file with Unix line endings ---
# paste0() builds ">label" + line break + sequence for each entry
fasta_text <- paste0(">", names(seqs), "\n", seqs)
out_path <- file.path(data_dir, "sashimi_sequences_v2.fasta")
con <- file(out_path, open = "wb") # "wb" = write binary, so R doesn't covert \n into \r\n
writeLines(fasta_text, con, sep ="\n")
close(con) # closing connection

# ----- Pairwise genetic diversity -----
# Step 1: Load alignment
library(ape) # read.dna(), dist.dna()
library(pegas) # nuc.div()

results_dir <- "C:/Users/nelin/Desktop/GWU Information/GWU Fall 2026/Princip_Bioinformatics/LAB1/results"

# Versions
R.version.string
packageVersion("ape")
packageVersion("pegas")

# Read the MAFFT alignment; because every aligned sequence has the same length,
# ape stores it as a matrix: row = sequences, columns = alignment positions
aln <- read.dna(file.path(results_dir, "sashimi_aligned.fasta"), format = "fasta")
dim(aln) # should print 6 1551
labels(aln) # six sequence names 

# Step 2: Analysis A - Full alignment
# p-distance: the proportion of compared sites that differ between two sequences
# pairwise.deletion = TRUE: for each pair, skip any column where either sequence has a gap
D_full <- dist.dna(aln, model = "raw", pairwise.deletion = TRUE)
round(as.matrix(D_full), 4)

# The same comparison as a RAW COUNT of differing sites (model = "N")
N_full <- dist.dna(aln, model = "N", pairwise.deletion = TRUE)
as.matrix(N_full)

# Helper function: how many sites were actually compared for each pair?
sites_compared <- function(dna) {
  m <- as.character(dna) # convert to a matrix of letters; gaps appear as "-"
  has_base <- m != "-"
  n <- nrow(m)
  out <- matrix(NA, n, n, dimnames = list(rownames(m), rownames(m)))
  for (i in 1:n) {
    for (j in 1:n) {
      out[i,j] <- sum(has_base[i, ] & has_base[j, ]) # sites where BOTH have a base
    }
  }
  out
}
sites_compared(aln)
      
# Step 3: Analysis B - Trim alignment to query region
query_row <- as.character(aln["Query_sashimi", ])
q_cols <- which(query_row != "-")
range(q_cols) # first and last query columns (expect 112 592)

# Keep only those columns for all 6 sequences
aln_trim <- aln[, min(q_cols):max(q_cols)]
dim(aln_trim) # should print 6 481

D_trim <- dist.dna(aln_trim, model = "raw", pairwise.deletion = TRUE)
round(as.matrix(D_trim), 4)

N_trim <- dist.dna(aln_trim, model = "N", pairwise.deletion = TRUE)
as.matrix(N_trim)

sites_compared(aln_trim) # every pair should now be 481

# Step 4: Nucleotide diversity and saving
# pi = average proportion of differing sites across all pairs of sequences
pi_trim <- nuc.div(aln_trim)
pi_trim

# Save matrices as CSV files
write.csv(round(as.matrix(D_trim), 4), file.path(results_dir, "pdist_trimmed.csv"))
write.csv(round(as.matrix(D_full), 4), file.path(results_dir, "pdist_full.csv"))
write.csv(as.matrix(N_full),           file.path(results_dir, "count_full.csv"))
write.csv(as.matrix(N_trim), file.path(results_dir, "count_trimmed.csv"))