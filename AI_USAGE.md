# AI USAGE LOG
- **Tool:** 
- **What I was doing:**
- **What I asked (verbatim):**
- **Exact output (verbatim):**
- **Verification:**
## Entry 1
- **Tool:** Claude (Sonnet 5)
- **What I was doing:** Trying to write a script for analysis to read the FASTA file.
- **What I asked (verbatim):** "Can you help me write a script in R to analyze this fasta file for its: count, check the first 22 and the last 22 nt's to check against the Primer3Plus output manually, <insert what the script does and looks at>? (pasted the query_sashimi.fasta from my computer)
- **Exact output (verbatim:)** 
```r
# Lab 1: DNA Barcoding Sashimi
setwd("C:/Users/nelin/Desktop/GWU Information/GWU Fall 2026/Princip_Bioinformatics/LAB1/scripts")

# Read FASTA
fasta_lines <- readLines("data/query_sashimi.fasta")

# Drop the header line, join the 7 sequence lines into one string, convert to uppercase
query_seq <- toupper(paste(fasta_lines[-1], collapse = ""))

# Count the nucleotides, should be 481
nchar(query_seq)

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

# Store each primer's length so the code adapts if the primers change (if rerun for any reason)
fwd_len <- nchar(fwd_primer)
rev_len <- nchar(rev_primer)
seq_len <- nchar(query_seq)

# Forward primer check:
# Identical to the start of the seq, take bases 1-19
first_bases <- substr(query_seq, 1, fwd_len)
first_bases # "CGGGCAGAACTAAGCCAGC"
identical(first_bases, fwd_primer) # Should print TRUE

# Reverse primer check:
# The reverse complement of the end of the sequence, take last 22 bases (460-481) and reverse-complement them
last_bases <-substr(query_seq, seq_len - rev_len + 1, seq_len)
last_bases # "ACCGCAGTACTACTCCTTCTAT"
revcomp(last_bases) # "ATAGAAGGAGTAGTACTGCGGT"
identical(revcomp(last_bases), rev_primer) # Should print TRUE
```
-**Verification:** I typed out manually and ran the code which printed true for `identical()` and then visually compared the primers to the 481 nucleotide sequence. They all matched. 

## Entry 2
- **Tool:** Claude (Sonnet 5)
- **What I was doing:** After collecting all FASTA documents from the BLAST searches for the sample, Hirame (P. Olivaceus), Sake (S. salar), Hamachi (S. quinqueradiata), Maguro (T. orientalis), and Nile Tilapia (O. niloticus - the top match), I needed to make a combined FASTA to align.
- **What I asked (verbatim):** "Can you write me a script that combines these FASTA into a combined FASTA to align?"
- **Exact output (verbatim):**
# ---- Build the combined FASTA of query + reference sequences ----
data_dir <- "C:/Users/nelin/Desktop/GWU Information/GWU Fall 2026/Princip_Bioinformatics/LAB1/data"

# Helper function: read a single-sequence FASTA and return the sequence as one uppercase string
read_one_seq <- function(path) {
  lines <- readLines(path)                      # each line of the file becomes one vector element
  toupper(paste(lines[-1], collapse = ""))      # drop the ">" header, join the sequence lines
}

# --- Extract COX1 from the Hamachi mitochondrial genome ---
hamachi_mito <- read_one_seq(file.path(data_dir, "s_quinqueradiata.fasta"))
nchar(hamachi_mito)                              # should print 16537

cox1_start <- ____    # start coordinate of COX1 from the GenBank FEATURES section
cox1_end   <- ____    # end coordinate of COX1

hamachi_cox1 <- substr(hamachi_mito, cox1_start, cox1_end)  # keep only the COX1 gene
nchar(hamachi_cox1)                              # should be about 1,550

# --- Collect all sequences, each with a short label ---
# Short labels matter: ape uses the header as the row/column name in the distance matrix,
# so the full GenBank headers would make the matrix unreadable
seqs <- c(
  Query_sashimi              = query_seq,
  Tilapia_Oniloticus_MN756471 = read_one_seq(file.path(data_dir, "o_niloticus.fasta")),
  Hirame_Polivaceus_MH032483  = read_one_seq(file.path(data_dir, "p_olivaceus.fasta")),
  Sake_Ssalar_MZ407797        = read_one_seq(file.path(data_dir, "s_salar.fasta")),
  Hamachi_Squinq_NC016868     = hamachi_cox1,
  Maguro_Talbacares_PQ812466  = read_one_seq(file.path(data_dir, "t_albacares.fasta"))
)

nchar(seqs)   # length of each sequence; note how different they are

# --- Write everything into one FASTA file ---
# paste0 builds ">label" + line break + sequence for each entry
fasta_text <- paste0(">", names(seqs), "\n", seqs)
writeLines(fasta_text, file.path(data_dir, "sashimi_sequences.fasta"))
- **Verification:** I ran the code but the output was not as expected. See below (collapsed)
```r

 # --- Extract COX1 from the Hamachi mitochondrial genome ---
> hamachi_mito <- read_one_seq(file.path(data_dir, "s_quinqueradiata.fasta"))
> nchar(hamachi_mito) #16537
  [1] 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71
 [44] 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71
 [87] 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71
[130] 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71
[173] 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71
[216] 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 71 18  1
> 
> cox1_start <- 5508 # Start coordinate of COX1 from the GenBank FEATURES section
> cox1_end <- 7058 # End coordinate of COX1
> 
> hamachi_cox1 <- substr(hamachi_mito, cox1_start, cox1_end) # Keep only the COX1 gene
> nchar(hamachi_cox1) # 1550
  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 [65] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
[129] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
[193] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
> 
> # --- Collect all sequences, each with short label ---
> # Short labels for ape - uses the header as the row/column name in the distance matrix
> # Full GenBank header would make the matrix unreadable
> seqs <- c(
+   Query_sashimi = query_seq,
+   Tilapia_Oniloticus_MN756471 = read_one_seq(file.path(data_dir, "o_niloticus.fasta")),
+   Hirame_Polivaceus_MH032483  = read_one_seq(file.path(data_dir, "p_olivaceus.fasta")),
+   Sake_Ssalar_MZ407797        = read_one_seq(file.path(data_dir, "s_salar.fasta")),
+   Hamachi_Squinq_NC016868     = hamachi_cox1,
+   Maguro_Talbacares_PQ812466  = read_one_seq(file.path(data_dir, "t_albacares.fasta"))
+ )
> nchar(seqs) # Length of each sequence
               Query_sashimi Tilapia_Oniloticus_MN7564711 Tilapia_Oniloticus_MN7564712 Tilapia_Oniloticus_MN7564713 
                         481                           71                           71                           71 
Tilapia_Oniloticus_MN7564714 Tilapia_Oniloticus_MN7564715 Tilapia_Oniloticus_MN7564716 Tilapia_Oniloticus_MN7564717 
                          71                           71                           71                           71 
Tilapia_Oniloticus_MN7564718 Tilapia_Oniloticus_MN7564719  Hirame_Polivaceus_MH0324831  Hirame_Polivaceus_MH0324832 
                           6                            1                           71                           71 
 Hirame_Polivaceus_MH0324833  Hirame_Polivaceus_MH0324834  Hirame_Polivaceus_MH0324835  Hirame_Polivaceus_MH0324836 
                          71                           71                           71                           71 
 Hirame_Polivaceus_MH0324837  Hirame_Polivaceus_MH0324838  Hirame_Polivaceus_MH0324839 Hirame_Polivaceus_MH03248310 
                          71                           71                           71                           71 
Hirame_Polivaceus_MH03248311 Hirame_Polivaceus_MH03248312 Hirame_Polivaceus_MH03248313 Hirame_Polivaceus_MH03248314 
                          71                           71                           71                           71 
Hirame_Polivaceus_MH03248315 Hirame_Polivaceus_MH03248316 Hirame_Polivaceus_MH03248317 Hirame_Polivaceus_MH03248318 
                          71                           71                           71                           71 
Hirame_Polivaceus_MH03248319 Hirame_Polivaceus_MH03248320 Hirame_Polivaceus_MH03248321 Hirame_Polivaceus_MH03248322 
                          71                           71                           71                           71 
Hirame_Polivaceus_MH03248323 Hirame_Polivaceus_MH03248324        Sake_Ssalar_MZ4077971        Sake_Ssalar_MZ4077972 
                          12                            1                           71                           71 
       Sake_Ssalar_MZ4077973        Sake_Ssalar_MZ4077974        Sake_Ssalar_MZ4077975        Sake_Ssalar_MZ4077976 
                          71                           71                           71                           71 
       Sake_Ssalar_MZ4077977        Sake_Ssalar_MZ4077978        Sake_Ssalar_MZ4077979       Sake_Ssalar_MZ40779710 
                          71                           71                           71                            1 
    Hamachi_Squinq_NC0168681     Hamachi_Squinq_NC0168682     Hamachi_Squinq_NC0168683     Hamachi_Squinq_NC0168684 
                           0                            0                            0                            0 
    Hamachi_Squinq_NC0168685     Hamachi_Squinq_NC0168686     Hamachi_Squinq_NC0168687     Hamachi_Squinq_NC0168688 
                           0                            0                            0                            0 
    Hamachi_Squinq_NC0168689    Hamachi_Squinq_NC01686810    Hamachi_Squinq_NC01686811    Hamachi_Squinq_NC01686812 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686813    Hamachi_Squinq_NC01686814    Hamachi_Squinq_NC01686815    Hamachi_Squinq_NC01686816 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686817    Hamachi_Squinq_NC01686818    Hamachi_Squinq_NC01686819    Hamachi_Squinq_NC01686820 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686821    Hamachi_Squinq_NC01686822    Hamachi_Squinq_NC01686823    Hamachi_Squinq_NC01686824 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686825    Hamachi_Squinq_NC01686826    Hamachi_Squinq_NC01686827    Hamachi_Squinq_NC01686828 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686829    Hamachi_Squinq_NC01686830    Hamachi_Squinq_NC01686831    Hamachi_Squinq_NC01686832 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686833    Hamachi_Squinq_NC01686834    Hamachi_Squinq_NC01686835    Hamachi_Squinq_NC01686836 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686837    Hamachi_Squinq_NC01686838    Hamachi_Squinq_NC01686839    Hamachi_Squinq_NC01686840 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686841    Hamachi_Squinq_NC01686842    Hamachi_Squinq_NC01686843    Hamachi_Squinq_NC01686844 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686845    Hamachi_Squinq_NC01686846    Hamachi_Squinq_NC01686847    Hamachi_Squinq_NC01686848 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686849    Hamachi_Squinq_NC01686850    Hamachi_Squinq_NC01686851    Hamachi_Squinq_NC01686852 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686853    Hamachi_Squinq_NC01686854    Hamachi_Squinq_NC01686855    Hamachi_Squinq_NC01686856 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686857    Hamachi_Squinq_NC01686858    Hamachi_Squinq_NC01686859    Hamachi_Squinq_NC01686860 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686861    Hamachi_Squinq_NC01686862    Hamachi_Squinq_NC01686863    Hamachi_Squinq_NC01686864 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686865    Hamachi_Squinq_NC01686866    Hamachi_Squinq_NC01686867    Hamachi_Squinq_NC01686868 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686869    Hamachi_Squinq_NC01686870    Hamachi_Squinq_NC01686871    Hamachi_Squinq_NC01686872 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686873    Hamachi_Squinq_NC01686874    Hamachi_Squinq_NC01686875    Hamachi_Squinq_NC01686876 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686877    Hamachi_Squinq_NC01686878    Hamachi_Squinq_NC01686879    Hamachi_Squinq_NC01686880 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686881    Hamachi_Squinq_NC01686882    Hamachi_Squinq_NC01686883    Hamachi_Squinq_NC01686884 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686885    Hamachi_Squinq_NC01686886    Hamachi_Squinq_NC01686887    Hamachi_Squinq_NC01686888 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686889    Hamachi_Squinq_NC01686890    Hamachi_Squinq_NC01686891    Hamachi_Squinq_NC01686892 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686893    Hamachi_Squinq_NC01686894    Hamachi_Squinq_NC01686895    Hamachi_Squinq_NC01686896 
                           0                            0                            0                            0 
   Hamachi_Squinq_NC01686897    Hamachi_Squinq_NC01686898    Hamachi_Squinq_NC01686899   Hamachi_Squinq_NC016868100 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868101   Hamachi_Squinq_NC016868102   Hamachi_Squinq_NC016868103   Hamachi_Squinq_NC016868104 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868105   Hamachi_Squinq_NC016868106   Hamachi_Squinq_NC016868107   Hamachi_Squinq_NC016868108 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868109   Hamachi_Squinq_NC016868110   Hamachi_Squinq_NC016868111   Hamachi_Squinq_NC016868112 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868113   Hamachi_Squinq_NC016868114   Hamachi_Squinq_NC016868115   Hamachi_Squinq_NC016868116 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868117   Hamachi_Squinq_NC016868118   Hamachi_Squinq_NC016868119   Hamachi_Squinq_NC016868120 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868121   Hamachi_Squinq_NC016868122   Hamachi_Squinq_NC016868123   Hamachi_Squinq_NC016868124 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868125   Hamachi_Squinq_NC016868126   Hamachi_Squinq_NC016868127   Hamachi_Squinq_NC016868128 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868129   Hamachi_Squinq_NC016868130   Hamachi_Squinq_NC016868131   Hamachi_Squinq_NC016868132 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868133   Hamachi_Squinq_NC016868134   Hamachi_Squinq_NC016868135   Hamachi_Squinq_NC016868136 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868137   Hamachi_Squinq_NC016868138   Hamachi_Squinq_NC016868139   Hamachi_Squinq_NC016868140 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868141   Hamachi_Squinq_NC016868142   Hamachi_Squinq_NC016868143   Hamachi_Squinq_NC016868144 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868145   Hamachi_Squinq_NC016868146   Hamachi_Squinq_NC016868147   Hamachi_Squinq_NC016868148 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868149   Hamachi_Squinq_NC016868150   Hamachi_Squinq_NC016868151   Hamachi_Squinq_NC016868152 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868153   Hamachi_Squinq_NC016868154   Hamachi_Squinq_NC016868155   Hamachi_Squinq_NC016868156 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868157   Hamachi_Squinq_NC016868158   Hamachi_Squinq_NC016868159   Hamachi_Squinq_NC016868160 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868161   Hamachi_Squinq_NC016868162   Hamachi_Squinq_NC016868163   Hamachi_Squinq_NC016868164 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868165   Hamachi_Squinq_NC016868166   Hamachi_Squinq_NC016868167   Hamachi_Squinq_NC016868168 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868169   Hamachi_Squinq_NC016868170   Hamachi_Squinq_NC016868171   Hamachi_Squinq_NC016868172 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868173   Hamachi_Squinq_NC016868174   Hamachi_Squinq_NC016868175   Hamachi_Squinq_NC016868176 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868177   Hamachi_Squinq_NC016868178   Hamachi_Squinq_NC016868179   Hamachi_Squinq_NC016868180 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868181   Hamachi_Squinq_NC016868182   Hamachi_Squinq_NC016868183   Hamachi_Squinq_NC016868184 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868185   Hamachi_Squinq_NC016868186   Hamachi_Squinq_NC016868187   Hamachi_Squinq_NC016868188 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868189   Hamachi_Squinq_NC016868190   Hamachi_Squinq_NC016868191   Hamachi_Squinq_NC016868192 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868193   Hamachi_Squinq_NC016868194   Hamachi_Squinq_NC016868195   Hamachi_Squinq_NC016868196 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868197   Hamachi_Squinq_NC016868198   Hamachi_Squinq_NC016868199   Hamachi_Squinq_NC016868200 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868201   Hamachi_Squinq_NC016868202   Hamachi_Squinq_NC016868203   Hamachi_Squinq_NC016868204 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868205   Hamachi_Squinq_NC016868206   Hamachi_Squinq_NC016868207   Hamachi_Squinq_NC016868208 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868209   Hamachi_Squinq_NC016868210   Hamachi_Squinq_NC016868211   Hamachi_Squinq_NC016868212 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868213   Hamachi_Squinq_NC016868214   Hamachi_Squinq_NC016868215   Hamachi_Squinq_NC016868216 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868217   Hamachi_Squinq_NC016868218   Hamachi_Squinq_NC016868219   Hamachi_Squinq_NC016868220 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868221   Hamachi_Squinq_NC016868222   Hamachi_Squinq_NC016868223   Hamachi_Squinq_NC016868224 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868225   Hamachi_Squinq_NC016868226   Hamachi_Squinq_NC016868227   Hamachi_Squinq_NC016868228 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868229   Hamachi_Squinq_NC016868230   Hamachi_Squinq_NC016868231   Hamachi_Squinq_NC016868232 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868233   Hamachi_Squinq_NC016868234   Hamachi_Squinq_NC016868235   Hamachi_Squinq_NC016868236 
                           0                            0                            0                            0 
  Hamachi_Squinq_NC016868237   Hamachi_Squinq_NC016868238  Maguro_Talbacares_PQ8124661  Maguro_Talbacares_PQ8124662 
                           0                            0                           71                           71 
 Maguro_Talbacares_PQ8124663  Maguro_Talbacares_PQ8124664  Maguro_Talbacares_PQ8124665  Maguro_Talbacares_PQ8124666 
                          71                           71                           71                           71 
 Maguro_Talbacares_PQ8124667  Maguro_Talbacares_PQ8124668  Maguro_Talbacares_PQ8124669 Maguro_Talbacares_PQ81246610 
                          71                           71                           71                           26 
Maguro_Talbacares_PQ81246611 
                           1 
> # Result:
> 
> # --- Write into one FASTA file ---
> # paste0() builds ">label" + line break + sequence for each entry
> fasta_text <- paste0(">", names(seqs), "\n", seqs)
> writeLines(fasta_text, file.path(data_dir, "sashimi_sequences.fasta"))
```
- **What I asked next (verbatim):** (Copy and pasted above output in RStudio console + added the `sashimi_sequences.fasta`) "This was the output of the code, it does not look as expected. It appears to have returned the file lines as separate pieces. 
- **Exact output (verbatim):** 
No, this isn't correct yet. The good news is that all three problems you're seeing trace back to **one small error in the `read_one_seq()` helper function**, so one fix solves everything.

## What went wrong

**Clue 1:** `nchar(hamachi_mito)` should print a single number (16537). Instead it printed **238 numbers**, mostly 71. So `read_one_seq()` returned the file's lines as **separate pieces** instead of joining them into one string.

**Clue 2:** The number 71 is informative. NCBI writes FASTA files with **70 bases per line**, so each piece has one extra character. I checked the `sashimi_sequences.fasta` file you uploaded, and every sequence line ends with a **space**.

That combination almost always means the `paste()` line inside the function is missing `collapse =`, and probably reads something like this:

```r
paste(lines[-1], "")             # what it likely says: pastes "" onto EACH line, separated by a space
paste(lines[-1], collapse = "")  # what it should say: joins ALL lines into one string
```

Without `collapse =`, R treats `""` as a second thing to attach to each line, using `paste()`'s default separator, a space. So every line got a space added at the end and nothing was joined. Compare the line in your helper with the second one above.

**Everything else follows from that:**

- **Hamachi = 0s:** `substr()` asked each 71-character piece for characters 5,508–7,058. None of the pieces is that long, so every result came back empty.
- **Labels with extra numbers** (e.g., `MN7564711`, `MN7564712`): `c()` was given 9 pieces labeled "Tilapia...", so R numbered them 1, 2, 3, and so on.
- **The written FASTA is broken:** it contains 293 "sequences" instead of 6, so don't use it.

Your query was fine (481) because you built it earlier with the correct `collapse = ""`.

## The fix

Replace your helper with this version, which also includes a check that stops with an error if the problem ever comes back:

```r
# Helper function: read a single-sequence FASTA and return the sequence as ONE uppercase string
read_one_seq <- function(path) {
  lines <- readLines(path)                          # each line of the file = one vector element
  s <- toupper(paste(lines[-1], collapse = ""))     # drop header; collapse = "" joins ALL lines into one
  stopifnot(length(s) == 1)                         # safety check: must be exactly one string
  s
}
```

**Re-run the function definition first**, then everything below it (the Hamachi extraction, `seqs`, and the writing step). R keeps using the old version of a function until you redefine it.

## One change to the writing step

The uploaded file also shows Windows line endings (`\r\n`). R and most tools handle these, but some command-line tools in WSL, possibly including the aligner we'll use next, can have trouble with them. Writing the file in binary mode (`"wb"`) keeps the plain `\n` line endings:

```r
# --- Write into one FASTA file with Unix line endings ---
fasta_text <- paste0(">", names(seqs), "\n", seqs)
out_path <- file.path(data_dir, "sashimi_sequences.fasta")
con <- file(out_path, open = "wb")     # "wb" = write binary, so R doesn't convert \n into \r\n
writeLines(fasta_text, con, sep = "\n")
close(con)                              # always close a connection you opened
```

## What correct output looks like

```
nchar(hamachi_mito)   # 16537        (one number)
nchar(hamachi_cox1)   # 1551
length(seqs)          # 6
nchar(seqs)
# Query_sashimi 481 | Tilapia 495 | Hirame 1551 | Sake 630 | Hamachi 1551 | Maguro 655
```

**About 1551 vs. the 1550 in your comment:** this is the same `+ 1` counting rule as with the reverse primer. `substr()` includes *both* endpoints, so positions 5,508 through 7,058 are 7,058 − 5,508 + 1 = **1,551** bases. That matches the complete flounder COI gene exactly, which is good evidence your coordinates are right. They also contain the region where your query matched (~5,620–6,100).

Update the comment to 1551, re-run, and send me the new `nchar(seqs)` output. If it matches the numbers above, we're ready to align.
- **Verification:** I ran the code this time with the adjustments and it printed the correct output and correct FASTA `sashimi_sequences_v2.fasta`. 
```r
> # --- Extract COX1 from the Hamachi mitochondrial genome ---
> hamachi_mito <- read_one_seq(file.path(data_dir, "s_quinqueradiata.fasta"))
> nchar(hamachi_mito) #16537
[1] 16537
> 
> cox1_start <- 5508 # Start coordinate of COX1 from the GenBank FEATURES section
> cox1_end <- 7058 # End coordinate of COX1
> 
> hamachi_cox1 <- substr(hamachi_mito, cox1_start, cox1_end) # Keep only the COX1 gene
> nchar(hamachi_cox1) # 1550
[1] 1551
> 
> # --- Collect all sequences, each with short label ---
> # Short labels for ape - uses the header as the row/column name in the distance matrix
> # Full GenBank header would make the matrix unreadable
> seqs <- c(
+   Query_sashimi = query_seq,
+   Tilapia_Oniloticus_MN756471 = read_one_seq(file.path(data_dir, "o_niloticus.fasta")),
+   Hirame_Polivaceus_MH032483  = read_one_seq(file.path(data_dir, "p_olivaceus.fasta")),
+   Sake_Ssalar_MZ407797        = read_one_seq(file.path(data_dir, "s_salar.fasta")),
+   Hamachi_Squinq_NC016868     = hamachi_cox1,
+   Maguro_Talbacares_PQ812466  = read_one_seq(file.path(data_dir, "t_albacares.fasta"))
+ )
> nchar(seqs) # Length of each sequence
              Query_sashimi Tilapia_Oniloticus_MN756471  Hirame_Polivaceus_MH032483        Sake_Ssalar_MZ407797 
                        481                         495                        1551                         630 
    Hamachi_Squinq_NC016868  Maguro_Talbacares_PQ812466 
                       1551                         655 
```

## Entry 3
- **Tool:** Claude (Sonnet 5)
- **What I was doing:** I was at the aligning step after combining the FASTAs.
- **What I asked (verbatim):** "How can I run MAFFT from my terminal to align `sashimi_sequences_v2.fasta`? (+ screenshot of my file paths and the `sashimi_sequences_v2.fasta`)
- **Exact output (verbatim):**
# Align the six sequences; ">" saves MAFFT's output to a file instead of printing it on screen
mafft --auto data/sashimi_sequences_v2.fasta > results/sashimi_aligned.fasta

# Confirm the aligned file still contains 6 sequences
grep -c ">" results/sashimi_aligned.fasta

What --auto does: MAFFT picks its alignment method based on how much data you give it. With only six sequences, it should choose its most accurate method, called L-INS-i. While it runs, MAFFT prints a line beginning with Strategy:. Copy that line into your notes.

What to expect

I ran the same alignment on your file so we have a reference to compare against. Your result should look like this:

Alignment length: 1,551 columns (the length of the complete COI gene).
Your query occupies columns 112–592, which is exactly 481 columns.
Every reference sequence covers that entire window with zero gaps inside it.

Why there are no internal gaps: COI encodes a protein, and it's highly conserved. An insertion or deletion of 1 or 2 bases would shift the reading frame and break the protein, so they're rare. The gaps sit only at the ends, where the shorter barcode fragments don't reach as far as the complete genes. The main thing MAFFT is doing here is lining up the starting points.

MAFFT writes the sequences in lowercase. That's normal, and R's ape package handles it.

- **Verification:** I ran the command it told me to and MAFFT v7.526 (2024/Apr/26) successfully ran on my FASTA with the L-INS-i strategy. I checked it with the grep function and 6 sequences came back as expected.

## Entry 4 
- **Tool:** Claude (Sonnet 5)
- **What I was doing:** With the aligning finished, I needed help creating a script using `ape` and `pegas` to calculate a pairwise distance matrix to estimate pairwise genetic diversity.
- **What I asked (verbatim):** "Can you help me write a script that estimates pairwise genetic diversity in R using `ape` and `pegas` to create a pairwise distance matrix?"
- **Exact output (verbatim):**
## Phase 5: Pairwise genetic diversity in R

We'll calculate the distances **two ways** and compare them.

### Step 1. Load the alignment

Add this to your R script. Run `install.packages(c("ape", "pegas"))` once first if you don't have those packages.

```r
# ===== Phase 5: Pairwise genetic diversity =====
library(ape)     # read.dna(), dist.dna()
library(pegas)   # nuc.div()

results_dir <- "C:/Users/nelin/Desktop/GWU Information/GWU Fall 2026/Princip_Bioinformatics/LAB1/results"

# Record versions for your Methods section
R.version.string
packageVersion("ape")
packageVersion("pegas")

# Read the MAFFT alignment; because every aligned sequence has the same length,
# ape stores it as a matrix: rows = sequences, columns = alignment positions
aln <- read.dna(file.path(results_dir, "sashimi_aligned.fasta"), format = "fasta")
dim(aln)       # should print 6 1551
labels(aln)    # the six sequence names
```

### Step 2. Analysis A: the full alignment

```r
# p-distance: the proportion of compared sites that differ between two sequences
# pairwise.deletion = TRUE: for each pair, skip any column where either sequence has a gap
D_full <- dist.dna(aln, model = "raw", pairwise.deletion = TRUE)
round(as.matrix(D_full), 4)

# The same comparison as a RAW COUNT of differing sites (model = "N")
N_full <- dist.dna(aln, model = "N", pairwise.deletion = TRUE)
as.matrix(N_full)

# Helper function: how many sites were actually compared for each pair?
sites_compared <- function(dna) {
  m <- as.character(dna)          # convert to a matrix of letters; gaps appear as "-"
  has_base <- m != "-"            # TRUE wherever a sequence has a nucleotide
  n <- nrow(m)
  out <- matrix(NA, n, n, dimnames = list(rownames(m), rownames(m)))
  for (i in 1:n) {
    for (j in 1:n) {
      out[i, j] <- sum(has_base[i, ] & has_base[j, ])   # sites where BOTH have a base
    }
  }
  out
}
sites_compared(aln)
```

### Step 3. Analysis B: trim the alignment to the query's region

```r
# Find the alignment columns where the query has a nucleotide (not a gap)
query_row <- as.character(aln["Query_sashimi", ])
q_cols <- which(query_row != "-")
range(q_cols)                  # first and last query columns (expect 112 592)

# Keep only those columns for all six sequences
aln_trim <- aln[, min(q_cols):max(q_cols)]
dim(aln_trim)                  # should print 6 481

D_trim <- dist.dna(aln_trim, model = "raw", pairwise.deletion = TRUE)
round(as.matrix(D_trim), 4)

N_trim <- dist.dna(aln_trim, model = "N", pairwise.deletion = TRUE)
as.matrix(N_trim)

sites_compared(aln_trim)       # every pair should now be 481
```

### Step 4. Nucleotide diversity (π) and saving the results

```r
# pi = average proportion of differing sites across all pairs of sequences
pi_trim <- nuc.div(aln_trim)
pi_trim

# Save the matrices as CSV files to use in your Results tables
write.csv(round(as.matrix(D_trim), 4), file.path(results_dir, "pdist_trimmed.csv"))
write.csv(round(as.matrix(D_full), 4), file.path(results_dir, "pdist_full.csv"))
write.csv(as.matrix(N_full),           file.path(results_dir, "count_full.csv"))
```

### Why two analyses

In the **full** alignment, each pair is compared over whatever length the two sequences share. The query is compared over 481 sites, but flounder vs. Hamachi, both complete genes, is compared over all 1,551. In the **trimmed** alignment, every pair is compared over the **same 481 sites**, so the numbers across the matrix measure the same thing.

**When you have your output, look for these:**

1. In the **Query row**, which menu species has the smallest p-distance? How does it compare with the distance to tilapia?
2. Compare the **raw counts** (`N_full`) for Hirame vs. Hamachi against Query vs. Hirame, then compare their **p-distances**. Does the raw count point the same way as the p-distance?
3. Which set of numbers would you report as your main result, and why?
- **Verification:** I typed out the code and manually/visually inspected the CSV files that were output and the console terminal output. Output files: `count_full.csv`, `count_trimmed.csv`, `pdist_full.csv`, and `pdist_trimmed.csv`
```r
> # ----- Pairwise genetic diversity -----
> # Step 1: Load alignment
> library(ape) # read.dna(), dist.dna()
> library(pegas) # nuc.div()

Attaching package: ‘pegas’

The following object is masked from ‘package:ape’:

    mst
> 
> results_dir <- "C:/Users/nelin/Desktop/GWU Information/GWU Fall 2026/Princip_Bioinformatics/LAB1/results"
> 
> # Versions
> R.version.string
[1] "R version 4.6.1 (2026-06-24 ucrt)"
> packageVersion("ape")
[1] ‘5.8.1’
> packageVersion("pegas")
[1] ‘1.4’
> 
> # Read the MAFFT alignment; because every aligned sequence has the same length,
> # ape stores it as a matrix: row = sequences, columns = alignment positions
> aln <- read.dna(file.path(results_dir, "sashimi_aligned.fasta"), format = "fasta")
> dim(aln) # should print 6 1551
[1]    6 1551
> labels(aln) # six sequence names 
[1] "Query_sashimi"               "Tilapia_Oniloticus_MN756471"
[3] "Hirame_Polivaceus_MH032483"  "Sake_Ssalar_MZ407797"       
[5] "Hamachi_Squinq_NC016868"     "Maguro_Talbacares_PQ812466" 
> 
> # Step 2: Analysis A - Full alignment
> # p-distance: the proportion of compared sites that differ between two sequences
> # pairwise.deletion = TRUE: for each pair, skip any column where either sequence has a gap
> D_full <- dist.dna(aln, model = "raw", pairwise.deletion = TRUE)
> round(as.matrix(D_full), 4)
                            Query_sashimi
Query_sashimi                      0.0000
Tilapia_Oniloticus_MN756471        0.0000
Hirame_Polivaceus_MH032483         0.2204
Sake_Ssalar_MZ407797               0.1933
Hamachi_Squinq_NC016868            0.2204
Maguro_Talbacares_PQ812466         0.1850
                            Tilapia_Oniloticus_MN756471
Query_sashimi                                    0.0000
Tilapia_Oniloticus_MN756471                      0.0000
Hirame_Polivaceus_MH032483                       0.2323
Sake_Ssalar_MZ407797                             0.2020
Hamachi_Squinq_NC016868                          0.2323
Maguro_Talbacares_PQ812466                       0.1899
                            Hirame_Polivaceus_MH032483
Query_sashimi                                   0.2204
Tilapia_Oniloticus_MN756471                     0.2323
Hirame_Polivaceus_MH032483                      0.0000
Sake_Ssalar_MZ407797                            0.2143
Hamachi_Squinq_NC016868                         0.1838
Maguro_Talbacares_PQ812466                      0.2107
                            Sake_Ssalar_MZ407797
Query_sashimi                             0.1933
Tilapia_Oniloticus_MN756471               0.2020
Hirame_Polivaceus_MH032483                0.2143
Sake_Ssalar_MZ407797                      0.0000
Hamachi_Squinq_NC016868                   0.2079
Maguro_Talbacares_PQ812466                0.2131
                            Hamachi_Squinq_NC016868
Query_sashimi                                0.2204
Tilapia_Oniloticus_MN756471                  0.2323
Hirame_Polivaceus_MH032483                   0.1838
Sake_Ssalar_MZ407797                         0.2079
Hamachi_Squinq_NC016868                      0.0000
Maguro_Talbacares_PQ812466                   0.2015
                            Maguro_Talbacares_PQ812466
Query_sashimi                                   0.1850
Tilapia_Oniloticus_MN756471                     0.1899
Hirame_Polivaceus_MH032483                      0.2107
Sake_Ssalar_MZ407797                            0.2131
Hamachi_Squinq_NC016868                         0.2015
Maguro_Talbacares_PQ812466                      0.0000
> 
> # The same comparison as a RAW COUNT of differing sites (model = "N")
> N_full <- dist.dna(aln, model = "N", pairwise.deletion = TRUE)
> as.matrix(N_full)
                            Query_sashimi
Query_sashimi                           0
Tilapia_Oniloticus_MN756471             0
Hirame_Polivaceus_MH032483            106
Sake_Ssalar_MZ407797                   93
Hamachi_Squinq_NC016868               106
Maguro_Talbacares_PQ812466             89
                            Tilapia_Oniloticus_MN756471
Query_sashimi                                         0
Tilapia_Oniloticus_MN756471                           0
Hirame_Polivaceus_MH032483                          115
Sake_Ssalar_MZ407797                                100
Hamachi_Squinq_NC016868                             115
Maguro_Talbacares_PQ812466                           94
                            Hirame_Polivaceus_MH032483
Query_sashimi                                      106
Tilapia_Oniloticus_MN756471                        115
Hirame_Polivaceus_MH032483                           0
Sake_Ssalar_MZ407797                               135
Hamachi_Squinq_NC016868                            285
Maguro_Talbacares_PQ812466                         138
                            Sake_Ssalar_MZ407797
Query_sashimi                                 93
Tilapia_Oniloticus_MN756471                  100
Hirame_Polivaceus_MH032483                   135
Sake_Ssalar_MZ407797                           0
Hamachi_Squinq_NC016868                      131
Maguro_Talbacares_PQ812466                   133
                            Hamachi_Squinq_NC016868
Query_sashimi                                   106
Tilapia_Oniloticus_MN756471                     115
Hirame_Polivaceus_MH032483                      285
Sake_Ssalar_MZ407797                            131
Hamachi_Squinq_NC016868                           0
Maguro_Talbacares_PQ812466                      132
                            Maguro_Talbacares_PQ812466
Query_sashimi                                       89
Tilapia_Oniloticus_MN756471                         94
Hirame_Polivaceus_MH032483                         138
Sake_Ssalar_MZ407797                               133
Hamachi_Squinq_NC016868                            132
Maguro_Talbacares_PQ812466                           0
> 
> # Helper function: how many sites were actually compared for each pair?
> sites_compared <- function(dna) {
+   m <- as.character(dna) # convert to a matrix of letters; gaps appear as "-"
+   has_base <- m != "-"
+   n <- nrow(m)
+   out <- matrix(NA, n, n, dimnames = list(rownames(m), rownames(m)))
+   for (i in 1:n) {
+     for (j in 1:n) {
+       out[i,j] <- sum(has_base[i, ] & has_base[j, ]) # sites where BOTH have a base
+     }
+   }
+   out
+ }
> sites_compared(aln)
                            Query_sashimi
Query_sashimi                         481
Tilapia_Oniloticus_MN756471           481
Hirame_Polivaceus_MH032483            481
Sake_Ssalar_MZ407797                  481
Hamachi_Squinq_NC016868               481
Maguro_Talbacares_PQ812466            481
                            Tilapia_Oniloticus_MN756471
Query_sashimi                                       481
Tilapia_Oniloticus_MN756471                         495
Hirame_Polivaceus_MH032483                          495
Sake_Ssalar_MZ407797                                495
Hamachi_Squinq_NC016868                             495
Maguro_Talbacares_PQ812466                          495
                            Hirame_Polivaceus_MH032483
Query_sashimi                                      481
Tilapia_Oniloticus_MN756471                        495
Hirame_Polivaceus_MH032483                        1551
Sake_Ssalar_MZ407797                               630
Hamachi_Squinq_NC016868                           1551
Maguro_Talbacares_PQ812466                         655
                            Sake_Ssalar_MZ407797
Query_sashimi                                481
Tilapia_Oniloticus_MN756471                  495
Hirame_Polivaceus_MH032483                   630
Sake_Ssalar_MZ407797                         630
Hamachi_Squinq_NC016868                      630
Maguro_Talbacares_PQ812466                   624
                            Hamachi_Squinq_NC016868
Query_sashimi                                   481
Tilapia_Oniloticus_MN756471                     495
Hirame_Polivaceus_MH032483                     1551
Sake_Ssalar_MZ407797                            630
Hamachi_Squinq_NC016868                        1551
Maguro_Talbacares_PQ812466                      655
                            Maguro_Talbacares_PQ812466
Query_sashimi                                      481
Tilapia_Oniloticus_MN756471                        495
Hirame_Polivaceus_MH032483                         655
Sake_Ssalar_MZ407797                               624
Hamachi_Squinq_NC016868                            655
Maguro_Talbacares_PQ812466                         655
>       
> # Step 3: Analysis B - Trim alignment to query region
> query_row <- as.character(aln["Query_sashimi", ])
> q_cols <- which(query_row != "-")
> range(q_cols) # first and last query columns (expect 112 592)
[1] 112 592
> 
> # Keep only those columns for all 6 sequences
> aln_trim <- aln[, min(q_cols):max(q_cols)]
> dim(aln_trim) # should print 6 481
[1]   6 481
> 
> D_trim <- dist.dna(aln_trim, model = "raw", pairwise.deletion = TRUE)
> round(as.matrix(D_trim), 4)
                            Query_sashimi
Query_sashimi                      0.0000
Tilapia_Oniloticus_MN756471        0.0000
Hirame_Polivaceus_MH032483         0.2204
Sake_Ssalar_MZ407797               0.1933
Hamachi_Squinq_NC016868            0.2204
Maguro_Talbacares_PQ812466         0.1850
                            Tilapia_Oniloticus_MN756471
Query_sashimi                                    0.0000
Tilapia_Oniloticus_MN756471                      0.0000
Hirame_Polivaceus_MH032483                       0.2204
Sake_Ssalar_MZ407797                             0.1933
Hamachi_Squinq_NC016868                          0.2204
Maguro_Talbacares_PQ812466                       0.1850
                            Hirame_Polivaceus_MH032483
Query_sashimi                                   0.2204
Tilapia_Oniloticus_MN756471                     0.2204
Hirame_Polivaceus_MH032483                      0.0000
Sake_Ssalar_MZ407797                            0.2266
Hamachi_Squinq_NC016868                         0.1913
Maguro_Talbacares_PQ812466                      0.2037
                            Sake_Ssalar_MZ407797
Query_sashimi                             0.1933
Tilapia_Oniloticus_MN756471               0.1933
Hirame_Polivaceus_MH032483                0.2266
Sake_Ssalar_MZ407797                      0.0000
Hamachi_Squinq_NC016868                   0.2037
Maguro_Talbacares_PQ812466                0.1996
                            Hamachi_Squinq_NC016868
Query_sashimi                                0.2204
Tilapia_Oniloticus_MN756471                  0.2204
Hirame_Polivaceus_MH032483                   0.1913
Sake_Ssalar_MZ407797                         0.2037
Hamachi_Squinq_NC016868                      0.0000
Maguro_Talbacares_PQ812466                   0.1892
                            Maguro_Talbacares_PQ812466
Query_sashimi                                   0.1850
Tilapia_Oniloticus_MN756471                     0.1850
Hirame_Polivaceus_MH032483                      0.2037
Sake_Ssalar_MZ407797                            0.1996
Hamachi_Squinq_NC016868                         0.1892
Maguro_Talbacares_PQ812466                      0.0000
> 
> N_trim <- dist.dna(aln_trim, model = "N", pairwise.deletion = TRUE)
> as.matrix(N_trim)
                            Query_sashimi
Query_sashimi                           0
Tilapia_Oniloticus_MN756471             0
Hirame_Polivaceus_MH032483            106
Sake_Ssalar_MZ407797                   93
Hamachi_Squinq_NC016868               106
Maguro_Talbacares_PQ812466             89
                            Tilapia_Oniloticus_MN756471
Query_sashimi                                         0
Tilapia_Oniloticus_MN756471                           0
Hirame_Polivaceus_MH032483                          106
Sake_Ssalar_MZ407797                                 93
Hamachi_Squinq_NC016868                             106
Maguro_Talbacares_PQ812466                           89
                            Hirame_Polivaceus_MH032483
Query_sashimi                                      106
Tilapia_Oniloticus_MN756471                        106
Hirame_Polivaceus_MH032483                           0
Sake_Ssalar_MZ407797                               109
Hamachi_Squinq_NC016868                             92
Maguro_Talbacares_PQ812466                          98
                            Sake_Ssalar_MZ407797
Query_sashimi                                 93
Tilapia_Oniloticus_MN756471                   93
Hirame_Polivaceus_MH032483                   109
Sake_Ssalar_MZ407797                           0
Hamachi_Squinq_NC016868                       98
Maguro_Talbacares_PQ812466                    96
                            Hamachi_Squinq_NC016868
Query_sashimi                                   106
Tilapia_Oniloticus_MN756471                     106
Hirame_Polivaceus_MH032483                       92
Sake_Ssalar_MZ407797                             98
Hamachi_Squinq_NC016868                           0
Maguro_Talbacares_PQ812466                       91
                            Maguro_Talbacares_PQ812466
Query_sashimi                                       89
Tilapia_Oniloticus_MN756471                         89
Hirame_Polivaceus_MH032483                          98
Sake_Ssalar_MZ407797                                96
Hamachi_Squinq_NC016868                             91
Maguro_Talbacares_PQ812466                           0
> 
> sites_compared(aln_trim) # every pair should now be 481
                            Query_sashimi
Query_sashimi                         481
Tilapia_Oniloticus_MN756471           481
Hirame_Polivaceus_MH032483            481
Sake_Ssalar_MZ407797                  481
Hamachi_Squinq_NC016868               481
Maguro_Talbacares_PQ812466            481
                            Tilapia_Oniloticus_MN756471
Query_sashimi                                       481
Tilapia_Oniloticus_MN756471                         481
Hirame_Polivaceus_MH032483                          481
Sake_Ssalar_MZ407797                                481
Hamachi_Squinq_NC016868                             481
Maguro_Talbacares_PQ812466                          481
                            Hirame_Polivaceus_MH032483
Query_sashimi                                      481
Tilapia_Oniloticus_MN756471                        481
Hirame_Polivaceus_MH032483                         481
Sake_Ssalar_MZ407797                               481
Hamachi_Squinq_NC016868                            481
Maguro_Talbacares_PQ812466                         481
                            Sake_Ssalar_MZ407797
Query_sashimi                                481
Tilapia_Oniloticus_MN756471                  481
Hirame_Polivaceus_MH032483                   481
Sake_Ssalar_MZ407797                         481
Hamachi_Squinq_NC016868                      481
Maguro_Talbacares_PQ812466                   481
                            Hamachi_Squinq_NC016868
Query_sashimi                                   481
Tilapia_Oniloticus_MN756471                     481
Hirame_Polivaceus_MH032483                      481
Sake_Ssalar_MZ407797                            481
Hamachi_Squinq_NC016868                         481
Maguro_Talbacares_PQ812466                      481
                            Maguro_Talbacares_PQ812466
Query_sashimi                                      481
Tilapia_Oniloticus_MN756471                        481
Hirame_Polivaceus_MH032483                         481
Sake_Ssalar_MZ407797                               481
Hamachi_Squinq_NC016868                            481
Maguro_Talbacares_PQ812466                         481
> 
> # Step 4: Nucleotide diversity and saving
> # pi = average proportion of differing sites across all pairs of sequences
> pi_trim <- nuc.div(aln_trim)
> pi_trim
[1] 0.1901594
```