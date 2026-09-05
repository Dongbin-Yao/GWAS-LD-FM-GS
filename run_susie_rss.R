#!/usr/bin/env Rscript

# =============================================================================
# run_susie_rss.R
# =============================================================================
#
# SuSiE-RSS fine-mapping using GWAS summary statistics and an LD correlation
# matrix.
#
# Workflow
# --------
#
# GWAS summary statistics
#          |
#          v
#      z-scores
#          +
#      LD matrix
#          |
#          v
#       SuSiE-RSS
#          |
#          +-------------------+
#          |                   |
#          v                   v
#      SNP-level PIP      Credible sets
#
# The script does not perform GWAS or LD calculation.
# It uses pre-computed GWAS summary statistics and an LD matrix.
#
# =============================================================================


# =============================================================================
# 1. Required packages
# =============================================================================

required_packages <- c(
  "susieR"
)

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    stop(
      paste0(
        "Required R package '",
        pkg,
        "' is not installed.\n",
        "Install it with:\n",
        "install.packages('",
        pkg,
        "')"
      )
    )
  }
}

library(susieR)


# =============================================================================
# 2. Command-line arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)


get_arg <- function(flag, default = NULL) {

  idx <- which(args == flag)

  if (length(idx) == 0) {
    return(default)
  }

  if (idx == length(args)) {
    stop(paste0("Missing value for argument: ", flag))
  }

  args[idx + 1]
}


gwas_file <- get_arg("--gwas")
ld_file <- get_arg("--ld")
output_prefix <- get_arg("--output")

L <- as.integer(
  get_arg("--L", default = "10")
)

coverage <- as.numeric(
  get_arg("--coverage", default = "0.95")
)

min_abs_corr <- as.numeric(
  get_arg("--min_abs_corr", default = "0")
)


# =============================================================================
# 3. Check arguments
# =============================================================================

if (is.null(gwas_file)) {
  stop("--gwas is required.")
}

if (is.null(ld_file)) {
  stop("--ld is required.")
}

if (is.null(output_prefix)) {
  stop("--output is required.")
}

if (!file.exists(gwas_file)) {
  stop("GWAS file does not exist: ", gwas_file)
}

if (!file.exists(ld_file)) {
  stop("LD matrix file does not exist: ", ld_file)
}


# =============================================================================
# 4. Read GWAS summary statistics
# =============================================================================

message("[1/6] Reading GWAS summary statistics...")

gwas <- read.table(
  gwas_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# -----------------------------------------------------------------------------
# Required columns
# -----------------------------------------------------------------------------

required_columns <- c(
  "SNP",
  "P"
)

missing_columns <- setdiff(
  required_columns,
  colnames(gwas)
)

if (length(missing_columns) > 0) {

  stop(
    paste0(
      "GWAS file is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# =============================================================================
# 5. Calculate z-scores
# =============================================================================

message("[2/6] Calculating GWAS z-scores...")


# -----------------------------------------------------------------------------
# If the GWAS file already contains Z, use it.
#
# Otherwise calculate:
#
#     z = qnorm(1 - P/2) * sign(BETA)
#
# If BETA is not available, the z-score is assumed to be positive.
# -----------------------------------------------------------------------------

if ("Z" %in% colnames(gwas)) {

  gwas$Z <- as.numeric(gwas$Z)

} else {

  if ("BETA" %in% colnames(gwas)) {

    gwas$BETA <- as.numeric(gwas$BETA)

    gwas$Z <- qnorm(
      1 - as.numeric(gwas$P) / 2
    ) * sign(gwas$BETA)

  } else {

    warning(
      paste0(
        "BETA/Z column not found. ",
        "Z-scores will be reconstructed from P-values ",
        "without effect-direction information."
      )
    )

    gwas$Z <- qnorm(
      1 - as.numeric(gwas$P) / 2
    )
  }
}


# -----------------------------------------------------------------------------
# Remove invalid records
# -----------------------------------------------------------------------------

gwas <- gwas[
  is.finite(gwas$Z) &
  is.finite(gwas$P) &
  gwas$P > 0 &
  gwas$P <= 1,
]


if (nrow(gwas) == 0) {
  stop("No valid GWAS records remain.")
}


# =============================================================================
# 6. Read LD matrix
# =============================================================================

message("[3/6] Reading LD correlation matrix...")


# -----------------------------------------------------------------------------
# Expected format:
#
#       SNP1 SNP2 SNP3 ...
# SNP1  1    r12  r13
# SNP2  r21  1    r23
# SNP3  r31  r32  1
#
# The row names and column names must correspond to SNP IDs.
# -----------------------------------------------------------------------------

R <- as.matrix(
  read.table(
    ld_file,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )
)

storage.mode(R) <- "numeric"


if (nrow(R) != ncol(R)) {
  stop("LD matrix must be square.")
}


if (is.null(rownames(R)) || is.null(colnames(R))) {
  stop("LD matrix must contain SNP IDs as row and column names.")
}


# =============================================================================
# 7. Match GWAS SNPs and LD matrix
# =============================================================================

message("[4/6] Matching GWAS SNPs with LD matrix...")


common_snps <- intersect(
  gwas$SNP,
  rownames(R)
)

if (length(common_snps) < 2) {

  stop(
    paste0(
      "Fewer than two SNPs were shared between GWAS and LD matrix."
    )
  )
}


# Preserve GWAS order
common_snps <- gwas$SNP[
  gwas$SNP %in% common_snps
]

common_snps <- unique(common_snps)


gwas_sub <- gwas[
  match(common_snps, gwas$SNP),
]


R_sub <- R[
  common_snps,
  common_snps,
  drop = FALSE
]


z <- gwas_sub$Z


# =============================================================================
# 8. Basic LD matrix checks
# =============================================================================

message("[5/6] Checking LD matrix...")


# -----------------------------------------------------------------------------
# Symmetrize numerical rounding differences
# -----------------------------------------------------------------------------

R_sub <- (R_sub + t(R_sub)) / 2


# -----------------------------------------------------------------------------
# Set diagonal to one
# -----------------------------------------------------------------------------

diag(R_sub) <- 1


# -----------------------------------------------------------------------------
# Optional correlation threshold
# -----------------------------------------------------------------------------

if (min_abs_corr > 0) {

  R_sub[
    abs(R_sub) < min_abs_corr
  ] <- 0

  diag(R_sub) <- 1
}


# -----------------------------------------------------------------------------
# Check positive semi-definiteness
# -----------------------------------------------------------------------------

eigen_values <- eigen(
  R_sub,
  symmetric = TRUE,
  only.values = TRUE
)$values


if (min(eigen_values) < -1e-6) {

  warning(
    paste0(
      "LD matrix is not positive semi-definite. ",
      "Minimum eigenvalue = ",
      signif(min(eigen_values), 4)
    )
  )
}


# =============================================================================
# 9. Run SuSiE-RSS
# =============================================================================

message("[6/6] Running SuSiE-RSS...")

message(
  paste0(
    "      Number of SNPs: ",
    length(z)
  )
)

message(
  paste0(
    "      Maximum number of effects (L): ",
    L
  )
)

message(
  paste0(
    "      Credible set coverage: ",
    coverage
  )
)


fit <- susie_rss(
  z = z,
  R = R_sub,
  L = L,
  coverage = coverage,
  check_prior = TRUE,
  estimate_residual_variance = FALSE
)


# =============================================================================
# 10. SNP-level PIP results
# =============================================================================

pip_results <- data.frame(
  SNP = common_snps,
  PIP = fit$pip
)


# Add GWAS statistics
pip_results <- merge(
  pip_results,
  gwas_sub,
  by = "SNP",
  sort = FALSE
)


# Sort by PIP
pip_results <- pip_results[
  order(
    -pip_results$PIP,
    pip_results$P
  ),
]


# Add PIP rank
pip_results$PIP_Rank <- seq_len(
  nrow(pip_results)
)


# Reorder columns
preferred_columns <- c(
  "PIP_Rank",
  "SNP",
  "PIP",
  "P",
  "Z"
)

remaining_columns <- setdiff(
  colnames(pip_results),
  preferred_columns
)

pip_results <- pip_results[
  c(preferred_columns, remaining_columns)
]


# =============================================================================
# 11. Save PIP results
# =============================================================================

pip_file <- paste0(
  output_prefix,
  "_PIP.txt"
)

write.table(
  pip_results,
  file = pip_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# =============================================================================
# 12. Extract credible sets
# =============================================================================

credible_sets <- data.frame()


if (!is.null(fit$sets$cs)) {

  cs_list <- fit$sets$cs

  for (cs_id in names(cs_list)) {

    snp_indices <- cs_list[[cs_id]]

    if (length(snp_indices) == 0) {
      next
    }

    cs_data <- data.frame(
      Credible_Set = cs_id,
      SNP = common_snps[snp_indices],
      PIP = fit$pip[snp_indices],
      stringsAsFactors = FALSE
    )

    credible_sets <- rbind(
      credible_sets,
      cs_data
    )
  }
}


# =============================================================================
# 13. Save credible sets
# =============================================================================

cs_file <- paste0(
  output_prefix,
  "_credible_sets.txt"
)

if (nrow(credible_sets) > 0) {

  write.table(
    credible_sets,
    file = cs_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

} else {

  write.table(
    data.frame(
      Credible_Set = character(),
      SNP = character(),
      PIP = numeric()
    ),
    file = cs_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


# =============================================================================
# 14. Save model summary
# =============================================================================

summary_file <- paste0(
  output_prefix,
  "_summary.txt"
)

sink(summary_file)

cat("SuSiE-RSS Fine-mapping Summary\n")
cat("==============================\n\n")

cat("Input GWAS file:\n")
cat(gwas_file, "\n\n")

cat("Input LD matrix:\n")
cat(ld_file, "\n\n")

cat("Number of SNPs:\n")
cat(length(z), "\n\n")

cat("Maximum number of effects (L):\n")
cat(L, "\n\n")

cat("Credible set coverage:\n")
cat(coverage, "\n\n")

cat("Number of credible sets:\n")

if (!is.null(fit$sets$cs)) {
  cat(length(fit$sets$cs), "\n\n")
} else {
  cat(0, "\n\n")
}

cat("Maximum PIP:\n")
cat(max(fit$pip), "\n\n")

cat("SuSiE model summary:\n\n")

print(fit)

sink()


# =============================================================================
# 15. Final message
# =============================================================================

message("")
message("SuSiE-RSS fine-mapping completed.")
message("")
message(
  paste0(
    "PIP results: ",
    pip_file
  )
)

message(
  paste0(
    "Credible sets: ",
    cs_file
  )
)

message(
  paste0(
    "Model summary: ",
    summary_file
  )
)
