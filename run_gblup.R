#!/usr/bin/env Rscript

# =============================================================================
# run_gblup.R
# =============================================================================
#
# PIP-prioritized GBLUP
#
# Workflow
# --------
#
# SuSiE-RSS PIP
#      |
#      v
# Select top-N SNPs
#      |
#      v
# Genotype matrix
#      |
#      v
# Genomic relationship matrix
#      |
#      v
# GBLUP
#      |
#      v
# 10-fold cross-validation
#
# The GBLUP model is implemented using rrBLUP::mixed.solve().
#
# =============================================================================


# =============================================================================
# 1. Required packages
# =============================================================================

required_packages <- c(
  "rrBLUP"
)

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    stop(
      paste0(
        "Package '", pkg,
        "' is required. Install it with:\n",
        "install.packages('", pkg, "')"
      )
    )
  }
}

library(rrBLUP)


# =============================================================================
# 2. Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)


get_arg <- function(flag, default = NULL) {

  idx <- which(args == flag)

  if (length(idx) == 0) {
    return(default)
  }

  if (idx == length(args)) {
    stop("Missing value for argument: ", flag)
  }

  args[idx + 1]
}


genotype_file <- get_arg("--genotype")
phenotype_file <- get_arg("--phenotype")
pip_file <- get_arg("--pip")
output_prefix <- get_arg("--output")

top_n <- as.integer(
  get_arg("--top_n", default = "0")
)

nfold <- as.integer(
  get_arg("--folds", default = "10")
)

seed <- as.integer(
  get_arg("--seed", default = "12345")
)


if (is.null(genotype_file)) stop("--genotype is required.")
if (is.null(phenotype_file)) stop("--phenotype is required.")
if (is.null(pip_file)) stop("--pip is required.")
if (is.null(output_prefix)) stop("--output is required.")


# =============================================================================
# 3. Read genotype
# =============================================================================

message("[1/7] Reading genotype matrix...")

geno <- read.table(
  genotype_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"ID" %in% colnames(geno)) {

  stop(
    "Genotype file must contain an 'ID' column."
  )
}

sample_ids <- geno$ID

geno$ID <- NULL

geno <- as.matrix(geno)

storage.mode(geno) <- "numeric"

rownames(geno) <- sample_ids


# =============================================================================
# 4. Read phenotype
# =============================================================================

message("[2/7] Reading phenotype...")

pheno <- read.table(
  phenotype_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("ID", "Phenotype") %in% colnames(pheno))) {

  stop(
    "Phenotype file must contain columns: ID and Phenotype"
  )
}

pheno$Phenotype <- as.numeric(
  pheno$Phenotype
)


# =============================================================================
# 5. Read PIP results
# =============================================================================

message("[3/7] Reading PIP results...")

pip <- read.table(
  pip_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("SNP", "PIP") %in% colnames(pip))) {

  stop(
    "PIP file must contain columns: SNP and PIP."
  )
}

pip <- pip[
  order(-pip$PIP),
]


# =============================================================================
# 6. Select SNPs
# =============================================================================

if (top_n <= 0) {

  selected_snps <- intersect(
    pip$SNP,
    colnames(geno)
  )

} else {

  selected_snps <- pip$SNP[
    seq_len(
      min(top_n, nrow(pip))
    )
  ]

  selected_snps <- intersect(
    selected_snps,
    colnames(geno)
  )
}


if (length(selected_snps) < 2) {

  stop(
    "Fewer than two selected SNPs were found in the genotype matrix."
  )
}


geno <- geno[
  ,
  selected_snps,
  drop = FALSE
]


message(
  "      Selected SNPs: ",
  ncol(geno)
)


# =============================================================================
# 7. Match individuals
# =============================================================================

message("[4/7] Matching individuals...")

common_ids <- intersect(
  pheno$ID,
  rownames(geno)
)

if (length(common_ids) < nfold) {

  stop(
    "Number of individuals is smaller than the requested number of folds."
  )
}


pheno <- pheno[
  match(common_ids, pheno$ID),
]

geno <- geno[
  common_ids,
  ,
  drop = FALSE
]


y <- pheno$Phenotype


# =============================================================================
# 8. 10-fold cross-validation
# =============================================================================

message("[5/7] Running cross-validation...")

set.seed(seed)

fold_id <- sample(
  rep(
    1:nfold,
    length.out = length(y)
  )
)


oof_prediction <- rep(
  NA_real_,
  length(y)
)


fold_results <- data.frame()


for (fold in 1:nfold) {

  message(
    "      Fold ",
    fold,
    "/",
    nfold
  )

  test_idx <- which(
    fold_id == fold
  )

  train_idx <- setdiff(
    seq_along(y),
    test_idx
  )


  X_train <- geno[
    train_idx,
    ,
    drop = FALSE
  ]

  X_test <- geno[
    test_idx,
    ,
    drop = FALSE
  ]

  y_train <- y[
    train_idx
  ]

  y_test <- y[
    test_idx
  ]


  # ---------------------------------------------------------------------------
  # Mean-impute genotype using training data only
  # ---------------------------------------------------------------------------

  marker_means <- colMeans(
    X_train,
    na.rm = TRUE
  )

  for (j in seq_len(ncol(X_train))) {

    missing_train <- is.na(
      X_train[, j]
    )

    missing_test <- is.na(
      X_test[, j]
    )

    if (any(missing_train)) {

      X_train[missing_train, j] <-
        marker_means[j]
    }

    if (any(missing_test)) {

      X_test[missing_test, j] <-
        marker_means[j]
    }
  }


  # ---------------------------------------------------------------------------
  # Construct genomic relationship matrix
  # ---------------------------------------------------------------------------

  p <- colMeans(
    X_train,
    na.rm = TRUE
  ) / 2

  Z_train <- sweep(
    X_train,
    2,
    2 * p,
    FUN = "-"
  )

  denominator <- sum(
    2 * p * (1 - p)
  )

  if (denominator <= 0) {

    stop(
      "Invalid GRM denominator."
    )
  }

  G_train <- tcrossprod(
    Z_train
  ) / denominator


  # ---------------------------------------------------------------------------
  # GBLUP
  # ---------------------------------------------------------------------------

  fit <- rrBLUP::mixed.solve(
    y = y_train,
    K = G_train
  )


  # ---------------------------------------------------------------------------
  # Prediction
  #
  # mixed.solve returns marker/random-effect solutions when Z is supplied.
  # Here we use the GBLUP equivalent through the genomic relationship matrix.
  # ---------------------------------------------------------------------------

  G_test_train <- {

    Z_test <- sweep(
      X_test,
      2,
      2 * p,
      FUN = "-"
    )

    tcrossprod(
      Z_test,
      Z_train
    ) / denominator
  }


  # Solve genomic prediction:
  #
  # u_hat = G (G + lambda I)^(-1) (y - mu)
  #

  lambda <- fit$Ve / fit$Vu

  centered_y <- y_train - fit$beta

  u_train <- solve(
    G_train + lambda * diag(nrow(G_train)),
    centered_y
  )

  pred <- as.numeric(
    fit$beta +
      G_test_train %*% u_train
  )


  oof_prediction[
    test_idx
  ] <- pred


  # ---------------------------------------------------------------------------
  # Fold metrics
  # ---------------------------------------------------------------------------

  pcc <- cor(
    y_test,
    pred,
    use = "complete.obs"
  )

  mse <- mean(
    (y_test - pred)^2
  )

  rmse <- sqrt(
    mse
  )


  fold_results <- rbind(
    fold_results,
    data.frame(
      Fold = fold,
      N_Test = length(test_idx),
      PCC = pcc,
      MSE = mse,
      RMSE = rmse
    )
  )
}


# =============================================================================
# 9. Pooled OOF metrics
# =============================================================================

message("[6/7] Calculating pooled OOF metrics...")

pooled_pcc <- cor(
  y,
  oof_prediction,
  use = "complete.obs"
)

pooled_mse <- mean(
  (y - oof_prediction)^2
)

pooled_rmse <- sqrt(
  pooled_mse
)


summary_results <- data.frame(
  Model = "GBLUP",
  Top_N_SNP = top_n,
  Selected_SNPs = ncol(geno),
  N = length(y),
  Fold_Mean_PCC = mean(fold_results$PCC),
  Fold_SD_PCC = sd(fold_results$PCC),
  Fold_Mean_MSE = mean(fold_results$MSE),
  Fold_SD_MSE = sd(fold_results$MSE),
  Fold_Mean_RMSE = mean(fold_results$RMSE),
  Fold_SD_RMSE = sd(fold_results$RMSE),
  Pooled_OOF_PCC = pooled_pcc,
  Pooled_OOF_MSE = pooled_mse,
  Pooled_OOF_RMSE = pooled_rmse
)


# =============================================================================
# 10. Save results
# =============================================================================

message("[7/7] Saving results...")

dir.create(
  dirname(output_prefix),
  recursive = TRUE,
  showWarnings = FALSE
)


write.table(
  fold_results,
  paste0(output_prefix, "_fold_results.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  summary_results,
  paste0(output_prefix, "_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


oof_results <- data.frame(
  ID = common_ids,
  Observed = y,
  Predicted = oof_prediction,
  Fold = fold_id
)


write.table(
  oof_results,
  paste0(output_prefix, "_OOF_predictions.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


message("")
message("GBLUP completed.")
message(
  "Pooled OOF PCC: ",
  round(pooled_pcc, 4)
)
message(
  "Pooled OOF MSE: ",
  round(pooled_mse, 4)
)
message(
  "Selected SNPs: ",
  ncol(geno)
)
