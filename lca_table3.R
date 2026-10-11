#!/usr/bin/env Rscript
#
# Latent class analysis of COVID-19 misinformation beliefs (4-class model)
# with posterior-weighted baseline characteristics and Wald tests.
#
# Fits a 4-class LCA via EM algorithm (100 random starts) to four common
# misinformation items from the EKOS COVID-19 Monitor, Wave 27:
#   A. COVID-19 deaths exaggerated
#   B. Vaccine deaths hidden by government
#   C. Vaccines cause infertility
#   D. Vaccines alter DNA
#
# Item coding (all are false statements):
#   0 = Correct (disagree/strongly disagree)
#   1 = Don't Know
#   2 = Incorrect (agree/strongly agree)
#
# Outputs a CSV of posterior-weighted proportions, means (SD), and
# Wald test p-values (chi-square for categorical, F for continuous).
#
# Usage:
#   Rscript lca_table3.R path/to/05021W27_client_data_file.sav [--output table3.csv]
#
# Dependencies: haven, dplyr

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})


# ---------------------------------------------------------------------------
# 1. Data loading and recoding
# ---------------------------------------------------------------------------

load_and_recode <- function(sav_path) {
  w27 <- read_sav(sav_path)
  cat(sprintf("Wave 27 raw: %d rows, %d columns\n", nrow(w27), ncol(w27)))

  recode_false <- function(x) {
    # Recode false-statement items: 1,2->Correct; 98->DK; 3,4->Incorrect
    ifelse(x %in% c(1, 2), 0L,
    ifelse(x == 98, 1L,
    ifelse(x %in% c(3, 4), 2L, NA_integer_)))
  }

  item_cols_raw <- c("Q34A", "Q34B", "Q34C", "Q34D")
  item_cols     <- paste0(item_cols_raw, "_r")

  for (i in seq_along(item_cols_raw)) {
    w27[[item_cols[i]]] <- recode_false(as.numeric(w27[[item_cols_raw[i]]]))
  }

  complete <- complete.cases(w27[, item_cols])
  w27 <- w27[complete, ]
  cat(sprintf("Complete on 4 items: N = %d\n", nrow(w27)))

  list(data = w27, item_cols = item_cols)
}


# ---------------------------------------------------------------------------
# 2. EM algorithm for K-class LCA with 3-category items
# ---------------------------------------------------------------------------

em_lca <- function(X, K, n_starts = 100, max_iter = 500, tol = 1e-8) {
  # X: integer matrix N x J  (values 0..C-1)
  # Returns: list(pi, rho, post, ll)

  N <- nrow(X)
  J <- ncol(X)
  C <- max(X) + 1L

  best_ll   <- -Inf
  best_pi   <- NULL
  best_rho  <- NULL

  for (start in seq_len(n_starts)) {
    set.seed(start * 42 + 7)

    # Random initialisation
    alpha <- rep(1, K)
    pi_k  <- as.numeric(MCMCpack::rdirichlet(1, alpha))
    if (!requireNamespace("MCMCpack", quietly = TRUE)) {
      # Fallback: simple Dirichlet via gamma
      g <- rgamma(K, 1)
      pi_k <- g / sum(g)
    }
    # Use gamma-based Dirichlet for robustness
    g <- rgamma(K, 1)
    pi_k <- g / sum(g)

    rho <- array(0, dim = c(K, J, C))
    for (k in seq_len(K)) {
      for (j in seq_len(J)) {
        g <- rgamma(C, 2)
        rho[k, j, ] <- g / sum(g)
      }
    }

    ll_prev <- -Inf

    for (iter in seq_len(max_iter)) {
      # E-step
      log_lik <- matrix(log(pi_k), nrow = N, ncol = K, byrow = TRUE)
      for (j in seq_len(J)) {
        for (k in seq_len(K)) {
          log_lik[, k] <- log_lik[, k] + log(rho[k, j, X[, j] + 1L] + 1e-300)
        }
      }

      max_ll_row <- apply(log_lik, 1, max)
      log_denom  <- max_ll_row + log(rowSums(exp(log_lik - max_ll_row)))
      post <- exp(log_lik - log_denom)
      ll   <- sum(log_denom)

      if (abs(ll - ll_prev) < tol) break
      ll_prev <- ll

      # M-step
      Nk   <- colSums(post)
      pi_k <- Nk / N
      for (k in seq_len(K)) {
        for (j in seq_len(J)) {
          for (cc in seq_len(C)) {
            rho[k, j, cc] <- sum(post[, k] * (X[, j] == (cc - 1L))) / (Nk[k] + 1e-300)
          }
        }
      }
    }

    if (ll > best_ll) {
      best_ll  <- ll
      best_pi  <- pi_k
      best_rho <- rho
    }
  }

  # Final E-step with best parameters
  log_lik <- matrix(log(best_pi), nrow = N, ncol = K, byrow = TRUE)
  for (j in seq_len(J)) {
    for (k in seq_len(K)) {
      log_lik[, k] <- log_lik[, k] + log(best_rho[k, j, X[, j] + 1L] + 1e-300)
    }
  }
  max_ll_row <- apply(log_lik, 1, max)
  log_denom  <- max_ll_row + log(rowSums(exp(log_lik - max_ll_row)))
  best_post  <- exp(log_lik - log_denom)

  list(pi = best_pi, rho = best_rho, post = best_post, ll = best_ll)
}


# ---------------------------------------------------------------------------
# 3. Auto-label classes
# ---------------------------------------------------------------------------

label_classes <- function(rho, K) {
  labels <- character(K)
  for (k in seq_len(K)) {
    avg_correct   <- mean(rho[k, , 1])   # column 1 = category 0 (Correct)
    avg_dk        <- mean(rho[k, , 2])   # column 2 = category 1 (DK)
    avg_incorrect <- mean(rho[k, , 3])   # column 3 = category 2 (Incorrect)

    if (avg_correct > 0.7) {
      labels[k] <- "Well-Informed"
    } else if (avg_incorrect > 0.4) {
      labels[k] <- "Misinformed"
    } else if (avg_dk > 0.3) {
      labels[k] <- "Uncertain/Disengaged"
    } else {
      labels[k] <- "Selectively Skeptical"
    }
  }

  # Resolve duplicates by rank-ordering on P(Correct)
  if (length(unique(labels)) < K) {
    avg_corr <- sapply(seq_len(K), function(k) mean(rho[k, , 1]))
    ord <- order(-avg_corr)
    manual <- c("Well-Informed", "Selectively Skeptical",
                "Uncertain/Disengaged", "Misinformed")
    labels <- character(K)
    for (rank in seq_along(ord)) {
      labels[ord[rank]] <- manual[rank]
    }
  }

  labels
}


# ---------------------------------------------------------------------------
# 4. Posterior-weighted Wald tests
# ---------------------------------------------------------------------------

pw_chisq_test <- function(binary_vec, post_matrix, labels, class_order) {
  # binary_vec: numeric 0/1 with NA
  # Returns list(W, p)
  K     <- length(class_order)
  valid <- !is.na(binary_vec)
  x     <- binary_vec[valid]
  pv    <- post_matrix[valid, , drop = FALSE]

  Nk      <- colSums(pv)
  N_total <- sum(Nk)

  pk <- sapply(class_order, function(cls) {
    k <- match(cls, labels)
    sum(pv[, k] * x) / Nk[k]
  })

  p_bar <- sum(rowSums(pv) * x) / N_total

  if (p_bar <= 0 || p_bar >= 1) return(list(W = NA_real_, p = NA_real_))

  Nk_ordered <- sapply(class_order, function(cls) Nk[match(cls, labels)])
  W     <- sum(Nk_ordered * (pk - p_bar)^2) / (p_bar * (1 - p_bar))
  p_val <- 1 - pchisq(W, df = K - 1)
  list(W = W, p = p_val)
}


pw_chisq_multi_test <- function(cat_vec, categories, post_matrix, labels, class_order) {
  # cat_vec: categorical variable with NA
  # categories: vector of valid category values
  K     <- length(class_order)
  valid <- !is.na(cat_vec)
  x     <- cat_vec[valid]
  pv    <- post_matrix[valid, , drop = FALSE]

  Nk      <- colSums(pv)
  N_total <- sum(Nk)
  C_cat   <- length(categories)

  W <- 0
  for (cat_val in categories) {
    indicator <- as.numeric(x == cat_val)
    pk <- sapply(class_order, function(cls) {
      k <- match(cls, labels)
      sum(pv[, k] * indicator) / Nk[k]
    })
    p_bar <- sum(rowSums(pv) * indicator) / N_total
    if (p_bar > 0) {
      Nk_ordered <- sapply(class_order, function(cls) Nk[match(cls, labels)])
      W <- W + sum(Nk_ordered * (pk - p_bar)^2) / p_bar
    }
  }

  df_val <- (K - 1) * (C_cat - 1)
  p_val  <- 1 - pchisq(W, df = df_val)
  list(W = W, p = p_val)
}


pw_F_test <- function(cont_vec, post_matrix, labels, class_order) {
  # cont_vec: numeric with NA
  K     <- length(class_order)
  valid <- !is.na(cont_vec)
  x     <- cont_vec[valid]
  pv    <- post_matrix[valid, , drop = FALSE]

  Nk      <- colSums(pv)
  N_total <- sum(Nk)

  mu_k <- sapply(class_order, function(cls) {
    k <- match(cls, labels)
    sum(pv[, k] * x) / Nk[k]
  })
  mu_bar <- sum(rowSums(pv) * x) / N_total

  Nk_ordered <- sapply(class_order, function(cls) Nk[match(cls, labels)])
  SSB <- sum(Nk_ordered * (mu_k - mu_bar)^2)
  dfB <- K - 1

  SSW <- 0
  for (i_c in seq_along(class_order)) {
    k <- match(class_order[i_c], labels)
    SSW <- SSW + sum(pv[, k] * (x - mu_k[i_c])^2)
  }
  dfW <- N_total - K

  if (dfW <= 0 || SSW <= 0) return(list(F_stat = NA_real_, p = NA_real_))

  F_stat <- (SSB / dfB) / (SSW / dfW)
  p_val  <- 1 - pf(F_stat, dfB, dfW)
  list(F_stat = F_stat, p = p_val)
}


# ---------------------------------------------------------------------------
# 5. Posterior-weighted descriptives (helpers)
# ---------------------------------------------------------------------------

pw_mean <- function(var_vec, weights) {
  valid <- !is.na(var_vec)
  if (sum(valid) == 0) return(NA_real_)
  weighted.mean(var_vec[valid], w = weights[valid])
}

pw_prop <- function(condition, weights) {
  valid <- !is.na(condition)
  if (sum(valid) == 0) return(NA_real_)
  weighted.mean(as.numeric(condition[valid]), w = weights[valid])
}

pw_sd <- function(var_vec, weights) {
  valid <- !is.na(var_vec)
  if (sum(valid) == 0) return(NA_real_)
  v <- var_vec[valid]
  w <- weights[valid]
  mu <- weighted.mean(v, w = w)
  sqrt(weighted.mean((v - mu)^2, w = w))
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  if (p < 0.01)  return(sprintf("%.3f", p))
  sprintf("%.2f", p)
}


# ---------------------------------------------------------------------------
# 6. Compute class-specific estimates and tests
# ---------------------------------------------------------------------------

compute_table <- function(w27, pi_est, rho, post, labels, class_order, K) {

  results <- list()
  p_values <- list()
  N <- nrow(w27)

  for (cls in class_order) {
    k_idx   <- match(cls, labels)
    weights <- post[, k_idx]
    row     <- list()
    row[["Class"]] <- cls
    row[["N"]]     <- sprintf("%.0f (%.1f%%)", pi_est[k_idx] * N, pi_est[k_idx] * 100)

    # Sex
    if ("QGENDR" %in% names(w27)) {
      qg <- as.numeric(w27$QGENDR)
      male   <- ifelse(!is.na(qg), qg == 1, NA)
      female <- ifelse(!is.na(qg), qg == 2, NA)
      row[["Male (%)"]]   <- sprintf("%.1f", pw_prop(male, weights) * 100)
      row[["Female (%)"]] <- sprintf("%.1f", pw_prop(female, weights) * 100)
    }

    # Age groups
    if ("QAGEMIX" %in% names(w27)) {
      qa <- as.numeric(w27$QAGEMIX)
      for (spec in list(
        list(vals = c(1, 2), lbl = "Age 18-34 (%)"),
        list(vals = c(3, 4), lbl = "Age 35-54 (%)"),
        list(vals = 5,       lbl = "Age 55-64 (%)"),
        list(vals = 6,       lbl = "Age >= 65 (%)")
      )) {
        cond <- ifelse(!is.na(qa), qa %in% spec$vals, NA)
        row[[spec$lbl]] <- sprintf("%.1f", pw_prop(cond, weights) * 100)
      }
    }

    # Education
    if ("QEDUC" %in% names(w27)) {
      qe <- as.numeric(w27$QEDUC)
      univ <- ifelse(!is.na(qe), qe %in% c(7, 8), NA)
      row[["Bachelor's or postgraduate (%)"]] <- sprintf("%.1f", pw_prop(univ, weights) * 100)
    }

    # Rural
    if ("QRURBAN" %in% names(w27)) {
      qr <- as.numeric(w27$QRURBAN)
      rural <- ifelse(!is.na(qr), qr == 1, NA)
      row[["Rural community (%)"]] <- sprintf("%.1f", pw_prop(rural, weights) * 100)
    }

    # Social class
    if ("SOCIALC1" %in% names(w27)) {
      sc <- as.numeric(w27$SOCIALC1)
      poor_wc <- ifelse(!is.na(sc), sc %in% c(1, 2), NA)
      row[["Poor or working social class (%)"]] <- sprintf("%.1f", pw_prop(poor_wc, weights) * 100)
    }

    # Vaccination
    if ("Q15" %in% names(w27)) {
      q15 <- as.numeric(w27$Q15)
      unvax <- ifelse(!is.na(q15), q15 == 4, NA)
      vax3  <- ifelse(!is.na(q15), q15 == 3, NA)
      row[["Unvaccinated (%)"]]              <- sprintf("%.1f", pw_prop(unvax, weights) * 100)
      row[["Three or more vaccine doses (%)"]] <- sprintf("%.1f", pw_prop(vax3, weights) * 100)
    }

    # Trust (1-7 scale)
    trust_map <- list(
      list(var = "Q3A", name = "Trust in doctors"),
      list(var = "Q3E", name = "Trust in scientists"),
      list(var = "Q3F", name = "Trust in public health experts"),
      list(var = "Q3H", name = "Trust in journalists")
    )
    for (tm in trust_map) {
      if (tm$var %in% names(w27)) {
        vals <- as.numeric(w27[[tm$var]])
        vals[vals == 99] <- NA
        m <- pw_mean(vals, weights)
        s <- pw_sd(vals, weights)
        row[[paste0(tm$name, " (1-7)")]] <- sprintf("%.2f (%.2f)", m, s)
      }
    }

    # Approval of PHAC
    if ("Q10D" %in% names(w27)) {
      vals <- as.numeric(w27$Q10D)
      vals[vals == 99] <- NA
      m <- pw_mean(vals, weights)
      s <- pw_sd(vals, weights)
      row[["Approval of PHAC (1-5)"]] <- sprintf("%.2f (%.2f)", m, s)
    }

    # Media reliance
    media_map <- list(
      list(var = "Q35F", name = "Reliance on health professionals"),
      list(var = "Q35C", name = "Reliance on public health agencies"),
      list(var = "Q35G", name = "Reliance on traditional media"),
      list(var = "Q35A", name = "Reliance on social media")
    )
    for (mm in media_map) {
      if (mm$var %in% names(w27)) {
        vals <- as.numeric(w27[[mm$var]])
        vals[vals == 99] <- NA
        m <- pw_mean(vals, weights)
        s <- pw_sd(vals, weights)
        row[[paste0(mm$name, " (1-5)")]] <- sprintf("%.2f (%.2f)", m, s)
      }
    }

    # Concern about disinformation
    if ("Q9" %in% names(w27)) {
      vals <- as.numeric(w27$Q9)
      vals[vals == 99] <- NA
      m <- pw_mean(vals, weights)
      s <- pw_sd(vals, weights)
      row[["Concern about disinformation (1-5)"]] <- sprintf("%.2f (%.2f)", m, s)
    }

    results[[cls]] <- row
  }

  # --- Statistical tests ---
  cat("\n=== Statistical tests (posterior-weighted Wald) ===\n\n")

  # Sex (binary)
  qg <- as.numeric(w27$QGENDR)
  male_bin <- ifelse(!is.na(qg), qg == 1, NA)
  res <- pw_chisq_test(as.numeric(male_bin), post, labels, class_order)
  p_values[["Sex"]] <- res$p
  cat(sprintf("Sex: W=%.2f, p=%s\n", res$W, fmt_p(res$p)))

  # Age (4-category)
  qa <- as.numeric(w27$QAGEMIX)
  age_grp <- ifelse(!is.na(qa),
    ifelse(qa %in% c(1, 2), 1,
    ifelse(qa %in% c(3, 4), 2,
    ifelse(qa == 5, 3,
    ifelse(qa == 6, 4, NA)))), NA)
  res <- pw_chisq_multi_test(age_grp, c(1, 2, 3, 4), post, labels, class_order)
  p_values[["Age"]] <- res$p
  cat(sprintf("Age group: W=%.2f, p=%s\n", res$W, fmt_p(res$p)))

  # Education (binary)
  qe <- as.numeric(w27$QEDUC)
  univ <- ifelse(!is.na(qe), as.numeric(qe %in% c(7, 8)), NA)
  res <- pw_chisq_test(univ, post, labels, class_order)
  p_values[["Education"]] <- res$p
  cat(sprintf("Education: W=%.2f, p=%s\n", res$W, fmt_p(res$p)))

  # Rural (binary)
  qr <- as.numeric(w27$QRURBAN)
  rural <- ifelse(!is.na(qr), as.numeric(qr == 1), NA)
  res <- pw_chisq_test(rural, post, labels, class_order)
  p_values[["Rural"]] <- res$p
  cat(sprintf("Rural: W=%.2f, p=%s\n", res$W, fmt_p(res$p)))

  # Social class (binary)
  sc <- as.numeric(w27$SOCIALC1)
  poor_wc <- ifelse(!is.na(sc), as.numeric(sc %in% c(1, 2)), NA)
  res <- pw_chisq_test(poor_wc, post, labels, class_order)
  p_values[["Social class"]] <- res$p
  cat(sprintf("Social class: W=%.2f, p=%s\n", res$W, fmt_p(res$p)))

  # Vaccination (3-category)
  q15 <- as.numeric(w27$Q15)
  vax_valid <- ifelse(q15 %in% c(2, 3, 4), q15, NA)
  res <- pw_chisq_multi_test(vax_valid, c(2, 3, 4), post, labels, class_order)
  p_values[["Vaccination"]] <- res$p
  cat(sprintf("Vaccination: W=%.2f, p=%s\n", res$W, fmt_p(res$p)))

  # Trust variables (continuous, F-test)
  trust_vars <- list(
    list(var = "Q3A", name = "Trust: doctors"),
    list(var = "Q3E", name = "Trust: scientists"),
    list(var = "Q3F", name = "Trust: PH experts"),
    list(var = "Q3H", name = "Trust: journalists")
  )
  for (tv in trust_vars) {
    vals <- as.numeric(w27[[tv$var]])
    vals[vals == 99] <- NA
    res <- pw_F_test(vals, post, labels, class_order)
    p_values[[tv$name]] <- res$p
    cat(sprintf("%s: F=%.2f, p=%s\n", tv$name, res$F_stat, fmt_p(res$p)))
  }

  # Approval of PHAC
  vals <- as.numeric(w27$Q10D)
  vals[vals == 99] <- NA
  res <- pw_F_test(vals, post, labels, class_order)
  p_values[["Approval PHAC"]] <- res$p
  cat(sprintf("Approval PHAC: F=%.2f, p=%s\n", res$F_stat, fmt_p(res$p)))

  # Media reliance (continuous, F-test)
  media_vars <- list(
    list(var = "Q35F", name = "Reliance: health prof"),
    list(var = "Q35C", name = "Reliance: PH agencies"),
    list(var = "Q35G", name = "Reliance: trad media"),
    list(var = "Q35A", name = "Reliance: social media")
  )
  for (mv in media_vars) {
    vals <- as.numeric(w27[[mv$var]])
    vals[vals == 99] <- NA
    res <- pw_F_test(vals, post, labels, class_order)
    p_values[[mv$name]] <- res$p
    cat(sprintf("%s: F=%.2f, p=%s\n", mv$name, res$F_stat, fmt_p(res$p)))
  }

  # Concern about disinformation
  vals <- as.numeric(w27$Q9)
  vals[vals == 99] <- NA
  res <- pw_F_test(vals, post, labels, class_order)
  p_values[["Concern disinfo"]] <- res$p
  cat(sprintf("Concern disinfo: F=%.2f, p=%s\n", res$F_stat, fmt_p(res$p)))

  # --- Assemble data frame ---
  # Collect all row names (from first class's result list)
  row_names <- names(results[[1]])
  row_names <- row_names[row_names != "Class"]

  df <- data.frame(row.names = row_names)
  for (cls in class_order) {
    df[[cls]] <- sapply(row_names, function(rn) {
      val <- results[[cls]][[rn]]
      if (is.null(val)) "" else val
    })
  }

  # P-value column
  pval_map <- list(
    "N"                                    = "",
    "Male (%)"                             = p_values[["Sex"]],
    "Female (%)"                           = "",
    "Age 18-34 (%)"                        = p_values[["Age"]],
    "Age 35-54 (%)"                        = "",
    "Age 55-64 (%)"                        = "",
    "Age >= 65 (%)"                        = "",
    "Bachelor's or postgraduate (%)"       = p_values[["Education"]],
    "Rural community (%)"                  = p_values[["Rural"]],
    "Poor or working social class (%)"     = p_values[["Social class"]],
    "Unvaccinated (%)"                     = p_values[["Vaccination"]],
    "Three or more vaccine doses (%)"      = "",
    "Trust in doctors (1-7)"               = p_values[["Trust: doctors"]],
    "Trust in scientists (1-7)"            = p_values[["Trust: scientists"]],
    "Trust in public health experts (1-7)" = p_values[["Trust: PH experts"]],
    "Trust in journalists (1-7)"           = p_values[["Trust: journalists"]],
    "Approval of PHAC (1-5)"               = p_values[["Approval PHAC"]],
    "Reliance on health professionals (1-5)" = p_values[["Reliance: health prof"]],
    "Reliance on public health agencies (1-5)" = p_values[["Reliance: PH agencies"]],
    "Reliance on traditional media (1-5)"  = p_values[["Reliance: trad media"]],
    "Reliance on social media (1-5)"       = p_values[["Reliance: social media"]],
    "Concern about disinformation (1-5)"   = p_values[["Concern disinfo"]]
  )

  pval_col <- sapply(row_names, function(rn) {
    val <- pval_map[[rn]]
    if (is.null(val) || identical(val, "")) return("")
    if (is.numeric(val) && !is.na(val)) return(fmt_p(val))
    ""
  })
  df[["P value"]] <- pval_col

  df
}


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    stop("Usage: Rscript lca_table3.R path/to/wave27.sav [--output table3.csv]")
  }

  sav_path <- args[1]
  output   <- "table3_baseline_by_class.csv"

  # Parse --output flag
  if ("--output" %in% args || "-o" %in% args) {
    flag_pos <- which(args %in% c("--output", "-o"))
    if (length(flag_pos) > 0 && flag_pos[1] < length(args)) {
      output <- args[flag_pos[1] + 1]
    }
  }

  K <- 4L
  CLASS_ORDER <- c("Well-Informed", "Selectively Skeptical",
                    "Uncertain/Disengaged", "Misinformed")

  # Load and recode
  loaded    <- load_and_recode(sav_path)
  w27       <- loaded$data
  item_cols <- loaded$item_cols
  X <- as.matrix(w27[, item_cols])
  storage.mode(X) <- "integer"

  # Fit LCA
  cat("\nRunning EM with 100 random starts...\n")
  fit <- em_lca(X, K)
  labels_vec <- label_classes(fit$rho, K)
  cat(sprintf("Labels: %s\n", paste(labels_vec, collapse = ", ")))
  cat(sprintf("Class proportions: %s\n",
              paste(round(fit$pi, 3), collapse = ", ")))
  cat(sprintf("Log-likelihood: %.2f\n", fit$ll))

  # Build table
  df <- compute_table(w27, fit$pi, fit$rho, fit$post, labels_vec, CLASS_ORDER, K)

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("TABLE 3: Baseline Characteristics by Latent Class (Wave 27)\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  print(df)

  write.csv(df, output)
  cat(sprintf("\nSaved to %s\n", output))
}

main()
