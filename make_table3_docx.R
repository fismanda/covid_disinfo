#!/usr/bin/env Rscript
#
# Generate Table 3 as a formatted Word document from the CSV output
# of lca_table3.R.
#
# Usage:
#   Rscript make_table3_docx.R [--input table3.csv] [--output Table_3.docx]
#
# Dependencies: flextable, officer

suppressPackageStartupMessages({
  library(flextable)
  library(officer)
})


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  input_file  <- "table3_baseline_by_class.csv"
  output_file <- "Table_3_baseline_by_class.docx"

  # Parse arguments
  if ("--input" %in% args || "-i" %in% args) {
    flag_pos <- which(args %in% c("--input", "-i"))
    if (length(flag_pos) > 0 && flag_pos[1] < length(args)) {
      input_file <- args[flag_pos[1] + 1]
    }
  }
  if ("--output" %in% args || "-o" %in% args) {
    flag_pos <- which(args %in% c("--output", "-o"))
    if (length(flag_pos) > 0 && flag_pos[1] < length(args)) {
      output_file <- args[flag_pos[1] + 1]
    }
  }

  # Read CSV
  df <- read.csv(input_file, row.names = 1, check.names = FALSE)

  col_order <- c("Well-Informed", "Selectively Skeptical",
                 "Uncertain/Disengaged", "Misinformed")
  df <- df[, c(col_order, "P value")]

  # Row definitions: list of (type, display_label, csv_key)
  rows_spec <- list(
    list("header", "Demographics", NA),
    list("data",   "N (% of sample)", "N"),
    list("data",   "Male (%)", "Male (%)"),
    list("data",   "Female (%)", "Female (%)"),
    list("data",   "Age 18-34 (%)", "Age 18-34 (%)"),
    list("data",   "Age 35-54 (%)", "Age 35-54 (%)"),
    list("data",   "Age 55-64 (%)", "Age 55-64 (%)"),
    list("data",   "Age ≥ 65 (%)", "Age >= 65 (%)"),
    list("data",   "Bachelor's or postgraduate education (%)",
                   "Bachelor's or postgraduate (%)"),
    list("data",   "Rural community (%)", "Rural community (%)"),
    list("data",   "Poor or working social class (%)",
                   "Poor or working social class (%)"),
    list("header", "Vaccination status", NA),
    list("data",   "Unvaccinated (%)", "Unvaccinated (%)"),
    list("data",   "Three or more doses (%)",
                   "Three or more vaccine doses (%)"),
    list("header", "Trust (1-7 scale), mean (SD)", NA),
    list("data",   "Doctors", "Trust in doctors (1-7)"),
    list("data",   "Scientists", "Trust in scientists (1-7)"),
    list("data",   "Public health experts",
                   "Trust in public health experts (1-7)"),
    list("data",   "Journalists", "Trust in journalists (1-7)"),
    list("header", "Information sources and attitudes", NA),
    list("data",   "Approval of PHAC (1-5), mean (SD)",
                   "Approval of PHAC (1-5)"),
    list("data",   "Reliance on health professionals (1-5), mean (SD)",
                   "Reliance on health professionals (1-5)"),
    list("data",   "Reliance on public health agencies (1-5), mean (SD)",
                   "Reliance on public health agencies (1-5)"),
    list("data",   "Reliance on traditional media (1-5), mean (SD)",
                   "Reliance on traditional media (1-5)"),
    list("data",   "Reliance on social media (1-5), mean (SD)",
                   "Reliance on social media (1-5)"),
    list("data",   "Concern about disinformation (1-5), mean (SD)",
                   "Concern about disinformation (1-5)")
  )

  # Build output data frame for the table
  n_rows <- length(rows_spec)
  out <- data.frame(
    Characteristic = character(n_rows),
    WI  = character(n_rows),
    SS  = character(n_rows),
    UD  = character(n_rows),
    MI  = character(n_rows),
    Pvalue = character(n_rows),
    stringsAsFactors = FALSE
  )

  header_rows <- integer(0)

  for (i in seq_along(rows_spec)) {
    spec <- rows_spec[[i]]
    if (spec[[1]] == "header") {
      out$Characteristic[i] <- spec[[2]]
      out$WI[i] <- ""
      out$SS[i] <- ""
      out$UD[i] <- ""
      out$MI[i] <- ""
      out$Pvalue[i] <- ""
      header_rows <- c(header_rows, i + 1)  # +1 because flextable adds header row
    } else {
      out$Characteristic[i] <- spec[[2]]
      csv_key <- spec[[3]]
      if (csv_key %in% rownames(df)) {
        out$WI[i] <- as.character(df[csv_key, "Well-Informed"])
        out$SS[i] <- as.character(df[csv_key, "Selectively Skeptical"])
        out$UD[i] <- as.character(df[csv_key, "Uncertain/Disengaged"])
        out$MI[i] <- as.character(df[csv_key, "Misinformed"])
        pv <- df[csv_key, "P value"]
        out$Pvalue[i] <- if (is.na(pv)) "" else as.character(pv)
      }
    }
  }

  # Create flextable
  ft <- flextable(out)
  ft <- set_header_labels(ft,
    Characteristic = "Characteristic",
    WI     = "Well-Informed",
    SS     = "Selectively Skeptical",
    UD     = "Uncertain/Disengaged",
    MI     = "Misinformed",
    Pvalue = "P value"
  )

  # Formatting
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- bold(ft, part = "header")
  ft <- italic(ft, j = "Pvalue", part = "header")
  ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- align(ft, j = 2:6, align = "center", part = "all")

  # Bold and shade header rows
  for (hr in header_rows) {
    body_row <- hr - 1  # convert back to body row index
    if (body_row >= 1 && body_row <= nrow(out)) {
      ft <- bold(ft, i = body_row, j = 1)
      ft <- bg(ft, i = body_row, bg = "#F2F2F2")
    }
  }

  # Column widths
  ft <- width(ft, j = 1, width = 2.4)
  ft <- width(ft, j = 2:5, width = 0.9)
  ft <- width(ft, j = 6, width = 0.7)

  # Create Word document
  doc <- read_docx()

  # Set page to US Letter
  doc <- body_set_default_section(doc,
    prop_section(
      page_size = page_size(width = 8.5, height = 11, orient = "portrait"),
      page_margins = page_mar(
        bottom = 0.75, top = 0.75, left = 0.75, right = 0.75,
        header = 0.3, footer = 0.3
      )
    )
  )

  # Title
  title_fp <- fp_text(font.family = "Times New Roman", font.size = 11)
  title_bold <- update(title_fp, bold = TRUE)
  doc <- body_add_fpar(doc,
    fpar(
      ftext("Table 3. ", prop = title_bold),
      ftext(paste0(
        "Baseline (Wave 27) characteristics by latent class, estimated using ",
        "posterior-weighted proportions and means from the four-class common-item ",
        "latent class model."
      ), prop = title_fp)
    )
  )

  # Add table
  doc <- body_add_flextable(doc, ft, align = "center")

  # Footnote
  footnote_text <- paste0(
    "Note: Estimates are posterior-weighted from the four-class latent class model ",
    "fitted to the four common misinformation items (COVID-19 deaths exaggerated, ",
    "vaccine deaths hidden, vaccines cause infertility, vaccines alter DNA). ",
    "N = 1,185 respondents with complete data on all four items. ",
    "P values are from posterior-weighted Wald tests: chi-square tests for ",
    "categorical variables and F-tests for continuous variables. ",
    "Trust scales range from 1 (no trust at all) to 7 (a great deal of trust). ",
    "Reliance and approval scales range from 1 to 5. ",
    "PHAC = Public Health Agency of Canada."
  )

  footnote_fp <- fp_text(font.family = "Times New Roman", font.size = 9, italic = TRUE)
  doc <- body_add_fpar(doc, fpar(ftext(footnote_text, prop = footnote_fp)))

  print(doc, target = output_file)
  cat(sprintf("Saved to %s\n", output_file))
}

main()
