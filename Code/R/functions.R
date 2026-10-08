# Shared helpers for the Wheeler et al. (2016) pipeline. Sourced by every step.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# ---- set-up -----------------------------------------------------------------

code_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) return(normalizePath("."))
  normalizePath(dirname(sub("^--file=", "", file_arg[1])))
}

load_config <- function() {
  cfg_file <- Sys.getenv("RHYH2_CONFIG", unset = file.path(code_dir(), "config.R"))
  if (!file.exists(cfg_file)) stop("Config not found: ", cfg_file)
  source(cfg_file, local = globalenv())
  message("Config: ", normalizePath(cfg_file))
  invisible(cfg_file)
}

# Parse "--key=value" arguments into a named list
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  kv <- str_match(args, "^--([^=]+)=(.*)$")
  ok <- !is.na(kv[, 1])
  setNames(as.list(kv[ok, 3]), kv[ok, 2])
}

run_label <- function(perm_seed = NA) {
  if (is.na(perm_seed)) "observed" else paste0("perm", perm_seed)
}

out_paths <- function(out_dir, label = "observed") {
  p <- list(
    geno    = file.path(out_dir, "genotypes"),
    grm     = file.path(out_dir, "grm"),
    genes   = file.path(out_dir, "genes"),
    pheno   = file.path(out_dir, "pheno", label),
    results = file.path(out_dir, "results", label),
    summary = file.path(out_dir, "summary", label),
    logs    = file.path(out_dir, "logs")
  )
  walk(p, dir.create, recursive = TRUE, showWarnings = FALSE)
  p
}

# ---- external programs ------------------------------------------------------

# Run a command; stop with its log tail if it fails (unless allow_fail = TRUE)
run_cmd <- function(exe, args, log_file = NULL, allow_fail = FALSE) {
  args <- as.character(args)
  out <- suppressWarnings(system2(exe, shQuote(args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status") %||% 0L
  if (!is.null(log_file)) write_lines(c(paste(exe, paste(args, collapse = " ")), out), log_file)
  if (status != 0 && !allow_fail) {
    stop(exe, " failed (status ", status, "):\n", paste(tail(out, 20), collapse = "\n"))
  }
  invisible(list(status = status, output = out))
}

# ---- file readers -------------------------------------------------------------

read_bim <- function(prefix) {
  read_tsv(paste0(prefix, ".bim"), col_names = c("chr", "id", "cm", "pos", "a1", "a2"),
           col_types = "cccicc", progress = FALSE)
}

read_fam <- function(prefix) {
  read_table(paste0(prefix, ".fam"), col_names = c("fid", "iid", "pat", "mat", "sex", "phen"),
             col_types = "cccccc", progress = FALSE)
}

# Genes from a GTF: one row per "gene" feature
read_gtf_genes <- function(gtf) {
  read_tsv(gtf, comment = "#", col_names = FALSE, col_types = "ccciicccc", progress = FALSE) |>
    filter(X3 == "gene") |>
    transmute(
      chr       = str_remove(X1, "^chr"),
      start     = X4,
      end       = X5,
      strand    = X7,
      gene_id   = str_match(X9, 'gene_id "([^"]+)"')[, 2],
      gene_name = str_match(X9, 'gene_name "([^"]+)"')[, 2],
      gene_type = str_match(X9, 'gene_type "([^"]+)"')[, 2]
    )
}

# GCTA .hsq -> one-row tibble. Handles 1 or 2 genetic components.
read_hsq <- function(file) {
  lines <- read_lines(file)
  fields <- str_split(lines, "\t")
  get <- function(key, i = 2) {
    hit <- keep(fields, ~ .x[1] == key)
    if (length(hit) == 0) NA_real_ else suppressWarnings(as.numeric(hit[[1]][i]))
  }
  tibble(
    h2_local     = get("V(G1)/Vp"), se_local  = get("V(G1)/Vp", 3),
    h2_distal    = get("V(G2)/Vp"), se_distal = get("V(G2)/Vp", 3),
    h2_single    = get("V(G)/Vp"),  se_single = get("V(G)/Vp", 3),
    vp           = get("Vp"),
    logL         = get("logL"),
    lrt          = get("LRT"),
    p_lrt        = get("Pval"),
    n            = get("n")
  ) |>
    # with a single GRM, GCTA labels the component V(G) instead of V(G1)
    mutate(h2_local = coalesce(h2_local, h2_single),
           se_local = coalesce(se_local, se_single)) |>
    select(-h2_single, -se_single)
}

# ---- misc -------------------------------------------------------------------

# Harmonise tissue labels across GTEx v8, Wheeler Table 1 and GenArchDB
normalise_tissue <- function(x) {
  x |>
    str_remove("_(TW|TS)$") |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]", "")
}

strip_version <- function(ensg) str_remove(ensg, "\\.\\d+$")
