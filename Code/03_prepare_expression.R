# Step 03 — expression phenotypes and covariates per tissue (GCTA format).
#   Expression: GTEx v8 normalised matrices (TMM + inverse-normal transform),
#   restricted to genes in genes.tsv and to donors with genotypes.
#   Covariates: N_FACTORS PEER factors (or expression PCs), optional genotype
#   PCs, and sex (dropped when constant, e.g. sex-specific tissues).
#
# Usage: Rscript 03_prepare_expression.R [--perm-seed=<n>]
#   With --perm-seed, donor labels are shuffled within tissue (expression and
#   covariates move together), breaking only the genotype-expression link.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
args <- parse_args()
perm_seed <- if (is.null(args$`perm-seed`)) NA else as.integer(args$`perm-seed`)
label <- run_label(perm_seed)
p <- out_paths(OUT_DIR, label)

genes <- read_tsv(file.path(p$genes, "genes.tsv"), col_types = cols(chr = "c"), progress = FALSE)
genotyped <- read_fam(file.path(p$geno, paste0("chr", AUTOSOMES[1])))$iid

expr_files <- list.files(GTEX_EXPR_DIR, pattern = "\\.v8\\.normalized_expression\\.bed\\.gz$", full.names = TRUE)
tissue_names <- str_remove(basename(expr_files), "\\.v8\\.normalized_expression\\.bed\\.gz$")
if (!is.null(TISSUES)) {
  keep <- tissue_names %in% TISSUES
  expr_files <- expr_files[keep]; tissue_names <- tissue_names[keep]
}
if (length(expr_files) == 0) stop("No expression files found in ", GTEX_EXPR_DIR)

prepare_tissue <- function(tissue, expr_file, tissue_index) {
  expr <- read_tsv(expr_file, col_types = cols(.default = "d", `#chr` = "c", gene_id = "c"),
                   progress = FALSE) |>
    mutate(gene_key = strip_version(gene_id)) |>
    semi_join(genes |> mutate(gene_key = strip_version(gene_id)), by = "gene_key")

  cov_raw <- read_tsv(file.path(GTEX_COV_DIR, paste0(tissue, ".v8.covariates.txt")),
                      col_types = cols(.default = "d", ID = "c"), progress = FALSE)

  donors <- intersect(intersect(names(expr), genotyped), names(cov_raw))
  if (length(donors) < MIN_SAMPLES) {
    message(tissue, ": n = ", length(donors), " < ", MIN_SAMPLES, ", skipped")
    return(tibble(tissue = tissue, n = length(donors), n_genes = 0L, status = "skipped_small_n"))
  }

  y <- as.matrix(expr[, donors])                      # genes x donors
  rownames(y) <- genes$gene_id[match(expr$gene_key, strip_version(genes$gene_id))]

  # ---- covariates ----
  cov_t <- cov_raw |> select(ID, all_of(donors))
  get_rows <- function(ids) {
    m <- cov_t |> filter(ID %in% ids)
    if (nrow(m) < length(ids)) stop(tissue, ": covariate rows missing: ",
                                    paste(setdiff(ids, m$ID), collapse = ", "))
    m <- m[match(ids, m$ID), ]
    t(as.matrix(m[, donors])) |> `colnames<-`(ids)
  }
  factors <- switch(COV_SOURCE,
    gtex_peer = get_rows(paste0("InferredCov", seq_len(N_FACTORS))),
    expression_pcs = {
      pcs <- prcomp(t(y), center = TRUE, scale. = TRUE, rank. = N_FACTORS)$x
      colnames(pcs) <- paste0("exprPC", seq_len(N_FACTORS)); pcs
    },
    stop("Unknown COV_SOURCE: ", COV_SOURCE)
  )
  qcov <- if (N_GENO_PCS > 0) cbind(factors, get_rows(paste0("PC", seq_len(N_GENO_PCS)))) else factors
  sex <- if (USE_SEX) get_rows("sex")[, 1] else NULL
  use_sex <- !is.null(sex) && length(unique(sex)) > 1

  # ---- donor labels (permuted for the negative control) ----
  geno_id <- donors
  if (!is.na(perm_seed)) {
    set.seed(perm_seed * 1000 + tissue_index)
    geno_id <- sample(donors)
  }
  id_cols <- tibble(FID = geno_id, IID = geno_id)

  stem <- file.path(p$pheno, tissue)
  bind_cols(id_cols, as_tibble(t(y), .name_repair = "minimal")) |>
    write_tsv(paste0(stem, ".phen"), col_names = FALSE, na = "NA")
  write_tsv(tibble(gene_id = rownames(y), mpheno = seq_len(nrow(y))), paste0(stem, ".genes.tsv"))
  bind_cols(id_cols, as_tibble(qcov)) |> write_tsv(paste0(stem, ".qcovar"), col_names = FALSE)
  if (use_sex) {
    bind_cols(id_cols, tibble(sex = sex)) |> write_tsv(paste0(stem, ".covar"), col_names = FALSE)
  } else {
    unlink(paste0(stem, ".covar"))
  }
  write_tsv(id_cols, paste0(stem, ".keep"), col_names = FALSE)

  message(tissue, ": n = ", length(donors), ", genes = ", nrow(y),
          ", quantitative covariates = ", ncol(qcov), ", sex = ", use_sex)
  tibble(tissue = tissue, n = length(donors), n_genes = nrow(y), status = "ok",
         n_qcovar = ncol(qcov), sex_covariate = use_sex, cov_source = COV_SOURCE)
}

tissues <- pmap(list(tissue_names, expr_files, seq_along(tissue_names)), prepare_tissue) |>
  list_rbind()
write_tsv(tissues, file.path(p$pheno, "tissues.tsv"))
print(tissues, n = Inf)
