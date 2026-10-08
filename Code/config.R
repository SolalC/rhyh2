# =============================================================================
# rhyh2 — Wheeler et al. (2016) local (cis) heritability pipeline: configuration
#
# Every script sources this file. To use another config (e.g. the toy test),
# set the environment variable RHYH2_CONFIG to its path.
# =============================================================================

# ---- PATHS TO EDIT (Bunya) --------------------------------------------------
# GTEx v8 WGS genotypes (dbGaP phs000424.v8), one multi-chromosome VCF (GRCh38)
GTEX_VCF      <- "/PATH/TO/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"
# Directory with <Tissue>.v8.normalized_expression.bed.gz (GTEx_Analysis_v8_eQTL_expression_matrices)
GTEX_EXPR_DIR <- "/PATH/TO/GTEx_Analysis_v8_eQTL_expression_matrices"
# Directory with <Tissue>.v8.covariates.txt (GTEx_Analysis_v8_eQTL_covariates)
GTEX_COV_DIR  <- "/PATH/TO/GTEx_Analysis_v8_eQTL_covariates"
# GTEx v8 collapsed gene model (GENCODE v26, GRCh38): gencode.v26.GRCh38.genes.gtf
GENE_GTF      <- "/PATH/TO/gencode.v26.GRCh38.genes.gtf"
# Optional BED (GRCh38, 0-based) of SNP positions to keep. Wheeler et al. kept
# HapMap Phase II SNPs only; supply a HapMap II list lifted to GRCh38 to mimic
# that. NA = keep all common biallelic SNPs (a documented deviation).
SNP_KEEP_BED  <- NA
# Optional file of donor IDs (one per line) to restrict the analysis to, e.g.
# European-ancestry donors. NA = all genotyped donors.
DONOR_KEEP    <- NA
# Output root (use /scratch on Bunya; nothing here is small)
OUT_DIR       <- "/PATH/TO/scratch/rhyh2/wheeler2016"
# Optional: GenArchDB SQLite file with Wheeler et al.'s per-gene estimates
# (see reference/get_genarchdb.sh). NA = skip the per-gene comparison.
GENARCH_DB    <- NA
# Executables (module names on Bunya may differ; full paths are safest)
PLINK2        <- "plink2"
GCTA          <- "gcta64"
# Node-local temporary directory for per-gene GRMs
TMP_ROOT      <- Sys.getenv("TMPDIR", unset = tempdir())

# ---- ANALYSIS PARAMETERS (defaults follow Wheeler et al. 2016) ---------------
TISSUES       <- NULL        # NULL = every tissue found in GTEX_EXPR_DIR
MIN_SAMPLES   <- 70          # Wheeler: tissues with n > 70
MAF_MIN       <- 0.05        # common SNPs only
DROP_AMBIGUOUS_SNPS <- TRUE  # drop A/T and C/G SNPs
AUTOSOMES     <- 1:22
GENE_TYPES    <- "protein_coding"
CIS_WINDOW_BP <- 1e6         # SNPs within 1 Mb of gene start and end
MODEL         <- "local_only"   # Wheeler's final model; "local_distal" adds a
                                # leave-one-chromosome-out GRM
REML_MODES    <- c("unconstrained", "constrained")  # means use unconstrained
# Covariates: Wheeler used 15 PEER factors + sex, no genotype PCs.
COV_SOURCE    <- "gtex_peer"    # "gtex_peer" = first N_FACTORS InferredCov
                                # columns of the GTEx v8 covariate file;
                                # "expression_pcs" = top N_FACTORS PCs of the
                                # normalised expression, computed here
N_FACTORS     <- 15
N_GENO_PCS    <- 0              # GTEx v8 file provides PC1-PC5
USE_SEX       <- TRUE
REML_MAXIT    <- 200
REML_ALGS     <- c(0, 1, 2)     # retry order: AI (GCTA default), Fisher
                                # scoring, EM. c(0) = GCTA default only
N_CHUNKS      <- 200            # SLURM array size for step 04
THREADS       <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
FDR_LEVEL     <- 0.1
SEED          <- 20261008
# Permutation (negative control): run steps 03-05 with --perm-seed=<n>.
# Genotype-to-sample labels are shuffled within tissue; expected mean h2 = 0.
PERM_GENE_FRACTION <- 0.1       # fraction of genes analysed in permuted runs
