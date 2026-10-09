# Toy-test configuration: the main config with paths pointing at simulated data.
# Set by run_toy_test.sh via RHYH2_CONFIG; TOY_DIR must be set in the environment.
source(file.path(Sys.getenv("RHYH2_CODE_DIR"), "config.R"), local = globalenv())

TOY_DIR       <- Sys.getenv("TOY_DIR")
GTEX_RELEASE  <- "v8"                # toy files are named .v8.
GTEX_VCF      <- file.path(TOY_DIR, "toy.vcf.gz")
GTEX_EXPR_DIR <- file.path(TOY_DIR, "expression")
GTEX_COV_DIR  <- file.path(TOY_DIR, "covariates")
GENE_GTF      <- file.path(TOY_DIR, "toy.genes.gtf")
OUT_DIR       <- file.path(TOY_DIR, "out")
PLINK2        <- Sys.getenv("PLINK2", "plink2")
GCTA          <- Sys.getenv("GCTA", "gcta64")
AUTOSOMES     <- 1:3
N_CHUNKS      <- 4
PERM_GENE_FRACTION <- 1
MODEL         <- Sys.getenv("TOY_MODEL", "local_only")
