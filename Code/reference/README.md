# Reference data

- `wheeler2016_table1.csv` — Table 1 of Wheeler et al. (2016) PLoS Genet 12:e1006423,
  doi:10.1371/journal.pgen.1006423, extracted programmatically from the PMC open-access
  JATS XML (PMC5106030, `pmc-oa-opendata` bucket) on 2026-10-08. Columns: tissue,
  n, mean unconstrained local h² across genes, SE "estimated from the h2 estimates
  obtained when the expression sample labels are shuffled 100 times", % and number
  of genes with FDR < 0.1, number of expressed genes (mean RPKM > 0.1).
- `get_genarchdb.sh` — downloads the authors' per-gene estimates (GenArchDB). These are
  constrained estimates from the preprint stage and do not reproduce Table 1 exactly
  (e.g. DGN mean 0.126 constrained vs 0.149 unconstrained in the paper).
