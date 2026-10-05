# Open items

## Needed from Griffin

- [ ] Ginsburg access (Slurm account; his scripts use `rent`)
- [ ] SRA/ENA run accessions for all 11 samples → `config/samples.tsv`
- [ ] Location of existing Cell Ranger outputs + BAMs, and SoupX folders
- [ ] Timepoint, inoculum dose and mouse sex per study
- [ ] Preferred long-term storage (AWS S3 or other) and who pays
- [ ] Final annotation gene lists (one version)

## Milestones

| # | Milestone | Status |
|---|---|---|
| M0 | Repo, config, environment, demo, tests | done |
| M1 | Raw data on Ginsburg, Cell Ranger for all 11 samples | waiting on access + accessions |
| M2 | Final preprocessing: reproduce Griffin (fixed mode), then compare MAD mode | code ready |
| M3 | Streamlit QC viewer | done (demo) |
| M4 | Analysis plan agreed with Prof. Prince (start with *S. aureus*) | |
| M5 | Analysis: Augur, Milo, gene-set scoring, monocyte pseudotime, Gram+ vs Gram− | |
| M6 | Report | |

## Planned pipeline additions

- CellBender rule (ambient RNA) and comparison with SoupX
- scVI integration option, benchmarked against Harmony
- velocyto / STARsolo rule for spliced/unspliced counts (RNA velocity), behind `velocity.enabled`
