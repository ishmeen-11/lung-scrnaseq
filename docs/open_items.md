# Open items

## Received from Griffin (2026-10-09)

- [x] SRA/GEO accessions for all 11 samples → `config/samples.tsv`
- [x] SoupX outputs (`SoupX_Outs.zip`, 2 GB, shared Drive folder) → unzip into `data/`
- [x] Doses, timepoints, mouse sex → `config/samples.tsv` and below
- [x] Final annotation marker lists → `config/config.yaml` (`annotation.markers`)
- [ ] Cell Ranger outs + BAMs (~430 GB total): waiting on a lab AWS account (Columbia enterprise agreement) or an external drive
- [ ] Ginsburg access

## Study design facts (from Griffin)

| Study | GEO | Inoculum | Dose (CFU in 50 µL) | Harvest |
|---|---|---|---|---|
| Wong Fok Lung 2022 | GSE190262 | *K. pneumoniae* MKP103 | 1 × 10⁸ | 48 hpi |
| Tomlinson 2023 (Pa) | GSE203352 | *P. aeruginosa* PAO1 or 17 clinical isolates (CS17) | 1 × 10⁶ | 24 hpi |
| Tomlinson 2023 (Sa) | GSE215195 | *S. aureus* USA300 LAC | 2 × 10⁷ | 24 hpi |

- Mice: C57BL/6NJ, WT or *Irg1*-/-; adults of similar age.
- Each sample pools cells from several mice, ~50:50 male:female. So the sex check should read "mixed"; a single-sex call would be a red flag.
- The Pa and Sa papers share the same PBS samples (PBS_WT_Tom, PBS_IRG1).

## Confounds to keep in mind

- **Kp differs in dose and timepoint** (10⁸ CFU, 48 h vs 24 h for the others). Kp vs. others mixes pathogen with infection stage.
- **Study = batch.** Each pathogen comes from one study; the two WT PBS controls (Wong vs Tom) are the only direct batch comparison.
- **One pooled sample per condition.** Differences between conditions are descriptive; no per-condition replicates for statistics.

## Milestones

| # | Milestone | Status |
|---|---|---|
| M0 | Repo, config, environment, demo, tests | done |
| M1 | Real data in hand (SoupX outputs) | done; raw/Cell Ranger pending AWS |
| M2 | Preprocessing on real data; compare with `currentObject.rds` | next |
| M3 | Streamlit QC viewer (+ cell-type composition) | done (demo) |
| M4 | Analysis plan agreed with Prof. Prince (start with *S. aureus*) | meeting Wed 10/21 3pm |
| M5 | Analysis: Augur, Milo, gene-set scoring, monocyte pseudotime, Gram+ vs Gram− | |
| M6 | Report | |
