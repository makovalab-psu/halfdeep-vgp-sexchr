# halfdeep-vgp-sexchr

Scripts used to run HalfDeep for VGP sex chromosome assembly (phase 1 freeze) validation.

This repository includes:
- modified HalfDeep execution script: `halfdeep_ko.sh`
- wrapper scripts for running HalfDeep with local/AWS or SRA read inputs
- master scripts for multi-species and per-species (mostly >20 Gbp genomes) parallel execution

`halfdeep_ko.sh` is adapted from the original `halfdeep.sh` to run HalfDeep downstream steps after manual mapping/depth generation.
