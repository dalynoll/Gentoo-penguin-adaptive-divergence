# UCE Phylogenomics Pipeline

This folder contains the pipeline used to extract ultraconserved elements (UCEs), align sequences, and reconstruct phylogenetic trees for the study:  
**"Adaptive divergence in Gentoo penguins: An integrative approach and future projections."**

---

## 1. UCE identification

We re-identified the set of ultraconserved elements (UCEs) shared across all genomes included in this study, including the South Georgia reference genome (GCA_030674165.1), to ensure consistency across samples.

### a) Download probe set
```bash
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta


### b) Match probes to genomes using LASTZ (via PHYLUCE)
