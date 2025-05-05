# Gentoo-penguin-adaptive-divergence

This repository contains the scripts, input files, and documentation used in the study:
"Adaptive divergence in Gentoo penguins: An integrative approach and future projections"

## Repository structure

- **1.read mapping, variant calling/**  
  Scripts and parameters for mapping Illumina reads to the reference genome, BAM filtering, and variant calling using GATK and bcftools.

- **2.phylogeography/**  
  Analyses of population structure, PCA, ADMIXTURE, heterozygosity estimates, F<sub>ST</sub>, D<sub>XY</sub>, and Treemix to assess historical connectivity and divergence among colonies.

- **3.phylogenomics/**  
    UCE alignment, ASTRAL species tree inference, and divergence time estimation using StarBEAST3.

- **4.selective sweep analysis/**  
  Detection of selective sweeps using RAiSD and nSL, including BED file generation and intersection analyses between methods and lineages.

- **5.Morphology/**  
  Raw data and statistical analyses (MANOVA, pairwise tests) of morphometric differences among gentoo penguin lineages using museum specimens.

- **README.md**  
  This file.

- **genbank_biosample_access.tsv**  
  Sample metadata used for GenBank and BioSample submission, including accession numbers and locality information.

---

## Overview
This study investigates the evolutionary processes driving divergence among four major gentoo penguin lineages. It integrates:

- Whole-genome resequencing
- Population genomic and phylogenomic analyses
- Selective sweep detection
- Genotype–environment associations
- Morphometric comparisons
- Species distribution modeling (provided in supplementary)

---

## Software and Dependencies

- **Read mapping & variant calling:** `bwa`, `samtools`, `picard`, `gatk`, `bcftools`
- **Population analyses:** `vcftools`, `plink`, `pixy`, `RAiSD`, `angsd`
- **Phylogenetics & divergence times:** `IQ-TREE`, `ASTRAL`, `BEAST2 / StarBEAST3`
- **Morphology:** `R` with `vegan`, `ggplot2`, `factoextra`, `mvtnorm`
- **Sweep detection:** `RAiSD`, `selscan (nSL)`


## Contact

For questions or collaborations:

**Daly Noll**  
dalynolll@gmail.com
GitHub: [@dalynoll](https://github.com/dalynoll)
