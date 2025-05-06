# UCE Phylogenomics Pipeline

This folder contains the pipeline used to extract ultraconserved elements (UCEs), align sequences, and reconstruct phylogenetic trees for the study:  
**"Adaptive divergence in Gentoo penguins: An integrative approach and future projections."**

---

## 1. UCE identification

We re-identified the set of ultraconserved elements (UCEs) shared across all genomes included in this study, including the South Georgia reference genome (GCA_030674165.1), to ensure consistency across samples.

### a) Download probe set
```bash
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta
```

### b) Match probes to genomes using LASTZ (via PHYLUCE)
faToTwoBit $MSK/${sample}_maskedPygo.fa $UCE/${sample}/${sample}.2bit
twoBitInfo $UCE/${sample}/${sample}.2bit stdout | sort -k2rn > $UCE/${sample}/sizes.tab

phyluce_probe_run_multiple_lastzs_sqlite \
    --db penguin.sqlite \
    --output output-lastz \
    --scaffoldlist <samples> \
    --genome-base-path ./ \
    --probefile uce-5k-probes.fasta \
    --cores 30

### c) Extract UCEs with 750 bp flanks
phyluce_probe_slice_sequence_from_genomes \
    --lastz $UCE/output-lastz \
    --conf $UCE/genomes.conf \
    --flank 750 \
    --name-pattern "uce-5k-probes.fasta_v_{}.lastz.clean" \
    --output $UCE/UCEs_all_fasta


## 2. UCE filtering criteria
Only UCE loci present in 100% of the following genomes were retained:
**Adelie** (adel)
**Chinstrap** (chins)
**South Georgia** (southgeorgia)
**Signy** (gsa10, gsa14)
**King George Island** (p0407, p0409)
**O'Higgins** (p0196, p0200)
**GGV** (p0048, p0052)
**Martillo Island** (pam10, pam1, pam3)
**Falkland Islands** (g5, p1884)
**Kerguelen** (15a2, 22a1, 24a1)
**Macquarie Island** (mq9209, mq9312)
**Crozet Islands** (cro07, cro08, p2136)
**Marion Island** (gtp06, gen04, gtp12)

This strict filtering ensured that all retained loci were homologous and aligned across all individuals used in the phylogenomic and divergence time analyses.

## 3. Alignment with MAFFT
mafft --thread $THREADS_PER_JOB \
      --genafpair \
      --maxiterate 1000 \
      --ep 0.123 \
      "$input_file" > "$output_file" 2> "$log_file"

## 4. Phylogenetic inference
iqtree2 -s $alignment \
        --prefix $GENE_TREES_DIR/$locus_name \
        -m MFP \
        -B 1000 \
        -nm 4000 \
        --scfl 100
