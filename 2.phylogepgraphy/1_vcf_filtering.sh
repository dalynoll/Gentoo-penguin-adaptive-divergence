### Outlier remotion
bedtools sort -i gentoo_all_outlier_windows.bed > gentoo_all_outlier_windows_sorted.bed

bedtools merge -i gentoo_all_outlier_windows_sorted.bed > gentoo_all_outlier_windows_sorted_merged.bed

vcftools --gzvcf $vcf/gentoo.vcf.gz --keep gentoo64_toVCFoutlier.txt \
--exclude-bed gentoo_all_outlier_windows_sorted_merged.bed --recode --out gentoo_no_outliers

# MAF 0.05
vcftools --gzvcf gentoo_no_outliers.vcf.gz \
--maf 0.05 --recode --stdout | bgzip > gentoo_no_outliers_maf05.vcf.gz

# LD pruning
plink --vcf gentoo_no_outliers_maf05.vcf.gz \
--indep-pairwise 50 10 0.1 \
--set-missing-var-ids @:# --double-id --allow-extra-chr \
--out gentoo_no_outliers_maf05_prun

plink --vcf gentoo_no_outliers_maf05.vcf.gz  \
--extract gentoo_no_outliers_prun.prune.in \
--set-missing-var-ids @:# --double-id \
--allow-extra-chr --make-bed --out gentoo_no_outliers_maf05_pruned

plink --vcf gentoo_no_outliers_maf05.vcf.gz  \
--extract gentoo_no_outliers_maf05_prun.prune.in \
--set-missing-var-ids @:# --double-id --allow-extra-chr \
--recode vcf --out gentoo_no_outliers_maf05_pruned
