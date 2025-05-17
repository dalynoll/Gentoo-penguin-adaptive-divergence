
#####################################	PCA	######################################

plink --bfile $vcf/gentoo_no_outliers_maf05_pruned \
--pca --out pca_maf05 --set-missing-var-ids @:# --double-id --allow-extra-chr


##################################	ADMIXTURE	##################################

plink --bfile $vcf/gentoo_no_outliers_maf05_pruned \
--recode 12 --out $vcf/gentoo_no_outliers_maf05_pruned_admix \
--allow-extra-chr --set-missing-var-ids @:# --double-id


################################## PAIRWISE FST ##################################
#!/bin/bash
# pop list
localities=("CROZ.txt" "MARI.txt" "KERG.txt" "FALK.txt" "MRTI.txt" "SIGN.txt" "OHIG.txt" "KING.txt" "GGV.txt")

# Iterate over all combinations of location pairs
for ((i=0; i<${#localidades[@]}; i++)); do
  for ((j=i+1; j<${#localidades[@]}; j++)); do
    loc1=${localidades[i]}
    loc2=${localidades[j]}
    echo "Calculating FST between $loc1 y $loc2"
    
    vcftools --gzvcf $vcf/gentoo_no_outliers_maf05_pruned.vcf.gz \
      --weir-fst-pop $loc1 \
      --weir-fst-pop $loc2 \
      --out ${loc1%.txt}_vs_${loc2%.txt}
  done
done
