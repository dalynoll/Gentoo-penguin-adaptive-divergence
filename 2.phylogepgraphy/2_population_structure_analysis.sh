




####### Pairwise Fst ####

#!/bin/bash
# pop list
localities=("CROZ.txt" "MARI.txt" "KERG.txt" "FALK.txt" "MRTI.txt" "SIGN.txt" "OHIG.txt" "KING.txt" "GGV.txt")

# Iterate over all combinations of location pairs
for ((i=0; i<${#localidades[@]}; i++)); do
  for ((j=i+1; j<${#localidades[@]}; j++)); do
    loc1=${localidades[i]}
    loc2=${localidades[j]}
    echo "Calculating FST between $loc1 y $loc2"
    
    vcftools --gzvcf $vcf \
      --weir-fst-pop $loc1 \
      --weir-fst-pop $loc2 \
      --out ${loc1%.txt}_vs_${loc2%.txt}
  done
done
