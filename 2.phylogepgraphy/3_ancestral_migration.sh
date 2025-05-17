### Create plink file including 62 gentoo penguin samples and one outgroup (Chinstrap penguin)
plink --vcf $vcf_gentoo_chin/gentoo62_chin_nomiss.vcf.gz --recode --out vcf_gentoo_chin --allow-no-sex --allow-extra-chr

### create bfile
plink --file vcf_gentoo_chin --make-bed --out vcf_gentoo_chin --allow-no-sex --allow-extra-chr

### calculate allele frequencies within sampling localities
plink --bfile vcf_gentoo_chin --freq --missing --within vcf_prun_sinmq_chins.clust --out vcf_gentoo_chin --allow-no-sex --allow-extra-chr

### compress the frequeny file
gzip vcf_gentoo_chin.frq.strat

### convert plink frequency file to treemix input file (https://github.com/thomnelson/tools/blob/master/plink2treemix.py)
python plink2treemix.py vcf_gentoo_chin.frq.strat.gz vcf_gentoo_chin.treemix.frq.gz

### Run treemix
for m in {1..5}
do
  for i in {1..10}
do
# Generate random seed
s=$RANDOM
echo "Random seed = ${s}"
  echo "Ejecutando TreeMix con m=${m}"
  treemix \
    -i vcf_gentoo_chin.treemix.frq.gz \
    -o run_bootstrap/vcf_gentoo_chin.${i}.${m} \
    -clust vcf_prun_sinmq.clust \
    -m ${m} \
    -global \
    -bootstrap \
    -k 500 \
    -root Chinstrap
    -seed ${s}
  done
done


### Treemix plotting was performing using OptM (Fitak et al 2021)



############################################################################

#!/bin/bash

### vcf without Macquarie
vcftools --gzvcf $vcf_gentoo_chin/gentoo62_chin_nomiss.vcf.gz --keep keep_dsuite.txt --recode --out dsuite_tree
### keep_dsuite.txt file: one sample for each locality
#chins, 13A2, CRO08, GEN04, GSA10, G4,P0167, P0402, P1884, PAM10


vcftools --gzvcf $vcf_gentoo_chin/gentoo62_chin_nomiss.vcf.gz --keep keep_dsuite.txt --recode --stdout |bcftools view -c 1 -Oz -o dsuite_tree.vcf.gz
plink --vcf dsuite_tree.vcf.gz --indep-pairwise 50 10 0.1 --set-missing-var-ids @:# --double-id --allow-extra-chr --out dsuite_tree_prun
plink --vcf dsuite_tree.vcf.gz --extract dsuite_tree_prun.prune.in --set-missing-var-ids @:# --double-id --allow-extra-chr --recode vcf --out dsuite_tree_pruned


### https://github.com/edgardomortiz/vcf2phylip.git
python $tovcf/vcf2phylip.py -i dsuite_tree_pruned.vcf -o chins
iqtree -s dsuite_tree.recode.min4.phy -m GTR+G -st DNA -bb 1000 -nt AUTO


# RUN DSUITE
dsuite=/Users/dalynoll/Documents/Programas/genetica/Dsuite/Build
$dsuite/./Dsuite Dtrios --KS-test-for-homoplasy $vcf_gentoo_chin/gentoo62_chin_nomiss.vcf.gz SETS_pop.txt -t dsuite_tree_n10.nwk

#!/bin/bash
$dsuite/Dsuite Fbranch dsuite_tree_n10.nwk  SETS_pop_tree.txt

python $dtools/dtools.py -n pops_gentoo --outgroup chins --use_distances --ladderize fbranch.txt dsuite_tree_n10.nwk

