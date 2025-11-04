#install.packages(c("vcfR", "adegenet", "poppr", "pegas", "hierfstat"))
library(vcfR)
library(adegenet)
library(poppr)
library(pegas)

vcf <- read.vcfR("/Users/dalynoll/Dropbox/Investigaciones/gentoo/PAPERS_TESIS_DOCTORAL/Adaptive_divergence/reanalisis/neutral/polarizado/input_treemix_samba/gentoo64_gt_nosexnomiss_polarizado_maf001_pruned_sinMQ_nooutlier.vcf.gz")

gl <- vcfR2genlight(vcf)
meta <- read.table("gentoo_popmap.txt", header = FALSE, col.names = c("ID", "Lineage", "Locality"))

indNames(gl) <- meta$ID

strata(gl) <- data.frame(Lineage = meta$Lineage, Locality = meta$Locality)
setPop(gl) <- ~Lineage/Locality  # estructura jerárquica: Lineage > Locality

amova_result <- poppr.amova(gl, ~Lineage/Locality)

print(amova_result)

set.seed(42)
amova_test <- randtest(amova_result, nrepet = 999)
plot(amova_test)
