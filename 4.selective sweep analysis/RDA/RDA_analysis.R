library(adegenet)
library(psych)     
library(vegan)     
library(dplyr)
library(sp)         
library(Cairo)  

gen <- read.PLINK('gentoo_outliers.raw', header=T)

summary(gen)
dim(gen)
gen$ind.names
# ver si hay missing en genind
sum(is.na(gen@other$genind))


gen_matriz <- tab(gen)
#gen_matriz
sum(is.na(gen_matriz))


gen.imp <- apply(gen_matriz, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs


env <- read.csv("marinas_gentoo.txt", header = TRUE, sep = "\t")
str(env) # Look at the structure of the data frame
env

identical(rownames(gen.imp), env[,1]) 

pairs.panels(env[,5:11], scale=F)
cor(env[,5:11])

pred <- subset(env, select=-c(o14_tempmi))
pred <- pred[,5:10]
pred

pairs.panels(pred, scale=F)
cor(pred)

gentoo.rda <- rda(gen.imp ~ ., data = pred, scale = TRUE)

gentoo.rda
RsquareAdj(gentoo.rda)
summary(eigenvals(gentoo.rda, model = "constrained"))
screeplot(gentoo.rda)

signif.full <- anova.cca(gentoo.rda, parallel = getOption("mc.cores")) # permutation=999 por defecto
signif.full

signif.axis <- anova.cca(gentoo.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis

vif.cca(gentoo.rda)


compression_quality <- 2

CairoPDF("RDA_marinas_segovia_sincoordenadas.pdf", 
         width = 6, height = 6, pointsize = 8, 
         compression = compression_quality)

par(plt=c(0.1, 0.9, 0.1, 0.9), 
    mar=c(4, 4, 2, 1))  # Ajusta plt y mar según tus preferencias

# Plot de la RDA parcial
plot(gentoo.rda, type="n", scaling=3)


# Puntos para SNPs
points(gentoo.rda, 
       display="species", 
       pch=20, cex=0.5, col="gray32", 
       scaling=3)           # the SNPs


levels(env$sites) <- c("aKerguelen","bCrozet","cMarion",
                       "eSigny","fOhiggins","gSPoint","hGGV",
                       "iFalkland","jMartillo")
eco <- env$sites
bg <- c("#33a02c","#1f78b4","#a6cee3","#E78E4A",
        "#ffff33","#EDE185","#AAAE01","#e31a1c",
        "#FF33F9" ) # 6 nice colors for our ecotypes


# Puntos para sitios (con leyenda)
points(gentoo.rda, display="sites", 
       pch=21, cex=3, col=NA, 
       scaling=3, bg=bg[as.factor(env$sites)]) # the wolves


# Variables predictoras
text(gentoo.rda, scaling=3, 
     display="bp", col="#0868ac", cex=1) 

#legend("bottomright", 
#       legend=levels(eco), bty="n", 
#       col="gray32", pch=21, cex=0.8, pt.bg=bg)

dev.off()
