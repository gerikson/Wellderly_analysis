### plot PCs in R ###
library(ggplot2)
dat <- read.table("test-PCA.eigenvec", F)
# first 519 are wellderly, remaining are inova #
COLOR <- c(rep("red", times=519), rep("blue", times=dim(dat)[1]-519))
p <- ggplot() +
    geom_point(aes(dat$V3, dat$V4), col=COLOR) +
    xlab("pc1") +
    ylab("pc2")
tiff("pca.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()
q()
'''

'''
##### ASSOCIATION TESTING #####

# create phenotype file in R #
# wellderly = "affected" #
dat <- read.table("test-MAF_LD_pruned.fam", F)
pheno <- c(rep(2, times=519), rep(1, times=dim(dat)[1]-519))
out <- data.frame(dat$V1, dat$V2, pheno)
write.table(out, "test.pheno", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
q()