### plot PCs in R ###
library(ggplot2)
dat <- read.table("test-PCA.eigenvec", F)
# first 519 are wellderly, remaining are inova #
COLOR <- c(rep("blue", times=dim(dat)[1]-511), rep("red", times=511))
p <- ggplot() +
    geom_point(aes(dat$V3, dat$V4), col=COLOR, size = 0.8) +
    xlab("Principal Component 1") +
    ylab("Principal Component 2") + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))
tiff("Figure 2A - PCA.tiff", width=3, height=3, res=300, units = 'in', compression="lzw")
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