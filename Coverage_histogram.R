
Total_coverage <- read.csv("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/MedianCovereage_wellderly/total_coverage.csv", header=FALSE)
cov <- Total_coverage$V2
Coverage <- cov[cov < 200]

tiff("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/MedianCovereage_wellderly/cov.tiff", width=2000, height=2000, res=300, compression="lzw")
hist(Coverage)
dev.off()
q()

whiteOnly_wg_and_exome <- read.delim("~/Desktop/Coverage/CG_coverage/whiteOnly_wg_and_exome.txt", header=FALSE)
hist(whiteOnly_wg_and_exome$V3, main = paste("Whole genome coverage"), xlab="Coverage")
hist(whiteOnly_wg_and_exome$V4, main = paste("Exome coverage"), xlab="Coverage")

#Whole Genome
coverage <- read.table("~/Desktop/Coverage/CG_coverage/total_cov_wg_and_exome.txt", header=FALSE)
cov <- coverage$V3
Coverage <- cov[cov < 200]

p <- ggplot() + 
        geom_histogram(aes(Coverage, y=..density..), color="black", fill="white", binwidth=5) +
        scale_x_continuous(expand=c(0,0), limits=c(0,125), breaks=(c(0,20,40,60,80,100,120))) +
        scale_y_continuous(expand=c(0,0)) +
        ylab("") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.ticks=element_line(color="black"), 
              axis.text=element_text(color="black"))

tiff("/Users/gerikson/Desktop/Coverage/CG_coverage/WG_coverage_wellderly.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()


#EXOME Coverage
coverage <- read.table("~/Desktop/Coverage/CG_coverage/total_cov_wg_and_exome.txt", header=FALSE)
cov <- coverage$V4
Coverage <- cov[cov < 200]

p <- ggplot() + 
        geom_histogram(aes(Coverage, y=..density..), color="black", fill="white", binwidth=0.001) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        ylab("") +
        xlab("Exome fraction where weightSumSequenceCoverage >= 10x") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.ticks=element_line(color="black"), 
              axis.text=element_text(color="black"))

tiff("/Users/gerikson/Desktop/Coverage/CG_coverage/Exome_coverage_wellderly.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()

#Inova
#Whole Genome
coverage <- read.csv("~/Desktop/Coverage/CG_coverage/Inova_coverage.txt", header=FALSE)
cov <- coverage$V2

p <- ggplot() + 
        geom_histogram(aes(cov, y=..density..), color="black", fill="white", binwidth=5) +
        scale_x_continuous(expand=c(0,0), limits=c(0,125), breaks=(c(0,20,40,60,80,100,120))) +
        scale_y_continuous(expand=c(0,0)) +
        ylab("") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.ticks=element_line(color="black"), 
              axis.text=element_text(color="black"))

tiff("/Users/gerikson/Desktop/Coverage/CG_coverage/WG_coverage_Inova.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()

#EXOME Coverage
coverage <- read.csv("~/Desktop/Coverage/CG_coverage/Inova_coverage.txt", header=FALSE)
cov <- coverage$V3
Coverage <- cov[cov < 200]

p <- ggplot() + 
        geom_histogram(aes(Coverage, y=..density..), color="black", fill="white", binwidth=0.001) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        ylab("") +
        xlab("Exome fraction where weightSumSequenceCoverage >= 10x") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.ticks=element_line(color="black"), 
              axis.text=element_text(color="black"))

tiff("/Users/gerikson/Desktop/Coverage/CG_coverage/Exome_coverage_Inova.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()
