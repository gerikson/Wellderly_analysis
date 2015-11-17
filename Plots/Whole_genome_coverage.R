library("gridExtra")
library("cowplot")
#Whole Genome
coverage <- read.table("~/Desktop/Coverage/CG_coverage/total_cov_wg_and_exome.txt", header=FALSE)
cov <- coverage$V3
Coverage3 <- cov[cov < 200]

p_wg <- ggplot() + 
        geom_histogram(aes(Coverage3, y=..density..), color="black", fill="white", binwidth=5) +
        scale_x_continuous(expand=c(0,0), limits=c(30,110), breaks=(c(30,50,70,90,110))) +
        scale_y_continuous(expand=c(0,0), limits=c(0,0.07), breaks=c(0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07)) +
        ylab("Count") +
        xlab("Wellderly WGS Coverage") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.ticks=element_line(color="black"), 
              axis.text=element_text(color="black"), text = element_text(size=8), axis.text.x=element_text(size = 8), axis.text.y=element_text(size = 8))


#Inova
#Whole Genome
coverage <- read.csv("~/Desktop/Coverage/CG_coverage/Inova_coverage.txt", header=FALSE)
cov4 <- coverage$V2

p_inova <- ggplot() + 
        geom_histogram(aes(cov4, y=..density..), color="black", fill="white", binwidth=5) +
        scale_x_continuous(expand=c(0,0), limits=c(30,110), breaks=(c(30,50,70,90,110))) +
        scale_y_continuous(expand=c(0,0), limits=c(0,0.07), breaks=c(0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07)) +
        ylab("Count") +
        xlab("ITMI WGS Coverage") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.ticks=element_line(color="black"), 
              axis.text=element_text(color="black"), text = element_text(size=8), axis.text.x=element_text(size = 8), axis.text.y=element_text(size = 8))


tiff("/Users/gerikson/Desktop/FinalPlots/Figure S1A – WellderlyWGS.Coverage.tiff",width= 4, height = 4, units = 'in', res=300)
par(mar=c(4,4,1,1)+0.1)
plot_grid(p_wg, p_inova, labels=c("A", "B"))
dev.off()

