library("gridExtra")
library("cowplot")

#EXOME Coverage
coverage <- read.table("~/Desktop/Coverage/CG_coverage/total_cov_wg_and_exome.txt", header=FALSE)
cov <- coverage$V4
Coverage <- cov[cov < 200]

p <- ggplot() + 
        geom_histogram(aes(Coverage, y=..density..), color="black", fill="white", binwidth=0.001) +
        scale_x_continuous(expand=c(0,0), limits=c(0.965,1), breaks=(c(0.97,0.98,0.99,1.00))) +
        scale_y_continuous(expand=c(0,0), limits=c(0,160), breaks=c(0,40,80,120,160)) +
        ylab("Percentage") +
        xlab("Wellderly % Exome > 10X Coverage") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.ticks=element_line(color="black"), 
              axis.text=element_text(color="black"), text = element_text(size=8), axis.text.x=element_text(size = 8), axis.text.y=element_text(size = 8))


#EXOME Coverage
coverage2 <- read.csv("~/Desktop/Coverage/CG_coverage/Inova_coverage.txt", header=FALSE)
cov2 <- coverage2$V3
Coverage2 <- cov2[cov2 < 200]

inova <- ggplot() + 
        geom_histogram(aes(Coverage2, y=..density..), color="black", fill="white", binwidth=0.001) +
        scale_x_continuous(expand=c(0,0), limits=c(0.965,1), breaks=(c(0.97,0.98,0.99,1.00))) +
        scale_y_continuous(expand=c(0,0), limits=c(0,160), breaks=c(0,40,80,120,160)) +
        ylab("Percentage") +
        xlab("ITMI % Exome > 10X Coverage") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              text = element_text(size=8), axis.text.x=element_text(size = 8), axis.text.y=element_text(size = 8))

tiff("/Users/gerikson/Desktop/FinalPlots/Figure S1B – WellderlyExome.Coverage.tiff",width= 4, height = 4, units = 'in', res=300)
par(mar=c(4,4,1,1)+0.1)
plot_grid(p, inova, labels=c("C", "D"))
dev.off()



