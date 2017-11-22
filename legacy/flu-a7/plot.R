library(ggplot2)
library(reshape2)
library(scales)
library(grid)

d = read.table(file=TABLE_P,
               sep="\t",
               comment.char="",
               na.strings=c("NA","-"),
               header=T
)

d <- d[d$readsize != 100, ]
d <- d[d$readsize != 250, ]


divr_labels <- c("0.001"="0.005 SNPs/base", "0.005"="0.025 SNPs/base", "0.01"="0.05 SNPs/base", "0.015"="0.075 SNPs/base", "0.02"="0.1 SNPs/base", "0.05"="0.25 SNPs/base", "0.1" = "0.5 SNPs/base")
read_l_labels <- c("50"="50 bp", "100"="100 bp", "150"="150 bp", "250"="250 bp", "500"="500 bp")

p <- ggplot(d, aes(x=factor(cov), fill=factor(cov), recovery)) + geom_boxplot() + facet_grid(readsize~., labeller=labeller(snpn = divr_labels, readsize=read_l_labels)) + theme_grey(base_size=30)

p <- p + ylab("Correctly Recovered SNPs (% Proportion)  -  Facets: Synthetic Read Length")
#p <- p + ggtitle("Gretel Recovery Rates: Flu-A7 (24, n=100)         09d1e503-3e17-455c-bd47-dbaea3f96b3e")

p <- p + theme(axis.title=element_text(size=30)) + theme(plot.title=element_text(size=20, face="bold"))
p <- p + scale_y_continuous(expand=c(0,0), limits=c(-1,100), breaks = c(0,10,20,30,40,50,60,70,80,90,95,100), minor_breaks=NULL)

#p <- p + theme(text=element_text(size=20)) + theme(axis.text.x=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
p <- p + scale_fill_discrete(guide = guide_legend(title = "Cov."))
#p <- p + theme(axis.title.x=element_blank())
p <- p + xlab("Per-Haplotype Read Depth")
p <- p + theme(panel.grid.major.x=element_blank())
p <- p + theme(panel.grid.major.y=element_line(color="gray"))
p <- p + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p <- p + theme(strip.text=element_text(size=30))
p <- p + theme(legend.position="none")


d$snptot <- d$snptot/277*100
p2 <- ggplot(d, aes(x=factor(cov), fill=factor(cov), snptot)) + geom_boxplot(outlier.size=3) + facet_grid(readsize~., labeller=labeller(readsize=read_l_labels)) + theme_grey(base_size=30)
p2 <- p2 + theme(panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
p2 <- p2 + scale_color_discrete(guide = guide_legend(title = "Hap."))
p2 <- p2 + theme(panel.grid.major.x=element_blank())
p2 <- p2 + theme(panel.grid.major.y=element_line(color="gray"))
p2 <- p2 + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p2 <- p2 + theme(strip.text=element_text(size=30))
p2 <- p2 + theme(legend.position="none")
p2 <- p2 + scale_y_continuous(expand=c(0,0), limits=c(-1,100), breaks = c(0,10,20,30,40,50,60,70,80,90,95,100), minor_breaks=NULL)
p2 <- p2 + ylab("SNPs found during Variant Calling (% Proportion)  -  Facets: Synthetic Read Length")
p2 <- p2 + xlab("Per-Haplotype Read Depth")




#p3 <- ggplot(d, aes(x=factor(cov), col=factor(cov), reads_propdrop)) + geom_point() + facet_grid(~readsize) + theme(strip.text= element_text(size=30)) + theme(axis.title=element_text(size=30)) + theme(axis.text.y=element_text(size=20))
#p3 <- p3 + theme(text=element_text(size=20)) + theme(panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
#p3 <- p3 + scale_color_discrete(guide = guide_legend(title = "Hap."))
#p3 <- p3 + theme(strip.text.x = element_blank())
#p3 <- p3 + xlab("Read Coverage   -   Facets: Synthetic Read Length")
#p3 <- p3 + ylab("%Dropped")
p3 <- ggplot(d, aes(x=factor(cov), fill=factor(cov), reads_propdrop)) + geom_boxplot() + facet_grid(readsize~in_hname, labeller=labeller(readsize=read_l_labels)) + theme_grey(base_size=30)
p3 <- p3 + theme(panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
p3 <- p3 + scale_color_discrete(guide = guide_legend(title = "Hap."))
p3 <- p3 + theme(panel.grid.major.x=element_blank())
p3 <- p3 + theme(panel.grid.major.y=element_line(color="gray"))
p3 <- p3 + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p3 <- p3 + theme(strip.text=element_text(size=30))
p3 <- p3 + theme(legend.position="none")
p3 <- p3 + scale_y_continuous(expand=c(0,0), limits=c(-1,100), breaks = c(0,10,20,30,40,50,60,70,80,90,95,100), minor_breaks=NULL)
p3 <- p3 + ylab("Reads Dropped by Alignment (% Proportion)  -  Facets: Synthetic Read Length")
p3 <- p3 + xlab("Per-Haplotype Read Depth   -   Facets: Input Haplotype (Reference Identity)")



# <3 http://stackoverflow.com/questions/18427455/multiple-ggplots-of-different-sizes
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(8, 3)))
print(p, vp = vplayout(1:8, 1))  # ROWS, COLS
#print(p2, vp = vplayout(1:8, 1))
#print(p3, vp = vplayout(1:8, 1:3))

ggsave("flu-main.png")
