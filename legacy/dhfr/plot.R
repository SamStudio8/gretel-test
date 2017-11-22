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

d <- d[d$readsize != 75, ]
d <- d[d$readsize != 90, ]
d <- d[d$readsize != 125, ]
d <- d[d$readsize != 250, ]
d <- d[d$readsize != 100, ]
d <- d[d$readsize != 500, ]

d <- d[d$cov != 500,]
d <- d[d$cov != 250,]
d <- d[d$cov != 1,]

read_l_labels <- c("50"="50 bp", "100"="100 bp", "125"="125 bp", "150"="150 bp")
d$in_hname_o = factor(d$in_hname, levels=c('BC070280', 'XR_634888', 'AK232978', 'M19237', 'XM_014960529'))

ho_labels <- c(BC070280 = "BC070280 (99.8%)", XR_634888 = "XR_634888 (97.3%)", AK232978="AK232978 (90.1%)", M19237="M19237 (83.5%)", XM_014960529="XM_014960529 (78.7%)")
ho_s_labels <- c(BC070280 = "BC", XR_634888 = "XR", AK232978="AK", M19237="M", XM_014960529="XM")

p <- ggplot(d, aes(x=factor(cov), fill=factor(cov), recovery)) + geom_boxplot(outlier.size=3) + facet_grid(readsize~in_hname_o, labeller=labeller(in_hname_o = ho_labels, readsize=read_l_labels)) + theme_grey(base_size=30)

p <- p + ylab("Correctly Recovered SNPs (% Proportion)  -  Facets: Synthetic Read Length")
p <- p + ggtitle("Gretel Recovery Rates: DHFR (Masked, n100, NS)         fbc8eada-a2a7-4ede-b1b5-757e3640fa93")
#p <- p + theme(axis.title=element_text(size=30)) + theme(plot.title=element_text(size=20, face="bold"))
p <- p + scale_y_continuous(expand=c(0,0), limits=c(-1,100), breaks = c(0,10,20,30,40,50,60,70,80,90,95,100), minor_breaks=NULL)

#p <- p + theme(text=element_text(size=20)) 
p <- p + theme(panel.grid.major.x=element_blank())
p <- p + scale_fill_discrete(guide = guide_legend(title = "Cov."))
#p <- p + theme(axis.title.x=element_blank())
p <- p + xlab("Per-Haplotype Read Depth   -   Facets: Input Haplotype (Reference Identity)")
p <- p + theme(panel.grid.major.x=element_blank())
p <- p + theme(panel.grid.major.y=element_line(color="gray"))
p <- p + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p <- p + theme(strip.text=element_text(size=30))
p <- p + theme(legend.position="none")

p3 <- ggplot(d, aes(x=factor(cov), fill=factor(cov), reads_propdrop)) + geom_boxplot() + facet_grid(readsize~in_hname_o, labeller=labeller(in_hname_o = ho_labels, readsize=read_l_labels)) + theme_grey(base_size=30)
p3 <- p3 + theme(panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
p3 <- p3 + scale_color_discrete(labels=ho_s_labels, guide = guide_legend(title = "Hap."))
p3 <- p3 + theme(panel.grid.major.x=element_blank())
p3 <- p3 + theme(panel.grid.major.y=element_line(color="gray"))
p3 <- p3 + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p3 <- p3 + theme(strip.text=element_text(size=30))
p3 <- p3 + theme(legend.position="none")
p3 <- p3 + scale_y_continuous(expand=c(0,0), limits=c(-1,100), breaks = c(0,10,20,30,40,50,60,70,80,90,95,100), minor_breaks=NULL)
p3 <- p3 + ylab("Reads Dropped by Alignment (% Proportion)  -  Facets: Synthetic Read Length")
p3 <- p3 + xlab("Per-Haplotype Read Depth   -   Facets: Input Haplotype (Reference Identity)")

d$snptot <- d$snptot/196*100
p2 <- ggplot(d, aes(x=factor(cov), fill=factor(cov), snptot)) + geom_boxplot(outlier.size=3) + facet_grid(readsize~., labeller=labeller(readsize=read_l_labels)) + theme_grey(base_size=30)
p2 <- p2 + theme(panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
p2 <- p2 + scale_color_discrete(labels=ho_s_labels, guide = guide_legend(title = "Hap."))
p2 <- p2 + theme(panel.grid.major.x=element_blank())
p2 <- p2 + theme(panel.grid.major.y=element_line(color="gray"))
p2 <- p2 + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p2 <- p2 + theme(strip.text=element_text(size=30))
p2 <- p2 + theme(legend.position="none")
p2 <- p2 + scale_y_continuous(expand=c(0,0), limits=c(-1,100), breaks = c(0,10,20,30,40,50,60,70,80,90,95,100), minor_breaks=NULL)
p2 <- p2 + ylab("SNPs found during Variant Calling (% Proportion)  -  Facets: Synthetic Read Length")
p2 <- p2 + xlab("Per-Haplotype Read Depth")


# <3 http://stackoverflow.com/questions/18427455/multiple-ggplots-of-different-sizes
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(8, 3)))
#print(p, vp = vplayout(1:8, 1:3))  # ROWS, COLS
#print(p2, vp = vplayout(1:8, 1))
print(p3, vp = vplayout(1:8, 1:3))

#ggsave("dhfr-main.png")
#ggsave("dhfr-snps.png")
ggsave("dhfr-dropped.png")
