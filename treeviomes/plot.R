library(ggplot2)
library(reshape2)
library(scales)
library(grid)


d = read.table(file="collate_wmuscle.txt",
               sep="\t",
               comment.char="",
               na.strings=c("NA","-"),
               header=T
)

d <- d[d$readsize != 75, ]
#d <- d[d$readsize != 250, ]


#d <- d[d$snpn != 0.001, ]
d <- d[d$snpn != 0.0025, ]
d <- d[d$snpn != 0.0075, ]
d <- d[d$snpn != 0.0125, ]
#d <- d[d$snpn != 0.015, ]
d <- d[d$snpn != 0.0175, ]
d <- d[d$snpn != 0.025, ]
d <- d[d$snpn != 0.25, ]

#d <- d[d$snpn != 0.1, ]

d <- d[d$cov != 500,]
d <- d[d$cov != 250,]
d <- d[d$cov != 75,]
d <- d[d$cov != 100,]


d <- d[d$readsize != 50, ]
#d <- d[d$readsize != 150, ]
d <- d[d$readsize != 500, ]


#divr_labels <- c("0.001"="0.005 SNPs/base", "0.005"="0.025 SNPs/base", "0.01"="0.05 SNPs/base", "0.015"="0.075 SNPs/base", "0.02"="0.1 SNPs/base", "0.05"="0.25 SNPs/base", "0.1" = "0.5 SNPs/base")
divr_labels <- c("0.001"="0.001 SNPs/hb", "0.005"="0.005 SNPs/hb", "0.01"="0.01 SNPs/hb", "0.015"="0.015 SNPs/hb", "0.02"="0.02 SNPs/hb", "0.05"="0.05 SNPs/hb", "0.1" = "0.1 SNPs/hb")

#read_l_labels <- c("50"="50 bp (1.67%)", "100"="100 bp (3.33%)", "125"="125 bp (4.17%)", "150"="150 bp (5%)", "250"="250 bp (8.33%)", "500"="500 bp (16.67%)")
read_l_labels <- c("50"="50 bp", "100"="100 bp", "125"="125 bp", "150"="150 bp", "250"="250 bp", "500"="500 bp")



#p <- ggplot(d, aes(x=factor(cov), fill=factor(cov), recovery)) + geom_boxplot(outlier.size=2) + facet_grid(readsize~snpn, labeller=labeller(snpn = divr_labels, readsize=read_l_labels)) + theme(strip.text= element_text(size=30)) + theme(axis.title=element_text(size=30)) + theme(axis.text.y=element_text(size=25), axis.text.x=element_text(size=25)) + theme(legend.text=element_text(size=25)) + theme_bw(base_size=30)
p <- ggplot(d, aes(x=factor(cov), fill=factor(cov), recovery)) + geom_boxplot(outlier.size=2) + facet_grid(readsize~snpn, labeller=labeller(snpn = divr_labels, readsize=read_l_labels)) + theme_grey(base_size=30)

p <- p + ylab("Correctly Recovered SNPs (% Proportion)  -  Facets: Synthetic Read Length")
#p <- p + ggtitle("Gretel Recovery Rates: Treeviomes (39K)         710dd289-621e-456c-9274-e080b4a5b5fd")
p <- p + theme(axis.title=element_text(size=30)) + theme(plot.title=element_text(size=20, face="bold"))
p <- p + scale_y_continuous(expand=c(0,0), limits=c(-1,100), breaks = c(0,10,20,30,40,50,60,70,80,90,95,100), minor_breaks=NULL)

##p <- p + theme(text=element_text(size=20)) + theme(axis.text.x=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
#p <- p + theme(text=element_text(size=20)) 
p <- p + scale_fill_discrete(guide = guide_legend(title = "Hap. Depth"))
#p <- p + theme(axis.title.x=element_blank())
p <- p + xlab("Average Per-Haplotype Read Depth   -   Facets: Simultated Haplotype Variation (#SNPs per haplotype base)")
p <- p + theme(panel.grid.major.x=element_blank())
p <- p + theme(panel.grid.major.y=element_line(color="gray"))
p <- p + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p <- p + theme(strip.text=element_text(size=30))
p <- p + theme(legend.position="none")









p2 <- ggplot(d, aes(x=factor(cov), group=factor(cov), col=factor(cov), snptot)) + geom_point() + facet_grid(hapn~snpn) + theme(strip.text= element_text(size=30)) + theme(axis.title=element_text(size=30)) + theme(axis.text.y=element_text(size=20))
p2 <- p2 + theme(text=element_text(size=20)) + theme(axis.text.x=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
p2 <- p2 + theme(axis.title.x=element_blank())
p2 <- p2 + theme(strip.text.x = element_blank())
p2 <- p2 + scale_color_discrete(guide = guide_legend(title = "Cov."))
#p2 <- p2 + geom_hline(yintercept=102)
p2 <- p2 + ylab("#SNPs")




##p3 <- ggplot(d, aes(x=factor(cov), col=factor(hapi), propdrop)) + geom_point() + facet_grid(hapn~snpn) + theme(strip.text= element_text(size=30)) + theme(axis.title=element_text(size=30)) + theme(axis.text.y=element_text(size=20))
##p3 <- p3 + theme(text=element_text(size=20)) + theme(panel.grid.major.x=element_blank(), panel.grid.major.x=element_blank(), axis.ticks.x=element_blank())
##p3 <- p3 + scale_color_discrete(guide = guide_legend(title = "Hap."))
##p3 <- p3 + theme(strip.text.x = element_blank())
##p3 <- p3 + xlab("Read Coverage   -   Facets: Synthetic Read Length")
##p3 <- p3 + ylab("%Dropped")



# <3 http://stackoverflow.com/questions/18427455/multiple-ggplots-of-different-sizes
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(7, 1)))
print(p, vp = vplayout(1:7, 1))  # ROWS, COLS
#print(p2, vp = vplayout(7, 1))

ggsave("tree-main.png")

library("plyr")
ddply(d, readsize~snpn, summarise, mean=mean(snptot), median=median(snptot))

