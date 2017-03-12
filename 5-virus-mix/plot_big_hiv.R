library(ggplot2)
library(reshape2)
library(scales)
library(grid)



d = read.table(file="HIVG-303/hiv.c.ss.tab",
               sep="\t",
               comment.char="",
               na.strings=c("NA","-"),
               header=T
)

d$bitscore_norm  <- d$bitscore
d[d$gene == "pol",]$bitscore_norm <- (d[d$gene == "pol",]$bitscore - min(d[d$gene == "pol",]$bitscore)) / (max(d[d$gene == "pol",]$bitscore) - min(d[d$gene == "pol",]$bitscore))
d[d$gene == "nef",]$bitscore_norm <- (d[d$gene == "nef",]$bitscore - min(d[d$gene == "nef",]$bitscore)) / (max(d[d$gene == "nef",]$bitscore) - min(d[d$gene == "nef",]$bitscore))
d[d$gene == "vif",]$bitscore_norm <- (d[d$gene == "vif",]$bitscore - min(d[d$gene == "vif",]$bitscore)) / (max(d[d$gene == "vif",]$bitscore) - min(d[d$gene == "vif",]$bitscore))
d[d$gene == "env",]$bitscore_norm <- (d[d$gene == "env",]$bitscore - min(d[d$gene == "env",]$bitscore)) / (max(d[d$gene == "env",]$bitscore) - min(d[d$gene == "env",]$bitscore))
d[d$gene == "gag",]$bitscore_norm <- (d[d$gene == "gag",]$bitscore - min(d[d$gene == "gag",]$bitscore)) / (max(d[d$gene == "gag",]$bitscore) - min(d[d$gene == "gag",]$bitscore))

d$llik_norm  <- d$llik
d[d$gene == "pol",]$llik_norm <- (d[d$gene == "pol",]$llik - min(d[d$gene == "pol",]$llik)) / (max(d[d$gene == "pol",]$llik) - min(d[d$gene == "pol",]$llik))
d[d$gene == "gag",]$llik_norm <- (d[d$gene == "gag",]$llik - min(d[d$gene == "gag",]$llik)) / (max(d[d$gene == "gag",]$llik) - min(d[d$gene == "gag",]$llik))
d[d$gene == "vif",]$llik_norm <- (d[d$gene == "vif",]$llik - min(d[d$gene == "vif",]$llik)) / (max(d[d$gene == "vif",]$llik) - min(d[d$gene == "vif",]$llik))
d[d$gene == "nef",]$llik_norm <- (d[d$gene == "nef",]$llik - min(d[d$gene == "nef",]$llik)) / (max(d[d$gene == "nef",]$llik) - min(d[d$gene == "nef",]$llik))
d[d$gene == "env",]$llik_norm <- (d[d$gene == "env",]$llik - min(d[d$gene == "env",]$llik)) / (max(d[d$gene == "env",]$llik) - min(d[d$gene == "env",]$llik))

p <- ggplot(d, aes(x=bitscore_norm, y=llik_norm, color=strain)) + geom_point(size=8) + facet_grid(strain~gene) + theme_grey(base_size=32)
p <- p + theme(panel.grid.major.y=element_line(color="gray"))
p <- p + theme(panel.grid.major.x=element_line(color="gray"))
p <- p + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p <- p + theme(strip.text=element_text(size=35))
p <- p + theme(legend.position="none")
p <- p + scale_y_continuous( breaks = c(0.0,0.2,0.4,0.6,0.8,1.0), minor_breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
p <- p + scale_x_reverse(expand=c(0,0), limits=c(1.025,-0.025), breaks = c(0,0.2,0.4,0.6,0.8,1.0), minor_breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())

p <- p + ylab("Normalised Gretel Likelihood (Best Decreasing)")
p <- p + xlab("Normalised BLAST Bitscore (Best Decreasing)")


p2 <- ggplot(d, aes(x=bitscore_norm, y=llik_norm, color=strain)) + geom_point(size=8) + theme_grey(base_size=32)
p2 <- p2 + theme(panel.grid.major.y=element_line(color="gray"))
p2 <- p2 + theme(panel.grid.major.x=element_line(color="gray"))
p2 <- p2 + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p2 <- p2 + theme(strip.text=element_text(size=35))
p2 <- p2 + theme(legend.position="none")
p2 <- p2 + scale_y_continuous( breaks = c(0.0,0.2,0.4,0.6,0.8,1.0), minor_breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
p2 <- p2 + scale_x_reverse(expand=c(0,0), limits=c(1.025,-0.025), breaks = c(0,0.2,0.4,0.6,0.8,1.0), minor_breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
p2 <- p2 + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p2 <- p2 + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())

p2 <- p2 + ylab("Normalised Gretel Likelihood (Best Decreasing)")
p2 <- p2 + xlab("Normalised BLAST Bitscore (Best Decreasing)")


# <3 http://stackoverflow.com/questions/18427455/multiple-ggplots-of-different-sizes
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p2, vp = vplayout(1, 1))  # ROWS, COLS
print(p, vp = vplayout(1, 2))

