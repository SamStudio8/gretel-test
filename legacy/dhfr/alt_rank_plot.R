library("ggplot2");
#d <- read.table(TABLE_P, header=T);
#d <- read.table("fbc_hamming_wpd_wmeta.txt", header=T);
d <- read.table("fbc_hamming_wpd_wmeta.withbestworstlikl.withabslikl.txt", header=T);


#d <- d[d$cov == 10,]

#p <- ggplot(d, aes(recovery, rank_w, colour=factor(in_hname))) + facet_grid(readsize~cov) + scale_y_reverse( lim=c(75,0)) + scale_x_reverse( lim=c(100,0)) + theme_grey(base_size=30) + geom_jitter(size=2)

d <- d[d$n_outhaps > 0, ]
d <- d[d$rank_w > -1, ]

d <- d[d$worstlikl < 112.14, ]

read_l_labels <- c("50"="50 bp", "100"="100 bp", "125"="125 bp", "150"="150 bp")
d$in_hname_o = factor(d$in_hname, levels=c('BC070280', 'XR_634888', 'AK232978', 'M19237', 'XM_014960529'))

ho_labels <- c(BC070280 = "BC070280 (99.8%)", XR_634888 = "XR_634888 (97.3%)", AK232978="AK232978 (90.1%)", M19237="M19237 (83.5%)", XM_014960529="XM_014960529 (78.7%)")

d$rank_w <- (d$rank_w/(d$n_outhaps-1))

d$liklscale <- (d$abslik - min(d$bestlikll)) / (max(d$worstlikl) - min(d$bestlikll))
#d$liklscale <- (d$abslik - (d$bestlikll)) / ((d$worstlikl) - (d$bestlikll))
#d$liklscale <- (d$abslik - d$bestlikll) / (d$worstlikl - d$bestlikll)
#d$liklscale <- (d$abslik - 0) / (d$worstlikl - 0)

p <- ggplot(d, aes(recovery, rank_w, colour=factor(cov))) + facet_grid(readsize~in_hname_o, labeller=labeller(in_hname_o = ho_labels, readsize=read_l_labels)) + scale_y_reverse(lim=c(1.01,-0.01)) + scale_x_reverse( lim=c(101,39)) + theme_grey(base_size=30) + geom_point(size=2.5, alpha=0.8)
#p <- p + theme(panel.grid.major.x=element_line(color="gray"))
#p <- p + theme(panel.grid.major.y=element_line(color="gray"))
p <- p + theme(legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines'))
p <- p + theme(strip.text=element_text(size=30))
p <- p + theme(legend.position="bottom")
p <- p + xlab("Correctly Recovered SNPs (% Proportion)   -   Facets: Input Haplotype (Reference Identity)")
p <- p + ylab("Haplotype Scaled Likelihood  -  Facets: Synthetic Read Length")
p <- p + guides(colour=guide_legend(title="Depth",nrow=1))

ggsave("dhfr-ranks.png")
