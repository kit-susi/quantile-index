require(ggplot2)
require(tikzDevice)
tikz('exp2.tex', width=5, height=5)

data = read.table("results.csv", sep=";", header=TRUE)
collections <- unique(data$collection)

data$bitsperinput = data$indexsize / data$inputsize * 8
d1 <- subset(data,  data$q == 64 & data$s == 16)
d1$index <- "s=16, q=64"

data <- rbind(data, d1)

p <- ggplot(data, aes(x=bitsperinput,y=avg, color=index, shape=index))
p <- p + theme(legend.position="top")
p <- p + facet_wrap(~collection)
p <- p + scale_colour_hue(l=50)
p <- p + geom_point(size=2)
p <- p + scale_y_log10(name="Average query time $[\\mu s]$", breaks=c(50,100, 200, 400, 800, 1600, 3200))
p <- p + scale_x_continuous(name="Space [Bits per input character]", limits = c(0,32), breaks=c(0,8,16,24,32))
plot(p)
dev.off()
