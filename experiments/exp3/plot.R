require(ggplot2)
require(reshape)
require(tikzDevice)
tikz('exp3.tex', width=5, height=5)

data = read.table("results.csv", sep=";", header=TRUE)

data <- melt(data, id=colnames(data)[c(-5,-6,-7,-8)])
data$value = data$value / data$inputsize * 8

# reorder from (grid, bv, doc, csa) components to (csa, bv, grid, doc)
data$variable <- factor(data$variable, levels(data$variable)[c(2,4,3,1)], ordered=T)

p <- ggplot(data, aes(index, value, fill=variable)) + geom_bar(position="stack",stat="identity")
p <- p + theme(legend.position="top")
p <- p + facet_wrap(~collection)
p <- p + scale_x_discrete(name="")
p <- p + scale_y_continuous(name="Space [Bits per input character]", limits = c(0,32), breaks=c(0,8,16,24,32))
plot(p)
dev.off()
