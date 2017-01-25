require(ggplot2)
require(tikzDevice)
tikz('exp2.tex', standAlone = TRUE, width=5, height=5)

data = read.table("results.csv", sep=";", header=TRUE)
collections <- unique(data$collection)

data$bitsperinput = data$indexsize / data$inputsize
for (col in collections) {
    coldata = data[data$collection==col,]
    p <- ggplot(coldata, aes(x=bitsperinput,y=avg, color=index))
    p <- p + scale_colour_hue(l=50)
    p <- p + geom_point(shape=1)
    p <- p + ggtitle(col) 
    p <- p + scale_x_continuous(limits = c(0,4))
    plot(p)
}

#plot(pp)
dev.off()
