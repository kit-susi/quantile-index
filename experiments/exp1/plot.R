require(ggplot2)
require(tikzDevice)
tikz('exp1.tex', standAlone = TRUE, width=5, height=5)

data = read.table("results.csv", sep=";", header=TRUE)
collections <- unique(data$collection)
for (col in collections) {
    max = data[data$collection==col & data$q == 1,]$G
    data[data$collection==col,]$G <- 100 * data[data$collection==col,]$G / max
}
yticks = c(1,2,5,10,25,50,100)
pp <- ggplot(data, aes(q, G,color=collection))
pp <- pp + geom_line() 
pp <- pp + geom_point()
pp <- pp + scale_y_log10()
#pp <- pp + scale_y_log10(breaks=yticks)
plot(pp)
dev.off()
