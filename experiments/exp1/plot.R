require(ggplot2)
require(tikzDevice)
tikz('exp1.tex', width=5, height=3)

data = read.table("results.csv", sep=";", header=TRUE)
collections <- unique(data$collection)
for (col in collections) {
    max = data[data$collection==col & data$q == 1,]$Gq
    data[data$collection==col,]$Gq <- 100 * data[data$collection==col,]$Gq / max
}
pp <- ggplot(data, aes(q, Gq,color=collection))
pp <- pp + geom_line() 
pp <- pp + geom_point()
pp <- pp + theme(legend.position="top")
#pp <- pp + scale_y_continuous()
pp <- pp + scale_y_log10(name="$|G_q| / |G_1|$ [\\%]", breaks=c(1,2,4,8,16,32,64,100))
pp <- pp + scale_x_log10(name="$q$", breaks=c(1,2,4,8,16,32,64,128,256))
pp <- pp + expand_limits(y=1)
plot(pp)
dev.off()
