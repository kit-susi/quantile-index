library(ggplot2)
library(grid)
library(scales)
library(tikzDevice)
library(stringr)
library(plyr)
library(dplyr)

#if (!exists( "tikzDeviceLoaded")) {
  #require(tikzDevice)
  #options(tikzLatexPackages =
    #c(getOption('tikzLatexPackages'), 
      #paste("\\input{",getwd(),"/macros.tex}",sep=""))
  #)
  #tikzDeviceLoaded = T
#}

#d <- read.csv2("fig-var-txt-size.csv",sep=";",header=TRUE)
d <- read.csv2("results_64_only.csv",sep=";",header=TRUE)
d$timenum <- as.numeric(as.character(d$time))
d$intervalsznum <- as.numeric(as.character(d$intervalsz))

roundit <- function(intervalsznum) { floor(log(intervalsznum)/log(10)) }
d <- cbind(d, x = mapply(roundit, d$intervalsznum))

#d$dataset <- revalue(d$dataset, c("KERNEL"="\\kernel", 
                                  #"CC"="\\commoncrawl",
                                  #"DNA"="\\dna"))
#d$method  <- revalue(d$method, c("QGRAM-RGXP"="\\qgramrgxp~~~",
                                 #"REGEXP"="\\regexp~~~~",
                                 #"WTDFS"="\\wtdfs~~~",
                                 #"SASCAN"="\\sascan~~~"))

#tikz("fig-var-txt-size.tex",width = 6.0, height = 2.8)

plot <- ggplot(d,aes(factor(x),timenum,fill=algo,color=algo))
#plot <- plot + geom_boxplot(outlier.size = 1)
plot <- plot + geom_boxplot(outlier.shape = NA)
#plot <- plot + geom_point()
plot <- plot + facet_grid(. ~ instance)
plot <- plot + scale_x_discrete(name = "SA Interval size")
plot <- plot + scale_y_log10(name = "Query time [µs]")
                             #limits = quantile(d$timenum, c(0.1,0.9)))
print(plot)

#dev.off()
