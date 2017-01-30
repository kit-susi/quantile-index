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
d <- read.csv2("results.csv",sep=";",header=TRUE)
d$timenum <- as.numeric(as.character(d$time))
d$knum <- as.numeric(as.character(d$k))
d$result_count_num <- as.numeric(as.character(d$result_count))

#d$dataset <- revalue(d$dataset, c("KERNEL"="\\kernel", 
                                  #"CC"="\\commoncrawl",
                                  #"DNA"="\\dna"))
#d$method  <- revalue(d$method, c("QGRAM-RGXP"="\\qgramrgxp~~~",
                                 #"REGEXP"="\\regexp~~~~",
                                 #"WTDFS"="\\wtdfs~~~",
                                 #"SASCAN"="\\sascan~~~"))

#tikz("fig-var-txt-size.tex",width = 6.0, height = 2.8)

plot <- ggplot(d,aes(factor(result_count_num),timenum,fill=algo,color=algo))
plot <- plot + geom_boxplot(outlier.size = 1)
#plot <- plot + geom_point()
plot <- plot + facet_grid(. ~ instance)
plot <- plot + scale_x_discrete(name = "Result count")
plot <- plot + scale_y_log10(name = "Query time [µs]") #,breaks=c(0.01,0.1,1,10,100,1000,10000),labels=c("0.01","0.1","1","10","100","1k","10k"))
print(plot)

#dev.off()
