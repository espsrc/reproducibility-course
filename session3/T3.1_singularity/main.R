temperature <- read.csv(url("https://raw.githubusercontent.com/datasets/global-temp/master/data/annual.csv"))

pdf("temerature.pdf") 
plot(x=temperature$Year,y=temperature$Mean)
dev.off()
