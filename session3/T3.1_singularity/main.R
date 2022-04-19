temperature <- read.csv(url("https://raw.githubusercontent.com/datasets/global-temp/master/data/annual.csv"))

jpeg(file="temperature.jpg")
plot(x=temperature$Year,y=temperature$Mean)
dev.off()
