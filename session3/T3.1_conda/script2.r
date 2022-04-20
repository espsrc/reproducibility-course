library(ggplot2) 

carsUSA <- read.csv("cars.csv",header= TRUE)

carsUSA$gear <- factor(carsUSA$gear,levels=c(3,4,5),
                       labels=c("3","4","5")) 
carsUSA$am <- factor(carsUSA$am,levels=c(0,1),
                     labels=c("0","1")) 
carsUSA$cyl <- factor(carsUSA$cyl,levels=c(4,6,8),
                      labels=c("4","6","8")) 

qplot(mpg, data=carsUSA, geom="density", fill=gear, alpha=I(.5), 
      main="Distribution of Gas Milage", xlab="Miles Per Gallon", 
      ylab="Density")

ggsave("mtcars.png")
