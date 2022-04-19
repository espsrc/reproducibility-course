install.packages("ggplot2")
library(ggplot2)

mydata<-economics
theme_set(theme_light())

ggplot(mydata, aes(x=date)) +
   geom_line(aes(y=unemploy), color = "#00AFBB") +
   labs(y="Unemployment", x="Year")

