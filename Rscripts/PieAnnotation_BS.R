library(ggplot2)
mydata <- read.table("/home/vimal/Desktop/pieDMR.txt",header=TRUE)
ggplot(mydata, aes(x=DMC/2, y = Per, fill = feature, width = DMC)) +
  geom_bar(position="fill", stat="identity") +
  facet_grid(Sample1 ~ Context) +
  coord_polar("y")
