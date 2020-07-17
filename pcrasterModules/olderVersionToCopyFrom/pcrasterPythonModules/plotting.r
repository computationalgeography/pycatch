a <- read.table("test.txt")
## transpose the table,remove the time value and plot the qq plot
qqnorm(t(a)[,2][2,100])

