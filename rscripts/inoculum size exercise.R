#As the volume of inoculum decreases, the variance increases. I'm worried that if we use very low volumes of inoculum, that's going to introduce a lot of variation that I don't want to deal with. I simulate how the range of values, mean and SE change as we increase the inoculum load. I start with a vector of 10 values, which I got by counting all the sporangia in a 3 ul drop. This corresponds to 3000-7667 sporangia/ml, with a mean of 5167. I assume a normal distribution of this sample and simulate additional samples and aggregate them according to inoclum size. I show that range and se decrease as size increases, but mean remains the same.

#vector of values after counting sporangia in 3 ul of spore solution. These are actual counts that I made in the lab by dropping 3 ul onto a slide and counting all the sporangia I saw.
v <- c(18, 14, 16, 23, 15, 21, 15, 11, 9, 13) 
#vector corresponding to sporangia/ml
v/3*1000
hist(v)
#I need to write a function that creates a normal distribution using the mean and sd of `v` above. From this normal distribution, I will sample some "amount" of points and then sum together every "y" points so that the total number of points is always the same. For example, mcmc(amount=10, y=1) will give me 10/1=10 sample points. mcmc(amount=20, 2) will give me 20 points, but every 2 points are summed together giving me a total of 10 points.
mcmc <- function(amount, y){
  x <- round(rnorm(amount, mean=mean(v), sd=sd(v)), 0)#create a random sample, modeled after the distribution of my real sample
  x2 <- tapply(x, rep(1:(length(x)/y), each=y), FUN = sum) #aggregate it
  #convert to sporangia/ml
  x3 <- x2/(3*y)*1000
return(x3)
}

#function for calculating standard error
se <- function(x){sd(x)/sqrt(length(x))}

#sample 10 points from the normal distribution using the `mcmc` function and replicate that 500 times.
#each line represents a different simulated inoculum volume. `r3` corresponds to 3 ul, `r6` to 6 ul, and so on.
r3 <- replicate(500, mcmc(amount=10, y=1))
r6 <- replicate(500, mcmc(20, 2))
r9 <- replicate(500, mcmc(30, 3))
r12 <- replicate(500, mcmc(40, 4))
r15 <- replicate(500, mcmc(50, 5))
r18 <- replicate(500, mcmc(60, 6))
r21 <- replicate(500, mcmc(70, 7))
r24 <- replicate(500, mcmc(80, 8))
r27 <- replicate(500, mcmc(90, 9))
r30 <- replicate(500, mcmc(100, 10))
r100 <- replicate(500, mcmc(330,33))

#find an efficient way to plot the results
l <- list(r3, r6, r9, r12, r15, r18, r21, r24, r27, r30, r100)
names(l) <- c("r3", "r6", "r9", "r12", "r15", "r18", "r21", "r24", "r27", "r30", "r100")
#plot the range of values. Can see that the whiskers on the box plots decrease with increasing volume. Seems like after 15 ul, the decrease in variance is mostly negligible. 
boxplot(l, main="range of values", ylab="#sporangia/ml")
#plot the mean of values
boxplot(lapply(l, function(x)apply(x, 2, mean)), main="mean of 3-30 ul of inoculum", ylab="#sporangia/ml")
#plot the SE of values
boxplot(lapply(l, function(x)apply(x, 2, se)), main="se of 3-30 ul of inoculum", ylab="#sporangia/ml")

