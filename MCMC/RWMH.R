##################  RWMH ################## 
oneD.RWMH <- function(lprob, sigma=1.0, theta1=0.0, n.iter=10000){
  output <- numeric(n.iter)
  output[1] <- theta1
  for(t in 2:n.iter){
    thetaStar <- output[t-1] + rnorm(1, sd=sigma)
    alpha <- exp(min(0.0, lprob(thetaStar) - lprob(output[t-1])))
    if(runif(1) < alpha && is.finite(alpha)){
      output[t] <- thetaStar
    } else {
      output[t] <- output[t-1]
    }
  }
  return(output)
}

lp_std_norm <- function(x) { return(dnorm(x, log=TRUE)) }

out <- oneD.RWMH(lprob = lp_std_norm, theta1 = 0.0, sigma = 2.4, n.iter = 10000)

# Plot the results
par(mfrow = c(2, 1))
plot(1:length(out), out, pch = 20, cex = 0.1, xlab = "MCMC iteration #", main = "Random Walk MCMC")
hist(out, main = "Histogram of MCMC Output", probability = TRUE)

# Check the mean and STD
print(paste("Mean of output:", mean(out)))
print(paste("Standard deviation of output:", sd(out)))
