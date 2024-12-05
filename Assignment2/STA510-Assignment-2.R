# Clear the environment
rm(list=ls())


################################################################
# Problem_1
################################################################
#a
################################################################
set.seed(123)  
n <- 10000 

f <- function(x) {
  exp(-x) / (1 + x^2)^2
}

a <- -2 - (sqrt(3))
b <- -1/sqrt(3)

random_samples <- runif(n, min = a, max = b)

# Estimate the integral
I_CMC <- (b - a) * mean(f(random_samples))
print(paste("Estimate of I using crude Monte Carlo (I_CMC):", I_CMC))

###################### cmc.fun1 ###################### 
cmc.fun1 <- function(Nsim=10000,a,b){
  x <- runif(Nsim,a,b)
  int <- (b-a)*mean(f(x))
  return(int)
}

#b
################################################################
x_vals <- seq(-5, 5, length.out = 1000)
f_vals <- f(x_vals)

plot(x_vals, f_vals, type = "l",xlim=c(-5,5), ylim=c(0,2), col = "blue", lwd = 2,
     main = "Plot of f(x) = exp(-x) / (1 + x^2)^2",
     xlab = "x", ylab = "f(x)")
abline(v=c(a,b), lty=2, col="blue")

#c
################################################################
n <-10000
n2 <- n/2

u <- runif(n2, min = a, max = b)
v <- b + a +(-u)  

I_AT <- (b - a) * mean((f(u) + f(v)) / 2)
print(paste("Estimate of I using antithetic variables (I_AT):", I_AT))

###################### amc.fun1 ###################### 
amc.fun1 <- function(Nsim=10000,a,b){
  int <- numeric(Nsim)
  x1 <- runif(Nsim/2,a,b)
  x2 <- a+b-x1
  x <- c(x1,x2)
  int <- (b-a)*mean(f(x))
  return(int)
}

#e
################################################################
f <- function(x) {
  exp(-x) / (1 + x^2)^2
}

g <- function(x) {
  c <- 6 / pi
  ifelse(x >= -2/sqrt(3) & x <= -1/sqrt(3), c / (1 + x^2), 0)
}

x_vals <- seq(-5, 5, length.out = 1000)

f_vals <- f(x_vals)
g_vals <- g(x_vals)
max_y <- max(f_vals, g_vals)


plot(x_vals, f_vals, type = "l", col = "blue", lwd = 2, ylim = c(0, max_y * 1.1),
     main = "Plot of f(x) and g(x)", xlab = "x", ylab = "Value")
lines(x_vals, g_vals, col = "red", lwd = 2)
legend("topright", legend = c("f(x)", "g(x)"), col = c("blue", "red"), lwd = 2)



#g
################################################################
set.seed(123)
nrep <- 10000 

###################### imc.fun1 ###################### 
imc.fun1 <- function(nrep) {
  x <- runif(nrep,0,1)
  rv <- tan( pi*(x -5/3)/4 )
  # Importance sampling
  wts <- f(rv)/f(rv)
  int <- mean(wts)
  return(int)
}
print(paste0("Estimated integral I by antithetic variables: ",
             round(imc.fun1(nrep),digits=3) ))

#h
################################################################
set.seed(123)  
nSim <- 1000
intCMC <- numeric(nrep)
intAMC <- numeric(nrep)
intIMC <- numeric(nrep)
for(i in 1:nrep){
  intCMC[i] <- cmc.fun1(nSim,a,b)
  intAMC[i] <- amc.fun1(nSim,a,b)
  intIMC[i] <- imc.fun1(nSim)
}
print( rbind( paste0("Estimated integral I by CMC: ", round(mean(intCMC),digits=3)),
              paste0("Estimated integral I by antithetic MC: ", round(mean(intAMC),digits=3)),
              paste0("Estimated integral I by importance sampling: ", round(mean(intIMC),digits=3)) ))

print( rbind (paste0("Estimated SD of the CMC: ", round(sd(intCMC),digits=5)),
              paste0("Estimated SD of the antithetic MC: ", round(sd(intAMC),digits=5)),
              paste0("Estimated SD of the importance sampling: ", round(sd(intIMC),digits=5))))

#i
################################################################
alpha <- 0.01
epsilon <- 0.001
ceiling(( (qnorm(1-alpha/2)^2)*((b-a)^2)*var(intCMC) ) / epsilon^2)


################################################################
# Problem_2
################################################################
#a
################################################################

# Implementation based on simulating times between the events
plotHPP <- function(lambda, stoptime){
  expectednumber <- stoptime * lambda  # Expected number of events until stoptime
  Nsim <- 3 * expectednumber  # Simulate more than the expected number to be certain to exceed stoptime
  timesbetween <- rexp(Nsim, lambda) # Simulate interarrival times
  timesto <- cumsum(timesbetween)   # Calculate arrival times
  timesto <- timesto[timesto < stoptime] # Discard the times larger than stoptime
  Nevents <- length(timesto) # Count the number of events
  plot(timesto, 1:Nevents, type="s", xlab="arrival time", 
       ylab="Event number", lwd=1.5, ylim=c(0, Nevents), 
       main=paste("Poisson Process N(t) with λ =", lambda))
  points(timesto, rep(0, Nevents), pch=21, bg="red") # Mark the events
}

lambda <- 3
stoptime <- 30

plotHPP(lambda, stoptime)
##### Alternative method #######
n <- 5000
x <- runif(n,1,4)
y <- runif(n,-2,2)
z <- runif(n,0,1)
int1 <- (3*4*1) * mean(1/(1+xˆ2+yˆ2+zˆ2))

# Second integral
w <- rexp(n, rate=1/4)
int2 <-4*mean(w<10)
# Estimate
int1*int2

#c
################################################################
# Function to plot multiple realizations
plot_multiple_HPP <- function(lambda, stoptime, n_realizations) {
  plot(0, 0, type="n", xlim=c(0, stoptime), ylim=c(0, stoptime), 
       xlab="Arrival Time", ylab="Event Number", 
       main=paste("40 Realizations of Poisson Process N(t) with λ =", lambda))
  
  for (i in 1:n_realizations) {
    expectednumber <- stoptime * lambda
    Nsim <- 3 * expectednumber
    timesbetween <- rexp(Nsim, lambda)
    timesto <- cumsum(timesbetween)
    timesto <- timesto[timesto < stoptime]
    Nevents <- length(timesto)
    
    lines(timesto, 1:Nevents, type="s", lwd=1.5, col=rainbow(n_realizations)[i])
    points(timesto, rep(0, Nevents), pch=21, bg=rainbow(n_realizations)[i])
  }
  
  #legend("topright", legend=paste("Realization", 1:n_realizations), col=rainbow(n_realizations), lwd=1.5, bty="n")
}

lambda <- 3
stoptime <- 30
n_realizations <- 40

plot_multiple_HPP(lambda, stoptime, n_realizations)


#d
################################################################
# Function to simulate the Poisson process
simulate_poisson_process <- function(lambda, stoptime, n_realizations) {
  event_counts <- numeric(n_realizations)
  
  for (i in 1:n_realizations) {
    # Simulate inter-arrival times
    timesbetween <- rexp(3 * lambda * stoptime, lambda)  # More than expected events
    timesto <- cumsum(timesbetween)  # Calculate arrival times
    timesto <- timesto[timesto < stoptime]  # Discard times larger than stoptime
    event_counts[i] <- length(timesto)  # Count the number of events
  }
  
  return(event_counts)
}

lambda <- 3
n_simulations <- 10000  

# i) Calculate the probability of observing at least 80 events in [0, 30]
event_counts_30 <- simulate_poisson_process(lambda, 30, n_simulations)
prob_at_least_80 <- mean(event_counts_30 >= 80)

# ii) Calculate the probability of observing less than 30 events in [0, 10]
event_counts_10 <- simulate_poisson_process(lambda, 10, n_simulations)
prob_less_than_30 <- mean(event_counts_10 < 30)

print(paste("Probability of observing at least 80 events in [0, 30]:", prob_at_least_80))
print(paste("Probability of observing less than 30 events in [0, 10]:", prob_less_than_30))


#e
################################################################
library(truncnorm)
library(stats)

f_truncated_normal <- function(x, mu, sigma, a, b) {
  numerator <- dnorm((x-mu)/sigma) / sigma
  denominator <- pnorm((b-mu)/sigma) - pnorm((a-mu)/sigma)
  return(numerator / denominator)
}

lambda <- function(t) {
  100 * f_truncated_normal(t, mu=12.5, sigma=1, a=8, b=17.5)
}

t_seq <- seq(7, 18, by=0.1)

# Calculate intensity values
lambda_values <- sapply(t_seq, lambda)

plot(t_seq, lambda_values, 
     type="l", 
     col="blue", 
     lwd=2,
     xlab="Time (hours)",
     ylab="λ(t)",
     main="Store Arrival Intensity Function")

grid()

abline(v=8, col="red", lty=2)
abline(v=17.5, col="red", lty=2)
#legend("topleft",legend=c("Intensity function", "Store hours"), col=c("blue", "red"), lty=c(1, 2), lwd=c(2, 1))


#f
################################################################

# Function for simulating arrival times for a NHPP between a and b using thinning
simtNHPP <- function(a,b,lambdamax,lambdafunc){
  if(max(lambdafunc(seq(a,b,length.out = 100)))>lambdamax)
    stop("lambdamax is smaller than max of the lambdafunction")
  expectednumber <- (b-a)*lambdamax  
  Nsim <- 3*expectednumber 
  timesbetween <- rexp(Nsim,lambdamax) 
  timesto <- a+cumsum(timesbetween)   
  timesto <- timesto[timesto<b] 
  Nevents <- length(timesto) 
  U <- runif(Nevents)
  timesto <- timesto[U<lambdafunc(timesto)/lambdamax]  
  timesto  
}
t_seq <- seq(7, 18, by=0.1)
lambdamax <- max(sapply(t_seq, lambda))  # Calculate maximum intensity
set.seed(123)  # For reproducibility

# Generate one realization of the NHPP
arrival_times <- simtNHPP(a=7, b=18, 
                          lambdamax=lambdamax, 
                          lambdafunc=lambda)

par(mar=c(4,4,2,1))

plot(arrival_times, 1:length(arrival_times), 
     type="s",
     xlab="Time (hours)", 
     ylab="Cumulative number of arrivals",
     main="Non-homogeneous Poisson Process Realization",
     lwd=1.5)

points(arrival_times, rep(0, length(arrival_times)), 
       pch=21, 
       bg="red",
       cex=0.8)
grid()

abline(v=c(8, 17.5), col="red", lty=2)
legend("topleft", 
       legend=c("Cumulative arrivals", "Store hours"),
       col=c("black", "red"), 
       lty=c(1, 2),
       lwd=c(1.5, 1))

print(paste("\nSimulation Summary:\n"))
print(paste("==================\n"))
print(paste(sprintf("Total number of arrivals: %d\n", length(arrival_times))))
print(paste(sprintf("First arrival time: %.2f\n", min(arrival_times))))
print(paste(sprintf("Last arrival time: %.2f\n", max(arrival_times))))
print(paste(sprintf("Average time between arrivals: %.3f hours\n", 
            mean(diff(arrival_times)))))
#g
################################################################

n_sim <- 10000
lambdamax <- max(sapply(seq(7, 18, by=0.1), lambda))
set.seed(123) 

arrivals_before_10 <- numeric(n_sim)
arrivals_11_to_13 <- numeric(n_sim)

# Function to count arrivals in a specific time interval
count_arrivals <- function(times, start, end) {
  sum(times >= start & times < end)
}

# Run simulations
for(i in 1:n_sim) {
  # Generate one realization
  times <- simtNHPP(a=7, b=18, lambdamax=lambdamax, lambdafunc=lambda)
  
  arrivals_before_10[i] <- count_arrivals(times, 7, 10)
  arrivals_11_to_13[i] <- count_arrivals(times, 11, 13)
}

# Calculate statistics
avg_before_10 <- mean(arrivals_before_10)
sd_before_10 <- sd(arrivals_before_10)
avg_11_to_13 <- mean(arrivals_11_to_13)
sd_11_to_13 <- sd(arrivals_11_to_13)

print(paste("\nResults from", n_sim, "replications:\n"))
print(paste("=====================================\n"))
print(paste(sprintf("Average arrivals before 10:00: %.2f (SD: %.2f)\n", 
            avg_before_10, sd_before_10)))
print(paste(sprintf("Average arrivals between 11:00-13:00: %.2f (SD: %.2f)\n", 
            avg_11_to_13, sd_11_to_13)))

par(mfrow=c(1,2))

# Histogram for arrivals before 10:00
hist(arrivals_before_10, 
     main="Arrivals before 10:00",
     xlab="Number of arrivals",
     col="lightblue",
     breaks=30)
abline(v=avg_before_10, col="red", lwd=2)

# Histogram for arrivals between 11:00-13:00
hist(arrivals_11_to_13, 
     main="Arrivals between 11:00-13:00",
     xlab="Number of arrivals",
     col="lightgreen",
     breaks=30)
abline(v=avg_11_to_13, col="red", lwd=2)

par(mfrow=c(1,1))

################################################################
# Problem_3
################################################################
#a
################################################################
set.seed(123)
n <- 5000

simulate_integral <- function() {
  x <- runif(1, 1, 2)      
  y <- runif(1, -2, 1)     
  z <- runif(1, 0, 1)     
  w <- rexp(1, 4)          
  
  if(w <= 10) { 
    integrand <- 1 / (1 + x^2 + y^2 + z^2) * 4 * exp(-w/4)  # 4 is the rate parameter
    vol_xyz <- (2 - 1) * (1 - (-2)) * (1 - 0)  # Volume for x,y,z
    return(integrand * vol_xyz)
  } else {
    return(0)  
  }
}

# Perform Monte Carlo integration
results <- replicate(n, simulate_integral())

# Calculate estimate and standard error
estimate <- mean(results)
std_error <- sd(results) / sqrt(n)

# Calculate 95% confidence interval
ci_lower <- estimate - 1.96 * std_error
ci_upper <- estimate + 1.96 * std_error

# Print results
print(paste("\nMonte Carlo Integration Results\n"))
print(paste("==============================\n"))
print(paste(sprintf("Estimate: %.6f\n", estimate)))
print(paste(sprintf("Standard Error: %.6f\n", std_error)))
print(paste(sprintf("95%% Confidence Interval: [%.6f, %.6f]\n", ci_lower, ci_upper)))


par(mfrow=c(1,2))
hist(results, 
     breaks=50, 
     main="Distribution of Integral Estimates",
     xlab="Estimate Value",
     col="lightblue")
abline(v=estimate, col="red", lwd=2)

running_mean <- cumsum(results)/(1:n)
plot(1:n, running_mean, 
     type="l", 
     main="Running Mean",
     xlab="Number of Simulations",
     ylab="Estimate")
abline(h=estimate, col="red", lwd=2)
par(mfrow=c(1,1))


################################################################
# Problem_4
################################################################
#a
################################################################

install.packages("shapes")
library(shapes)


data(panf.dat)
data(panm.dat)

plotshapes(panf.dat)  
plotshapes(panm.dat)  

shape1 <- panf.dat
shape2 <- panm.dat

dim(panf.dat)  
dim(panm.dat)  
#b
################################################################

cs1 <- centroid.size(shape1)
cs2 <- centroid.size(shape2)
summary(cs1)

summary(cs2)

# Bootstrap estimate for standard deviation and bias
B <- 5000
n1 <- length(cs1)
n2 <- length(cs2)
meancentroid1 <- numeric(B)
meancentroid2 <- numeric(B)
for(i in 1:B){
  meancentroid1[i] <- mean(sample(cs1,size=n1,replace=TRUE))
  meancentroid2[i] <- mean(sample(cs2,size=n2,replace=TRUE))
}
hist(meancentroid1,prob=TRUE,nclass=sqrt(B)/2,
     main="Female chimpanzee")
abline(v=mean(cs1),col="blue",lwd=2)

hist(meancentroid2,prob=TRUE,nclass=sqrt(B)/2,
     main="Male chimpanzee")
abline(v=mean(cs2),col="blue",lwd=2)

#c
################################################################
# Mean centroid
mu1 <- mean(cs1)
mu2 <- mean(cs2)

# Standard error
sd(meancentroid1)
sd(meancentroid2)

mean(meancentroid1)-mu1
mean(meancentroid2)-mu2

# Percentile CI interval for centroid size
lowerP1 <- quantile(meancentroid1,0.025)
upperP1 <- quantile(meancentroid1,0.975)
lowerP2 <- quantile(meancentroid2,0.025)
upperP2 <- quantile(meancentroid2,0.975)
method <- c("Female","Male")
lowerCI <- round(as.vector(c(lowerP1,lowerP2)),2)
upperCI <- round(as.vector(c(upperP1,upperP2)),2)
rbind(method,lowerCI,upperCI)

