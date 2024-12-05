rm(list=ls())
##############
# Problem_1
#e)
##############

set.seed(123)
n_simulations <- 10000
         
############### Simulations a ##############

throws <- replicate(n_simulations, sum(sample(1:6, 3, replace = TRUE)))
sim_prob_sum_9 <- mean(throws == 9)
sim_prob_sum_10 <- mean(throws == 10)

print(paste("P(sum 9):", sim_prob_sum_9))
print(paste("P(sum 10):", sim_prob_sum_10))


############### Simulations b ##############


throws_4 <- replicate(n_simulations, sum(sample(1:6, 4, replace = TRUE) == 6))
prob_at_least_one_six_4 <- mean(throws_4 >= 1)

throws_24 <- replicate(n_simulations, sum(sample(1:6, 24, replace = TRUE) == 6))
prob_at_least_two_sixes_24 <- mean(throws_24 >= 2)

print(paste("P(at least one six in 4 throws):", prob_at_least_one_six_4))
print(paste("P(at least two sixes in 24 throws):", prob_at_least_two_sixes_24))

############### Simulations c ##############

throws_6 <- replicate(n_simulations, sum(sample(1:6, 6, replace = TRUE) == 6))
prob_at_least_one_six_6 <- mean(throws_6 >= 1)

throws_12 <- replicate(n_simulations, sum(sample(1:6, 12, replace = TRUE) == 6))
prob_at_least_two_sixes_12 <- mean(throws_12 >= 2)

throws_18 <- replicate(n_simulations, sum(sample(1:6, 18, replace = TRUE) == 6))
prob_at_least_three_sixes_18 <- mean(throws_18 >= 3)

print(paste("P(at least one six in 6 throws):", prob_at_least_one_six_6))
print(paste("P(at least two sixes in 12 throws):", prob_at_least_two_sixes_12))
print(paste("P(at least three sixes in 18 throws):", prob_at_least_three_sixes_18))

############### Simulations d ##############
success_in_8_throws <- replicate(n_simulations, sum(sample(1:6, 8, replace = TRUE) == 6) == 3)
prob_exactly_8_throws <- mean(success_in_8_throws)

print(paste("P(exactly 8 throws to get 3 sixes):", prob_exactly_8_throws))

##############
# Problem_2
#e
##############

# Inverse CDF function 
inverse_Fa <- function(U, a) {
  return(a - log(-log(U)))
}

# generate n random numbers(uniform distribution)
generate_random_gumbel <- function(n, a) {
  U <- runif(n) 
  return(inverse_Fa(U, a))
}

n <- 1000 
a <- 5     

random_samples <- generate_random_gumbel(n, a)

# Plot KDE
hist(random_samples, prob = TRUE, main = "Histogram of Random Gumbel Samples (n = 1000)", 
     xlab = "Random Samples", col = "lightblue", breaks = 30)

lines(density(random_samples), col = "red", lwd = 2)




#second way for KDE
#random_samples <- generate_random_gumbel(n = 1000, a = 5)
#plot(density(random_samples), main = "Gumbel KDE and True Gumbel Density", xlab = "x")
#curve(gumbel_density(x, a, b), add = TRUE, col = "red", lwd = 2)

##############
#f
##############
gamma <- 0.5772
sample_mean <- mean(random_samples)
theoretical_mean <- a + gamma

# Compare the sample and the theoretical mean
print(paste("Sample Mean:", sample_mean))
print(paste("Theoretical Mean (E(T)):", theoretical_mean))

print(paste("Difference:", abs(sample_mean - theoretical_mean)))

##############
#g
##############
sample_variance <- var(random_samples)
#sample_sd <- sqrt(sample_variance)

precision <- 0.01

# Number of sample - desired precision
required_n <- round(4*(sample_variance) / (precision^2))

print(paste("Sample Variance:", sample_variance))
print(paste("Required number of samples for precision 0.01:", ceiling(required_n)))

##############
#h
##############
Fa <- function(t, a) {
  return(exp(-exp(-(t - a))))
}

# Theoretical CDF
curve(Fa(x, a), from = min(random_samples), to = max(random_samples), 
      col = 'blue', lwd = 8, ylab = 'CDF', xlab = 'T', 
      main = "Theoretical CDF vs. ECDF")

# Add ECDF 
lines(ecdf(random_samples), col = 'red', lwd = 3)


legend("bottomright", legend = c("Theoretical CDF", "Empirical CDF"), 
       col = c("blue", "red"), lwd = 2)


# The theoretical CDF is a smooth curve, while the ECDF is a step function.
# As the sample size increases, the ECDF should converge to the theoretical CDF, illustrating that the generated random samples accurately reflect the underlying probability distribution.

##############
#Problem_3
#b
##############

# Compute MLE 
sample_data <- c(1.3, -0.6, 0.2, 0.4, -0.8)
lambda_hat <- mean(abs(sample_data))
print(paste("The MLE for lambda is:", lambda_hat))

##############
# c
##############

# Plot density - lambda = 0.6
lambda <- 0.6
f_laplace <- function(x, lambda) {
  return((1 / (2 * lambda)) * exp(-abs(x) / lambda))
}

x_vals <- seq(-5, 5, by = 0.01)
pdf_vals <- f_laplace(x_vals, lambda)

plot(x_vals, pdf_vals, type = "l", col = "blue", lwd = 2, 
     main = expression(paste("PDF of Laplace Distribution with ", lambda, "=0.6")),
     ylab = "Density", xlab = "x")

##############
# d
##############

# Generate random numbers - accept-reject method
gen_laplace <- function(n, lambda) {
  samples <- numeric(n)
  for (i in 1:n) {
    while (TRUE) {
      samples[i] <- runif(1, -5, 5) 
      aprob <- (6/5)*f_laplace(samples[i],lambda = 0.6)
      u <- runif(1)
      if (u <= aprob) break
      }
    }
  }
  return(samples)
}

set.seed(123)
n <- 500
samples <- gen_laplace(n, lambda)

hist(samples, freq = FALSE, main="Histogram of Generated Laplace Samples", 
     ylim=c(0, max(pdf_vals)*1.1)) 

lines(density(samples), col = "red", lwd = 2)
curve(f_lambda(x, lambda), add = TRUE, col = "blue", lwd = 3)

legend("topright", legend = c(" KDE", " PDF"), 
       col = c("red", "blue"), lwd = c(2, 3))

##############
# e
##############

# Investigate different bin-width rules and bandwidth rules for KDE

# Histogram 
hist_sturges <- hist(random_data, breaks = "Sturges", plot = FALSE)
hist_scott <- hist(random_data, breaks = "Scott", plot = FALSE)
hist_fd <- hist(random_data, breaks = "FD", plot = FALSE)

# KDE
kde_silverman <- density(random_data, bw = "nrd0")
kde_ucv <- density(random_data, bw = "ucv")
kde_bcv <- density(random_data, bw = "bcv")
kde_sj <- density(random_data, bw = "SJ")

# Different bin-widths
par(mfrow = c(1, 3))
plot(hist_sturges, col = "lightblue", main = "Sturges' formula", xlim = c(-5, 5))
plot(hist_scott, col = "lightgreen", main = "Scott's formula", xlim = c(-5, 5))
plot(hist_fd, col = "lightcoral", main = "Freedman-Diaconis formula", xlim = c(-5, 5))

# Different bandwidths
par(mfrow = c(1, 4))
plot(kde_silverman, main = "KDE with Silverman's Rule", col = "blue", xlim = c(-5, 5))
plot(kde_ucv, main = "KDE with UCV's Rule", col = "red", xlim = c(-5, 5))
plot(kde_bcv, main = "KDE with BCV's Rule", col = "green", xlim = c(-5, 5))
plot(kde_sj, main = "KDE with BCV's Rule", col = "purple", xlim = c(-5, 5))


#Histogram bin-width rules:
  
#1.Sturges' rule tends to use fewer bins, which can oversmooth the data.
#2.Scott's rule and Freedman-Diaconis (F-D) rule often provide more bins, revealing more detail in the distribution.
#3.For this Laplace distribution, F-D rule might provide the best balance between detail and smoothness.


#KDE bandwidth rules:
  
#1.Silverman's rule of thumb is simple but can oversmooth for non-normal distributions.
#2.Unbiased and Biased Cross-Validation methods can be more adaptive but may be sensitive to outliers.
#3.Sheather & Jones method often provides a good balance and is generally recommended for a wide range of distributions.

##############
# f
##############
F_lambda <- function(x, lambda) {
  0.5 + 0.5 * sign(x) * (1 - exp(-abs(x) / lambda))
}

x_vals_cdf <- seq(-5, 5, length.out = 100)
cdf_vals <- F_lambda(x_vals_cdf, lambda)

# ECDF
ecdf_vals <- ecdf(random_data)

# Plot CDF and ECDF
plot(x_vals_cdf, cdf_vals, type = "l", col = "blue", lwd = 10, main = "CDF vs ECDF", xlim = c(-5, 5))
plot(ecdf_vals, col = "red", add = TRUE, lwd=1)
legend("topleft", legend = c("CDF", "ECDF"), col = c("blue", "red"), lty = 1, lwd = 1)


##############
# Problem_4
# a
##############

# Return the number of trials
rps_trials <- function(k) {
  lives1 <- k
  lives2 <- k
  trials <- 0 
  
  while(lives1 > 0 & lives2 > 0) {
    trials <- trials + 1
    # Simulate a round(1 for player1 win, -1 for player2 win, 0 for tie)
    result <- sample(c(1, -1, 0), 1, prob = c(1/3, 1/3, 1/3))
    if (result == 1) {
      # Player2 loses
      lives2 <- lives2 - 1
    } else if (result == -1) {
      # Player1 loses 
      lives1 <- lives1 - 1
    }
  }
  return(trials)
}

##############
#b
##############
#Compare the result with 10^5 simulation in geometric distribution with prob(2/3)
n_sim <- 10^5
trials_k1 <-numeric(n_sim)
p <- 2/3
mean_geom <- 1 / p

par(mfrow = c(1, 1))
mean_simulated <- mean(trials_k1)
print(paste("Simulated mean for k = 1:", mean_simulated))
print(paste("Theoretical mean for k = 1:", mean_geom))


hist(trials_k1, probability = TRUE, breaks = max(trials_k1), main = "Distribution of Trials for k = 1")
curve(dgeom(0:9, p), col = "red", add = TRUE, lwd=2)

##############
#c
##############

# Simulate for k = 3
trials_k3 <- replicate(n_sim, rps_trials(3))
mean_k3 <- mean(trials_k3)
prob_8_or_more_k3 <- mean(trials_k3 >= 8)

# Simulate for k = 5
trials_k5 <- replicate(n_sim, rps_trials(5))
mean_k5 <- mean(trials_k5)
prob_8_or_more_k5 <- mean(trials_k5 >= 8)

print(paste("Expected number of trials for k = 3:", mean_k3))
print(paste("Probability of 8 or more trials for k = 3:", prob_8_or_more_k3))

print(paste("Expected number of trials for k = 5:", mean_k5))
print(paste("Probability of 8 or more trials for k = 5:", prob_8_or_more_k5))

hist(trials_k3, probability = TRUE, breaks = max(trials_k3), main = "Distribution of Trials for k = 3")
hist(trials_k5, probability = TRUE, breaks = max(trials_k5), main = "Distribution of Trials for k = 5")
