#################### Gibbs Sampler ##############################
gibbs_sampler_bivariate_normal <- function(n.iter = 10000, rho = 0.95, theta_init = c(3, -2)) {
  theta <- theta_init 
  res <- matrix(0.0, 2 * n.iter + 1, 2)  
  res[1,] <- theta
  k <- 2
  
  for (i in 1:n.iter) {
    # Update first block 
    theta[1] <- rnorm(1, mean = rho * theta[2], sd = sqrt(1 - rho^2))
    res[k,] <- theta
    k <- k + 1
    
    # Update second block 
    theta[2] <- rnorm(1, mean = rho * theta[1], sd = sqrt(1 - rho^2))
    res[k,] <- theta
    k <- k + 1
  }
  
  return(res)
}

set.seed(123)  
n.iter <- 10000  
rho <- 0.95  
theta_init <- c(3, -2)  

result <- gibbs_sampler_bivariate_normal(n.iter, rho, theta_init)

# Plot the results
plot(result[,1], result[,2], type = "l", main = "Gibbs Sampler for Bivariate Normal",
     xlab = "Theta[1]", ylab = "Theta[2]")

cov(result)
colMeans(result)
