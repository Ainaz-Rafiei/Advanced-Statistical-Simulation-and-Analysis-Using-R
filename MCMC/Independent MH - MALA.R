################## Independent MH ##########################
independent_MH <- function(target_prob, n.iter = 10000, sigma = 0.7) {
  res <- numeric(n.iter)
  res[1] <- 0.0
  wt.old <- target_prob(res[1])/dnorm(res[1], sd = sigma)
  Nacc <- 0
  for (i in 2:n.iter) {
    # Proposal
    thetaStar <- rnorm(1, sd = sigma)
    wt.star <- target_prob(thetaStar)/dnorm(thetaStar, sd = sigma)
    alpha <- min(1.0, wt.star/wt.old)
    # Accept/reject
    if (runif(1) < alpha) {
      res[i] <- thetaStar
      wt.old <- wt.star
      Nacc <- Nacc + 1
    } else {
      res[i] <- res[i - 1]
    }
  }
  print(paste0("accept rate : ", Nacc / n.iter))
  return(res)
}

target_prob <- function(theta) { return(dunif(theta, min = -1, max = 1)) }

n.iter <- 10000
sigma <- 0.7
res <- independent_MH(target_prob, n.iter, sigma)

# Plot the results
plot(res, pch = 20, cex = 0.1, xlab = "MCMC iteration #")

################## MALA function ########################## 
oneDMALA <- function(lp_grad, theta1 = 0.0, Delta = 0.5, n.iter = 10000) {
  out <- numeric(n.iter)
  out[1] <- theta1
  old.lpg <- lp_grad(theta1)
  Nacc <- 0
  
  for (i in 2:n.iter) {
    # Proposal
    thetaStar <- out[i - 1] + 0.5 * Delta * old.lpg$grad + sqrt(Delta) * rnorm(1)
    new.lpg <- lp_grad(thetaStar)
    
    lqf <- dnorm(thetaStar, mean = out[i - 1] + 0.5 * Delta * old.lpg$grad, sd = sqrt(Delta), log = TRUE)
    lqb <- dnorm(out[i - 1], mean = thetaStar + 0.5 * Delta * new.lpg$grad, sd = sqrt(Delta), log = TRUE)
    
    alpha <- exp(min(0.0, new.lpg$lp + lqb - old.lpg$lp - lqf))
    
    if (runif(1) < alpha && is.finite(alpha)) {
      out[i] <- thetaStar
      old.lpg <- new.lpg
      Nacc <- Nacc + 1
    } else {
      out[i] <- out[i - 1]
    }
  }
  
  print(paste0("MALA done, accept rate: ", Nacc / (n.iter - 1)))
  return(out)
}

mu <- 1.0
sigma <- 0.5
lp_gauss <- function(theta) {
  return(list(
    lp = -(theta - mu)^2 / (2 * sigma^2), 
    grad = -(theta - mu) / (sigma^2)  # log-g gradient
  ))
}

out.mala <- oneDMALA(lp_gauss, theta1 = mu, Delta = 0.7, n.iter = 50000)

# Plot the results
par(mfrow = c(2, 1))
ts.plot(out.mala)
hist(out.mala, probability = TRUE, breaks = 30)
xg <- seq(from = mu - 5 * sigma, to = mu + 5 * sigma, length.out = 1000)
lines(xg, dnorm(xg, mean = mu, sd = sigma), col = "red")

