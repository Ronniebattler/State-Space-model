# Simulate some data
set.seed(123)
n <- 100
phi_true <- 0.5
sigma_true <- 1
y <- arima.sim(n=n, list(ar=phi_true, order=c(1,0,0), sd=sigma_true))

# Initialize parameters for Gibbs sampler
niter <- 10000
phi_samples <- numeric(niter)
sigma2_samples <- numeric(niter)
phi_samples[1] <- 0  # starting value for phi
sigma2_samples[1] <- 1  # starting value for sigma^2

s2 <- 10  # variance for phi prior
alpha <- 0.01
beta <- 0.01

for (i in 2:niter) {
  # Update phi using Metropolis-Hastings
  phi_candidate <- rnorm(1, phi_samples[i-1], 0.1)
  acceptance_ratio <- exp((sum(dnorm(y[-1], mean=phi_candidate*y[-n], sd=sqrt(sigma2_samples[i-1]), log=TRUE)) + dnorm(phi_candidate, 0, sqrt(s2), log=TRUE)) - (sum(dnorm(y[-1], mean=phi_samples[i-1]*y[-n], sd=sqrt(sigma2_samples[i-1]), log=TRUE)) + dnorm(phi_samples[i-1], 0, sqrt(s2), log=TRUE)))
  if (runif(1) < acceptance_ratio) {
    phi_samples[i] <- phi_candidate
  } else {
    phi_samples[i] <- phi_samples[i-1]
  }
  
  # Update sigma^2 using inverse gamma distribution
  residuals <- y[-1] - phi_samples[i] * y[-n]
  alpha_star <- alpha + n/2
  beta_star <- beta + 0.5 * sum(residuals^2)
  sigma2_samples[i] <- 1 / rgamma(1, alpha_star, beta_star)
}

mean(phi_samples[2000:niter])
mean(sigma2_samples[2000:niter])

# Plot results
plot(phi_samples, type="l")
plot(sigma2_samples, type="l")
