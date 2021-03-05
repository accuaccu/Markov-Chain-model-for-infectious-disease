# Markov susceptible-exposed-infected-recovered-death-immunity (SEIRDM) model, with randomized transition matrix

# Machine for Random transition matrix P:
rand.P <- function(x,k){
  
 s_e <- runif(1, min = x[1], max = x[2]) # (user defined)
 s_s <- 1 - s_e
 
 e_i <- runif(1, min = x[3], max = x[4]) # (user defined)
 e_s <- 1 - e_i
 
 i_r <- k[1] # Recovery-rate (user defined)
 i_d <- 1 - i_r # Death-rate
 
 r_s <- k[2] # Back-to-susceptible-rate (user defined)
 r_m <- 1 - r_s # Immunity-rate
 
 d_d <- 1
 m_m <- 1
 
 P <- matrix(c(s_s,s_e,0,0,0,0,
               e_s,0,e_i,0,0,0,
               0,0,0,i_r,i_d,0,
               r_s,0,0,0,0,r_m,
               0,0,0,0,1,0,
               0,0,0,0,0,1),
             ncol = 6, byrow = T)
 
  if(sum(rowSums(P) == 1) == 6) {
    if(sum(P >= 0 & P <= 1) == 36) {
    return(P)
    }
    else {
      warning('Rowsums of transition matrix P are not equal to 1 or some of the elements are below 0 or above 1')
    }
  }
 
}


# Monte Carlo machine to repeat rand.P:
mc.mc <- function(n, pii, x, k, P, random.P) {
  
  Z <- matrix(NA,ncol = 6, nrow = n)
  Z[1,] <- pii
  if(random.P == T) {
    for(i in 1:(n-1)) {
      Z[i+1,] <- Z[i,] %*% rand.P(x,k)
    }
  }
  if(random.P == F)
    for(i in 1:(n-1)) {
      Z[i+1,] <- Z[i,] %*% P
    }
  
    
  return(Z)
  
}

# Example:
# Number of cycles:
b <- 1000
# Back-to-susceptible rate, i.e., when recovered, the probability of not getting immunity.
# If this value is high, the disease will stay longer in the population whereas if this value is low (more people gets immunity)
# the shorter the disease will stay in the population:
srate <- 0.5
# Define random transition limits (uni[a,b]):
rand_limits <- c(0.2,0.3,0.3,0.4)
# Define probability to recover from the disease:
constant_rates <- c(0.97,srate)
# Save MCMC with random transition matrix into new matrix:
pai <- mc.mc(n = b, pii = c(1,0,0,0,0,0), x = rand_limits, k = constant_rates, NA, T)
# Plot, for example:
plot(0,0,pch=NA,ylim = c(0,1), xlim = c(0,b), yaxt = 'n', xlab = 'Cycle number', ylab = 'Proportion');axis(2, las = 2)
lines(pai[,1], col = 1)
lines(pai[,2], col = 2)
lines(pai[,3], col = 3)
lines(pai[,4], col = 4)
lines(pai[,5], col = 5)
lines(pai[,6], col = 6)
legend('topright', legend = c('susceptible','exposed','infected','recovered','dead','immune'), pch = 175, col = 1:6)
title(adj = 0, paste0('Random Transition Matrix Markov Chain equilibrium state\n',
                      'with immunity-rate: ', 100*(1-srate),'% of the population'))


