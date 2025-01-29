install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)

set.seed(1998)

setwd('/work/KatFred/hier')

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

# NB! mod(ntrials, nstruct) (aka. ntrials %% nstruct) must be 0
ntrials <- 100 # total number of trials in our payoff structure
nstruct <- 20 # size of our subdivisions for pseudorandomization (one block)
freq <- 0.5 # probability of our frequent losses (we have losses half of the time)
infreq <- 0.1 # probability of our infrequent losses (we have losses 1/10th of the time)
bad_r <- 100 # "bad" winnings
bad_freq_l <- -250 # "bad" frequent loss
bad_infreq_l <- -1250 # "bad" infrequent loss
good_r <- 50 # "good" winnings
good_freq_l <- -50 # "good" frequent loss
good_infreq_l <- -250 # "good" infrequent loss

# Bad frequent
A_R <- rep(bad_r, nstruct) # we win on every trial
A_L <- c(rep(bad_freq_l, nstruct*freq), rep(0, nstruct*(1-freq))) # we have losses half of the time

# Bad infrequent
B_R <- rep(bad_r, nstruct)
B_L <- c(rep(bad_infreq_l, nstruct*infreq), rep(0, nstruct*(1-infreq))) # we have losses 1/10th of the time

# Good frequent
C_R <- rep(good_r, nstruct)
C_L <- c(rep(good_freq_l, nstruct*freq), rep(0, nstruct*(1-freq)))

# Good infrequent
D_R <- rep(good_r, nstruct)
D_L <- c(rep(good_infreq_l, nstruct*infreq), rep(0, nstruct*(1-infreq)))

# create the pseudorandomized full payoff structure
A <- array(NA, ntrials) # setting up an empty array to be filled
B <- array(NA, ntrials)
C <- array(NA, ntrials)
D <- array(NA, ntrials)
for (i in 1:(ntrials/nstruct)) {
  A[(1+(i-1)*nstruct):(i*nstruct)] <- (A_R + sample(A_L)) # randomly shuffling the loss-array for every block (and adding those losses to the winnings)
  B[(1+(i-1)*nstruct):(i*nstruct)] <- (B_R + sample(B_L))
  C[(1+(i-1)*nstruct):(i*nstruct)] <- (C_R + sample(C_L))
  D[(1+(i-1)*nstruct):(i*nstruct)] <- (D_R + sample(D_L))
}

payoff <- cbind(A, B, C, D)/100 # NB! SCALING! combining all four decks as columns with each 100 trials - dividing our payoffs by 100 to make the numbers a bit easier to work with

# let's look at the payoff
colSums(payoff) # the two bad decks should sum to -25 (i.e. -2500), and the two good ones to 25 (i.e. 2500)
#----------------------------------------------------

avg_block_payoff <- colSums(payoff[1:20,])

##################### Parameters #####################################
nblocks_all <- rep(5, 55)
nblocks <- 5
ndecks <- 4            
blocksize <- 20
nsubs <- 55
niterations <- 50

# mu
true_mu_a_rew <- array(NA,c(niterations))
true_mu_a_pun <- array(NA,c(niterations))
true_mu_theta <- array(NA,c(niterations))
true_mu_omega_f <- array(NA,c(niterations))
true_mu_T_weight <- array(NA,c(niterations))

infer_mu_a_rew <- array(NA,c(niterations))
infer_mu_a_pun <- array(NA,c(niterations))
infer_mu_theta <- array(NA,c(niterations))
infer_mu_omega_f <- array(NA,c(niterations))
infer_mu_T_weight <- array(NA,c(niterations))

# sigma (SD for R) / lambda (precision for JAGS)
true_lambda_a_rew <- array(NA,c(niterations))
true_lambda_a_pun <- array(NA,c(niterations))
true_lambda_theta <- array(NA,c(niterations))
true_lambda_omega_f <- array(NA,c(niterations))
true_lambda_T_weight <- array(NA,c(niterations))

infer_lambda_a_rew <- array(NA,c(niterations))
infer_lambda_a_pun <- array(NA,c(niterations))
infer_lambda_theta <- array(NA,c(niterations))
infer_lambda_omega_f <- array(NA,c(niterations))
infer_lambda_T_weight <- array(NA,c(niterations))


start_time = Sys.time()
for (i in 1:niterations) {
  #nblocks <- nblocks_all
  
  
  # let's see how robust the model is. Does it recover all sorts of values?
  mu_a_rew <- runif(1,0,1)
  mu_a_pun <- runif(1,0,1)
  mu_theta <- runif(1,.2,2) # could also just be a set value (e.g. 1) to simplify the model a bit
  mu_omega_f <- runif(1,-2,2)
  mu_T_weight <- runif(1,0,1)
  
  sigma_a_rew <- runif(1,0,0.1)
  sigma_a_pun <- runif(1,0,0.1)
  sigma_theta <- runif(1,0,0.2) # if theta is just a set value (e.g. 1), then this isn't relevant anymore
  sigma_omega_f <- runif(1,0,0.4)
  sigma_T_weight <- runif(1,0,0.3)
  
  # sigma_a_rew <- runif(1,0,.5)
  # sigma_a_pun <- runif(1,0,.5)
  # sigma_K <- runif(1,0,.5)
  # sigma_theta <- runif(1,0,.5)
  # sigma_omega_f <- runif(1,0,.5)
  # sigma_omega_p <- runif(1,0,.5)
  
  source('test.R')
  ORL_sims <- hier_ORL_sim(avg_block_payoff, ntrials, blocksize, nblocks_all,
                           mu_a_rew,mu_a_pun,mu_theta,mu_omega_f,
                           sigma_a_rew,sigma_a_pun,sigma_theta,
                           sigma_omega_f, mu_T_weight, sigma_T_weight)
  
  x <- ORL_sims$x
  NetScore <- ORL_sims$NetScore
  score_weight <- ORL_sims$score_weight
  T_weight <- ORL_sims$T_weight
  
  data <- list("x","NetScore","score_weight", "nblocks", "nsubs", 
               "blocksize", "avg_block_payoff", "T_weight") 
  
  params <- c("mu_a_rew","mu_a_pun","mu_theta","mu_omega_f",
              "mu_T_weight", "lambda_a_rew","lambda_a_pun","lambda_theta",
              "lambda_omega_f", "lambda_T_weight")

  samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="F_hier.txt", n.chains=3, 
                           n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)
  
  # mu
  true_mu_a_rew[i] <- mu_a_rew
  true_mu_a_pun[i] <- mu_a_pun
  true_mu_theta[i] <- mu_theta
  true_mu_omega_f[i] <- mu_omega_f
  true_mu_T_weight[i] <- mu_T_weight
  
  # find maximum a posteriori
  Y <- samples$BUGSoutput$sims.list
  infer_mu_a_rew[i] <- MPD(Y$mu_a_rew)
  infer_mu_a_pun[i] <- MPD(Y$mu_a_pun)
  infer_mu_theta[i] <- MPD(Y$mu_theta)
  infer_mu_omega_f[i] <- MPD(Y$mu_omega_f)
  infer_mu_T_weight[i] <- MPD(Y$mu_T_weight)
  
  # lambda
  true_lambda_a_rew[i] <- sigma_a_rew
  true_lambda_a_pun[i] <- sigma_a_pun
  true_lambda_theta[i] <- sigma_theta
  true_lambda_omega_f[i] <- sigma_omega_f
  true_lambda_T_weight[i] <- sigma_T_weight
  
  # find maximum a posteriori
  infer_lambda_a_rew[i] <- MPD(Y$lambda_a_rew)
  infer_lambda_a_pun[i] <- MPD(Y$lambda_a_pun)
  infer_lambda_theta[i] <- MPD(Y$lambda_theta)
  infer_lambda_omega_f[i] <- MPD(Y$lambda_omega_f)
  infer_lambda_T_weight[i] <- MPD(Y$lambda_T_weight)
  
  print(i)
  
}

end_time = Sys.time()
end_time - start_time

# let's look at some scatter plots
# plotting code courtesy of Lasse
source('recov_plot.R')
pl1 <- recov_plot(true_mu_a_rew, infer_mu_a_rew, c("true mu_a_rew", "infer mu_a_rew"), 'smoothed linear fit')
pl2 <- recov_plot(true_mu_a_pun, infer_mu_a_pun, c("true mu_a_pun", "infer mu_a_pun"), 'smoothed linear fit')
pl4 <- recov_plot(true_mu_theta, infer_mu_theta, c("true mu_theta", "infer mu_theta"), 'smoothed linear fit')
pl5 <- recov_plot(true_mu_omega_f, infer_mu_omega_f, c("true mu_omega_f", "infer mu_omega_f"), 'smoothed linear fit')
pl6 <- recov_plot(true_mu_T_weight, infer_mu_T_weight, c("true mu_T_weight", "infer mu_T_weight"), 'smoothed linear fit')

mu_plots <- ggarrange(pl1, pl2, pl4, pl5, pl6, ncol = 1, nrow = 5)

ggsave("mu_recovery.png", mu_plots, width = 8, height = 20)
ggarrange(pl1, pl2, pl4, pl5, pl6)

pl1 <- recov_plot(true_lambda_a_rew, infer_lambda_a_rew, c("true lambda_a_rew", "infer lambda_a_rew"), 'smoothed linear fit')
pl2 <- recov_plot(true_lambda_a_pun, infer_lambda_a_pun, c("true lambda_a_pun", "infer lambda_a_pun"), 'smoothed linear fit')
pl4 <- recov_plot(true_lambda_theta, infer_lambda_theta, c("true lambda_theta", "infer lambda_theta"), 'smoothed linear fit')
pl5 <- recov_plot(true_lambda_omega_f, infer_lambda_omega_f, c("true lambda_omega_f", "infer lambda_omega_f"), 'smoothed linear fit')
pl6 <- recov_plot(true_lambda_T_weight, infer_lambda_T_weight, c("true lambda_T_weight", "infer lambda_T_weight"), 'smoothed linear fit')

lambda_plots <- ggarrange(pl1, pl2, pl4, pl5, pl6, ncol = 1, nrow = 5)
ggsave("lambda_recovery.png", lambda_plots, width = 8, height = 20)




trace1 <- traceplot(samples, c(2, 2))
dev.off()
traceplot(samples, c(2, 3))

trace2
dev.off()



