hier_ORL_sim <- function(avg_block_payoff, ntrials, blocksize, nblocks,
                         mu_a_rew, mu_a_pun, mu_theta, mu_omega_f,
                         sigma_a_rew, sigma_a_pun, sigma_theta,
                         sigma_omega_f, mu_T_weight, sigma_T_weight) {
  
  # --------- Block-level variables -----------------
  # Arrays to populate for simulation
  x <- array(NA, c(nsubs, nblocks[1], blocksize))
  NetScore <- array(NA, c(nsubs, nblocks[1]))
  Ev <- array(NA, c(nsubs, nblocks[1], 4))
  score_weight <- array(NA, c(nsubs, nblocks[1], 4))
  
  for (s in 1:nsubs) {
    
    # Free parameters - sampled from a normal distribution with group mean and SD
    a_rew <- rtruncnorm(1, 0, , mu_a_rew, sigma_a_rew)
    a_pun <- rtruncnorm(1, 0, , mu_a_pun, sigma_a_pun)
    theta <- rtruncnorm(1, 0, , mu_theta, sigma_theta)
    omega_f <- rtruncnorm(1, 0, , mu_omega_f, sigma_omega_f)
    T_weight <- rtruncnorm(1, 0, mu_T_weight, sigma_T_weight)
    
    # Arrays to populate for simulation 
    Ev_update <- array(NA, c(nblocks[s], 4))
    signX <- array(NA, c(nblocks[s]))
    Ef_cho <- array(NA, c(nblocks[s], 4))
    Ef_not <- array(NA, c(nblocks[s], 4))
    Ef <- array(NA, c(nblocks[s], 4))
    V <- array(NA, c(nblocks[s], 4))
    exp_p <- array(NA, c(nblocks[s], 4))
    p <- array(NA, c(nblocks[s], 4))
    
    # ----------- Initialize for the first block -----------------
    # Bad and good decks, redundant but better for readability
    T_update1 <- 1 + T_weight
    T_update2 <- 1 - T_weight
    
    # Initialize probabilities
    p[1, 1] <- 0.25 * T_update1
    p[1, 2] <- 0.25 * T_update1
    p[1, 3] <- 0.25 * T_update2
    p[1, 4] <- 0.25 * T_update2
    
    # Number of choices per deck
    x[s, 1, ] <- rcat(blocksize, p[1, ])
    
    # Initialize Ev and Ef
    Ev[s, 1, 1] <- 0.25 * T_update1
    Ev[s, 1, 2] <- 0.25 * T_update1
    Ev[s, 1, 3] <- 0.25 * T_update2
    Ev[s, 1, 4] <- 0.25 * T_update2
    
    Ef[1, 1] <- 0.25 * T_update1
    Ef[1, 2] <- 0.25 * T_update1
    Ef[1, 3] <- 0.25 * T_update2
    Ef[1, 4] <- 0.25 * T_update2
    
    # Calculate NetScore for the first block
    NetScore[s, 1] <- sum(x[s, 1, ] > 2) / blocksize * avg_block_payoff[3] + 
      (blocksize - sum(x[s, 1, ] > 2)) / blocksize * avg_block_payoff[2]
    
    # ------------------- Block-specific updates --------------------
    for (b in 2:nblocks[s]) {
      signX[b] <- ifelse(NetScore[s, b - 1] < 0, -1, 1)
      
      for (d in 1:4) {
        
        # Update expected values
        score_weight[s, b, d] <- sum(x[s, b - 1, ] == d) / blocksize#
        
        Ev_update[b, d] <- ifelse(NetScore[s, b - 1] >= 0,
                                  Ev[s, b - 1, d] + a_rew * (NetScore[s, b - 1] * score_weight[s, b, d] - Ev[s, b - 1, d]),
                                  Ev[s, b - 1, d] + a_pun * (NetScore[s, b - 1] * score_weight[s, b, d] - Ev[s, b - 1, d]))
        
        Ev[s, b, d] <- ifelse(sum(d == x[s, b - 1, ]) > 0, 
                              Ev_update[b, d],
                              Ev[s, b - 1, d])
        
        # Update expected frequencies
        Ef_cho[b, d] <- ifelse(NetScore[s, b - 1] >= 0,
                               Ef[b - 1, d] + a_rew * (signX[b] * score_weight[s, b, d] - Ef[b - 1, d]),
                               Ef[b - 1, d] + a_pun * (signX[b] * score_weight[s, b, d] - Ef[b - 1, d]))
        
        Ef_not[b, d] <- ifelse(NetScore[s, b - 1] >= 0,
                               Ef[b - 1, d] + a_pun * (-(signX[b] / 3) * score_weight[s, b, d] - Ef[b - 1, d]),
                               Ef[b - 1, d] + a_rew * (-(signX[b] / 3) * score_weight[s, b, d] - Ef[b - 1, d]))
        
        Ef[b, d] <- ifelse(sum(d == x[s, b - 1, ]) > 0, 
                           Ef_cho[b, d],
                           Ef_not[b, d])
        
        # Valence model
        V[b, d] <- Ev[s, b, d] + Ef[b, d] * omega_f
        
        # Softmax part 1
        exp_p[b, d] <- exp(theta * V[b, d])
      }
      
      # Softmax part 2
      for (d in 1:4) {
        p[b, d] <- exp_p[b, d] / sum(exp_p[b, ])
      }
      
      # Make a choice for the block
      x[s, b, ] <- rcat(blocksize, p[b, ])
      
      # Calculate NetScore for the block
      NetScore[s, b] <- sum(x[s, b, ] > 2) / blocksize * avg_block_payoff[3] + 
        (blocksize - sum(x[s, b, ] > 2)) / blocksize * avg_block_payoff[2]
    }
  }
  
  # ------------------- Return results --------------------
  result <- list(
    x = x,
    NetScore = NetScore,
    score_weight = score_weight,
    Ev = Ev,
    Ef = Ef,
    V = V
  )
  
  return(result)
}

