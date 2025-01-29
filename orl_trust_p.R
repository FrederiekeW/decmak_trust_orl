ORL <- function(avg_block_payoff, blocksize, nblocks, a_rew, a_pun, theta, omega_f, T_weight) {
  
  # --------- Block-level variables -----------------
  x <- array(NA, c(nblocks, blocksize))
  NetScore <- array(0, nblocks)      # NetScore for each block
  
  Ev_update <- array(NA, c(nblocks, 4))
  Ev <- array(NA, c(nblocks, 4))
  
  signX <- rep(NA, nblocks)
  score_weight <- array(NA, c(nblocks, 4))
  Ef_cho <- array(NA, c(nblocks, 4))
  Ef_not <- array(NA, c(nblocks, 4))
  Ef <- array(NA, c(nblocks, 4))
  
  V <- array(NA, c(nblocks, 4))
  
  exp_p <- array(NA, c(nblocks, 4))
  p <- array(NA, c(nblocks, 4))
  
  # ----------- Initialize for the first block -----------------
  # bad and good decks, redundant but better for brain
  T_update1 <- 1 + T_weight
  T_update2 <- 1 - T_weight
  
  # Initialize probs
  p[1, 1] <- 0.25 * T_update1
  p[1, 2] <- 0.25 * T_update1
  p[1, 3] <- 0.25 * T_update2
  p[1, 4] <- 0.25 * T_update2
  
  # number of choices per deck
  x[1, ] <- rcat(blocksize, p[1, ])
  
  #Ev and Ef
  Ev[1, 1] <- 0.25 * T_update1
  Ev[1, 2] <- 0.25 * T_update1
  Ev[1, 3] <- 0.25 * T_update2
  Ev[1, 4] <- 0.25 * T_update2
  
  Ef[1, 1] <- 0.25 * T_update1
  Ef[1, 2] <- 0.25 * T_update1
  Ef[1, 3] <- 0.25 * T_update2
  Ef[1, 4] <- 0.25 * T_update2
  
  
  # Calculate NetScore for the first block
  NetScore[1] <- sum(x[1, ] > 2) / blocksize * avg_block_payoff[3] + 
    (blocksize - sum(x[1, ] > 2)) / blocksize * avg_block_payoff[2]

  
  # ------------------- Block stufff --------------------
  for (b in 2:nblocks) {
    signX[b] <- ifelse(NetScore[b - 1] < 0, -1, 1)
    
    for (d in 1:4) {
      
      # ------------------- Update expected values ----------------------
      score_weight[b, d] <- sum(x[b-1,] == d) / blocksize
      Ev_update[b, d] <- ifelse(NetScore[b - 1] >= 0,
                                Ev[b - 1, d] + a_rew * (NetScore[b - 1] * score_weight[b, d] - Ev[b - 1, d]),
                                Ev[b - 1, d] + a_pun * (NetScore[b - 1] * score_weight[b, d] - Ev[b - 1, d]))
      
      Ev[b, d] <- ifelse(sum(d == x[b - 1, ]) > 0, 
                         Ev_update[b, d],
                         Ev[b - 1, d])
      
      # ------------------- Update expected frequencies -------------------
      Ef_cho[b, d] <- ifelse(NetScore[b - 1] >= 0,
                             Ef[b - 1, d] + a_rew * (signX[b] * score_weight[b, d] - Ef[b - 1, d]),
                             Ef[b - 1, d] + a_pun * (signX[b] * score_weight[b, d] - Ef[b - 1, d]))
      
      Ef_not[b, d] <- ifelse(NetScore[b - 1] >= 0,
                             Ef[b - 1, d] + a_pun * (-(signX[b] / 3) * score_weight[b, d] - Ef[b - 1, d]),
                             Ef[b - 1, d] + a_rew * (-(signX[b] / 3) * score_weight[b, d] - Ef[b - 1, d]))
      
      Ef[b, d] <- ifelse(sum(d == x[b - 1, ]) > 0, 
                         Ef_cho[b, d],
                         Ef_not[b, d])
      
      # ---------------------- Valence model -------------------------------
      V[b, d] <- Ev[b, d] + Ef[b, d] * omega_f
      
      # ----------------------- Softmax part 1 ---------------------------------
      exp_p[b, d] <- exp(theta * V[b, d])
    }
    
    # ----------------------------- Softmax part 2 -------------------------------
    for (d in 1:4) {
      p[b, d] <- exp_p[b, d] / sum(exp_p[b, ])
    }
    
    # Make a choice for the block
    x[b, ] <- rcat(blocksize, p[b, ])
    
    # Calculate NetScore for the block
    NetScore[b] <- sum(x[b,] > 2) / blocksize * avg_block_payoff[3] + 
      (blocksize - sum(x[b,] > 2)) / blocksize * avg_block_payoff[2]
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
