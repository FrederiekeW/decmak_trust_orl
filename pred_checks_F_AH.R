install.packages("pacman")
pacman::p_load(R2jags, parallel, ggplot2)

set.seed(3006)


setwd('/work/KatFred/trust')

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#load young data
#trial_data <- read_csv("/work/KatFred/Andreas/y_data_trial.csv")
agg_data <- read_csv("/work/KatFred/group/young_IC.csv")



#----------prepare data for jags models - want trial x subject arrays for choice and gain & loss ----
# identify and count unique subject IDs
ID <- unique(agg_data$ID)
nsubs <- length(ID)
nblocks <- length(unique(agg_data$Block))
ndecks <- 4
ntrials_max <- 100
deck_cols <- c("Acount", "Bcount", "Ccount", "Dcount")

# all outcomes (NetScore)
NetScore_scaled <- agg_data$NetScore/4 #scaled


#--- assign choices and outcomes in trial x sub matrix

# empty arrays to fill
x_all <- array(NA,c(nsubs,nblocks,blocksize))
counts_all <- array(NA,c(nsubs,nblocks,ndecks))
NetScore_all <- array(NA,c(nsubs,nblocks))

for (s in 1:nsubs) {
  
  #record n trials for subject s
  for (d in 1:ndecks) {
    counts_all[s,,d] <- unlist(agg_data[agg_data$ID==ID[s], deck_cols[d]], use.names = FALSE)
  }
  
  for (b in 1:nblocks) {
    onset <- 0
    offset <- 0
    for (d in 1:ndecks) {
      if (counts_all[s,b,d] != 0) {
        offset <- offset + counts_all[s,b,d]
        x_all[s,b, (1+onset):offset] <- d
        onset <- onset + counts_all[s,b,d]
      }
    }
  }
  
  NetScore_all[s,] <- NetScore_raw[agg_data$ID==ID[s]] 
  
}

# ---------------------- set up jags ----------------------------------
#----------testing our data curation by running JAGS on one subject

# Now we'll fit one subject just to make sure everything works

x <- x_all[1,,]
NetScore <- NetScore_all[1,]

# set up jags and run jags model
# What's our "data"
data <- list("x","NetScore","score_weight", "nblocks", "blocksize") 

# What's our "parameters"
params <- c("a_rew", "a_pun", "theta", "omega_f", "T_weight", "p")


# How to call jags?
samples <- jags.parallel(data, inits=NULL, params,
                         model.file ="jags_trust.txt", n.chains=3, 
                         n.iter=25000, n.burnin=5000, n.thin=1, n.cluster=3)

# let's look at the posteriors for the parameters
par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$a_rew))
plot(density(samples$BUGSoutput$sims.list$a_pun))
plot(density(samples$BUGSoutput$sims.list$theta))
plot(density(samples$BUGSoutput$sims.list$omega_f))
plot(density(samples$BUGSoutput$sims.list$T_weight))


#----------Posterior predictive checks of descriptive accuracy

# Posterior prediction - start by looking at posteriors for p parameter

p_post <- samples$BUGSoutput$sims.list$p # probabilities as the outcome from softmax

#plot probability of each deck on trial 84
par(mfrow=c(2,2))
plot(density(p_post[,4,1]))
plot(density(p_post[,4,2]))
plot(density(p_post[,4,3]))
plot(density(p_post[,4,4]))

# which option will be chosen?
counts_all[1,4,]
# is this a good prediction?

# let's write a loop that loops and see how the model goes at predicting responses for all trials 
count_predict <- array(nblocks)

for (t in 1:nblocks) {
  
  p_predict <- c(
    MPD(p_post[,t,1]),
    MPD(p_post[,t,2]),
    MPD(p_post[,t,3]),
    MPD(p_post[,t,4])
  )
  
  count_predict[t] <- which.max(p_predict)
}
# how well did our model do?
for (b in 1:nblocks) {
  print(count_predict[b]==which(max(counts_all[1,b,])==counts_all[1,b,]))
}



# till here... 
for (b in 1:dim(p_post_mean)[1]) {
  
  # Predicted rank
  pred_rank <- rank(-p_post_mean[b, ], ties.method = "min")  # Use "min" for ties
  
  observed_rank <- rank(-score_weight[b, ], ties.method = "min")  # Use "min" for ties
  # Compare ranks
  matches <- sum(pred_rank == observed_rank)  # n matching
  block_accuracy[b] <- matches / length(pred_rank)  # Prop of correct
  
  # Counting successes (was model right?)
  if (all(pred_rank == observed_rank)) {
    rank_success <- rank_success + 1
  }
}

# Proportion of correctly ranked blocks --> for plot?
rank_prediction_success <- rank_success / dim(p_post_mean)[1]

# Print results
print(rank_prediction_success)
print(block_accuracy)





