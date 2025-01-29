install.packages("pacman")
pacman::p_load(R2jags, parallel)

set.seed(1998)

setwd('/work/KatFred/hier')


# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#----------getting th

#load young data
y_data <- read_csv("/work/KatFred/Andreas/data/young_IC.csv")
trial_data <- read_csv("/work/KatFred/Andreas/y_data_trial.csv")





##################################################################
#------------------------ DATA Preparation GOES HERE ----------
##################################################################
ID <- unique(y_data$ID)
nsubs <- length(ID)
nblocks <- length(unique(y_data$Block))
ndecks <- 4
ntrials_max <- 100
deck_cols <- c("Acount", "Bcount", "Ccount", "Dcount")

# all outcomes (NetScore)
NetScore_scaled <- y_data$NetScore/4 #scaled


#--- assign choices and outcomes in trial x sub matrix

# empty arrays to fill
x_all <- array(NA,c(nsubs,nblocks,blocksize))
counts_all <- array(NA,c(nsubs,nblocks,ndecks))
NetScore_all <- array(NA,c(nsubs,nblocks))

for (s in 1:nsubs) {
  
  #record n trials for subject s
  for (d in 1:ndecks) {
    counts_all[s,,d] <- unlist(y_data[y_data$ID==ID[s], deck_cols[d]], use.names = FALSE)
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
  
  NetScore_all[s,] <- NetScore_scaled[y_data$ID==ID[s]] 
  
}

# ---------------- score weight array ----------------------------
nblocks <- 5
ndecks <- 4
participant_ids <- unique(y_data$ID)
nsubs <- length(participant_ids)

#participants x blocks x decks
score_weight_all <- array(NA, dim = c(nsubs, nblocks, ndecks))

# Fill the score_weight array with proportions
for (i in seq_along(participant_ids)) {
  participant_id <- participant_ids[i]
  participant_data <- y_data[y_data$ID == participant_id, ]
  
  for (block in 1:nblocks) {
    block_data <- participant_data[participant_data$Block == block, ]
    counts <- as.numeric(block_data[c("Acount", "Bcount", "Ccount", "Dcount")])
    total_choices <- sum(counts)
    proportions <- counts / total_choices  # Calculate proportions
    score_weight_all[i, block, ] <- proportions
  }
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
                         model.file ="/work/KatFred/trust/jags_trust.txt", n.chains=3, 
                         n.iter=25000, n.burnin=5000, n.thin=1, n.cluster=3)

# let's look at the posteriors for the parameters
par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$a_rew))
plot(density(samples$BUGSoutput$sims.list$a_pun))
plot(density(samples$BUGSoutput$sims.list$theta))
plot(density(samples$BUGSoutput$sims.list$omega_f))
plot(density(samples$BUGSoutput$sims.list$T_weight))

#########################################################################
#---------- run the hierarchical model on the young population --------
########################################################################


nsubs <- 55
x <- x_all
NetScore <- NetScore_all
score_weight <- score_weight_all

nblocks <- 5

# set up jags and run jags model
data <- list("x","NetScore","score_weight", "nblocks", "blocksize", "nsubs") 

params<-c("mu_a_rew","mu_a_pun","mu_theta","mu_omega_f", "mu_T_weight") 

samples <- jags.parallel(data, inits=NULL, params,
                         model.file ="/work/KatFred/hier/F_hier.txt",
                         n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)

par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$mu_a_rew))
plot(density(samples$BUGSoutput$sims.list$mu_a_pun))
plot(density(samples$BUGSoutput$sims.list$mu_theta))
plot(density(samples$BUGSoutput$sims.list$mu_omega_f))
plot(density(samples$BUGSoutput$sims.list$mu_T_weight))

# xlim scaled to Haines et al.
par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$mu_a_rew), xlim=c(0,0.4))
plot(density(samples$BUGSoutput$sims.list$mu_a_pun), xlim=c(0,0.08))
plot(density(samples$BUGSoutput$sims.list$mu_theta), xlim=c(0,5))
plot(density(samples$BUGSoutput$sims.list$mu_omega_f), xlim=c(0,5))























