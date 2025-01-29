install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline)

set.seed(1998)

setwd('/work/KatFred/hier')

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#load control data
old_data <- read_csv("data/old_IC.csv")
old_trial <- 

young_data <- read_csv("data/young_IC.csv")
young_trial


############################## young data prep ###############################################
ID_old <- unique(old_data$ID)
nsubs <- length(ID_old)
nblocks <- length(unique(old_data$Block))
ndecks <- 4
ntrials_max <- 100
deck_cols <- c("Acount", "Bcount", "Ccount", "Dcount")

# all outcomes (NetScore)
NetScore_scaled <- old_data$NetScore/4 #scaled


#--- assign choices and outcomes in trial x sub matrix

# empty arrays to fill
x_all_old <- array(NA,c(nsubs,nblocks,blocksize))
counts_all_old <- array(NA,c(nsubs,nblocks,ndecks))
NetScore_all_old <- array(NA,c(nsubs,nblocks))

for (s in 1:nsubs) {
  
  #record n trials for subject s
  for (d in 1:ndecks) {
    counts_all_old[s,,d] <- unlist(old_data[old_data$ID==ID[s], deck_cols[d]], use.names = FALSE)
  }
  
  for (b in 1:nblocks) {
    onset <- 0
    offset <- 0
    for (d in 1:ndecks) {
      if (counts_all_old[s,b,d] != 0) {
        offset <- offset + counts_all_old[s,b,d]
        x_all_old[s,b, (1+onset):offset] <- d
        onset <- onset + counts_all_old[s,b,d]
      }
    }
  }
  
  NetScore_all_old[s,] <- NetScore_raw_old[old_data$ID==ID[s]] 
  
}




#----------prepare data for jags models - want trial x subject arrays for choice, gain, and loss ----
# identify and count unique subject IDs
subIDs_old<- unique(old_data$ID)
nsubs_old <- length(subIDs_old)

subIDs_young<- unique(young_data$ID)
nsubs_young <- length(subIDs_young)

ntrials_max <- 100

#------- the netscore is already scaled down -------------




