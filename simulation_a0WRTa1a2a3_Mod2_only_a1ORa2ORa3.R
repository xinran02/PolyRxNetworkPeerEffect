# This script performs simulation for network peer-effect models.
# This script simulates the scenario when network peer-effect only transmits in inbound, outbound, or mutual dyads.
# And it compares the respective estimated peer-effects to a single peer-effect estimated using the influence matrix constructed from the undirected adjacency matrix 
# There are four sets of simulations:
# Network size n = 100, network density = 0.1 and 0.2
# Network size n = 500, network density = 0.05 and 0.1
# Assume Erdos-Renyi network

# This script only demonstrates when network peer-effect only transmits through mutual dyads (when only a3 is present)
# To simulate the scenarios when only a1 or a2 is present, may replace the corresponding parameters 

# True model: Yit = b0 + b1*Yi(t-1) + a3*(W_MUT(t-1)*Y(t-1))i + sigma
# Or true model: Yit = b0 + b1*Yi(t-1) + a2*(W_MUT(t-1)*Y(t-1))i + sigma
# Or true model: Yit = b0 + b1*Yi(t-1) + a1*(W_MUT(t-1)*Y(t-1))i + sigma
# Alternative model: Yit = b0 + b1*Yi(t-1) + a0*(W_UND(t-1)*Y(t-1))i + gamma

library(igraph)

# input network size and density
arg <- commandArgs(trailingOnly = TRUE) # used when running a Bash script in linux
nnodes <- as.numeric(arg[1]) # or manually specify
den <- as.numeric(arg[2])

n <- 1000 # number of runs

b0 <- 1
b1 <- 0.2

# setting a1 = 0, a2 = 0, and only a3 is present 
# can also set 1) only a1 is present, a2 = 0, a3 = 0; or 2) only a2 is present, a1 = 0, a3 = 0
a1 <- 0
a2 <- 0
a3 <- seq(-1,1,0.05)
values <- a3
nn <- length(values)

sigma <- sqrt(2)

# compare estimated peer-effect using an undirected influence matrix with the true peer-effect
# i.e., compare estimated a0 with a3
bias_mod2 <- rep(0, nn)
MSE_mod2 <- rep(0, nn)
cover_mod2 <- rep(0, nn)

for (r in 1:nn) {
  
  print(r)
  
  a3 <- values[r]
  
  bias_mod2_n <- rep(0, n)
  sqErr_mod2_n <- rep(0, n)
  coef_mod2_n <- rep(0, n)
  stdErr_mod2_n <- rep(0, n)
  
  for (q in 1:n) {
    # 1 - simulate random network
    net <- erdos.renyi.game(nnodes, den, type = "gnp", directed = T)
    adj <- as.matrix(as_adjacency_matrix(net))
    
    # 2 - make influence matrices
    adj_out <- matrix(rep(0, nnodes*nnodes), nrow = nnodes, ncol = nnodes)
    adj_in <- matrix(rep(0, nnodes*nnodes), nrow = nnodes, ncol = nnodes)
    adj_mut <- matrix(rep(0, nnodes*nnodes), nrow = nnodes, ncol = nnodes)
    
    for (i in 1:nnodes) {
      for (j in 1:nnodes) {
        adj_out[i,j] <- ifelse((adj[i,j] > 0 & adj[j,i] == 0), adj[i,j], 0)
        adj_in[i,j] <- ifelse((adj[i,j] == 0 & adj[j,i] > 0), adj[j,i], 0)
        adj_mut[i,j] <- ifelse((adj[i,j] > 0 & adj[j,i] > 0), adj[i,j] + adj[j,i], 0) 
      }
    }
    
    # turn influence matrices into row stochastic
    row_stoch <- function(mat, n) {
      mat <- t(apply(mat, 1, function(i) i/sum(i)))
      mat[is.nan(mat)] <- 1/(n-1)
      diag(mat) <- 0
      return(mat)
    }
    
    w_out <- row_stoch(adj_out, nnodes)
    w_in <- row_stoch(adj_in, nnodes)
    w_mut <- row_stoch(adj_mut, nnodes)
    
    # construct influence matrix if assuming network is undirected
    adj_und <- adj + t(adj)
    w_und <- t(apply(adj_und, 2, function(i) i/sum(i)))
    w_und[is.nan(w_und)]<- 1/(nnodes-1)
    diag(w_und) <- 0
    
    # 3 - construct peer variables and fit the models 
    lagY <- rnorm(nnodes, 0, 1)
    
    pvar_in <- w_in %*% lagY
    pvar_out <- w_out %*% lagY
    pvar_mut <- w_mut %*% lagY
    pvar_und <- w_und %*% lagY
    
    Y <- rnorm(nnodes, mean = (b0 + b1*lagY + a1*pvar_in + a2*pvar_out + a3*pvar_mut), sd=sigma)
    
    mod2 <- lm(Y ~ lagY + pvar_und) # fit using comparative model
    
    # estimated a0 compared to a3
    bias_mod2_n[q] <- coef(mod2)[3] - a3
    sqErr_mod2_n[q] <- (coef(mod2)[3] - a3)^2
    coef_mod2_n[q] <- coef(mod2)[3]
    stdErr_mod2_n[q] <- sqrt(diag(vcov(mod2)))[3]
  }
  
  # compute mean bias, MSE, and coverage
  bias_mod2[r] <- mean(bias_mod2_n)
  MSE_mod2[r] <- mean(sqErr_mod2_n)
  cover_mod2[r] <- mean(coef_mod2_n-1.96*stdErr_mod2_n<=a3 & coef_mod2_n+1.96*stdErr_mod2_n>=a3)
}

# plot the results
outdir <- "Path_Of_Out_Directory"

x <- seq(1:nn)

plotRes <- function(res, fn) {
  png(paste0(outdir, "netw_", as.character(nnodes), "_", as.character(den), "_", fn, ".png"), height=6,width=8,units="in",res=300)
  plot(x,res, ylab= fn, 
       xlab="Mutual Peer Variable Regression Coefficient")
  abline(lm(res ~ x))
  dev.off()
}

plotRes(bias_mod2, "a0WRTa3_bias_mod2_onlya3")
plotRes(MSE_mod2, "a0WRTa3_MSE_mod2_onlya3")
plotRes(cover_mod2, "a0WRTa3_cover_mod2_onlya3")







