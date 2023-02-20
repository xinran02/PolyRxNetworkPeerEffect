# This script performs simulation for network peer-effect models using directed influence matrices and an undirected influence matrix.
# This includes the following four sets of simulations:
# Network size n = 100, network density = (0.05, 0.1, 0.2)
# Network size n = 500, network density = (0.05, 0.1, 0.2)
# Assume Erdos-Renyi network

# True model: Yit = b0 + b1*Yi(t-1) + a1*(W_IN(t-1)*Y(t-1))i + a2*(W_OUT(t-1)*Y(t-1))i + a3*(W_MUT(t-1)*Y(t-1))i + sigma
# Alternative model: Yit = b0 + b1*Yi(t-1) + a0*(W_UND(t-1)*Y(t-1))i + gamma

library(igraph)

# input network size and density
arg <- commandArgs(trailingOnly = TRUE) # used when running a Bash script in linux
nnodes <- as.numeric(arg[1]) # or manually specify
den <- as.numeric(arg[2])

print(c(nnodes, den))

n <- 100 # number of runs

b0 <- 1
b1 <- 0.2

# varying a1, a2, and a3
a1 <- seq(-1,1,0.5)
a2 <- seq(-1,1,0.5)
a3 <- seq(-1,1,0.5)
values <- expand.grid(a1,a2,a3)
nn <- dim(values)[1] # number of unique values of (a1,a2,a3)

sigma <- sqrt(2)

a1_bias_mod1 <- rep(0, nn)
a1_MSE_mod1 <- rep(0, nn)
a1_cover_mod1 <- rep(0, nn)

a2_bias_mod1 <- rep(0, nn)
a2_MSE_mod1 <- rep(0, nn)
a2_cover_mod1 <- rep(0, nn)

a3_bias_mod1 <- rep(0, nn)
a3_MSE_mod1 <- rep(0, nn)
a3_cover_mod1 <- rep(0, nn)

bias_mod2_wrta1 <- rep(0, nn)
MSE_mod2_wrta1 <- rep(0, nn)
cover_mod2_wrta1 <- rep(0, nn)

bias_mod2_wrta2 <- rep(0, nn)
MSE_mod2_wrta2 <- rep(0, nn)
cover_mod2_wrta2 <- rep(0, nn)

bias_mod2_wrta3 <- rep(0, nn)
MSE_mod2_wrta3 <- rep(0, nn)
cover_mod2_wrta3 <- rep(0, nn)

# proportion of inbound, outbound, and mutual dyads
prop_dyads <- matrix(rep(0, 3*nn), nrow = nn, ncol = 3) # across unique values of (a1, a2, a3)
prop_dyads_n <- matrix(rep(0, 3*n), nrow = n, ncol = 3) # across simulated networks

for (r in 1:nn) {
  
  print(r)
  
  a1 <- values[r,1]
  a2 <- values[r,2]
  a3 <- values[r,3]
  
  a1_bias_mod1_n <- rep(0, n)
  a1_sqErr_mod1_n <- rep(0, n)
  a1_coef_mod1_n <- rep(0, n)
  a1_stdErr_mod1_n <- rep(0, n)
  
  a2_bias_mod1_n <- rep(0, n)
  a2_sqErr_mod1_n <- rep(0, n)
  a2_coef_mod1_n <- rep(0, n)
  a2_stdErr_mod1_n <- rep(0, n)
  
  a3_bias_mod1_n <- rep(0, n)
  a3_sqErr_mod1_n <- rep(0, n)
  a3_coef_mod1_n <- rep(0, n)
  a3_stdErr_mod1_n <- rep(0, n)
  
  bias_mod2_n_wrta1 <- rep(0, n)
  bias_mod2_n_wrta2 <- rep(0, n)
  bias_mod2_n_wrta3 <- rep(0, n)
  
  sqErr_mod2_n_wrta1 <- rep(0, n)
  sqErr_mod2_n_wrta2 <- rep(0, n)
  sqErr_mod2_n_wrta3 <- rep(0, n)
  
  coef_mod2_n <- rep(0, n)
  stdErr_mod2_n <- rep(0, n)
  
  for (q in 1:n) {
    # 1 - simulate random network
    net <- erdos.renyi.game(nnodes, den, type = "gnp", directed = T)
    adj <- as.matrix(as_adjacency_matrix(net))
    
    # 2 - make influence matrices
    adj_in <- matrix(rep(0, nnodes*nnodes), nrow = nnodes, ncol = nnodes)
    adj_out <- matrix(rep(0, nnodes*nnodes), nrow = nnodes, ncol = nnodes)
    adj_mut <- matrix(rep(0, nnodes*nnodes), nrow = nnodes, ncol = nnodes)
    
    for (i in 1:nnodes) {
      for (j in 1:nnodes) {
        adj_in[i,j] <- ifelse((adj[i,j] == 0 & adj[j,i] > 0), adj[j,i], 0)
        adj_out[i,j] <- ifelse((adj[i,j] > 0 & adj[j,i] == 0), adj[i,j], 0)
        adj_mut[i,j] <- ifelse((adj[i,j] > 0 & adj[j,i] > 0), adj[i,j] + adj[j,i], 0) 
      }
    }
    
    # prop. of dyads in simulated network
    n_in <- sum(apply(adj_in, c(1,2), function(x) x!=0))/2
    n_out <- sum(apply(adj_out, c(1,2), function(x) x!=0))/2
    n_mut <- sum(apply(adj_mut, c(1,2), function(x) x!=0))/2
    prop_dyads_n[q, c(1,2,3)] <- c(n_in,n_out,n_mut)/(sum(n_in,n_out,n_mut))
    
    # turn influence matrices into row stochastic
    row_stoch <- function(mat, n) {
      mat <- t(apply(mat, 1, function(i) i/sum(i)))
      mat[is.nan(mat)] <- 1/(n-1)
      diag(mat) <- 0
      return(mat)
    }
    
    w_in <- row_stoch(adj_in, nnodes)
    w_out <- row_stoch(adj_out, nnodes)
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
    
    mod1 <- lm(Y ~ lagY + pvar_in + pvar_out + pvar_mut) # true model
    mod2 <- lm(Y ~ lagY + pvar_und) # alternative model
    
    a1_bias_mod1_n[q] <- coef(mod1)[3] - a1
    a1_sqErr_mod1_n[q] <- (coef(mod1)[3] - a1)^2
    a1_coef_mod1_n[q] <- coef(mod1)[3]
    a1_stdErr_mod1_n[q] <- sqrt(diag(vcov(mod1)))[3]
    
    a2_bias_mod1_n[q] <- coef(mod1)[4] - a2
    a2_sqErr_mod1_n[q] <- (coef(mod1)[4] - a2)^2
    a2_coef_mod1_n[q] <- coef(mod1)[4]
    a2_stdErr_mod1_n[q] <- sqrt(diag(vcov(mod1)))[4]
    
    # in case network is sparse without mutual dyads
    a3_bias_mod1_n[q] <- coef(mod1)[5] - a3
    if (is.na(coef(mod1)[5] - a3)) {
      print(summary(mod1))
    }
    
    a3_sqErr_mod1_n[q] <- (coef(mod1)[5] - a3)^2
    a3_coef_mod1_n[q] <- coef(mod1)[5]
    a3_stdErr_mod1_n[q] <- sqrt(diag(vcov(mod1)))[5]
    
    # estimated a0 compared to each of a1, a2, and a3
    bias_mod2_n_wrta1[q] <- coef(mod2)[3] - a1
    sqErr_mod2_n_wrta1[q] <- (coef(mod2)[3] - a1)^2
    
    bias_mod2_n_wrta2[q] <- coef(mod2)[3] - a2
    sqErr_mod2_n_wrta2[q] <- (coef(mod2)[3] - a2)^2
    
    bias_mod2_n_wrta3[q] <- coef(mod2)[3] - a3
    sqErr_mod2_n_wrta3[q] <- (coef(mod2)[3] - a3)^2
    
    coef_mod2_n[q] <- coef(mod2)[3]
    stdErr_mod2_n[q] <- sqrt(diag(vcov(mod2)))[3]
  }
  
  # compute mean bias, MSE, and coverage
  a1_bias_mod1[r] <- mean(a1_bias_mod1_n)
  a1_MSE_mod1[r] <- mean(a1_sqErr_mod1_n)
  a1_cover_mod1[r] <- mean(a1_coef_mod1_n-1.96*a1_stdErr_mod1_n<=a1 & a1_coef_mod1_n+1.96*a1_stdErr_mod1_n>=a1)
  
  a2_bias_mod1[r] <- mean(a2_bias_mod1_n)
  a2_MSE_mod1[r] <- mean(a2_sqErr_mod1_n)
  a2_cover_mod1[r] <- mean(a2_coef_mod1_n-1.96*a2_stdErr_mod1_n<=a2 & a2_coef_mod1_n+1.96*a2_stdErr_mod1_n>=a2)
  
  a3_bias_mod1[r] <- mean(a3_bias_mod1_n)
  a3_MSE_mod1[r] <- mean(a3_sqErr_mod1_n)
  a3_cover_mod1[r] <- mean(a3_coef_mod1_n-1.96*a3_stdErr_mod1_n<=a3 & a3_coef_mod1_n+1.96*a3_stdErr_mod1_n>=a3)
  
  # a0
  bias_mod2_wrta1[r] <- mean(bias_mod2_n_wrta1)
  MSE_mod2_wrta1[r] <- mean(sqErr_mod2_n_wrta1)
  cover_mod2_wrta1[r] <- mean(coef_mod2_n-1.96*stdErr_mod2_n<=a1 & coef_mod2_n+1.96*stdErr_mod2_n>=a1)
  
  bias_mod2_wrta2[r] <- mean(bias_mod2_n_wrta2)
  MSE_mod2_wrta2[r] <- mean(sqErr_mod2_n_wrta2)
  cover_mod2_wrta2[r] <- mean(coef_mod2_n-1.96*stdErr_mod2_n<=a2 & coef_mod2_n+1.96*stdErr_mod2_n>=a2)
  
  bias_mod2_wrta3[r] <- mean(bias_mod2_n_wrta3)
  MSE_mod2_wrta3[r] <- mean(sqErr_mod2_n_wrta3)
  cover_mod2_wrta3[r] <- mean(coef_mod2_n-1.96*stdErr_mod2_n<=a3 & coef_mod2_n+1.96*stdErr_mod2_n>=a3)
  
  # prop. of dyads across unique values of (a1, a2, a3)
  prop_dyads[r, c(1,2,3)] <- apply(prop_dyads_n,2, function(x) mean(x))
}

# plot the results
# outdir <- "Path_Of_Out_Directory/"
outdir <- getwd()

x <- seq(1:nn)

# plotRes <- function(res, fn) {
#   png(paste0(outdir, "/netw_", as.character(nnodes), "_", as.character(den), "_", fn, ".png"), height=6,width=8,units="in",res=300)
#   plot(x,res, ylab= fn, 
#        xlab="Combinations of Directed and Mutual Peer Variables Regression Coefficients")
#   abline(lm(res ~ x))
#   dev.off()
# }

print(mean(a1_bias_mod1))
print(mean(a1_MSE_mod1))
print(mean(a1_cover_mod1))

print(mean(a2_bias_mod1))
print(mean(a2_MSE_mod1))
print(mean(a2_cover_mod1))

print(mean(a3_bias_mod1,na.rm=T))
print(mean(a3_MSE_mod1,na.rm=T))
print(mean(a3_cover_mod1,na.rm=T))


# plotRes(a1_bias_mod1, "a1_bias_mod1")
# plotRes(a1_MSE_mod1, "a1_MSE_mod1")
# plotRes(a1_cover_mod1, "a1_cover_mod1")
# 
# plotRes(a2_bias_mod1, "a2_bias_mod1")
# plotRes(a2_MSE_mod1, "a2_MSE_mod1")
# plotRes(a2_cover_mod1, "a2_cover_mod1")
# 
# plotRes(a3_bias_mod1, "a3_bias_mod1")
# plotRes(a3_MSE_mod1, "a3_MSE_mod1")
# plotRes(a3_cover_mod1, "a3_cover_mod1")

# plotRes(bias_mod2_wrta1, "a0_WRT_a1_bias_mod2")
# plotRes(MSE_mod2_wrta1, "a0_WRT_a1_MSE_mod2")
# plotRes(cover_mod2_wrta1, "a0_WRT_a1_cover_mod2")
# 
# plotRes(bias_mod2_wrta2, "a0_WRT_a2_bias_mod2")
# plotRes(MSE_mod2_wrta2, "a0_WRT_a2_MSE_mod2")
# plotRes(cover_mod2_wrta2, "a0_WRT_a2_cover_mod2")
# 
# plotRes(bias_mod2_wrta3, "a0_WRT_a3_bias_mod2")
# plotRes(MSE_mod2_wrta3, "a0_WRT_a3_MSE_mod2")
# plotRes(cover_mod2_wrta3, "a0_WRT_a3_cover_mod2")

# save results
writeFile <- function(res, fn) {
  write.csv(res, paste0(outdir, "/netw_", as.character(nnodes), "_", as.character(den), "_", fn, ".csv"), row.names = F)
}

writeFile(prop_dyads_n, "prop_dyads_n")
writeFile(prop_dyads, "prop_dyads")

writeFile(a1_bias_mod1, "a1_bias_mod1")
writeFile(a1_MSE_mod1, "a1_MSE_mod1")
writeFile(a1_cover_mod1, "a1_cover_mod1")

writeFile(a2_bias_mod1, "a2_bias_mod1")
writeFile(a2_MSE_mod1, "a2_MSE_mod1")
writeFile(a2_cover_mod1, "a2_cover_mod1")

writeFile(a3_bias_mod1, "a3_bias_mod1")
writeFile(a3_MSE_mod1, "a3_MSE_mod1")
writeFile(a3_cover_mod1, "a3_cover_mod1")

writeFile(bias_mod2_wrta1, "a0_WRT_a1_bias_mod2")
writeFile(MSE_mod2_wrta1, "a0_WRT_a1_MSE_mod2")
writeFile(cover_mod2_wrta1, "a0_WRT_a1_cover_mod2")

writeFile(bias_mod2_wrta2, "a0_WRT_a2_bias_mod2")
writeFile(MSE_mod2_wrta2, "a0_WRT_a2_MSE_mod2")
writeFile(cover_mod2_wrta2, "a0_WRT_a2_cover_mod2")

writeFile(bias_mod2_wrta3, "a0_WRT_a3_bias_mod2")
writeFile(MSE_mod2_wrta3, "a0_WRT_a3_MSE_mod2")
writeFile(cover_mod2_wrta3, "a0_WRT_a3_cover_mod2")


