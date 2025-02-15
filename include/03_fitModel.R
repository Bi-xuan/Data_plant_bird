################################################################################
################################################################################
# MODEL FITTING
################################################################################
if(! file.exists("data/jagsModelFit1.RData")) {
  cat(paste("Model is running for ", nIter + nAdpt, " iterations ...",
            "\n\nburn-in: ", nAdpt, " iterations",
            "\n\nposterior sampling: ", nIter, " iterations", sep = ""))
  ##############################################################################
  # 1.) Make cluster
  # 2.) Export the necessary R objects & functions
  # 3.) Run model
  ##############################################################################
  cluster1 <- makeCluster(detectCores())
  clusterExport(cl = cluster1, 
                varlist = c("jagsData1", "nAdpt", "nIter", "nThin", "initsFun"))
  postSamp1 <- 
    parLapply(cluster1, 1:nRuns,
              function(x) {
                ################################################################
                # 1.) load rjags
                # 2.) create initial values for chains
                # 3.) initialize model
                # 4.) sample from posterior distribution
                # 5.) stop the cluster
                ################################################################
                library(rjags)
                InitList1 <- 
                  initsFun(chain = x,
                           model = "include/jagsModel1.txt",
                           data = jagsData1)
                Init1 <- 
                  with(InitList1, {
                    jags.model(file = model, data = data, inits = inits,
                               n.chains = 1, n.adapt = nAdpt)
                  })
                postSamp1 <- 
                  coda.samples(model = Init1, 
                               variable.names = 
                                 c("mB", "mP", "rSq", "betaT", "beta",
                                   "tau", "typeRE", "siteRE", "lambda", "eta"),
                               n.iter = nIter,
                               thin = nThin)
                postSamp1 <- mcmc(postSamp1[[1]])
                return(postSamp1)
              }
    )
  stopCluster(cluster1)
  postSamp1 <- mcmc.list(postSamp1)
  ##############################################################################
  # save the results
  ##############################################################################
  save(list = c("postSamp1", "jagsData1"), file = "data/jagsModelFit1.RData")
} else {
  ##############################################################################
  # alternatively load, if already exists
  ##############################################################################
  load(file = "data/jagsModelFit1.RData")
}
################################################################################
# END OF SCRIPT
################################################################################
################################################################################

