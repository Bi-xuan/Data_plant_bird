################################################################################
################################################################################
# JAGS CODE OF THE HIERARCHICAL STRUCTURAL EQUATION MODEL
################################################################################
if(! file.exists("include/jagsModel1.txt")) {
  jagsModel({
    for (k in 1:nP) {# k = structural equation models [1:nP]
      ############################################################################
      # Indicator model selection (Kuo & Mallick 1998) with random effect
      # - global shrinkage parameter tauB is estimated by the model
      # - the prior on the inclusion probability (gamma) is fixed at 0.5
      ############################################################################
      # Priors 
      sigB[k] ~ dunif(0, 20)
      tauB[k] <- pow(sigB[k], -2)
      for (j in 1:nY) {# j = response variables (Y) [1:nY]
        ##### fixed effects
        # intercepts
        betaT[1, j, k] ~ dnorm(0, 1e-4)
        gamma[1, j, k] <- 1
        beta[1, j, k] <- gamma[1, j, k] * betaT[1, j, k]
        # coefficients for main predictors
        for (x in 2:(nX[j] + 1)) {# x = predictor variables (X) [1:nX]
          gamma[x, j, k] ~ dbern(0.5)
          betaT[x, j, k] ~ dnorm(0, tauB[k])
          beta[x, j, k] <- gamma[x, j, k] * betaT[x, j, k]
        }
        ##### random effects
        for (r in 1:nRE) {# r = random effects (RE) [1:nRE]
          tau[r, j, k] <- pow(sig[r, j, k], -2)
          sig[r, j, k] ~ dunif(0, 20)
        }
        for (s in 1:nSite) {# s = random effect levels for site (siteRE) [1:nSite]
          siteRE[s, j, k] ~ dnorm(0, tau[1, j, k])
        }
        for (t in 1:nType) {# t = random effect levels for type (typeRE) [1:nType]
          typeRE[t, j, k] ~ dnorm(0, tau[2, j, k])
        }
        ##### Residual variance
        tau[3, j, k] <- 1/psiStar[j, j, k]
        ##### Explained variance
        # 1.) Variance explained by fixed effects
        # 2.) Variance explained by random effects
        rSq[1, j, k] <- pow(sd(mum[, j, k]), 2) /
          (pow(sd(mum[, j, k]), 2) + sum(pow(tau[, j, k], -1)))
        rSq[2, j, k] <- sum(pow(tau[1:nRE, j, k], -1)) /
          (pow(sd(mum[, j, k]), 2) + sum(pow(tau[, j, k], -1)))
        for (i in 1:nObs) {# i = observations [1:nObs]
          # Linear regressions
          mu[i, j, k] <- inprod(X[i, 1:(nX[j] + 1), k], beta[1:(nX[j] + 1), j, k]) +
            lambda[j, k] * eta[i, nE[j], k] + 
            siteRE[RE[i, 1], j, k] + typeRE[RE[i, 2], j, k]
          # likelihood
          Y[i, j, k] ~ dnorm(mu[i, j, k], tau[3, j, k])
          # marginal predictions based on fixed effects
          mum[i, j, k] <- inprod(X[i, 1:(nX[j] + 1), k], beta[1:(nX[j] + 1), j, k]) +
            lambda[j, k] * eta[i, nE[j], k]
        }
      }
      ##### Latent variables for covariance terms
      for (i in 1:nObs) {
        eta[i, 1, k] ~ dnorm(0, 1)
        eta[i, 2, k] ~ dnorm(0, 1)
      }
      ############################################################################
      # Parameter assignments & equality constraints for covariance terms
      # Code for implementing covariance terms is taken from R package 'blavaan'
      # Edgar C. Merkle and Yves Rosseel (2016). blavaan: Bayesian structural
      # equation models via parameter expansion.
      # arXiv 1511.05604. URL https://arxiv.org/abs/1511.05604
      ############################################################################
      lvRho[1, 2, k] <- -1 + 2 * covPars[5, k]
      lvRho[3, 4, k] <- -1 + 2 * covPars[6, k]
      psiStar[1, 1, k] <- psi[1, 1, k] - (sqrt(abs(lvRho[1, 2, k]) * psi[1, 1, k]))^2
      psiStar[2, 2, k] <- psi[2, 2, k] - ((-1 + 2 * step(lvRho[1, 2, k])) *
                                            sqrt(abs(lvRho[1, 2, k]) * psi[2, 2, k]))^2
      psiStar[3, 3, k] <- psi[3, 3, k] - (sqrt(abs(lvRho[3, 4, k]) * psi[3, 3, k]))^2
      psiStar[4, 4, k] <- psi[4, 4, k] - ((-1 + 2 * step(lvRho[3, 4, k])) *
                                            sqrt(abs(lvRho[3, 4, k]) * psi[4, 4, k]))^2
      lambda[1, k] <- sqrt(abs(lvRho[1, 2, k]) * psi[1, 1, k])
      lambda[2, k] <- (-1 + 2 * step(lvRho[1, 2, k])) *
        sqrt(abs(lvRho[1, 2, k]) * psi[2, 2, k])
      lambda[3, k] <- sqrt(abs(lvRho[3, 4, k]) * psi[3, 3, k])
      lambda[4, k] <- (-1 + 2 * step(lvRho[3, 4, k])) *
        sqrt(abs(lvRho[3, 4, k]) * psi[4, 4, k])
      # Inferential covariances
      psi[1, 2, k] <- lambda[1, k] * lambda[2, k]
      psi[3, 4, k] <- lambda[3, k] * lambda[4, k]
      # Priors
      for (j in 1:nY) {
        psi[j, j, k] <- pow(covPars[j, k], -1)
        covPars[j, k] ~ dgamma(1, .5)
      }
      covPars[5, k] ~ dbeta(1, 1)
      covPars[6, k] ~ dbeta(1, 1)
      ############################################################################
      # output for pathmodel
      ############################################################################
      # matrix with coefficients
      for (j in 1:nY) {
        mB[1:nX[j], j + 3, k] <- beta[2:(nX[j] + 1), j, k]
      }
      mB[5, 4, k] <- psi[1, 2, k]
      mB[4, 5, k] <- mB[5, 4, k]
      mB[6, 7, k] <- psi[3, 4, k]
      mB[7, 6, k] <- mB[6, 7, k]
      # matrix with selection probabilities
      for (j in 1:nY) {
        mP[1:nX[j], j + 3, k] <- gamma[2:(nX[j] + 1), j, k]
      }
      mP[5, 4, k] <- 1
      mP[4, 5, k] <- mP[5, 4, k]
      mP[6, 7, k] <- 1
      mP[7, 6, k] <- mP[6, 7, k]
    }
  }, file = "include/jagsModel1.txt")
}
################################################################################
# END OF SCRIPT
################################################################################
################################################################################
