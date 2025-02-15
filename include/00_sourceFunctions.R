################################################################################
################################################################################
# SOURCE FILE WITH ALL IMPORTANT FUNCTIONS USED IN THE ANALYSIS
################################################################################
# Function to calculate shannon diversity
################################################################################
diversityFun <- 
  function(x, MARGIN = 1) {
    if (length(dim(x)) > 1) {
      total <- apply(x, MARGIN = MARGIN, FUN = sum)
      x <- sweep(x, MARGIN = MARGIN, total, "/")
    } else {
      x <- x/(total <- sum(x))
    }
    x <- -x * log(x)
    if (length(dim(x)) > 1) {
      H <- apply(x, MARGIN = MARGIN, FUN = sum, na.rm = TRUE)
    } else {
      H <- sum(x, na.rm = TRUE)
    }
    return(H)
  }
################################################################################
# Function to calculate functional dispersion
################################################################################
fdFun <- 
  function(x, which, commStand = NULL, traitScale = NULL, disMethod = "gower") {
    comm <- switch(which,
                     p = xtabs(frequency ~ plot_date_type + plant_species, data = x$N),
                     a = xtabs(frequency ~ plot_date_type + animal_species, data = x$N))
    traits <- switch(which,
                     p = sqrt(x$pT[, 2:4]),
                     a = sqrt(x$aT[, 2:4]))
    comm <- comm[, rownames(traits)]
    if (!is.null(commStand)) {
      comm <- vegan::decostand(comm, MARGIN = 1, method = commStand)
    }
    if (!is.null(traitScale)) {
      traits <- vegan::decostand(traits, MARGIN = 2, method = traitScale)
    }
    select <- list("All" = 1:3, "1" = 1, "2" = 2, "3" = 3)
    traitList <- lapply(1:4, FUN = function(x) subset(traits, select = select[[x]]))
    disList <- lapply(traitList, FUN = function(x) vegdist(x, method = disMethod))
    FDisList <- lapply(disList, FUN = function(x) fdisp(d = x, a = comm)$FDis)
    out <- do.call(cbind, FDisList)
    colnames(out) <- paste("FDis", names(select), toupper(which), sep = "")
    return(out)
  }
################################################################################
# function to generate initial values for JAGS models
################################################################################
initsFun <- function(chain, model, data) {
  load.module("lecuyer")
  list(model = model,
       data = data,
       inits = with(data, {
         list(sigB = runif(nP, 0, 1),
              covPars = array(runif(6 * nP, 0, 1), dim = c(6, nP)),
              sig = array(runif(nRE * nY * nP, 0, 1), dim = c(nRE, nY, nP)),
              eta = array(runif(nObs * 2 * nP, -2, 2), dim = c(nObs, 2, nP)),
              betaT = {
                m <- array(runif((nX[length(nX)] + 1) * nY * nP, -2, 2),
                           dim = c((nX[length(nX)] + 1), nY, nP))
                m[5:6, 1:2, 1:nP] <- NA
                m
              },
              .RNG.name = "lecuyer::RngStream",
              .RNG.state = parallel.seeds("lecuyer::RngStream", 1)[[1]]$.RNG.state)
       }))
}
################################################################################
# function to set up model for input in jags
################################################################################
jagsModel <-
  function(model, file) {
    cl <- match.call()
    model <- cl$model
    cat("", file = file, append = FALSE)
    cat("model ", file = file, append = TRUE)
    cat(deparse(model, width.cutoff = 200L), sep = "\n", file = file, append = TRUE)
    invisible()
  }
################################################################################
# function to calculate  of detecting K of N tested relationships at
# a given significance level alpha by chance
# Moran, M. D. (2003). Arguments for rejecting the sequential Bonferroni in
# ecological studies. Oikos, 100, 403–405.
# https://doi.org/10.1034/j.1600-0706.2003.12010.x
################################################################################
moranFun <-
  function(K, N, alpha) {
    (factorial(N)/(factorial(N - K) * factorial(K))) * alpha^K * (1 - alpha)^(N - K)
  }
################################################################################
# Function to extract trait-specific results from RLQ-fourth-corner analysis
################################################################################
rlqExtract <- 
  function(data) {
    rlqData <- 
      rbind(data.frame(type = rep(c("Bird-fruit", "Bird-flower", "Insect-flower"), each = 6),
                       trophic_level = "plants",
                       trait_type = rep(c("Matching", "Energy", "Foraging"), 6),
                       trait = sapply(data, FUN = function(x) 
                         unlist(strsplit(x$rAxes$tabD[6]$names,
                                         split = " / ", 
                                         fixed = TRUE)))[which(1:36 %% 2 == 1)],
                       axis = substr(sapply(data, FUN = function(x) 
                         unlist(strsplit(x$rAxes$tabD[6]$names, 
                                         split = " / ", 
                                         fixed = TRUE)))[which(1:36 %% 2 != 1)], 5, 5),
                       r = c(sapply(data, FUN = function(x) x$rAxes$tabD[1]$obs)),
                       n = rep(sapply(data, FUN = function(x) length(unique(x$N$plant_species))), each = 6),
                       p = c(sapply(data, FUN = function(x) x$rAxes$tabD[7]$pvalue)),
                       col = rep(c("gray", "deepskyblue1", "firebrick1"), each = 6)),
            data.frame(type = rep(c("Bird-fruit", "Bird-flower", "Insect-flower"), each = 6),
                       trophic_level = "animals",
                       trait_type = rep(rep(c("Matching", "Energy", "Foraging"), each = 2), 3),
                       trait = sapply(data, FUN = function(x) 
                         unlist(strsplit(x$qAxes$tabD[6]$names,
                                         split = " / ",
                                         fixed = TRUE)))[which(1:36 %% 2 != 1)],
                       axis = substr(sapply(data, FUN = function(x) 
                         unlist(strsplit(x$qAxes$tabD[6]$names,
                                         split = " / ",
                                         fixed = TRUE)))[which(1:36 %% 2 == 1)], 5, 5),
                       r = c(sapply(data, FUN = function(x) x$qAxes$tabD[1]$obs)),
                       n = rep(sapply(data, FUN = function(x) length(unique(x$N$animal_species))), each = 6),
                       p = c(sapply(data, FUN = function(x) x$qAxes$tabD[7]$pvalue)),
                       col = rep(c("gray", "deepskyblue1", "firebrick1"), each = 6)))
    rlqData$trait_type <- factor(rlqData$trait_type, levels = c("Matching", "Energy", "Foraging"))
    rlqData[, c(6, 8)] <- signif(rlqData[, c(6, 8)], 2)
    
    rlqDataMean <- aggregate(p < 0.05 ~ trait_type + axis, sum, data = rlqData)
    names(rlqDataMean)[3] <- "k"
    rlqDataMean$n <- aggregate(abs(r) ~ trait_type + axis, length, data = rlqData)[, 3]
    rlqDataMean$'P-value' <- with(rlqDataMean, moranFun(k, n, 0.05))
    rlqDataMean$abs.r <- aggregate(abs(r) ~ trait_type + axis, mean, data = rlqData)[, 3]
    rlqDataMean$sd <- aggregate(abs(r) ~ trait_type + axis, sd, data = rlqData)[, 3]
    rlqDataMean$se <- with(rlqDataMean, sd/sqrt(n))
    rlqDataMean[, 5:8] <- signif(rlqDataMean[, 5:8], 2)
    
    rlqGlobal <- sapply(data, FUN = function(x) signif(c(x$trRLQ$obs, x$trRLQ$pvalue), 2))
    dimnames(rlqGlobal) <- list(stat = c("sum(eigenvalues)", "P-value"),
                                type = c("Bird-fruit", "Bird-flower", "Insect-flower"))
    
    rlqAxes <- lapply(dataSet1, FUN = function(x) {
      with(x$fcRLQ$tabD, {
        m <- signif(cbind(obs[c(1, 4)], pvalue[c(1, 4)]), 2)
        dimnames(m) <- list(axes = names[c(1, 4)], stat = c("r", "P-value"))
        return(m)
      })
    })
    names(rlqAxes) <- c("Bird-fruit", "Bird-flower", "Insect-flower")
    
    rlqEigenvalues <- sapply(dataSet1, FUN = function(x) signif(x$RLQ$eig/sum(x$RLQ$eig), 2))
    dimnames(rlqEigenvalues) <- 
      list(Axis = 1:3, type = c("Bird-fruit", "Bird-flower", "Insect-flower"))
    
    return(list(rlqGlobal = rlqGlobal, rlqAxes = rlqAxes,
                rlqEigenvalues = rlqEigenvalues,
                rlqData = rlqData, rlqDataMean = rlqDataMean))  
  }
################################################################################
# Function to extract parameter estimates of path analysis
################################################################################
pathExtract <- 
  function(mcmcList) {
    mcmcPooled <- mcmc(do.call(rbind, mcmcList))
    b <- apply(mcmcPooled[, grep("mB", colnames(mcmcPooled))],
               MARGIN = 2, FUN = function(x) quantile(x, prob = c(0.5, 0.025, 0.975)))
    p <- apply(mcmcPooled[, grep("mP", colnames(mcmcPooled))],
               MARGIN = 2, FUN = mean)
    psrf <- gelman.diag(mcmcList[, grep("mB[", colnames(mcmcPooled), fixed = TRUE)],
                        transform = TRUE, autoburnin = FALSE, multivariate = FALSE)$psrf[, 2]
    effn <- effectiveSize(mcmcList[, grep("mB[", colnames(mcmcPooled), fixed = TRUE)])
    rSq <- array(apply(mcmcPooled[, grep("rSq", colnames(mcmcPooled))],
                       MARGIN = 2, FUN = median),
                 dim = c(2, 4, 8),
                 dimnames = list(type = c("m", "c"),
                                 response = c("FDp", "FDa", "IDp", "IDa"),
                                 model_index = c("FDisAll_eH", "FDisAll_d", "FDis1_eH", "FDis1_d",
                                                 "FDis2_eH", "FDis2_d", "FDis3_eH", "FDis3_d")))
    sig <- array(apply(mcmcPooled[, grep("tau", colnames(mcmcPooled))],
                       MARGIN = 2, FUN = function(x) median(1 / x^2)),
                 dim = c(3, 4, 8),
                 dimnames = list(type = c("site", "type", "residual"),
                                 response = c("FDp", "FDa", "IDp", "IDa"),
                                 model_index = c("FDisAll_eH", "FDisAll_d", "FDis1_eH", "FDis1_d",
                                                 "FDis2_eH", "FDis2_d", "FDis3_eH", "FDis3_d")))
    mb <- ml <- mu <- mp <- mpsrf <- meffn <- 
      array(NA, dim = c(7, 7, 8),
            dimnames = list(a = c("LU", "MAP", "MAT", "FDp", "FDa", "IDp", "IDa"),
                            b = c("LU", "MAP", "MAT", "FDp", "FDa", "IDp", "IDa"),
                            c = c("FDisAll_eH", "FDisAll_d", "FDis1_eH", "FDis1_d",
                                  "FDis2_eH", "FDis2_d", "FDis3_eH", "FDis3_d")))
    
    eval(parse(text = paste("mb", substr(colnames(b), 3, 9), " <- ", "b[1, ", 1:ncol(b), "]", sep = "")))
    eval(parse(text = paste("ml", substr(colnames(b), 3, 9), " <- ", "b[2, ", 1:ncol(b), "]", sep = "")))
    eval(parse(text = paste("mu", substr(colnames(b), 3, 9), " <- ", "b[3, ", 1:ncol(b), "]", sep = "")))
    eval(parse(text = paste("mp", substr(names(p), 3, 9), " <- ", "p[", 1:length(p), "]", sep = "")))
    eval(parse(text = paste("mpsrf", substr(names(psrf), 3, 9), " <- ", "psrf[", 1:length(psrf), "]", sep = "")))
    eval(parse(text = paste("meffn", substr(names(effn), 3, 9), " <- ", "effn[", 1:length(effn), "]", sep = "")))
    
    mbf <- 2 * log(1 / ((1 - mp) / mp))
    
    mTable <- lapply(1:8, FUN = function(x) {
      mTable <- expand.grid(rownames(mb[, , x]), colnames(mb[, , x]))
      names(mTable) <- c("x", "y")
      mTable$effect <- paste(signif(mb[, , x], 2), 
                             " (", signif(ml[, , x], 2),", ", 
                             signif(mu[, , x], 2), ")", 
                             sep = "")
      mTable$pSel <- c(signif(mp[, , x], 2))
      mTable$bf <- c(signif(mbf[, , x], 2))
      mTable$psrf <- c(signif(mpsrf[, , x], 4))
      mTable$effn <- c(signif(meffn[, , x], 4))
      mTable <- mTable[complete.cases(mTable), ]
      select <- c(NA, 1:3, 5:7, 9:13, 15:19, NA, 8, 20)
      mTable <- mTable[select, ]
      return(mTable)
    })
    names(mTable) <- dimnames(mb)[[3]]
    
    mc <- mp[, , 1]
    mc[c(22:24, 26, 29:31, 32, 36:40, 42, 43:47, 48)] <- 
      c(rep("deepskyblue1", 3), gray(0.5),
        rep("firebrick1", 3), rep(gray(0.5), 4),
        "deepskyblue1", "firebrick1", rep(gray(0.5), 4),
        "deepskyblue1", "firebrick1", gray(0.5))
    
    vars <- list(list("LU", "MAP", "MAT", expression(italic(FD[p])), expression(italic(FD[a])),
                      expression(italic(e**H[p])), expression(italic(e**H[a]))),
                 list("LU", "MAP", "MAT", expression(italic(FD[p])), expression(italic(FD[a])),
                      expression(italic(d*"'"[p])), expression(italic(d*"'"[a]))))
    vars <- c(vars, vars, vars, vars)
    cols <- c(rep("black", 3), "deepskyblue1", "firebrick1", "deepskyblue1", "firebrick1")
    
    nPos <- matrix(c(rep(1:3, each = 3), rep(1:3,  3)), ncol = 2)
    nPos <- nPos[-c(5, 8), ]
    nPos[1, 2] <- nPos[1, 2] + 0.25 
    nPos[3, 2] <- nPos[3, 2] - 0.25
    nPos[6, 2] <- nPos[6, 2] + 0.5 
    nPos[7, 2] <- nPos[7, 2] - 0.5
    rPos1 <- matrix(c(0.2, -1, 0.2, 1.15, 1.1, -0.75, 1.1, 0.9), 4, 2, byrow = TRUE)
    rPos2 <- matrix(c(0.2, -1.15, 0.2, 1, 1.1, -0.9, 1.1, 0.75), 4, 2, byrow = TRUE)
    mb[is.na(mb)] <- 0
    
    varsFull <- list("LU", "MAP", "MAT", expression(italic(FD[p])), expression(italic(FD[a])),
                     expression(italic(ID[p])), expression(italic(ID[a])))
    mFull <- {
      mFull <- matrix(NA, 7, 7)
      mFull[upper.tri(mFull)] <- 1
      mFull[1:2, 2:3] <- NA
      mFull[5, 4] <- 1
      mFull[7, 6] <- 1
      mFull
    }
    out <- list(mb = mb, ml = ml, mu = mu, mp = mp, mbf = mbf, mFull = mFull,
                mTable = mTable,
                rSq = rSq, rSqCum = apply(rSq, c(2, 3), cumsum), sig = sig,
                mc = mc, vars = vars, varsFull = varsFull, cols = cols,
                nPos = nPos, rPos1 = rPos1, rPos2 = rPos2)
    return(out)
  }
################################################################################
# Function to extract partial residuals from path analysis
################################################################################
residExtract <- 
  function(mcmcList, data) {
    mcmcPooled <- mcmc(do.call(rbind, mcmcList))
    mCoef <- array(NA, dim = c(9, 4, 8),
                   dimnames = list(coef = c("Int", "LU", "MAP", "MAT", "FDp", "FDa", "lambda", "siteRE", "typeRE"),
                                   response = c("FDp", "FDa", "IDp", "IDa"),
                                   model_index = c("FDisAll_eH", "FDisAll_d", "FDis1_eH", "FDis1_d",
                                                   "FDis2_eH", "FDis2_d", "FDis3_eH", "FDis3_d")))
    mBeta <- apply(mcmcPooled[, grep("beta", colnames(mcmcPooled))[!(grep("beta", colnames(mcmcPooled))
                                                                     %in% grep("betaT", colnames(mcmcPooled)))]],
                   MARGIN = 2, FUN = median)
    eval(parse(text = paste("mCoef", substr(names(mBeta), 5, 11),
                            " <- ", "mBeta[", 1:length(mBeta), "]", sep = "")))
    mCoef[7, , ] <- apply(mcmcPooled[, grep("lambda", colnames(mcmcPooled))],
                          MARGIN = 2, FUN = median)
    mCoef[8:9, , ] <- 1
    mCoef[5:6, 1:2, ] <- 0
    
    eta <- matrix(apply(mcmcPooled[, grep("eta", colnames(mcmcPooled))[!(grep("eta", colnames(mcmcPooled))
                                                                         %in% grep("beta", colnames(mcmcPooled)))]],
                        MARGIN = 2, FUN = median),
                  ncol = 2 * 8, byrow = FALSE)
    siteRE <- matrix(apply(mcmcPooled[, grep("siteRE", colnames(mcmcPooled))],
                           MARGIN = 2, FUN = median),
                     ncol = 4 * 8, byrow = FALSE)
    typeRE <- matrix(apply(mcmcPooled[, grep("typeRE", colnames(mcmcPooled))],
                           MARGIN = 2, FUN = median),
                     ncol = 4 * 8, byrow = FALSE)
    site <- siteRE[data$RE[, 1], ]
    type <- typeRE[data$RE[, 2], ]
    
    X1 <- as.matrix(cbind(data$X[, , 1], eta[, 1], site[, 1], type[, 1]))
    X2 <- as.matrix(cbind(data$X[, , 1], eta[, 1], site[, 2], type[, 2]))
    X3 <- as.matrix(cbind(data$X[, , 1], eta[, 2], site[, 3], type[, 3]))
    X4 <- as.matrix(cbind(data$X[, , 1], eta[, 2], site[, 4], type[, 4]))
    X5 <- as.matrix(cbind(data$X[, , 1], eta[, 3], site[, 5], type[, 5]))
    X6 <- as.matrix(cbind(data$X[, , 1], eta[, 3], site[, 6], type[, 6]))
    X7 <- as.matrix(cbind(data$X[, , 1], eta[, 4], site[, 7], type[, 7]))
    X8 <- as.matrix(cbind(data$X[, , 1], eta[, 4], site[, 8], type[, 8]))
    Xall <- sapply(mget(paste("X", 1:8, sep = "")), identity, simplify = "array")
    dimnames(Xall)[[2]][7:9] <- c("eta", "site", "type")
    
    rFDisP <- data$Y[, 1, 1] - Xall[, -3, 1] %*% mCoef[-3, 1, 1]
    rFDisA <- data$Y[, 2, 1] - Xall[, -4, 2] %*% mCoef[-4, 2, 1]
    reHP <- data$Y[, 3, 1] - Xall[, -6, 3] %*% mCoef[-6, 3, 1]
    reHA <- data$Y[, 4, 1] - Xall[, -5, 4] %*% mCoef[-5, 4, 1]
    rdP <- data$Y[, 3, 2] - Xall[, -6, 7] %*% mCoef[-6, 3, 2]
    rdA <- data$Y[, 4, 2] - Xall[, -5, 8] %*% mCoef[-5, 4, 2]
    pred <- with(data, cbind(as.data.frame(X[, , 1]), RE, rFDisP, rFDisA, reHP, reHA, rdP, rdA))
    pred$col <- c("darkgray", "deepskyblue1", "firebrick1")[pred$type]
    out <- list(pred = pred, mCoef = mCoef)
    return(out)
  }
################################################################################
# END OF SCRIPT
################################################################################
################################################################################
