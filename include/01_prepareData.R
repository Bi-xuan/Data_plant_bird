################################################################################
################################################################################
# LOAD AND PREPARE DATA
################################################################################
if(! file.exists("data/jagsData1.RData")) {
  load("data/dataSet1.RData")
  ##############################################################################
  # combined RLQ-fourth-corner analysis
  ##############################################################################
  dataSet1 <- lapply(1:3, FUN = function(x) {
    ############################################################################
    # setup data matrices
    ############################################################################
    dataSet1[[x]]$pTRLQ <- sqrt(dataSet1[[x]]$pT[, 2:4])
    dataSet1[[x]]$aTRLQ <- sqrt(dataSet1[[x]]$aT[, 2:4])
    dataSet1[[x]]$nRLQ <- xtabs(frequency ~ plant_species + animal_species,
                                data = dataSet1[[x]]$N)
    dataSet1[[x]]$nRLQ <- with(dataSet1[[x]],
                               nRLQ[rownames(pTRLQ), rownames(aTRLQ)])
    dataSet1[[x]]$nRLQ <- with(dataSet1[[x]],
                               decostand(nRLQ, "pa"))
    ############################################################################
    # perform correspondence and principal components analyses
    ############################################################################
    dataSet1[[x]]$coa <- with(dataSet1[[x]],
                              dudi.coa(nRLQ, scannf = FALSE, nf = 2))
    dataSet1[[x]]$duP <- with(dataSet1[[x]],
                              dudi.pca(pTRLQ, scannf = FALSE,
                                       center = TRUE, scale = TRUE,
                                       nf = 2, row.w = coa$lw))
    dataSet1[[x]]$duA <- with(dataSet1[[x]],
                              dudi.pca(aTRLQ, scannf = FALSE,
                                       center = TRUE, scale = TRUE,
                                       nf = 2, row.w = coa$cw))
    ############################################################################
    # perform RLQ analysis
    ############################################################################
    dataSet1[[x]]$RLQ <- with(dataSet1[[x]],
                              rlq(duP, coa, duA, scannf = FALSE, nf = 2))
    ############################################################################
    # perform permutation tests
    ############################################################################
    dataSet1[[x]]$trRLQ <- 
      with(dataSet1[[x]],
           fourthcorner2(as.data.frame(pTRLQ),
                         as.data.frame(nRLQ),
                         as.data.frame(aTRLQ),
                         p.adjust.method.G = "fdr",
                         modeltype = 6, nrepet = 9999)$trRLQ)
    attach(dataSet1[[x]])
    dataSet1[[x]]$fcRLQ <- 
      fourthcorner.rlq(RLQ, type = "axes",
                       p.adjust.method.G = "fdr",
                       p.adjust.method.D = "fdr",
                       p.adjust.D = "levels",
                       modeltype = 6, nrepet = 9999)
    dataSet1[[x]]$qAxes <- 
      fourthcorner.rlq(RLQ, type = "Q.axes",
                       p.adjust.method.G = "fdr",
                       p.adjust.method.D = "fdr",
                       p.adjust.D = "levels",
                       modeltype = 6, nrepet = 9999)
    dataSet1[[x]]$rAxes <- 
      fourthcorner.rlq(RLQ, type = "R.axes",
                       p.adjust.method.G = "fdr",
                       p.adjust.method.D = "fdr",
                       p.adjust.D = "levels",
                       modeltype = 6, nrepet = 9999)
    detach(dataSet1[[x]])
    return(dataSet1[[x]])
  })
  ##############################################################################
  # plot-level data compilation for analysis
  ##############################################################################
  dataSet1 <- lapply(1:3, FUN = function(x) {
    ############################################################################
    # tabulate interaction networks for extraction of network metrics 
    ############################################################################
    dataSet1[[x]]$networks <- 
      apply(xtabs(frequency ~ plant_species + animal_species + plot_date_type,
                  data = dataSet1[[x]]$N),
            MARGIN = 3,
            FUN = function(y) {
              y[rowSums(y) > 0, colSums(y) > 0]
            })
    return(dataSet1[[x]])
  })
  ##############################################################################
  # create data frame with information about origin of interaction networks
  # study site ID, date of sampling, mutualism type
  ##############################################################################
  data1 <- data.frame(id = unlist(sapply(dataSet1, FUN = function(x) names(xtabs(~ plot_date_type, data = x$N)))))
  data1$plot <- sapply(strsplit(as.character(data1$id), "_", fixed = TRUE), FUN = function(x) x[1])
  data1$date <- as.Date(sapply(strsplit(as.character(data1$id), "_", fixed = TRUE), FUN = function(x) x[2]))
  data1$type <- factor(sapply(strsplit(as.character(data1$id), "_", fixed = TRUE), FUN = function(x) x[3]),
                       levels = c("Bird-fruit", "Bird-flower", "Insect-flower"))
  ##############################################################################
  # compute network metrics
  ##############################################################################
  data1 <- cbind(data1, do.call(rbind, lapply(dataSet1, FUN = function(x) {
    m <- lapply(x$networks, FUN = function(y) {
      data.frame(I = sum(y),
                 L = sum(y > 0),
                 rP = nrow(y),
                 rA = ncol(y),
                 eHP = mean(exp(diversityFun(y, MARGIN = 1))),
                 eHA = mean(exp(diversityFun(y, MARGIN = 2))),
                 dP = mean(dfun(y)$dprime),
                 dA = mean(dfun(t(y))$dprime))
    })
    do.call(rbind, m)
  })))
  ##############################################################################
  # compute functional diversity metrics
  ##############################################################################
  data1 <- cbind(data1, do.call(rbind, lapply(dataSet1, fdFun, which = "p", disMethod = "gower")))
  data1 <- cbind(data1, do.call(rbind, lapply(dataSet1, fdFun, which = "a", disMethod = "gower")))
  ##############################################################################
  # assign plot-level data to data for path analysis
  ##############################################################################
  data1$habitat <- plotData$habitat[match(as.character(data1$plot), as.character(plotData$plot))]
  data1$lu <- plotData$lu[match(as.character(data1$plot), as.character(plotData$plot))]
  data1$elev <- plotData$elev[match(as.character(data1$plot), as.character(plotData$plot))]
  data1$mat <- plotData$mat[match(as.character(data1$plot), as.character(plotData$plot))]
  data1$map <- plotData$map[match(as.character(data1$plot), as.character(plotData$plot))]
  ##############################################################################
  # take subset of data (remove NAs if present)
  ##############################################################################
  data1Sub <- subset(data1, 
                     subset = !(is.na(lu) | is.na(map) | is.na(mat) | is.na(dA) | is.na(dP)),
                     select = c(plot, type, habitat, lu, elev, mat, map,
                                FDisAllP, FDis1P, FDis2P, FDis3P,
                                FDisAllA, FDis1A, FDis2A, FDis3A,
                                I, L, rP, rA, eHP, eHA, dP, dA))
  ##############################################################################
  # data preparation for JAGS model
  ##############################################################################
  # array X with predictor variables for each of the 8 structural equation models
  ##############################################################################
  data1Sub$lu <- as.numeric(data1Sub$lu) - 1
  X1 <- X2 <- cbind(1, as.matrix(subset(data1Sub, select = c(lu, map, mat, FDisAllP, FDisAllA))))
  X3 <- X4 <- cbind(1, as.matrix(subset(data1Sub, select = c(lu, map, mat, FDis1P, FDis1A))))
  X5 <- X6 <- cbind(1, as.matrix(subset(data1Sub, select = c(lu, map, mat, FDis2P, FDis2A))))
  X7 <- X8 <- cbind(1, as.matrix(subset(data1Sub, select = c(lu, map, mat, FDis3P, FDis3A))))
  X <- sapply(mget(paste("X", 1:8, sep = "")), identity, simplify = "array")
  names(dimnames(X)) <- c("id", "X", "Model_index")
  dimnames(X)[[2]] <- c("Intercept", "lu", "map", "mat", "FDisP", "FDisA")
  dimnames(X)[[3]] <- c("FDisAll_eH", "FDisAll_d", "FDis1_eH", "FDis1_d",
                        "FDis2_eH", "FDis2_d", "FDis3_eH", "FDis3_d")
  X[, -1, ] <- apply(X[, -1, ], c(2, 3), scale)
  ##############################################################################
  # array Y with response variables for each of the 8 structural equation models
  ##############################################################################
  Y1 <- as.matrix(subset(data1Sub, select = c(FDisAllP, FDisAllA, eHP, eHA)))
  Y2 <- as.matrix(subset(data1Sub, select = c(FDisAllP, FDisAllA, dP, dA)))
  Y3 <- as.matrix(subset(data1Sub, select = c(FDis1P, FDis1A, eHP, eHA)))
  Y4 <- as.matrix(subset(data1Sub, select = c(FDis1P, FDis1A, dP, dA)))
  Y5 <- as.matrix(subset(data1Sub, select = c(FDis2P, FDis2A, eHP, eHA)))
  Y6 <- as.matrix(subset(data1Sub, select = c(FDis2P, FDis2A, dP, dA)))
  Y7 <- as.matrix(subset(data1Sub, select = c(FDis3P, FDis3A, eHP, eHA)))
  Y8 <- as.matrix(subset(data1Sub, select = c(FDis3P, FDis3A, dP, dA)))
  Y <- sapply(mget(paste("Y", 1:8, sep = "")), identity, simplify = "array")
  names(dimnames(Y)) <- c("id", "Y", "Model_index")
  dimnames(Y)[[2]] <- c("FDp", "FDa", "IDp", "IDa")
  dimnames(Y)[[3]] <- c("FDisAll_eH", "FDisAll_d", "FDis1_eH", "FDis1_d",
                        "FDis2_eH", "FDis2_d", "FDis3_eH", "FDis3_d")
  Y[, 3:4, c(1, 3, 5, 7)] <- log(Y[, 3:4, c(1, 3, 5, 7)]) # log-transform eH metrics
  Y <- ifelse(is.na(Y), NA, Y)
  Y[, , ] <- apply(Y[, , ], c(2, 3), scale)
  ##############################################################################
  # matrix RE with random effects
  ##############################################################################
  site <- as.numeric(factor(data1Sub$plot))
  type <- as.numeric(factor(data1Sub$type))
  RE <- cbind(site, type)
  ##############################################################################
  # indicator variables for the model
  ##############################################################################
  nX <- c(3, 3, 5, 5) # index for predictor variables
  names(nX) <- c("FDp", "FDa", "IDp", "IDa")
  nY <- dim(Y)[2] # index for response variables
  nRE <- ncol(RE) # index for random effects
  nP <- dim(Y)[3] # index for pathmodels
  nObs <- dim(X)[1] # index for observations
  nSite <- length(unique(site)) # index for study sites
  nType <- length(unique(type)) # index for mutualisms
  nE <- c(1, 1, 2, 2) # index for latent variables for covariance terms
  names(nE) <- c("FDp", "FDa", "IDp", "IDa")
  ##############################################################################
  # bundle data
  ##############################################################################
  jagsData1 <- list(X = X, nX = nX,
                    Y = Y, nY = nY, nP = nP, nE = nE,
                    RE = RE, nSite = nSite, nType = nType,
                    nRE = nRE, nObs = nObs)
  ##############################################################################
  # save the results
  ##############################################################################
  save(list = c("jagsData1", "data1", "data1Sub", "dataSet1"),
       file = "data/jagsData1.RData")
} else {
  ##############################################################################
  # alternatively load, if already exists
  ##############################################################################
  load(file = "data/jagsData1.RData")
}
################################################################################
# END OF SCRIPT
################################################################################
################################################################################
