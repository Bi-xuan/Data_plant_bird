################################################################################
################################################################################
# MODEL INFERENCE
################################################################################
# EXTRACT INFORMAION FOR TABLES AND FIGURES
################################################################################
# RLQ ANALYSIS
# 1.) extract summary statistics from RLQ analysis
# 2.) write results of RLQ analysis to tables
################################################################################
rlqData <- rlqExtract(dataSet1)
rlqData$rlqDataMean
with(rlqData, {
  write.csv(subset(rlqData, select = -col), file = "output/02_rlqTable.csv", row.names = FALSE)
  write.csv(rlqDataMean, file = "output/03_rlqMeanTable.csv", row.names = FALSE)
})
################################################################################
# BAYESIAN HIERARCHICAL STRUCTURAL EQUATION MODEL
# 1.) extract summary statistics from model
# 2.) extract information for partial residual plots
# 3.) write model results to tables
################################################################################
pathData <- pathExtract(postSamp1)
resData <- residExtract(postSamp1, jagsData1)
with(pathData, {
  lapply(1:2, function(x) {
    write.csv(mTable[[x]], file = paste("output/0", x + 3, "_pathTable", x, ".csv", sep = ""), row.names = FALSE)
  })
})
################################################################################
################################################################################
# Figure 2. Plot of RLQ-fourth-corner analysis
################################################################################
pdf("output/figure_02.pdf", width = 0.394 * 12, height = 0.394 * 18)
m <- matrix(rep(c(1, 2, 4, 5, 7, 8), each = 4), 3, 8, byrow = TRUE)
m <- apply(m, 2, rep, each = 4)
m[c(40, 44, 48)] <- c(3, 6, 9)
layout(m)
par(cex = 0.7, oma = c(0, 1, 1, 0))
lapply(1:3, function(x) {
  par(mar = c(0.5, 0.5, 2, 0.5) + 0.1)
  RLQ <- dataSet1[[x]]$RLQ
  p <- RLQ$mR
  a <- RLQ$mQ
  N <- dataSet1[[x]]$nRLQ
  p <- cbind(p, (rowSums(N > 0)/max(rowSums(N > 0)) + 0.25) * 2)
  a <- cbind(a, (colSums(N > 0)/max(colSums(N > 0)) + 0.25) * 2)
  p <- p[order(p[, 3], decreasing = TRUE), ]
  a <- a[order(a[, 3], decreasing = TRUE), ]
  Ntab <- as.data.frame.table(N > 0)
  Ntab$pX <- p[as.character(Ntab$plant_species), 1]
  Ntab$pY <- p[as.character(Ntab$plant_species), 2]
  Ntab$aX <- a[as.character(Ntab$animal_species), 1]
  Ntab$aY <- a[as.character(Ntab$animal_species), 2]
  Ntab <- subset(Ntab, Freq > 0)
  plot.new()
  plot.window(xlim = c(-1, 1) * 1.1, ylim = c(-1, 1) * 1.1, asp = 1)
  abline(h = 0, v = 0, col = "gray")
  arrows(0, 0, RLQ$c1[, 1], RLQ$c1[, 2], col = "firebrick1",
         code = 2, angle = 30, length = 0.05, lwd = 1, lty = c(1, 5, 2))
  arrows(0, 0, RLQ$l1[, 1], RLQ$l1[, 2], col = "deepskyblue1",
         code = 2, angle = 30, length = 0.05, lwd = 1, lty = c(1, 5, 2))
  text(RLQ$c1 * 1.05, gsub("_", " ", rownames(RLQ$c1), fixed = TRUE),
       col = "firebrick1", font = 1, cex = 1, xpd = TRUE)
  text(RLQ$l1 * 1.05, gsub("_", " ", rownames(RLQ$l1), fixed = TRUE),
       col = "deepskyblue1", font = 1, cex = 1, xpd = TRUE)
  mtext(side = 3, line = 0, adj = -0.1, cex = 1, letters[1:6 %% 2 != 0][x], font = 2)
  mtext(side = 3, line = -0.5, adj = 0, cex = 1,
        c("Bird-fruit", "Bird-flower", "Insect-flower")[x], font = 2)
  if(x == 1) {
    legend("topright", bty = "n", lty = c(1, 5, 2), title = "Trait type",
           text.col = c("black", "black", "black"),
           legend = c("Matching", "Energy", "Foraging"))
  }
  plot.new()
  plot.window(xlim = range(c(p[, 1], a[, 1])), ylim = range(c(p[, 2], a[, 2])), asp = 1)
  with(Ntab, arrows(pX, pY, aX, aY, code = 0, lwd = 0.5, col = gray(0.75)))
  points(p[, 1:2], pch = 21, cex = p[, 3], lwd = 0.5,
         col = "black", bg = "deepskyblue1")
  points(a[, 1:2], pch = 21, cex = a[, 3], lwd = 0.5,
         col = "black", bg = "firebrick1")
  mtext(side = 3, line = 0, adj = 0, cex = 1, letters[1:6 %% 2 == 0][x], font = 2)
  if(x == 1) {
    legend("topleft", bty = "n",
           text.col = c("deepskyblue1", "firebrick1", "gray"),
           legend = c("Plants", "Animals", "Interactions"))
  }
  par(mar = c(0.5, 0.5, 0.5, 0.5) + 0.1)
  b <- barplot(RLQ$eig/sum(RLQ$eig), ylim = c(0, 1), axes = FALSE)
  mtext(side = 3, line = 0, adj = 0.5, cex = 0.5, "Eigenvalues")
  axis(2, at = seq(0, 1, 1), las = 1, mgp = c(1.5, 0.75, 0))
})
dev.off()
################################################################################
################################################################################
# Supplementary Figure 1. Trait-specific correlations with RLQ axes
# in the fourth-corner analysis
################################################################################
pdf("output/supplementary_figure_01.pdf", width = 0.394 * 10, height = 0.394 * 15)
m <- matrix(c(1:3), nrow = 3, byrow = TRUE)
layout(m, heights = c(1, 1, 2))
with(rlqData, {
  with(subset(rlqData, axis == 1), {
    par(cex = 0.7, mgp = c(2.25, 0.75, 0), mar = c(1, 5.5, 2, 1) + 0.1)
    plot.new()
    plot.window(xlim = c(0, 0.275), ylim = c(0.5, 3.5))
    points(abs(r), jitter(as.numeric(trait_type), factor = 0.5),
           cex = 1.5, lwd = 1, col = as.character(col),
           pch = 21)
    with(subset(rlqDataMean, axis == 1), {
      arrows(abs.r, as.numeric(trait_type) - 0.25,
             abs.r, as.numeric(trait_type) + 0.25, lwd = 2, code = 0)
      arrows(abs.r - se, as.numeric(trait_type),
             abs.r + se, as.numeric(trait_type), angle = 90, length = 0.05, code = 3)
    })
    axis(1)
    axis(2, at = 1:3, levels(trait_type), las = 1)
    box()
    mtext(side = 2, line = 4.5, cex = 0.7, las = 0, "Trait type")
    mtext(side = 3, line = -1, adj = 0.95, las = 0, "RLQ axis 1", cex = 0.7, font = 1)
    mtext(side = 3, line = 0.25, adj = -0.05, las = 0, "a", cex = 1, font = 2)
  })
  with(subset(rlqData, axis == 2), {
    plot.new()
    plot.window(xlim = c(0, 0.275), ylim = c(0.5, 3.5))
    points(abs(r), jitter(as.numeric(trait_type), factor = 0.5),
           cex = 1.5, lwd = 1, col = as.character(col),
           pch = 21)
    with(subset(rlqDataMean, axis == 2), {
      arrows(abs.r, as.numeric(trait_type) - 0.25,
             abs.r, as.numeric(trait_type) + 0.25, lwd = 2, code = 0)
      arrows(abs.r - se, as.numeric(trait_type),
             abs.r + se, as.numeric(trait_type), angle = 90, length = 0.05, code = 3)
    })
    axis(1)
    axis(2, at = 1:3, levels(trait_type), las = 1)
    box()
    mtext(side = 2, line = 4.5, cex = 0.7, las = 0, "Trait type")
    mtext(side = 3, line = -1, adj = 0.95, las = 0, "RLQ axis 2", cex = 0.7, font = 1)
    mtext(side = 3, line = 0.25, adj = -0.05, las = 0, "b", cex = 1, font = 2)
  })
  with(rlqData, {
    par(cex = 0.7, las = 1, mgp = c(2.25, 0.75, 0), mar = c(4, 5.5, 2, 1) + 0.1)
    plot.new()
    plot.window(xlim = c(0, 0.275), ylim = range(p), log = "y")
    points(abs(r), p,
           cex = 1.5, lwd = 1, col = as.character(col),
           pch = 21)
    axis(1)
    axis(2, at = c(1e0, 1e-1, 1e-2, 1e-3),
         labels = c(expression(10^0), expression(10^-1),
                    expression(10^-2), expression(10^-3)))
    axis(2, at = c(seq(0.2, 0.9, 0.1),
                   seq(0.02, 0.09, 0.01),
                   seq(0.002, 0.009, 0.001),
                   seq(0.0002, 0.0009, 0.0001)),
         labels = FALSE, tcl = -0.2)
    box()
    abline(h = 0.05, lty = 3)
    mtext(side = 2, line = 3.5, cex = 0.7, las = 0, expression(italic(P)-value))
    mtext(side = 3, line = 0.25, adj = -0.05, las = 0, "c", cex = 1, font = 2)
    mtext(side = 1, line = 2.75, cex = 0.7, las = 0,
          expression("Absolute correlation coefficient"~~italic(r)))
    legend("bottomleft", title = "Mutualism", title.col = "black",
           text.col = c("gray", "deepskyblue1", "firebrick1"), bty = "n",
           legend = c("Bird-fruit", "Bird-flower", "Insect-flower"))
  })
})
dev.off()
#################################################################################
#################################################################################
# Figure 3. pathmodel
#################################################################################
pdf("output/figure_03.pdf", width = 0.394 * 7, height = 0.394 * 17)
par(mfrow = c(3, 1), las = 1, mar = c(4, 4, 3, 5) + 0.1, cex = 0.7)
envData <- aggregate(cbind(elev, map, mat, lu) ~ plot + habitat,
                     FUN = mean, data = data1Sub)
envData <- droplevels(envData)
plot.new()
with(envData[order(envData$lu), ], {
  plot.window(ylim = c(0, 1), xlim = range(elev))
  m1 <- loess((mat - min(mat)) / (max(mat) - min(mat)) ~ elev, span = 0.5)
  x <- seq(min(elev), max(elev), 1)
  lines(x, predict(m1, data.frame(elev = x)), col = "firebrick1", lwd = 1)
  m2 <- loess((map - min(map)) / (max(map) - min(map)) ~ elev, span = 0.5)
  lines(x, predict(m2, data.frame(elev = x)), col = "deepskyblue1", lwd = 1)
  points(I((mat - min(mat)) / (max(mat) - min(mat))) ~ elev, pch = 21,
         col = "firebrick4", bg = ifelse(lu < 1, "firebrick1", "transparent"))
  points(I((map - min(map)) / (max(map) - min(map))) ~ elev, pch = 21,
         col = "deepskyblue4", bg = ifelse(lu < 1, "deepskyblue1", "transparent"))
  legend("topleft", xpd = TRUE, inset = c(-0.2, -0.4),
         ncol = 1, lty = 1, col = c("firebrick1", "deepskyblue1"),
         legend = c("Temp. (MAT)", "Prec. (MAP)"), bty = "n")
  legend("topright", xpd = TRUE, inset = c(-0.2, -0.4),
         ncol = 1, pch = 21, pt.bg = c("black", "transparent"),
         legend = c("Near-natural", "Anthropogenic"), bty = "n")
  axis(1, at = seq(1000, 4000, 1000))
  axis(2, at = (seq(5, 25, 5) - min(mat)) / (max(mat) - min(mat)),
       labels = seq(5, 25, 5))
  axis(4, at = (seq(500, 2500, 500) - min(map)) / (max(map) - min(map)),
       labels = seq(500, 2500, 500))
  mtext(side = 1, line = 2.5, cex = 0.7, text = "Elevation (m a.s.l.)")
  mtext(side = 2, line = 3, las = 0, cex = 0.7, text = "Temperature (°C)")
  mtext(side = 4, line = 3.5, las = 0, cex = 0.7, text = "Precipitation (mm yr-1)")
  mtext(side = 3, line = 2, adj = -0.325, font = 2, cex = 1, text = letters[1])
  box()
})
with(pathData, {
  sapply(1:2, function(x) {
    qgraph(mb[, , x], layout = nPos, mar = c(4, 4, 4, 6), asize = 5, arrowAngle = 0.3,
           borders = TRUE, maximum = max(mb[, , x], na.rm = TRUE), normalize = FALSE, vsize = 6, 
           labels = vars[[x]], label.cex = rep(1, 7), label.scale = FALSE, label.color = cols,
           bidirectional = TRUE, esize = 6, edge.label.position = 0.5, edge.label.cex = rep(1, 18),
           edge.labels = ifelse(mbf[, , x] > 2, signif(mb[, , x], 2), ""), edge.label.bg = "white",
           fade = FALSE, edge.color = ifelse(mbf[, , x] > 2, mc, "transparent"))
    mtext(side = 3, line = 1, adj = 0.05, font = 2, cex = 1,
          text = letters[x + 1])
    sapply(1:4, function(y) {
      text(rPos1[y, 1], rPos1[y, 2], adj = 0, cex = 0.7,
           labels = bquote(italic(r[m])^2 == .(round(rSqCum[1, y, x], 2))), xpd = TRUE)
      text(rPos2[y, 1], rPos2[y, 2], adj = 0, cex = 0.7,
           labels = bquote(italic(r[c])^2 == .(round(rSqCum[2, y, x], 2))), xpd = TRUE)
    })
  })
})
dev.off()
#################################################################################
#################################################################################
# Supplementary Figure 2. full pathmodel
#################################################################################
pdf("output/supplementary_figure_02.pdf", width = 0.394 * 6, height = 0.394 * 5.5)
par(oma = c(0, 0, 0, 1))
with(pathData, {
  qgraph(!is.na(mFull), layout = nPos, mar = c(2, 2, 2, 2), asize = 2, arrowAngle = 0.3,
       borders = TRUE, maximum = 1, normalize = FALSE, vsize = 4, 
       labels = varsFull, label.cex = rep(0.7, 7), label.scale = FALSE, label.color = cols,
       bidirectional = TRUE, esize = 1, edge.label.position = 0.5, edge.label.cex = rep(0.7, 18),
       fade = FALSE, edge.color = gray(0.5))
})
dev.off()
#################################################################################
#################################################################################
# Supplementary Figure 3. partial relationships from the pathmodel
#################################################################################
pdf("output/supplementary_figure_03.pdf", width = 0.394 * 16, height = 0.394 * 10)
m <- matrix(c(1:6), nrow = 2, byrow = TRUE); m <- cbind(m, c(7, 7))
layout(m, widths = c(2, 2, 2, 1))
par(cex = 0.7, las = 1, mgp = c(2.25, 0.75, 0), mar = c(3.5, 3.5, 2, 1) + 0.1)
with(resData, {
  with(pred, {
    plot(rFDisP ~ map, col = col, ylim = c(-1.5, 3),
         xlab = "MAP (mm/year)",
         ylab = expression(italic(FD[p])))
    curve(mCoef[3, 1, 1] * x, min(map), max(map), 
          add = TRUE, col = "blue", lwd = 2)
    mtext(side = 3, line = 0, adj = -0.1, font = 2, text = "a")
    
    plot(reHA ~ FDisP, col = col, ylim = c(-1.5, 3),
         xlab = expression(italic(FD[p])),
         ylab = expression(italic(e**H[a])))
    curve(mCoef[5, 4, 1] * x, min(FDisP), max(FDisP), 
          add = TRUE, col = "blue", lwd = 2)
    mtext(side = 3, line = 0, adj = -0.1, font = 2, text = "b")
    
    plot(rdA ~ FDisP, col = col, ylim = c(-1.5, 3),
         xlab = expression(italic(FD[p])),
         ylab = expression(italic(d*"'"[a])))
    curve(mCoef[5, 4, 2] * x, min(FDisP), max(FDisP), 
          add = TRUE, col = "blue", lwd = 2)
    mtext(side = 3, line = 0, adj = -0.1, font = 2, text = "c")
    mtext(side = 4, line = 0.5, adj = 0.5, las = 0, cex = 0.7, font = 2,
          text = "Bottom-up effects")
    
    plot(rFDisA ~ mat, col = col, ylim = c(-1.5, 3),
         xlab = "MAT (°C)",
         ylab = expression(italic(FD[a])))
    curve(mCoef[4, 2, 1] * x, min(mat), max(mat), 
          add = TRUE, col = "blue", lwd = 2)
    mtext(side = 3, line = 0, adj = -0.1, font = 2, text = "d")
    
    plot(reHP ~ FDisA, col = col, ylim = c(-1.5, 3),
         xlab = expression(italic(FD[a])),
         ylab = expression(italic(e**H[p])))
    curve(mCoef[6, 3, 1] * x, min(FDisA), max(FDisA), 
          add = TRUE, col = "blue", lwd = 2)
    mtext(side = 3, line = 0, adj = -0.1, font = 2, text = "e")
    
    plot(rdP ~ FDisA, col = col, ylim = c(-1.5, 3),
         xlab = expression(italic(FD[a])),
         ylab = expression(italic(d*"'"[p])))
    curve(mCoef[6, 3, 2] * x, min(FDisA), max(FDisA), 
          add = TRUE, col = "blue", lwd = 2)
    mtext(side = 3, line = 0, adj = -0.1, font = 2, text = "f")
    mtext(side = 4, line = 0.5, adj = 0.5, las = 0, cex = 0.7, font = 2,
          text = "Top-down effects")
    par(mar = c(0, 0, 0, 0) + 0.1)
    plot.new()
    legend("left", inset = c(-0.2, 0), title = "Mutualism", title.col = "black",
           text.col = unique(col), bty = "n",
           legend = c("Bird-fruit", "Bird-flower", "Insect-flower"))
  })
})
dev.off()
#################################################################################
#################################################################################
# Supplementary Figure 4. Pathmodel for trait-specific FD metrics
#################################################################################
pdf("output/supplementary_figure_04.pdf", width = 0.394 * 14, height = 0.394 * 17)
par(mfrow = c(3, 2), cex = 0.7, las = 1, mar = c(4, 4, 3, 5) + 0.1, oma = c(0, 1, 1, 0))
with(pathData, {
  sapply(1:6, function(x) {
    qgraph(mb[, , x + 2], layout = nPos, mar = c(4, 4.5, 4, 5.5),
           asize = 5, arrowAngle = 0.3, borders = TRUE,
           maximum = max(mb[, , x + 2], na.rm = TRUE), normalize = FALSE, vsize = 6, 
           labels = vars[[x]], label.cex = rep(1, 7), label.scale = FALSE, label.color = cols,
           bidirectional = TRUE, esize = 6, edge.label.position = 0.5, edge.label.cex = rep(1, 18),
           edge.labels = ifelse(mbf[, , x + 2] > 2, signif(mb[, , x + 2], 2), ""),
           edge.label.bg = "white",
           fade = FALSE, edge.color = ifelse(mbf[, , x + 2] > 2, mc, "transparent"))
    mtext(side = 3, line = 1, adj = 0.05, font = 2, cex = 1,
          text = letters[x])
    mtext(side = 2, line = 3.75, adj = 0.5, las = 0, font = 1, cex = 1,
          text = c("Matching traits", "", "Energy traits", "", "Foraging traits", "")[x])
    if(x %in% 1:2) {
      mtext(side = 3, line = 2.75, adj = 0.5, font = 1, cex = 1,
            text = c("Niche breadth", "Niche partitioning")[x])
    }
    sapply(1:4, function(y) {
      text(rPos1[y, 1], rPos1[y, 2], adj = 0, cex = 0.7,
           labels = bquote(italic(r[m])^2 == .(round(rSqCum[1, y, x + 2], 2))), xpd = TRUE)
      text(rPos2[y, 1], rPos2[y, 2], adj = 0, cex = 0.7,
           labels = bquote(italic(r[c])^2 == .(round(rSqCum[2, y, x + 2], 2))), xpd = TRUE)
    })
  })
})
dev.off()
################################################################################
# create descriptive summary of the results
################################################################################
capture.output({
  ##############################################################################
  # number of species and interactions in the mutualisms
  # correlations between species in the networks and FD metrics
  ##############################################################################
  cat("\n##############################\n")
  cat("SUMMARY OF DESCRIPTIVE RESULTS")
  cat("\n##############################\n")
  
  cat("\n##############################\n")
  cat("NUMBER OF SPECIES ACROSS MUTUALISMS:\n",
      "\nPlant species:", length(unique(unlist(sapply(dataSet1, function(x) x$N$plant_species)))),
      "\nBird species:", length(unique(unlist(sapply(dataSet1[1:2], function(x) x$N$animal_species)))),
      "\nInsect species:", length(unique(unlist(sapply(dataSet1[3], function(x) x$N$animal_species)))),
      "\n")
  
  cat("\n##############################\n")
  cat("SAMPLE SIZE:\n",
      "\nN - Observations:", nrow(data1Sub),
      "\nN - Study sites:", length(unique(data1Sub$plot)),
      "\nN - Mutualisms:", length(unique(data1Sub$type)),
      "\n")
  
  netSum <- lapply(dataSet1, function(x) x$summary)
  names(netSum) <- c("Bird-fruit", "Bird-flower", "Insect-flower")
  cat("\n##############################\n")
  cat("NUMBER OF SPECIES AND INTERACTIONS IN THE MUTUALISMS:\n",
      "\nbefore and after merging the network data with the trait data\n\n")
  print(netSum)
  ##############################################################################
  # RLQ-fourth-corner analysis
  # extraction of
  # 1.) global test statistic
  # 2.) relative eigenvalue of each axes
  # 3.) statistics for individual axes
  # 4.) contribution of traits to individual axes
  ##############################################################################
  cat("\n##############################\n")
  cat("GLOBAL TEST STATISTIC FROM RLQ ANALYSIS:\n\n")
  print(rlqData$rlqGlobal, digits = 2)
  
  cat("\n##############################\n")
  cat("RELATIVE EIGENVALUE OF EACH AXIS FROM RLQ ANALYSIS:\n\n")
  print(rlqData$rlqEigenvalues, digits = 2)
  
  cat("\n##############################\n")
  cat("AXES-SPECIFIC TEST STATISTIC FROM RLQ ANALYSIS:\n\n")
  print(rlqData$rlqAxes, digits = 2)
  
  cat("\n##############################\n")
  cat("ASSOCIATIONS OF DIFFERENT TRAIT TYPES WITH RLQ AXES:\n",
      "\nVARIABLE DESCRIPTION:",
      "\ntrait_type, trait type",
      "\naxis, RLQ axis",
      "\nk, number of significant tests",
      "\nn, number of statistical tests",
      "\nP-value, P-value of Moran's test",
      "\nabs.r, mean absolute correlation coefficient",
      "\nsd, standard deviation",
      "\nse, standard error",
      "\n\n")
  print(rlqData$rlqDataMean, digits = 2)
  ##############################################################################
  # power analysis
  ##############################################################################
  m1 <- with(rlqData$rlqData, lm(log(p/(1-p)) ~ abs(r) * log(n), na.action = "na.fail"))
  cat("\n##############################\n")
  cat("POWER ANALYSIS\n",
      "\nTesting for the effects of absolute effect size [abs(r)]",
      "\nand sample size (n) on probability (p)",
      "\nof accepting the null hypothesis:\n\n")
  print(summary(m1), digits = 3)
  
  cat("\n##############################\n")
  cat("SEM 1 (MULTIVARIATE FUNCTIONAL DIVERSITY & NICHE BREADTH)",
      "\nPath coefficients and selection probabilities:\n",
      "\nVARIABLE DESCRIPTION:",
      "\nx, predictor variable",
      "\ny, response variable",
      "\neffect, median (95% Credible interval)",
      "\npSel, selection probability",
      "\nbf, 2log[e](Bayes factor)",
      "\npsrf, potential scale reduction factor",
      "\neffn, effective sample size",
      "\n\n")
  print(pathData$mTable[[1]])
  
  cat("\n##############################\n")
  cat("SEM 2 (MULTIVARIATE FUNCTIONAL DIVERSITY & NICHE PARTITIONING)",
      "\nPath coefficients and selection probabilities:\n",
      "\nVARIABLE DESCRIPTION:",
      "\nx, predictor variable",
      "\ny, response variable",
      "\neffect, median (95% Credible interval)",
      "\npSel, selection probability",
      "\nbf, 2log[e](Bayes factor)",
      "\npsrf, potential scale reduction factor",
      "\neffn, effective sample size",
      "\n\n")
  print(pathData$mTable[[2]])
  
  cat("\n##############################\n")
  cat("MAXIMUM CORRELATION BETWEEN SPECIES RICHNESS AND FD METRICS:\n",
      "\nPlants: r =",
      signif(max(cor(subset(data1Sub, select = c(rP, FDisAllP, FDis1P, FDis2P, FDis3P)))[1, -1]), 2),
      "\nAnimals: r =",
      signif(max(cor(subset(data1Sub, select = c(rA, FDisAllA, FDis1A, FDis2A, FDis3A)))[1, -1]), 2), "\n")
}, file = "output/01_descriptive_summary.txt")
cat("Results have been written into the subfolder .../data_package/output/")
################################################################################
# END OF SCRIPT
################################################################################
################################################################################
