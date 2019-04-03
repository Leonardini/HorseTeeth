curDir = getwd()
if (!(substr(curDir, nchar(curDir) - nchar("HorseTeeth") + 1, nchar(curDir)) == "HorseTeeth")) {
  setwd("HorseTeeth/")
}
source("AlternativeDataPrep.R")

### Useful packages
listOfPackages2 = c("stats", "rgl", "ape", "phangorn", "RColorBrewer", "logistf", "brglm", "rrlda", "randomForest")
requirePackages(listOfPackages2)

### Useful constants
numBootstraps = 1000
first = FALSE
tol = 1e-8
exploratory = FALSE

### Various auxiliary functions
getNJTree = function(PMatrix, outgroup = "Paleolithic") {
  tree = nj(as.matrix(dist(PMatrix)))
  treeR = root(tree, outgroup = outgroup, resolve.root = TRUE)
  treeR = rotate(treeR, nrow(PMatrix) + 1)
  treeR
}

computeVars = function(PCresult, plot = FALSE) {
  SDs = PCresult$sdev
  Vars = SDs^2
  VarNorm = Vars/sum(Vars)
  if (plot) {
    plot(VarNorm, type = "l")
    dev.off()
  }
  VarNorm
}

checkFullRank = function(myMatrix) {
  myRank = qr(myMatrix)$rank
  return(myRank == min(dim(myMatrix)))
}

findDependentColumns = function(myMatrix) {
  badCols = c()
  curMat = matrix(NA, nrow(myMatrix), 0)
  for (ind in 1:ncol(myMatrix)) {
    curMat = cbind(curMat, rawData[,ind])
    if (!checkFullRank(curMat)) {
      badCols = c(badCols, ind)
      curMat = curMat[,-ncol(curMat)]
    }
  }
  badCols
}

### Exploring the differences in PCA functions
explorePCA = function(rawData) {
  PCA0 = prcomp(rawData, scale = TRUE, center = TRUE)   ### uses correlation
  PCA1 = princomp(rawData, scores = TRUE)               ### uses covariances
  PCA2 = princomp(rawData, cor = TRUE, scores = TRUE)   ### uses correlation
  Vars0 = computeVars(PCA0)
  Vars1 = computeVars(PCA1)
  Vars2 = computeVars(PCA2)
  print(paste("The difference between the two correlation-based PCA functions is", sum(abs(Vars0 - Vars2))))
  print(paste("The difference between the correlation and the covariance-based PCA functions is", sum(abs(Vars0 - Vars1))))
  print(paste("The first 3 PCs account for a fraction of", sum(Vars1[1:3]), "of the variance"))
}

### Sanity check on the data
sanityCheck = function(rawData) {
  rown = rownames(rawData)
  rown = sapply(rown, function(x) {tail(unlist(strsplit(x, "_")),1)})
  rown = sapply(rown, function(x) {gsub("DSCN","",x)})
  rown = sapply(rown, function(x) {as.integer(head(unlist(strsplit(x,"\\.")),1))})
  redLabels = Labels[,1]
  redLabels = sapply(redLabels, function(x) {tail(unlist(strsplit(x, "_")),1)})
  redLabels = sapply(redLabels, function(x) {as.integer(head(unlist(strsplit(x,"\\.")),1))})
  stopifnot(all(rown == redLabels))                     ### check that everything matches up!
}

### First plot: 3-D PCA by domestication status
make3DPCAStatus = function(Labels, PCA) {
  Labels1 = Labels[,"Domestication.Status"]
  Labels1[Labels1 == "Mixed"] = "Botai"
  Labels1 = as.factor(Labels1)
  cols1 = rainbow(length(levels(Labels1)))
  cols1[1:2] = cols1[2:1] ### switch the colors for Botai and Domestic
  open3d()
  par3d(windowRect = c(100, 100, 900, 900))
  plot3d(PCA$scores[,1:3], xlab = "PC1", ylab = "PC2", zlab = "PC3", col = cols1[as.integer(Labels1)], size = 1.5, type='s')
  bgplot3d({
    plot.new()
    title(main = "Principal components by domestication status", line = 1, cex.main = 2)
    legend("topright", legend = levels(Labels1), col = cols1, pch = 16, cex = 2, lwd = 0, lty = "solid", inset = 0.02, bty = "n")
  })
  snapshot3d(filename = 'DomStatusPCA3D.png', fmt = 'png')
  rgl.close()
}

### Second plot: 3-D PCA by time period
make3DPCATime = function(Labels, PCA) {
  Labels2 = as.factor(Labels[,"Time.Period"])
  cols2 = c("red", "magenta", "forestgreen", "deepskyblue")
  myOrder = c("Modern","LBA", "Eneolithic", "Paleolithic")
  myCols = cols2[match(Labels2, myOrder)]
  open3d()
  par3d(windowRect = c(100, 100, 900, 900))
  plot3d(PCA1$scores[,1:3], xlab = "PC1", ylab = "PC2", zlab = "PC3", col = myCols, size = 1.5, type='s')
  bgplot3d({
    plot.new()
    title(main = "Principal components by time period", line = 1, cex.main = 2)
    legend("topright", legend = myOrder, col = cols2, pch = 16, cex = 1.5, lwd = 0, lty = "solid", bty = "n")
  })
  snapshot3d(filename = 'TimePeriodPCA3D.png', fmt = 'png')
  rgl.close()
}

### First NJ tree: by time period
makeNJTreeTime = function(rawData, Labels2, refined = FALSE, Labels = NULL) {
  if (refined) {
    Labels2 = as.vector(Labels2)
    Labels2[Labels[,"Site"] == "Botai"] = "Botai"
    Labels2[Labels2 == "Eneolithic"] = "Eneolithic excluding Botai"
  }
  means2 = t(sapply(split(rawData, Labels2), function(x) {colMeans(matrix(x, ncol = ncol(rawData)))}))
  tree2R = getNJTree(means2)
  bootstrap2 = boot.phylo(tree2R, means2, getNJTree, B = numBootstraps, trees = TRUE)
  pdf(paste0("NJTreeByTimePeriod", rep("Refined", refined), ".pdf"))
  plot(tree2R, use.edge.length = TRUE)
  nodelabels(round(bootstrap2$BP/numBootstraps * 100), adj = -0.2, col = "red", frame = "none")
  dev.off()
}

### Third NJ tree: by domestication status
if (first) {
  Labels4 = Labels3
  Labels4[Labels4 %in% c("Eneolithic excluding Botai", "Paleolithic")] = "Prehistoric, wild"
  Labels4[Labels4 == "LBA"] = "Prehistoric, domestic"
  Labels4[Labels4 == "Modern"] = "Modern, domestic"
  means4 = t(sapply(split(rawData, Labels4), function(x) {colMeans(matrix(x, ncol = ncol(rawData)))}))
  tree4R = getNJTree(means4, outgroup = "Prehistoric, wild")
  bootstrap4 = boot.phylo(tree4R, means4, function(x) {getNJTree(x, outgroup = "Prehistoric, wild")}, B = numBootstraps, trees = TRUE)
  pdf("NJTreeByDomesticationStatus.pdf")
  plot(tree4R, use.edge.length = TRUE)
  nodelabels(round(bootstrap4$BP/numBootstraps * 100), adj = -0.2, col = "red", frame = "none")
  dev.off()
}

### Coordinate-based t-test and pairwise correlations
computeCoordinateCorrelations = function(Dataset, oneVsAll = FALSE) {
  allLabels = unique(Dataset$Labels)
  numLabels = length(allLabels)
  output = vector("list", ifelse(oneVsAll, numLabels, choose(numLabels, 2)))
  index = 1
  if (!oneVsAll) {
    allPairs = t(combn(1:numLabels, 2))
  }
  else {
    allPairs = cbind(1:numLabels, NULL)
  }
  for (ind1 in 1:nrow(allPairs)) {
    print(index)
    curInds = allPairs[ind1,]
    curLow = curInds[1]
    if (!oneVsAll) {
      curHigh = curInds[2]
      tempDataset = dichtomizeLabels(Dataset, lowLabel = allLabels[curLow], highLabel = allLabels[curHigh])
    }
    else {
      tempDataset = dichtomizeLabels(Dataset, lowLabel = allLabels[curLow])
    }
    trainData = tempDataset$Data
    Labels = tempDataset$Labels
    allRes = vector("list", ncol(trainData))
    for (ind in 1:ncol(trainData)) {
      curCol = trainData[,ind]
      cur0 = curCol[Labels == 0]
      cur1 = curCol[Labels == 1] 
      allRes[[ind]] = t.test(cur0, cur1)
      if (min(cur0) >= max(cur1) || min(cur1) >= max(cur0)) {
        print(paste("Detected perfect separation by", ind))
      }
    }
    pVals = sapply(allRes, function(x) {x$p.value})
    bestCoord = which.min(pVals)
    Cors = cor(trainData)
    diag(Cors) = 0
    print(paste("The range of the correlations is", range(Cors)))
    output[[index]] = list(Cors = Cors, pVals = pVals, bestCoord = bestCoord)
    index = index + 1
  }
  output
}

transformDatasetToDataFrame = function(Dataset, labelName = "status") {
  DF = cbind(Dataset$Data, Dataset$Labels)
  colnames(DF)[ncol(DF)] = labelName
  DF
}

### Penalized logistic regression
penalizedLR = function(trainData, testData, Predictions = list(), dichotomize = FALSE) {
  modelLRF = logistf(status ~., data = as.data.frame(trainData))
  catPreds = modelLRF$predict
  if (dichotomize) {
    catPreds = dichotomize(catPreds)
    trainAcc = getAccuracy(catPreds, trainData[,"status"])
  }
  else {
    trainAcc = getAUC(catPreds, trainData[,"status"], TRUE, filename = "TrainAUCLRFirth.pdf")
  }
  betas = coef(modelLRF)
  catPredsNew = 1 / (1 + exp(-cbind(1, testData) %*% betas))
  if (dichotomize) {
    catPredsNewBin = dichotomize(catPredsNew)
    testAcc = getAccuracy(catPredsNewBin, testData[,"status"])
  }
  else {
    testAcc = getAUC(catPredsNewBin, testData[,"status"], TRUE, filename = "TestAUCLRFirth.pdf")
  }
  Predictions[["LRFirth"]] = list(predictions = catPredsNew, model = modelLRF, trainAcc = trainAcc, testAcc = testAcc)
  Predictions
}
 
### Alternative with the brglm package (bias reduction)
penalizedLRAlternative = function(trainData, Predictions = list()) {
  modelLRB = brglm(status ~., family = binomial(link = 'logit'), data = as.data.frame(trainData))
  PredsBR = predict(modelLRB, newdata = as.data.frame(trainData), type = "response")
  getAccuracy(dichotomize(PredsBR), trainData[,"status"])
  PredsBRNew = predict(modelLRB, newdata = as.data.frame(testData), type = "response")
  table(dichotomize(PredsBRNew))
  getAccuracy(dichotomize(PredsBRNew), dichotomize(catPredsNew))
  Predictions[["LRRedBias"]] = list(PredsBRNew, modelLRB)
  Predictions
}

### From here on, only those functions that are exploratory are listed. TODO: these need to be refactored eventually.
if (exploratory) {
  Predictions = vector("list", 20)

### Detecting linear dependence relations
# alias(modelLRB)
# badComps = which(PCA1$sdev < tol)
# relations = PCA1$loadings[, badComps, drop=FALSE]
# relations = - round(t(t(relations) / apply(abs(relations), 2, max)), digits = 2)
# badCols = findDependentColumns(rawData[,-numCol])
# allTuples = vector("list", numCol)
# for (ind in 3:(numCol - 3)) {
#   allTuples[[ind]] = combn(numCol - 1 - 2 * (ind >= 9), ind, function(x) {!checkFullRank(rawData[,x])})
# }
# lincomb = t(Null(t(rawData[,-(17:19)])))
# lincombN = lincomb/min(abs(lincomb))
# Predictions[["LinearCombos"]] = list(goodCombos[,1], goodCombos[,2], lincombN)

### Part of the logistic regression postprocessing
# rankMatrix(rawData)                                   ### the model is not full rank; removing some columns!
# PredsLRCat = dichotomize(PredsLR, "Paleolithic", "LBA.Modern")
# with(modelLR, null.deviance - deviance)
# with(modelLR, df.null - df.residual)
# with(modelLR, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
# vif(modelLR)

### Part of the lasso/ridge regression postprocesisng
# predStats = prediction(model.pred, trainData[,numCol])
# predAccuracy = (predStats@tp[[1]] + predStats@tn[[1]])/nrow(trainData)
# bestAccuracy = max(predAccuracy)
# plot(predStats@cutoffs[[1]], predAccuracy)
# bestAlpha = predStats@cutoffs[[1]][which.max(predAccuracy)]
# roc.perf = performance(predStats, measure = "tpr", x.measure = "fpr")
# plot(roc.perf)
# curPred = dichotomize(model.pred.new, threshold = bestAlpha)
# table(curPred)

### Robust LDA
robustLDA = function(trainData, numCol = ncol(trainData)) {
  RLDModel = rrlda(x = as.data.frame(trainData[,-numCol]), grouping = as.factor(trainData[,numCol]))
  RLDtest = predict(RLDModel, x = as.data.frame(trainData[,-numCol]), type = class)
  Tab = table(RLDtest$class, trainData[,numCol])
  Tab
}

### QDA model without LOOCV
QDA = function(trainData, numCol = ncol(trainData)) {
  QDModel = qda(x = as.data.frame(trainData[,-numCol]), grouping = as.factor(trainData[,numCol]), CV = FALSE)
  QDModel
}

### LDA with LOOCV
LDCV = function(trainData, numCol = ncol(trainData)) {
  LDCV = lda(x = as.data.frame(trainData[,1:14]), grouping = as.factor(trainData[,numCol]), CV = TRUE)
  Tab = table(LDCV$class, trainData[,numCol])
  Tab
}

### Decision trees
DecTree = rpart(status ~ ., method = "class", data = as.data.frame(trainData), control = rpart.control(maxsurrogate = 3))
printcp(DecTree)
plotcp(DecTree)
summary(DecTree)
plot(DecTree, uniform = TRUE, main="Classification Tree for Horse Teeth")

### Random forests
RF = randomForest(x = as.data.frame(trainData[,-numCol]), y = as.factor(trainData[,numCol]))
RF$confusion

### From ProcessingNew.R: Mahalanobis distance-based probabilities for the Eneolithic samples
# train1Class0 = reduceToFullRank(restrictLabels(dataset1, c("LBA/Modern")) , fromEnd = TRUE, oneMore = TRUE)
# train1Class1 = reduceToFullRank(restrictLabels(dataset1, c("Paleolithic")), fromEnd = TRUE, oneMore = TRUE)
# d00 = mahalanobis(train1Class0$Data, colMeans(train1Class0$Data), cov(train1Class0$Data))
# d11 = mahalanobis(train1Class1$Data, colMeans(train1Class1$Data), cov(train1Class1$Data))
# dists1Class0 = mahalanobis(test1$Data[,1:ncol(train1Class0$Data)], colMeans(train1Class0$Data), cov(train1Class0$Data))
# dists1Class1 = mahalanobis(test1$Data[,1:ncol(train1Class1$Data)], colMeans(train1Class1$Data), cov(train1Class1$Data))
# P1C0 = 1 - pchisq(dists1Class0, df = ncol(train1Class0$Data))
# P1C1 = 1 - pchisq(dists1Class1, df = ncol(train1Class1$Data))
# train2Class0 = reduceToFullRank(restrictLabels(dataset2, c("Iron Age/LBA/Modern")) , fromEnd = TRUE, oneMore = TRUE)
# train2Class1 = reduceToFullRank(restrictLabels(dataset2, c("Paleolithic")), fromEnd = TRUE, oneMore = TRUE)
# dists2Class0 = mahalanobis(test2$Data[,1:ncol(train2Class0$Data)], colMeans(train2Class0$Data), cov(train2Class0$Data))
# dists2Class1 = mahalanobis(test2$Data[,1:ncol(train2Class1$Data)], colMeans(train2Class1$Data), cov(train2Class1$Data))
# P2C0 = 1 - pchisq(dists2Class0, df = ncol(train2Class0$Data))
# P2C1 = 1 - pchisq(dists2Class1, df = ncol(train2Class1$Data))
}