curDir = getwd()
if (!(substr(curDir, nchar(curDir) - nchar("HorseTeeth") + 1, nchar(curDir)) == "HorseTeeth")) {
  setwd("HorseTeeth/")
}
source("AlternativeDataPrep.R")

### Useful packages
listOfPackages1 = c("pls", "caret", "c060", "pamr", "dplyr", "kernlab", "ROCR", "pROC", "glmnet", "plsgenomics")
listOfPackages1 = c(listOfPackages1, c("car", "MASS", "klaR", "glmnet", "e1071", "MatrixCorrelation"))
requirePackages(listOfPackages1)

### Useful constants
RANGE = -5:5
SVMVars = 0
otherVars = 0
P2 = FALSE            ### TRUE if we are analyzing the P2 teeth, FALSE if we are analyzing the M3's
ways = 3              ### equals 3 if we are doing a 3-way classification (Wild vs Botai vs modern)
mySeed = 1492         ### The seed value used for the random number generator
exploratory = FALSE   ### TRUE if all the exploratory functions should be applied

### Various auxiliary functions
dichotomize = function(x, lowValue = 0, highValue = 1, threshold = 0.5) {
  y = rep(NA, length(x))
  y[x < threshold] = lowValue
  y[x > threshold] = highValue
  y
}

getAccuracy = function(predicted, true) {
  L = length(true)
  stopifnot(length(predicted) == L)
  accuracy = sum(predicted == true)/L
  accuracy
}

### This function computes the area under the receiver operating curve (AUROC) for a binary or real-valued prediction 
getAUC = function(predicted, true, plot = FALSE, filename = NULL) {
  AUC = auc(true, predicted)
  if (plot) {
    predob = prediction(predicted, true)
    perf = performance(predob, "tpr", "fpr")
    pdf(filename)
    plot(perf)
    area = format(round(AUC, 4), nsmall = 4)
    text(x = 0.8, y = 0.1, labels = paste("AUC =", area))
    # the reference x=y line
    segments(x0=0, y0=0, x1=1, y1=1, col="gray", lty=2)
    dev.off()
  }
  AUC
}

### Data preparation - obsolete now
### P2-Clean.txt is from: www.dropbox.com/home/GMM%20Horse%20Teeth/Data%20goes%20here!?preview=P2-Clean.txt
prepareData = function() {
  if (P2) {
    rawData = read.csv("P2-Clean.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    rawData = read.csv("All M3s.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  badID = which(duplicated(rawData[,1]))                ### see if any identifiers are duplicated
  if (P2) {
    dupRows = which(rawData[,1] == rawData[badID,1])      ### see which rows have duplicate identifiers
    sum(rawData[dupRows[1],-1] - rawData[dupRows[2],-1])  ### check if the two duplicate rows are identical
  } else {
    for (ind in badID) {
      numChar = nchar(rawData[ind - 1, 2])
      stopifnot(substr(rawData[ind, 2], 1, numChar) == rawData[ind - 1, 2])
    }
    # badID = badID - 1 ### uncomment this if you want to use the rc version instead of original
  }
  rawData = rawData[-badID,]                            ### remove the duplicate row(s)
  rownames(rawData) = rawData[,1]                       ### use the identifiers as row names
  rawData = rawData[,-1]                                ### remove the identifiers from the data
  rawData = as.matrix(rawData)
  rawData
}

### Function for adding metadata, now obsolete
addMetadata = function(rawData) {
  ### Metadata edited.csv is from: www.dropbox.com/home/GMM%20Horse%20Teeth/Data%20goes%20here!?preview=Metadata+edited.xlsx
  ### The first sheet only was converted to csv format using Microsoft Excel's "Save As" command
  if (P2) {
    Labels = read.csv("Metadata edited.csv", header = TRUE, stringsAsFactors = FALSE)
    Labels = Labels[-nrow(Labels),]                       ### removing the empty line
  } else {
    Labels =  rawData[,1:3]
    rawData = matrix(as.numeric(rawData[,-(1:3)]), nrow = nrow(rawData))
  }
  output = list(rawData, Labels)
  output
}

if (exploratory) {
### Data preparation
if (!P2) {
  TimeMap = c("Eneolithic","LBA","LBA","Paleolithic","Paleolithic","Iron Age","Modern","Eneolithic")
  names(TimeMap) = c("Botai","Rogolik","Kent","31607","33128","Arzhan","Modern","Pavlodar") 
  Labels = cbind(Labels, Time.Period = TimeMap[Labels[,"Site"]])
}
Labels5 = Labels[,"Time.Period"]
Labels5[Labels[,"Site"] == "Botai"] = "Botai"
Labels5[Labels5 == "Eneolithic"] = "Eneolithic excluding Botai"
Labels5[Labels5 %in% c("LBA", "Modern")] = "LBA/Modern"

if (ways == 3) {
  trainRows = which(Labels5 %in% c("Paleolithic", "Botai", "LBA/Modern"))
  testRows  = which(Labels5 == "Eneolithic excluding Botai")
} else {
  trainRows = which(Labels5 %in% c("Paleolithic", "LBA/Modern"))
  testRows  = which(Labels5 == "Botai")
}
trainData = cbind(rawData[trainRows,], status = ways - as.numeric(as.factor(Labels5[trainRows])))
testData  = rawData[testRows ,]
Predictions = vector("list", 30)
numRow = nrow(trainData)
numCol = ncol(trainData)
}

### Function for creating folds and their transformed versions
### balanced = TRUE if we are forcing the folds to balance (previously set to FALSE)
createFoldsSpecial = function(trainData = NULL, numFolds = 5, numReps = 5, balanced = TRUE, mySeed = 1492) {
  set.seed(mySeed)
  numRow = nrow(trainData$Data)
  fullFolds = vector("list", numReps)
  if (balanced) {
    for (r in 1:numReps) {
      curFold = balancedFolds(trainData$Labels, numFolds)
      names(curFold) = 1:numRow
      fullFolds[[r]] = curFold
    }
    Segments = unlist(lapply(fullFolds, function(x) {lapply(split(x, x), function(y) {
      as.integer(names(y))})}), recursive = FALSE)
    Folds = lapply(Segments, function(x) {setdiff(1:numRow, x)})
  } else {
    Folds = createMultiFolds(trainData$Labels, k = numFolds, times = numReps)
    Segments = lapply(Folds, function(x) {setdiff(1:numRow, x)})
    for (r in 1:numReps) {
      curFold = rep(NA, numRow)
      for (k in 1:numFolds) {
        curFold[Segments[[(r - 1) * numReps + k]]] = k
      }
      fullFolds[[r]] = curFold
    }
  }
  output = list(fullFolds = fullFolds, miniFolds = Folds, Segments = Segments)
  output
}

### Support vector machine
trainSVMs = function(trainData, numCol = ncol(trainData), testData = NULL, Folds = NULL) {
  print("SVM")
  Predictions = vector("list", 3)
  ctrl = trainControl(method = "repeatedcv", repeats = numReps, index = Folds)
  radGrid   = expand.grid(sigma = 2^RANGE, C = 2^RANGE)
  linGrid  = expand.grid(C = 2^RANGE)
  polyGrid = expand.grid(degree = 1:5, scale = 1, C = 2^RANGE)
  Grids = list("Radial" = radGrid, "Linear" = linGrid, "Poly" = polyGrid)
  for (method in names(Grids)) {
    print(method)
    SVMVars = SVMVars + 1 # removed preProc = c("center","scale") from the following line!
    svm.tune = train(x = trainData[,-numCol], y = factor(make.names(trainData[,numCol])),
                     method = paste0("svm",method), tuneGrid = Grids[[method]], trControl = ctrl)
    svmAccuracy = svm.tune$results[as.numeric(rownames(svm.tune$bestTune)), svm.tune$metric]
    PredsSVM = predict(svm.tune, newdata = as.data.frame(testData))
    Predictions[[paste0("SVM", method)]] = list(PredsSVM, svm.tune, svmAccuracy, 1 - error(svm.tune$finalModel))
  }
  Predictions
}

if (exploratory) {
### Partial least-squares analysis
print("PLS")
otherVars = otherVars + 1
for (r in 1:numReps) {
  Segs = Segments[numFolds * (r-1) + 1:numFolds]
  PLSModel = plsr(status ~ ., ncomp = numComps, data = as.data.frame(trainData), validation = "CV", segments = Segs)
  PLSPreds = round(predict(PLSModel, ncomp = numComps, newdata = as.data.frame(testData), type = "scores"))
  Predictions[[paste0("PLS", r)]] = list(PLSPreds, PLSModel, 1 - min(RMSEP(PLSModel)$val[1,1,]))
}
}

if (exploratory) {
### Logistic regression
if (ways == 2) {
  print("Logistic Regression")
  otherVars = otherVars + 1
  modelLRFull = glm(status ~., family = binomial(link = 'logit'), data = as.data.frame(trainData[, c(1:(numCol-7), numCol)]), control = list(maxit = 50))
  PredsLR = predict(modelLRFull, newdata = as.data.frame(testData), type = "response")
  accuracyLR = rep(NA, numReps * numFolds)
  for (r in 1:numReps) {
    for (k in 1:numFolds) {
      curIndex = numFolds * (r-1) + k
      curFold = Folds[[curIndex]]
      curRest = Segments[[curIndex]]
      modelLR = glm(status ~., family = binomial(link = 'logit'), data = as.data.frame(trainData[curFold, c(1:(numCol-7), numCol)]), control = list(maxit = 50))
      testLR = dichotomize(predict(modelLR, newdata = as.data.frame(trainData[curRest, 1:(numCol-7)]), type = "response"))
      accuracyLR[curIndex] = getAccuracy(testLR, trainData[curRest, "status"])
    }
    Predictions[[paste0("LogReg",r)]] = list(PredsLR, modelLRFull, mean(accuracyLR[numFolds*(r-1) + (1:numFolds)]))
  }
}

### Ridge and lasso regression
print("Ridge and Lasso Regression")
otherVars = otherVars + 2
grid = 2^RANGE
for (alpha in c(0, 1)) {
  model = ifelse(alpha == 0, "Ridge", "Lasso")
  fam = ifelse(ways == 2, "binomial", "multinomial")
  Type = ifelse(ways == 2, "response", "class")
  out = glmnet(x = trainData[,-numCol], y = trainData[,numCol], alpha = alpha, lambda = grid, family = fam)
  for (r in 1:numReps) {
    cv.out = cv.glmnet(x = trainData[,-numCol], y = trainData[,numCol], alpha = alpha, lambda = grid, family = fam, type.measure = "class", foldid = fullFolds[[r]])
    bestLambda = cv.out$lambda.min
    predAccuracy = 1 - cv.out$cvm[cv.out$lambda == bestLambda]
    model.pred = as.numeric(predict(out, type = Type, s = bestLambda, newx = trainData[,-numCol]))
    if (ways == 2) {
      fullAccuracy = getAccuracy(dichotomize(model.pred[,1]), trainData[,numCol])
    } else {
      fullAccuracy = getAccuracy(model.pred, trainData[,numCol])
    }
    model.pred.new = predict(out, type = Type, s = bestLambda, newx = testData)
    Predictions[[paste0(model, r)]] = list(out, model.pred.new, predAccuracy, fullAccuracy, bestLambda)
  }
}

### PLS combined with LDA
print("PLS-LDA")
otherVars = otherVars + 1
for (r in 1:numReps) {
  CombAccuracies = rep(NA, numFolds)
  AllComps = rep(NA, numFolds)
  for (k in 1:numFolds) {
    curIndex = numFolds * (r-1) + k
    curFold = Folds[[curIndex]]
    curRest = Segments[[curIndex]]
    curPreds = pls.lda(trainData[curFold,-numCol], trainData[curFold,numCol], trainData[curRest,-numCol], ncomp = numComps)
    AllComps[k] = curPreds$ncomp
    CombAccuracies[k] = getAccuracy(curPreds$predclass, trainData[curRest, numCol])
  }
  bestNumComps = round(weighted.mean(AllComps, CombAccuracies))
  CombPreds = pls.lda(trainData[,-numCol], trainData[,numCol], testData, ncomp = bestNumComps)
  Predictions[[paste0("PLS-LDA", r)]] = list(bestNumComps, CombPreds, mean(CombAccuracies))
}

### Pure LDA - leave out the last 4 columns (2 landmarks) to avoid collinearity
print("Pure LDA")
otherVars = otherVars + 1
LDModel = lda(x = as.data.frame(trainData[,1:(numCol-5)]), grouping = as.factor(trainData[,numCol]))
LDtest = predict(LDModel, newdata = as.data.frame(trainData[,1:(numCol-5)]))
fullLDAccuracy = getAccuracy(LDtest$class, trainData[,numCol])
LDPreds = predict(LDModel, newdata = as.data.frame(testData[,1:(numCol-5)]))
for (r in 1:numReps) {
  LDAccuracies = rep(NA, numFolds)
  for (k in 1:numFolds) {
    curIndex = numFolds * (r-1) + k
    curFold = Folds[[curIndex]]
    curRest = Segments[[curIndex]]
    curModel = lda(x = as.data.frame(trainData[curFold,1:(numCol-5)]), grouping = as.factor(trainData[curFold,numCol]))
    curPred = predict(curModel, newdata = trainData[curRest, 1:(numCol-5)])
    LDAccuracies[k] = getAccuracy(curPred$class, trainData[curRest, numCol])
  }
  Predictions[[paste0("LDA", r)]] = list(LDModel, LDPreds, mean(LDAccuracies), fullLDAccuracy)
}

### Naive Bayes classifier
print("Naive Bayes")
otherVars = otherVars + 1
NBmodel = NaiveBayes(x = as.data.frame(trainData[,-numCol]), grouping = as.factor(trainData[,numCol]))
NBtest = predict(NBmodel, newdata = as.data.frame(trainData[,-numCol]), type = class)
fullNBAccuracy = getAccuracy(NBtest$class, trainData[,numCol])
NBPreds = predict(NBmodel, newdata = as.data.frame(testData))
for (r in 1:numReps) {
  NBAccuracies = rep(NA, numFolds)
  for (k in 1:numFolds) {
    curIndex = numFolds * (r-1) + k
    curFold = Folds[[curIndex]]
    curRest = Segments[[curIndex]]
    curModel = NaiveBayes(x = as.data.frame(trainData[curFold, -numCol]), grouping = as.factor(trainData[curFold,numCol]))
    curPred = predict(curModel, newdata = trainData[curRest,-numCol])
    NBAccuracies[k] = getAccuracy(curPred$class, trainData[curRest, numCol])
  }
  Predictions[[paste0("NaiveBayes",r)]] = list(NBmodel, NBPreds, mean(NBAccuracies), fullNBAccuracy)
}

### Cleaning up, inspecting and saving the predictions
Predictions = Predictions[sapply(Predictions, length) > 0]
Performance = sapply(Predictions, function(x) {x[[3]]})
splitP = sapply(split(Performance, c(1:SVMVars, rep(SVMVars + (1:otherVars), each = numReps))), mean)
names(splitP) = gsub("1", "", names(Predictions)[c(1:SVMVars, SVMVars + 1 + (0:(otherVars-1)) * numReps)])
bestMethod = names(splitP)[which.max(splitP)]
print(paste("The best model is", bestMethod))
if (bestMethod %in% c("PLS", "Ridge", "Lasso", "PLS-LDA", "LDA", "NaiveBayes")) {
  bestInd = which.max(Performance[paste0(bestMethod, 1:5)])
  bestModel = Predictions[[paste0(bestMethod, bestInd)]]
} else {
  bestModel = Predictions[[bestMethod]]
}
bestPredictions = bestModel[[1]]
print(table(bestPredictions))
filename = paste0("All", rep("M3", !P2), "Preds", numFolds, "Folds", numReps, "Reps", 
                  ways, "Way", rep("Balanced", balanced), format(Sys.Date()), ".RData")
save(Predictions, file = filename)
}

### Optional part: hierarchical clustering, to see where the Botai samples fit
# allData = rbind(trainData[,-numCol], testData)
# allD = dist(allData)
# allH = hclust(allD)
# pdf(paste0("All", ifelse(P2, "P2", "M3"), "SamplesDendrogram.pdf"), 12, 4)
# plot(allH, labels = c(trainData[,numCol], rep("?", nrow(testData))), main = "Dendrogram for all samples", xlab = "Samples by class")
# legend(70, 0.25, "0: Prehistoric\n 1: LBA/Modern\n ?: Botai", bty = "n")
# dev.off()

### Optional part: multidimensional scaling, to see where the Eneolithic samples fit
# goodPos = which(Labels5 != "Iron Age")
# Labels5 = Labels5[goodPos]
# allData = rawData[goodPos,]
# allD = dist(allData)
# Coords = cmdscale(allD)
# pdf(paste0("All", ifelse(P2, "P2", "M3"), "MDSPlot.pdf"), 12, 6)
# labs = as.factor(Labels5)
# cols = rainbow(length(unique(Labels5)))
# plot(Coords[,1], Coords[,2], col = cols[as.numeric(labs)], main = "Multi-dimensional scaling for all samples", xlab = "Dimension 1", ylab = "Dimension 2")
# legend(0.08, ifelse(P2, 0.1, 0.2), legend = levels(labs), lty = c(2,2), col = cols)
# dev.off()