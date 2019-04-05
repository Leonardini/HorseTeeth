curDir = getwd()
if (!(substr(curDir, nchar(curDir) - nchar("HorseTeeth") + 1, nchar(curDir)) == "HorseTeeth")) {
  setwd("HorseTeeth/")
}
source("http://www.phaget4.org/R/myImagePlot.R")
source("AlternativeDataPrep.R")

### Useful packages
listOfPackages = c("pls", "c060", "pamr", "dplyr", "kernlab", "ROCR", "pROC", "car", "MASS", "e1071", "klaR", "useful")
listOfPackages = c(listOfPackages, c("glmnet", "plsgenomics", "FactoMineR", "Morpho", "fossil", "geomorph", "caret"))
listOfPackages = c(listOfPackages, "dendextend", "party", "rpart")
requirePackages(listOfPackages)

### Useful constants
numReps = 5
numFolds = 5

### Prepares raw data based on a directory
readFilesInDir = function(directory = "Raw data/P2s/", save = FALSE) {
  initDir = getwd()
  setwd(directory)
  allFiles = list.files()
  n = length(allFiles)
  allData = vector("list", n)
  names(allData) = allFiles
  for (curFile in allFiles) {
    curData = readland.tps(curFile, specID = "imageID")
    allData[[curFile]] = curData[, , 1]
  }
  setwd(initDir)
  allDim = sapply(allData, dim)
  p = allDim[1,1]
  stopifnot(all(allDim[1,] == p))
  k = allDim[2,1]
  stopifnot(all(allDim[2,] == k))
  fullData = array(NA, dim = c(p, k, n), dimnames = list(c(), c(), allFiles))
  for (ind in 1:n) {
    fullData[, , ind] = allData[[ind]]
  }
  if (save) {
    save(fullData, file = paste0(gsub("RawData", "", gsub("\\/", "", directory)), ".RData")) 
  }
  fullData
}

### This function prepares all raw data files, puts them in the left orientation for generalized Procrustes analysis
prepareRawData = function() {
  P2s = readFilesInDir("Raw data/P2s/", save = TRUE)
  leftTemplateP2  = normalizeSample(P2s[,,"2LP.tps"])
  rightTemplateP2 = normalizeSample(P2s[,,"2RP.tps"])
  orientedP2s = apply(P2s, 3, function(x) {testLeftRight(x, leftTemplateP2, rightTemplateP2, normTemplates = TRUE)})
  orientationsP2 = sapply(orientedP2s, function(x) {x[[2]]})
  stopifnot(orientationsP2["2LP.tps"] == "left")
  stopifnot(orientationsP2["2RP.tps"] == "right")
  orientedP2s = t(sapply(orientedP2s, function(x) {x[[1]]}))
  M3s = readFilesInDir("Raw data/M3s/", save = TRUE)
  leftTemplateM3  = normalizeSample(M3s[,,"2LM.tps"])
  rightTemplateM3 = normalizeSample(M3s[,,"2RM.tps"])
  orientedM3s = apply(M3s, 3, function(x) {testLeftRight(x, leftTemplateM3, rightTemplateM3, normTemplates = TRUE)})
  orientationsM3 = sapply(orientedM3s, function(x) {x[[2]]})
  stopifnot(orientationsM3["2LM.tps"] == "left")
  stopifnot(orientationsM3["2RM.tps"] == "right")
  orientedM3s = t(sapply(orientedM3s, function(x) {x[[1]]}))
  output = list(P2s = orientedP2s, M3s = orientedM3s)
  output
}

### This function uses the metadata to filter out any unannotated samples from the analysis and to annotate the rest
filterData = function(inputList) {
  P2s = inputList[["P2s"]]
  M3s = inputList[["M3s"]]
  metadataP2 = readDataset("Metadata edited.csv", labelCol = "From.data", dataCols = 2:4)
  metadataP2 = removeEmptyLabels(metadataP2)
  labelsP2 = toupper(metadataP2$Labels)
  stopifnot(all(labelsP2 %in% rownames(P2s)))
  originalP2s = P2s[labelsP2, , drop = FALSE]
  extraP2s = P2s[c(grep("LP", rownames(P2s)), grep("RP", rownames(P2s))), , drop = FALSE]
  P2s = rbind(originalP2s, extraP2s)
  PRZMatrix = matrix(c("Mixed", "PRZ", "PRZ"), nrow = nrow(extraP2s), ncol = 3, byrow = TRUE)
  colnames(PRZMatrix) = colnames(metadataP2$Data)
  metadataP2 = rbind(metadataP2$Data, PRZMatrix)
  ### Changed according to Alan's request; remove if the original classification, Eneolithic (ie unknown), is desired
  metadataP2[metadataP2[,"Site"] == "Shiderty", "Time.Period"] = "Paleolithic"
  P2FullDataset = list(Data = P2s, Labels = metadataP2)
  P2RedDataset = list(Data = P2s, Labels = metadataP2[,"Time.Period"])
  P2RedDataset = mergeLabels(P2RedDataset, c("LBA", "Modern"))
  metadataM3 = readDataset("All M3s.txt", labelCol = 1, dataCols = 2:4, sep = "\t")
  metadataM3 = removeDuplicates(metadataM3, labelsOnly = TRUE)
  labelsM3 = metadataM3$Data[,1]
  originalM3s = M3s[labelsM3, , drop = FALSE]
  extraM3s = M3s[c(grep("LM", rownames(M3s)), grep("RM", rownames(M3s))), , drop = FALSE]
  M3s = rbind(originalM3s, extraM3s)
  TimeMap = c("Eneolithic", "LBA", "LBA", "Paleolithic", "Paleolithic", "Iron Age", "Modern", "Eneolithic")
  names(TimeMap) = c("Botai", "Rogolik", "Kent", "31607", "33128", "Arzhan", "Modern", "Pavlodar")
  metadataM3$Data = cbind(metadataM3$Data[,c(3,2)], Time.Period = TimeMap[metadataM3$Data$Site])
  metadataM3 = rbind(metadataM3$Data, PRZMatrix)
  M3FullDataset = list(Data = M3s, Labels = metadataM3)
  M3RedDataset = list(Data = M3s, Labels = as.vector(metadataM3[,"Time.Period"]))
  M3RedDataset = mergeLabels(M3RedDataset, c("Iron Age", "LBA", "Modern"))
  output = list(P2s = P2FullDataset, M3s = M3FullDataset, redP2s = P2RedDataset, redM3s = M3RedDataset)
  output
}

### This function transforms an n x kp matrix (n = sample size, p = number of landmarks, k = number of dimensions)
### to a p x k x n array suitable for use with the gpagen function
reshapeMatrixToArray = function(Matrix, k = 2) {
  Dims = dim(Matrix)
  n = Dims[1]
  p = Dims[2]/k
  Array = array(NA, dim = c(p, k, n), dimnames = list(NULL, NULL, rownames(Matrix)))
  for (ind in 1:n) {
    Array[,,ind] = matrix(Matrix[ind,], ncol = k)
  }
  Array
}

### This function is the inverse of the one above
reshapeArrayToMatrix = function(Array) {
  Dims = dim(Array)
  p = Dims[1]
  k = Dims[2]
  n = Dims[3]
  Matrix = matrix(NA, nrow = n, ncol = p * k, dimnames = list(dimnames(Array)[[3]], NULL))
  for (ind in 1:n) {
    Matrix[ind,] = as.vector(Array[,,ind])
  }
  Matrix
}

procrustesAnalysis = function(inputData) {
  if (length(inputData) == 4) {
    inputData = inputData[3:4]
  }
  P2s = inputData[[1]]
  M3s = inputData[[2]]
  P2s$Data = reshapeMatrixToArray(P2s$Data)
  M3s$Data = reshapeMatrixToArray(M3s$Data)
  gpaP2s = gpagen(P2s$Data)
  print(paste("The final convergence criterion value for P2 was", gpaP2s$Q, "after", gpaP2s$iter, "iterations"))
  gpaM3s = gpagen(M3s$Data)
  print(paste("The final convergence criterion value for M3 was", gpaM3s$Q, "after", gpaM3s$iter, "iterations"))
  P2s$Data = reshapeArrayToMatrix(gpaP2s$coords)
  M3s$Data = reshapeArrayToMatrix(gpaM3s$coords)
  output = list(P2s = P2s, M3s = M3s)
  output
}

### This is the new driver function (for the dataset including PRZ)
updatedProcess = function() {
  output1 = prepareRawData()
  output2 = filterData(output1)
  output3 = procrustesAnalysis(output2)
  P2s = output3[[1]]
  M3s = output3[[2]]
  ### unsupervised methods: exploring the data
  ### optional step: compute pairwise RV coefficients
  pairsRV1M = allPairsRV(P2s, diagVsAll = FALSE)
  pairsRV2M = allPairsRV(M3s, diagVsAll = FALSE)
  ### optional step: plot the pairwise RV coefficients
  pdf("RVforP2Merged.pdf"); myImagePlot(pairsRV1M); dev.off()
  pdf("RVforM3Merged.pdf"); myImagePlot(pairsRV2M); dev.off()
  ### optional step: produce a 3-D PCA for both sets of teeth
  make3DPCA(P2s)
  make3DPCA(M3s)
  ### optional step: compute K-means clustering for both
  P2Kmeans = performKMeans(P2s)
  M3KMeans = performKMeans(M3s)
  ### optional step: produce a hierarchical clustering for both
  performHClust(P2s, filename = "P2HierarchicalClustering.pdf")
  performHClust(M3s, filename = "M3HierarchicalClustering.pdf")
  ### creating folds and segments to be used in cross-validation
  trainLabelsP2s = c("LBA/Modern", "Paleolithic", "PRZ")
  trainP2s = restrictLabels(P2s, trainLabelsP2s)
  testP2s = restrictLabels(P2s, c("Eneolithic"))
  foldsAndSegsP2 = createFoldsSpecial(trainP2s, numFolds = numFolds, numReps = numReps)
  foldsP2 = foldsAndSegsP2$fullFolds
  miniFoldsP2 = foldsAndSegsP2$miniFolds
  segsP2  = foldsAndSegsP2$Segments
  trainLabelsM3s = c("Iron Age/LBA/Modern", "Paleolithic", "PRZ")
  trainM3s = restrictLabels(M3s, trainLabelsM3s)
  testM3s = restrictLabels(M3s, c("Eneolithic"))
  foldsAndSegsM3 = createFoldsSpecial(trainM3s, numFolds = numFolds, numReps = numReps)
  foldsM3 = foldsAndSegsM3$fullFolds
  miniFoldsM3 = foldsAndSegsM3$miniFolds
  segsM3  = foldsAndSegsM3$Segments
  ### also creating data frames to be used with statistical formulas
  DFP2 = transformDatasetToDataFrame(P2s)
  DFM3 = transformDatasetToDataFrame(M3s)
  trainDFP2 = transformDatasetToDataFrame(trainP2s)
  trainDFM3 = transformDatasetToDataFrame(trainM3s)
  testDFP2  = transformDatasetToDataFrame(testP2s)
  testDFM3  = transformDatasetToDataFrame(testM3s)
  ### supervised methods: classification with statistical methods; starting with model selection by cross-validation
  ### failed attempts at penalized logistic regression
  statPredictions = list()
  fullList = 1:3
  for (index in 1:numReps) {
    curFoldM3 = foldsM3[[index]]
    for (ind in 1:numFolds) {
      curTestFold = which(curFoldM3 == ind)
      curTrainFold = which(curFoldM3 != ind)
      curTrainM3s = restrictDataset(trainM3s, curTrainFold)
      curTestM3s  = restrictDataset(trainM3s, curTestFold)
      for (i in 1:3) {
        print(paste(index, ind, i))
        redList = setdiff(fullList, i)
        redTrainM3s = dichotomizeLabels(restrictLabels(curTrainM3s, trainLabelsM3s[redList]), redList[1])
        redTestM3s  = dichotomizeLabels(restrictLabels(curTestM3s,  trainLabelsM3s[redList]), redList[1])
        redTrainDFM3 = transformDatasetToDataFrame(redTrainM3s, number = FALSE)
        redTestDFM3 = transformDatasetToDataFrame(redTestM3s, number = FALSE)
        # statPredictions = penalizedLRAlt(redTrainDFM3, redTestDFM3, statPredictions, dichotomize = TRUE, extras = paste("M3", index, ind, i))
      }
    }
  }
  ### LDA
  P2LDA = LDCV(trainDFP2)
  accLDA1 = computeAccuracy(P2LDA)
  print(paste("The cross-validation accuracy of LDA on P2 is", accLDA1))
  M3LDA = LDCV(trainDFM3)
  accLDA2 = computeAccuracy(M3LDA)
  print(paste("The cross-validation accuracy of LDA on M3 is", accLDA2))
  ### Robust LDA
  P2RLDA = robustLDA(trainDFP2)
  accRLDA1 = computeAccuracy(P2RLDA)
  print(paste("The cross-validation accuracy of robust LDA on P2 is", accRLDA1))
  M3RLDA = robustLDA(trainDFM3)
  accRLDA2 = computeAccuracy(M3RLDA)
  print(paste("The cross-validation accuracy of robust LDA on M3 is", accRLDA2))
  ### QDA
  redDFP2 = trainDFP2[trainDFP2[,'status'] != 3,] ### excludes PRZ - this still fails!
  P2QDA = QDCV(redDFP2)
  accQDA1 = computeAccuracy(P2QDA)
  print(paste("The cross-validation accuracy of QDA on P2 is", accQDA1))
  redDFM3 = trainDFM3[trainDFM3[,'status'] != 3,] ### excludes PRZ - this still fails!
  M3QDA = QDCV(redDFM3)
  accQDA2 = computeAccuracy(M3QDA)
  print(paste("The cross-validation accuracy of QDA on M3 is", accQDA2))
  ### supervised methods: classification with machine learning methods
  ### decision trees
  DT1 = decisionTree(trainDFP2, title = "Classification Tree for Horse Second Premolar Teeth")
  printcp(DT1)
  DT2 = decisionTree(trainDFM3, title = "Classification Tree for Horse Third Molar Teeth")
  printcp(DT2)
  ### random forests
  RF1 = rForest(trainDFP2)
  accRF1 = computeAccuracy(RF1[,-ncol(RF1)])
  print(paste("The accuracy of random forests on P2 is", accRF1))
  RF2 = rForest(trainDFM3)
  accRF2 = computeAccuracy(RF2[,-ncol(RF2)])
  print(paste("The accuracy of random forests on M3 is", accRF2))
  ### Support vector machines (SVMs)
  U1 = trainSVMs(trainData = trainDFP2, testData = testDFP2, Folds = miniFoldsP2)
  U2 = trainSVMs(trainData = trainDFM3, testData = testDFM3, Folds = miniFoldsM3)
}

computeAccuracy = function(Table) {
  stopifnot(all(rownames(Table) == colnames(Table)))
  accuracy = sum(diag(Table))/sum(Table)
  accuracy
}

make3DPCA = function(Dataset, labelName = "time period", tooth = "P2s") {
  Labels = Dataset$Labels
  Data = Dataset$Data
  PCA = princomp(Data, scores = TRUE)
  Labels = as.factor(Labels)
  cols = rainbow(length(levels(Labels)))
  open3d()
  par3d(windowRect = c(100, 100, 900, 900))
  plot3d(PCA$scores[,1:3], xlab = "PC1", ylab = "PC2", zlab = "PC3", col = cols[as.integer(Labels)], size = 1.5, type='s')
  ### the user may want to perform some manual rotation here!
  bgplot3d({
    plot.new()
    title(main = paste("Principal components by", labelName, "for", tooth), line = 1, cex.main = 2)
    legend("topright", legend = levels(Labels), col = cols, pch = 16, cex = 2, lwd = 0, lty = "solid", inset = 0.02, bty = "n")
  })
  snapshot3d(filename = paste0("PCA3DFor", tooth, ".png"), fmt = 'png')
  rgl.close()
}

### This function computes the angles between each triplet of consecutive points in 2-D, by using inverse cosines
computeAngles = function(xCoords, yCoords) {
  stopifnot(length(xCoords) == length(yCoords))
  allPoints = cbind(xCoords, yCoords)
  N = nrow(allPoints)
  stopifnot(N >= 3)
  shiftInds = c(2:N, 1)
  allVectors = allPoints - allPoints[shiftInds, ]
  allNorms = sqrt(rowSums(allVectors^2))
  allVectorsN = allVectors/allNorms
  allDots = rowSums(allVectorsN * allVectorsN[shiftInds,])
  allAngles = acos(allDots) * (180/pi)
  allAngles
}

### This function computes the area of the polygon whose vertices are the given points in 2-D, via the shoelace algo
computeArea = function(xCoords, yCoords) {
  stopifnot(length(xCoords) == length(yCoords))
  N = length(xCoords)
  stopifnot(N >= 3)
  shiftInds = c(2:N, 1)
  forwardLace = sum(xCoords * yCoords[shiftInds])
  reverseLace = sum(xCoords[shiftInds] * yCoords)
  area = abs(forwardLace - reverseLace)/2
  area
}

compareFinalNumbers = function(Vector1, Vector2) {
  stopifnot(length(Vector1) == length(Vector2))
  jointVector = c(Vector1, Vector2)
  names(jointVector) = NULL
  match = regexpr("[0-9]+[\\.]*[A-Z]*$", jointVector)
  numbers = as.numeric(gsub("[A-Z\\.]","", substr(jointVector, match, match + attr(match, "match.length") - 1)))
  half1 = numbers[1:(length(numbers)/2)]
  half2 = numbers[-(1:(length(numbers)/2))]
  stopifnot(half1 == half2)
}

oneVsAllRV = function(Dataset) {
  labels = sort(unique(Dataset$Labels))
  L = length(labels)
  oneVsAll = rep(NA, L)
  names(oneVsAll) = labels
  for (ind in 1:L) {
    tempDataset = dichtomizeLabels(Dataset, labels[ind])
    # curCov = cov(cbind(tempDataset$Data, as.numeric(tempDataset$Labels)))
    # oneVsAll[ind] = RV(curCov[, -ncol(curCov), drop = FALSE], curCov[, ncol(curCov), drop = FALSE])
    oneVsAll[ind] = coeffRV(tempDataset$Data, matrix(as.numeric(tempDataset$Labels), ncol = 1))$rv
  }
  oneVsAll
}

### This function normalizes a sample, compares it to two templates (left and right) to determine the best orientation
### and returns it in the left orientation, along with the original optimal orientation itself ("left" or "right").
### If normTemplates is TRUE, the templates are assumed (but not checked) to be normalized to mean 0 and unit size (!).
testLeftRight = function(sample, leftTemplate, rightTemplate, normTemplates = FALSE) {
  sample = normalizeSample(sample)
  if (normTemplates) {
    leftTemplate = normalizeSample(leftTemplate)
    rightTemplate = normalizeSample(rightTemplate)
  }
  refSample = reflectSample(sample, axis = "x")
  score1L = findOptimalRotation(sample, leftTemplate)$score
  score1R = findOptimalRotation(sample, rightTemplate)$score
  score2L = findOptimalRotation(refSample, leftTemplate)$score
  score2R = findOptimalRotation(refSample, rightTemplate)$score
  sign1 = sign(score1L - score1R)
  sign2 = sign(score2L - score2R)
  if (sign1 == 0 || sign2 == 0) {
    stop("Error: perfect match detected!")
  }
  if (sign1 == sign2) {
    stop("Error: inconsistent sample orientation or perfect symmetry detected!")
  }
  relScore1 = min(score1L, score1R) / max(score1L, score1R)
  relScore2 = min(score2L, score2R) / max(score2L, score2R)
  if (max(relScore1, relScore2) >= 0.9) {
    print("Warning: the Procrustes distances differ by less than 10%, which may suggest an ambiguous orientation!")
  }
  orientation = ifelse(sign1 == -1, "left", "right")
  if (orientation == "right") {
    output = list(sample = refSample, orientation = orientation)
  }
  else {
    output = list(sample = sample, orientation = orientation)
  }
  output
}

### This function multiplies the coordinates of a given sample by -1 along the specified axis, i.e. "x", "y" or "z"
reflectSample = function(sample, axis) {
  Dim = match(axis, c("x", "y", "z"))
  if (is.na(Dim)) {
    stop("Error: please pick an axis among x, y, z")
  }
  if (ncol(sample) >= Dim) {
    sample[,Dim] = -sample[,Dim]
  }
  else {
    stop("Error: invalid axis for the sample")
  }
  sample
}

### This function finds the optimal angle of rotation in two dimensions to match the sample with the template.
### It is assumed (without checking) that both the sample and the template are normalized to mean 0 and unit size (!).
### It returns a list with two named values: angle, the optimal angle (in radians), and score, the resulting distance.
findOptimalRotation = function(sample, template, Dim = 2) {
  Xs = template[,1]
  Ys = template[,2]
  Ws = sample[,1]
  Zs = sample[,2] 
  XW = sum(Xs * Ws)
  XZ = sum(Xs * Zs)
  YW = sum(Ys * Ws)
  YZ = sum(Ys * Zs)
  numerator = YW - XZ
  denominator = XW + YZ
  constant = sum(template^2) + sum(sample^2)
  angle = atan2(denominator, numerator)
  score = sqrt(constant - sqrt(numerator^2 + denominator^2))
  output = list(angle = angle, score = score)
  output
}

### This function normalizes the input sample to mean 0 and size (measured by the root mean square distance) 1.
normalizeSample = function(sample) {
  K = nrow(sample)
  M = ncol(sample)
  means = colMeans(sample)
  sample = sample - matrix(means, nrow = K, ncol = M, byrow = TRUE)
  size = sqrt(sum(sample^2)/K)
  sample = sample/size
  sample
}

allPairsRV = function(Dataset, diagVsAll = TRUE) {
  labels = sort(unique(Dataset$Labels))
  L = length(labels)
  allPairs = matrix(0, L, L, dimnames = list(labels, labels))
  if (diagVsAll) {
    diag(allPairs) = oneVsAllRV(Dataset)
  }
  else {
    diag(allPairs) = 1
  }
  for (ind1 in 1:(L - 1)) {
    for (ind2 in (ind1 + 1):L) {
      tempDataset = dichtomizeLabels(Dataset, lowLabel = labels[ind1], highLabel = labels[ind2])
      # curCov = cov(cbind(tempDataset$Data, as.numeric(tempDataset$Labels)))
      # curRV = RV(curCov[, -ncol(curCov), drop = FALSE], curCov[, ncol(curCov), drop = FALSE])
      curRV = coeffRV(tempDataset$Data, matrix(as.numeric(tempDataset$Labels), ncol = 1))$rv
      allPairs[ind1, ind2] = curRV
      allPairs[ind2, ind1] = curRV
    }
  }
  allPairs
}

allLandmarksRV = function(Dataset) {
  uLabels = unique(Dataset$Labels)
  stopifnot(length(uLabels) == 2)
  Dataset = dichtomizeLabels(Dataset, highLabel = uLabels[1], lowLabel = uLabels[2])
  DatasetR = reformatData(Dataset)
  numLandmarks = dim(DatasetR$Data)[1]
  allRVs = rep(NA, numLandmarks)
  for (ind in 1:numLandmarks) {
    cur = t(DatasetR$Data[ind, ,])
    # curCov = cov(cbind(cur, as.numeric(Dataset$Labels)))
    # allRVs[ind] = RV(curCov[, -ncol(curCov), drop = FALSE], curCov[, ncol(curCov), drop = FALSE])
    allRVs[ind] = coeffRV(cur, matrix(as.numeric(Dataset$Labels), ncol = 1))$rv
  }
  allRVs
}

### Obsolete function
prepareDatasetsHT = function() {
  Dataset1 = readDataset("P2-Clean.txt", labelCol = 1, dataCols = 2:19, sep = "\t")
  Dataset1 = removeDuplicates(Dataset1)
  Dataset1Aux = readDataset("Metadata edited.csv", labelCol = "Time.Period", dataCols = 1)
  Dataset1Aux = removeEmptyLabels(Dataset1Aux)
  compareFinalNumbers(unlist(Dataset1Aux$Data), Dataset1$Labels)
  Dataset1$Labels = Dataset1Aux$Labels
  Dataset2 = readDataset("All M3s.txt", labelCol = 1, dataCols = 5:20, sep = "\t")
  Dataset2Aux = readDataset("All M3s.txt", labelCol = "Site", dataCols = 1, sep = "\t")
  Dataset2 = removeDuplicates(Dataset2, labelsOnly = TRUE)
  Dataset2Aux = removeDuplicates(Dataset2Aux, rowsOnly = TRUE)
  Dataset2$Labels = Dataset2Aux$Labels
  TimeMap = c("Eneolithic","LBA","LBA","Paleolithic","Paleolithic","Iron Age","Modern","Eneolithic")
  names(TimeMap) = c("Botai","Rogolik","Kent","31607","33128","Arzhan","Modern","Pavlodar") 
  Dataset2 = renameLabels(Dataset2, TimeMap)
  labels2 = sort(unique(Dataset2$Labels))
  pairsRV1 = allPairsRV(Dataset1, diagVsAll = TRUE)
  pairsRV2 = allPairsRV(Dataset2, diagVsAll = TRUE)
  pdf("RVforP2.pdf"); myImagePlot(pairsRV1); dev.off()
  pdf("RVforM3.pdf"); myImagePlot(pairsRV2); dev.off()
  Dataset1M = mergeLabels(Dataset1, c("LBA", "Modern"))
  Dataset2M = mergeLabels(Dataset2, c("Iron Age", "LBA", "Modern"))
  pairsRV1M = allPairsRV(Dataset1M, diagVsAll = TRUE)
  pairsRV2M = allPairsRV(Dataset2M, diagVsAll = TRUE)
  pdf("RVforP2Merged.pdf"); myImagePlot(pairsRV1M); dev.off()
  pdf("RVforM3Merged.pdf"); myImagePlot(pairsRV2M); dev.off()
  output = list(Dataset1M, Dataset2M)
  output
}

reformatData = function(Dataset, Dim = 2) {
  Data = Dataset$Data
  m = nrow(Data)
  n = ncol(Data)
  stopifnot(n %% Dim == 0)
  nRed = n / Dim
  Array = array(NA, dim = c(nRed, Dim, m))
  for (ind in 1:Dim) {
    Array[ , ind, ] = t(Data[, Dim * (0:(nRed - 1)) + ind])
  }
  Dataset$Data = Array
  Dataset
}

performCVA = function(Dataset) {
  reformattedDataset = reformatData(Dataset)
  result = CVA(dataarray = reformattedDataset$Data, groups = reformattedDataset$Labels, cv = TRUE)
  result
}

plotCVAresults = function(CVAresults, filename = "Plot.pdf") {
  Scores = CVAresults$CVscore
  stopifnot(ncol(Scores) == 2)
  Labels = rownames(Scores)
  uLabels = sort(unique(Labels))
  matches = match(Labels, uLabels)
  L = length(uLabels)
  cols = rainbow(L)
  pdf(filename)
  plot(Scores[,1], Scores[,2], col = cols[matches], xlab = "CV1", ylab = "CV2")
  legend("topright", legend = uLabels, col = cols, lty = rep(1,L), xjust = 1, yjust = 1)
  dev.off()
}

generateSplits = function(Labels, balanced = FALSE, numReps = 5, numFolds = 5, seed = 1492) {
  set.seed(seed)
  numRow = length(Labels)
  fullFolds = vector("list", numReps)
  if (balanced) {
    for (r in 1:numReps) {
      curFold = balancedFolds(Labels, numFolds)
      names(curFold) = 1:numRow
      fullFolds[[r]] = curFold
    }
    Segments = unlist(lapply(fullFolds,
                             function(x) { lapply(split(x, x), function(y) {as.integer(names(y))})}), recursive = FALSE)
    Folds = lapply(Segments, function(x) {setdiff(1:numRow, x)})
  } 
  else {
    Folds = createMultiFolds(Labels, k = numFolds, times = numReps)
    Segments = lapply(Folds, function(x) {setdiff(1:numRow, x)})
    for (r in 1:numReps) {
      curFold = rep(NA, numRow)
      for (k in 1:numFolds) {
        curFold[Segments[[(r - 1) * numReps + k]]] = k
      }
      fullFolds[[r]] = curFold
    }
  }
  output = list(Folds = Folds, fullFolds = fullFolds, Segments = Segments)
  output
}

### This function applies a number of machine learning classification methods to a dataset
processDataset = function(trainDataset, testDataset, Method, Splits, loocv = FALSE, pr = FALSE) {
  if (Method == "svm") {
    Predictions = list()
    RANGE = -5 + 2 * (0:5)
    if (loocv) {
      ctrl = trainControl(method = "loocv", classProbs = pr, savePredictions = "final")
    }
    else {
    ctrl = trainControl(method = "repeatedcv", number = max(Splits$fullFolds[[1]]), 
        repeats = length(Splits$fullFolds), index = Splits$Folds, classProbs = pr, savePredictions = "final")
    }
    radGrid  = expand.grid(sigma = 2^RANGE, C = 2^RANGE)
    linGrid  = expand.grid(C = 2^RANGE)
    polyGrid = expand.grid(degree = 1:5, scale = 1, C = 2^RANGE)
    Grids = list("Radial" = radGrid, "Linear" = linGrid, "Poly" = polyGrid)
    for (method in names(Grids)) {
      print(method)
      svm.tune = train(x = trainDataset$Data, y = factor(make.names(trainDataset$Labels)),
                       method = paste0("svm",method), trControl = ctrl, tuneGrid = Grids[[method]])
      svmAccuracy = svm.tune$results[as.numeric(rownames(svm.tune$bestTune)), svm.tune$metric]
      # output = list(accuracy = svmAccuracy, error = error(svm.tune$finalModel))
      PredsSVM = predict(svm.tune, newdata = as.data.frame(testDataset$Data), type = ifelse(pr,"prob","raw"))
      Predictions[[paste0("SVM", method)]] = list(PredsSVM, svm.tune, svmAccuracy, 1 - error(svm.tune$finalModel))
    }
  }
  Predictions
}

transformViaPCA = function(Dataset, cor = FALSE, minVar = 0.95, numComp = NULL) {
  PCA = princomp(Dataset$Data, cor = cor, scores = TRUE)
  if (!is.null(minVar)) {
    eigs = PCA$sdev^2
    percVar = cumsum(eigs)/sum(eigs)
    numComp = which(percVar > minVar)[1]
    print(paste('Using the first', numComp, 'principal components to explain', percVar[numComp], 'variance'))
  }
  newDataset = list(Data = PCA$scores[,1:numComp], Labels = Dataset$Labels)
  newDataset
}

computeLDA = function(Dataset, skipColumns = 4, fromEnd = TRUE) {
  if (fromEnd) {
    colRange = 1:(ncol(Dataset$Data) - skipColumns)
  }
  else {
    if (skipColumns == 0) {
      colRange = 1:ncol(Datset$Data)
    }
    else {
      colRange = -(1:skipColumns)
    }
  }
  LDA = lda(x = Dataset$Data[, colRange, drop = FALSE], grouping = Dataset$Labels, CV = TRUE)
  LDAacc1 = mean(LDA$class == Dataset$Labels)
  LDAacc2 = mean(LDA$posterior[cbind(1:nrow(LDA$posterior), as.factor(Dataset$Labels))])
  output = list(LDA, LDAacc1, LDAacc2)
  output
}

performKMeans = function(Dataset, PCAfirst = FALSE, plot = FALSE, filename = NULL) {
  numLabels = length(unique(Dataset$Labels))
  if (PCAfirst) {
    Dataset$Data = princomp(Dataset$Data, scores = TRUE)$scores[,1:2]
  }
  Cluster = kmeans(Dataset$Data, centers = numLabels, nstart = 10)
  if (plot) {
    pdf(filename)
    plot(Cluster, data = Dataset$Data, class = Dataset$Labels)
    dev.off()
  }
  Tab = table(Cluster$cluster, Dataset$Labels)
  Ind = adj.rand.index(Cluster$cluster, as.numeric(as.factor(Dataset$Labels)))
  output = list(Tab, Ind)
  output
}

performHClust = function(Dataset, filename = "Plot.pdf") {
  Dendro = hclust(dist(Dataset$Data))
  Dendro = as.dendrogram(Dendro)
  uLabels = unique(Dataset$Labels)
  L = length(uLabels)
  coloursToUse = as.numeric(as.factor(Dataset$Labels))
  coloursToUse = coloursToUse[order.dendrogram(Dendro)]
  pdf(filename)
  Dendro %>% set("labels", rep("", length(Dataset$Labels))) %>%
    set("leaves_pch", 19) %>%
    set("leaves_cex", 0.6) %>% 
    set("leaves_col", coloursToUse) %>% 
    plot
  legend("topright", legend = uLabels, col = rainbow(L), lty = rep(1,L), xjust = 1, yjust = 1)
  dev.off()
}

fullAnalysis = function() {
  datasets = prepareDatasetsHT()
  dataset1 = datasets[[1]]
  dataset2 = datasets[[2]]
  CVA1 = performCVA(dataset1)
  typProbs1 = typprobClass(CVA1$CVscores, groups = as.factor(dataset1$Labels))
  # plotCVAresults(CVA1, "CVforP2.pdf")
  CVA2 = performCVA(dataset2)
  typProbs2 = typprobClass(CVA2$CVscores, groups = as.factor(dataset2$Labels))
  # plotCVAresults(CVA2, "CVforM3.pdf")
  train1 = restrictLabels(dataset1, c("LBA/Modern", "Paleolithic"))
  test1  = restrictLabels(dataset1, c("Eneolithic"))
  splits1 = generateSplits(train1$Labels)
  # result1SVM = processDataset(train1, test1, "svm", splits1)
  # hist(result1SVM$SVMPoly[[1]][,1])
  train2 = restrictLabels(dataset2, c("Iron Age/LBA/Modern", "Paleolithic"))
  test2  = restrictLabels(dataset2, c("Eneolithic"))
  splits2 = generateSplits(train2$Labels)
  # result2SVM = processDataset(train2, test2, "svm", splits2)
  # hist(result2SVM$SVMRadial[[1]][,1])
  CVA1R = performCVA(train1)
  typProbs1R = typprobClass(CVA1R$CVscores, groups = as.factor(train1$Labels))
  typProbs1S = typprobClass(CVA1R$CVscores, groups = as.factor(train1$Labels), small = TRUE, method = "wilson")
  # plotCVAresults(CVA1R, "CVforP2Train.pdf")
  CVA2R = performCVA(train2)
  typProbs2R = typprobClass(CVA2R$CVscores, groups = as.factor(train2$Labels))
  typProbs2S = typprobClass(CVA2R$CVscores, groups = as.factor(train2$Labels), small = TRUE, method = "wilson")
  # plotCVAresults(CVA2R, "CVforM3Train.pdf")
  ### Pure LDA - dropping the last 2 landmarks to avoid collinearity; accuracy below 60%
  result1LDA = computeLDA(train1, skipColumns = 4, fromEnd = TRUE)
  result2LDA = computeLDA(train2, skipColumns = 4, fromEnd = TRUE)
  ### PCA followed by LDA - only use enough PCs to explain 95% of the variance; accuracy below 95%
  PCA1 = transformViaPCA(train1, cor = TRUE, minVar = 0.95)
  result1PCALDA = computeLDA(PCA1, skipColumns = 0)
  result1PCASVM = processDataset(PCA1, test1, "svm", splits1, loocv = TRUE)
  PCA2 = transformViaPCA(train2, cor = TRUE, minVar = 0.95)
  result2PCALDA = computeLDA(PCA2, skipColumns = 0)
  result2PCASVM = processDataset(PCA2, test2, "svm", splits2, loocv = TRUE)
  ### RV coefficients with one landmark at a time, followed by LDA and SVM
  RV1 = allLandmarksRV(train1)
  RV2 = allLandmarksRV(train2)
  # train1R = selectVariables(train1, RV1, min = 0.95, pairs = TRUE)
  # train2R = selectVariables(train2, RV2, min = 0.95, pairs = TRUE)
  # result1RVLDA = computeLDA(train1R, skipColumns = 0)
  # result2RVLDA = computeLDA(train2R, skipColumns = 0)
  # result1RSVM = processDataset(train1R, test1, "svm", splits1)
  # result2RSVM = processDataset(train2R, test2, "svm", splits2)
  ### Cluster analysis with K-means; assign points to 2/3 clusters
  miniCluster1 = performKMeans(train1)
  miniCluster2 = performKMeans(train2)
  Cluster1 = performKMeans(dataset1)
  Cluster2 = performKMeans(dataset2)
  ### Repeating this with hierarchical clustering instead
  performHClust(train1, "HCforP2Train.pdf")
  performHClust(train2, "HCforM3Train.pdf")
  performHClust(dataset1, "HCforP2.pdf")
  performHClust(dataset2, "HCforM3.pdf")
  performHClust(PCA1, "HCforP2PCA.pdf")
  performHClust(PCA2, "HCforM3PCA.pdf")
}