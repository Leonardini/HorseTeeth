### Format of the data: list with 2 elements, Data and Labels, that have the same length (number of rows)
### Format of a map: a vector whose names are the elements to be mapped and entries are the ones to map to

restrictDataset = function(Dataset, goodPos) {
  Dataset$Data = Dataset$Data[goodPos, , drop = FALSE]
  Dataset$Labels = Dataset$Labels[goodPos]
  Dataset
}

restrictLabels = function(Dataset, labelList) {
  goodOnes = which(Dataset$Labels %in% labelList)
  Dataset = restrictDataset(Dataset, goodOnes)
  Dataset
}

trimLabels = function(Dataset, minSize) {
  Tab = table(Dataset$Labels)
  goodLabels = names(Tab)[Tab >= minSize]
  Dataset = restrictLabels(Dataset, goodLabels)
  Dataset
}

removeEmptyLabels = function(Dataset) {
  goodOnes = which(Dataset$Labels != "")
  Dataset = restrictDataset(Dataset, goodOnes)
  Dataset
}

removeIncomplete = function(Dataset) {
  goodOnes = complete.cases(Dataset$Data)
  Dataset = restrictDataset(Dataset, goodOnes)
  Dataset
}

removeDuplicates = function(Dataset, labelsOnly = FALSE, rowsOnly = FALSE) {
  badLabels = which(duplicated(Dataset$Labels))
  if (labelsOnly) {
    badPos = badLabels
  }
  else {
    badRows = which(duplicated(Dataset$Data))
    if (rowsOnly) {
      badPos = badRows
    }
    else {
      badPos = intersect(badLabels, badRows)
    }
  }
  goodPos = setdiff(1:length(Dataset$Labels), badPos)
  Dataset = restrictDataset(Dataset, goodPos)
  Dataset
}

renameLabels = function(Dataset, labelMap, defaultOld = TRUE) {
  oldLabels = Dataset$Labels
  newLabels = labelMap[oldLabels]
  if (defaultOld) {
    newLabels[is.na(newLabels)] = oldLabels[is.na(newLabels)]
  }
  Dataset$Labels = newLabels
  Dataset
}

mergeLabels = function(Dataset, labelList, sep = "/") {
  newLabelMap = rep(paste(labelList, collapse = sep), length(labelList))
  names(newLabelMap) = labelList
  Dataset = renameLabels(Dataset, newLabelMap)
  Dataset
}

dichtomizeLabels = function(Dataset, lowLabel, highLabel = NULL) {
  if (is.null(highLabel)) {
    highLabel = paste0("non-", lowLabel)
    otherLabels = setdiff(unique(Dataset$Labels), lowLabel)
    otherLabelMap = rep(highLabel, length(otherLabels))
    names(otherLabelMap) = otherLabels
    redDataset = renameLabels(Dataset, otherLabelMap)
  }
  else {
    redDataset = restrictLabels(Dataset, c(lowLabel, highLabel))
  }
  labelMap = c(0,1)
  names(labelMap) = c(lowLabel, highLabel)
  redDataset = renameLabels(redDataset, labelMap)
  redDataset
}

reorderVariables = function(Dataset, otherDataset) {
  Dataset$Data = Dataset$Data[ , colnames(otherDataset$Data), drop = FALSE]
  Dataset
}
  
renameVariables = function(Dataset, nameMap, defaultOld = TRUE) {
  oldNames = colnames(Dataset$Data)
  newNames = nameMap[oldNames]
  if (defaultOld) {
    newNames[is.na(newNames)] = oldNames[is.na(newNames)]
  }
  colnames(Dataset$Data) = newNames
  Dataset
}

selectVariables = function(Dataset, scores, minScore, pairs = FALSE) {
  stopifnot(ncol(Dataset$Data) == length(scores) * (1 + pairs))
  goodCols = which(scores > minScore)
  if (pairs) {
    goodCols = 2 * goodCols - (1:0)
  }
  Dataset$Data = Dataset$Data[, goodCols, drop = FALSE]
  Dataset
}

reduceToFullRank = function(Dataset, fromEnd = TRUE, oneMore = TRUE) {
  rk = qr(Dataset$Data)$rank
  while (rk < ncol(Dataset$Data)) {
    if (fromEnd) {
      Dataset$Data = Dataset$Data[,-ncol(Dataset$Data)]
    }
    else {
      Dataset$Data = Dataset$Data[,-1]
    }
    rk = qr(Dataset$Data)$rank
  }
  if (oneMore) {
    if (fromEnd) {
      Dataset$Data = Dataset$Data[,-ncol(Dataset$Data)]
    }
    else {
      Dataset$Data = Dataset$Data[,-1]
    }
  }
  Dataset
}

convertEntries = function(Dataset, converter) {
  oldData = Dataset$Data
  m = nrow(oldData)
  n = ncol(oldData)
  newData = matrix(NA, m, n, dimnames = dimnames(oldData))
  for (ind1 in 1:m) {
    for (ind2 in 1:n) {
      newData[ind1, ind2] = converter(oldData[ind1, ind2])
    }
  }
  Dataset$Data = newData
  Dataset
}

mergeDatasets = function(Dataset1, Dataset2) {
  stopifnot(all(colnames(Dataset1$Data) == colnames(Dataset2$Data)))
  curData = rbind(Dataset1$Data, Dataset2$Data)
  curLabels = c(Dataset1$Labels, Dataset2$Labels)
  mergedDataset = list(Data = curData, Labels = curLabels)
  mergedDataset
}

readDataset = function(filename, dataCols, labelCol, sep = ",") {
  Tab = read.csv(filename, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  Dataset = list(Data = Tab[, dataCols, drop = FALSE], Labels = Tab[, labelCol])
  Dataset
}

writeDataset = function(Dataset, filename, allCaps = FALSE, noUnderscore = FALSE, label = TRUE, semi = FALSE) {
  ### The default settings correspond to the regular way of writing a file
  coln = colnames(Dataset$Data)
  if (allCaps) {
    coln = toupper(coln)
  }
  if (noUnderscore) {
    match = regexpr("\\_", coln)
    goodOnes = which(match > 0)
    nextLetters = substr(coln[goodOnes], match[goodOnes] + 1, match[goodOnes] + 1)
    delete = goodOnes[nextLetters %in% c(letters, LETTERS)]
    merge  = goodOnes[!nextLetters %in% c(letters, LETTERS)]
    coln[delete] = substr(coln[delete], 1, match[delete] - 1)
    coln[merge]  = gsub("\\_", "", coln[merge])
  }
  colnames(Dataset$Data) = coln
  writeFun = ifelse(semi, write.csv2, write.csv)
  if (label) {
    writeFun(cbind(Label = Dataset$Label, Dataset$Data), file = filename, row.names = FALSE)
  }
  else {
    writeFun(Dataset$Data, file = filename, row.names = FALSE)
  }
}

### Checks if packages specified in the input list are available, installs them if needed, and then requires them all
requirePackages = function(listOfPackages) {
  for (package in listOfPackages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE, repos = "https://cloud.r-project.org")
    }
    require(package, character.only = TRUE)
  }
}