# HorseTeeth

This repository contains various scripts, data and analysis results for the Horse Teeth project. 

Explanation of the files:

# Scripts

LegacyCode.R - contains some of the code that is now obsolete, though some remains useful

Processing.R - contains the code for both statistical and machine learning models

ProcessingNew.R - contains the main code for the processing we are currently doing

AlternativeDataPrep.R - contains the data manipulation functions used for machine learning

# Data files

Raw dataM3s.RData - contains the data files for all the M3 teeth as a geomorph data frame

Raw dataP2s.RData - contains the data files for all the P2 teeth as a geomorph data frame

Raw data folder - contains the landmarked PRZ teeth in two appropriately named subfolders

Metadata edited.csv - contains metadata for P2 teeth

AllM3.txt - only used for metadata (for M3 teeth) at this point

# Figures

PCA3D[M3/P2].png - the 3-D PCA plots for the corresponding teeth

[M3/P2]HierarchicalClustering.pdf - hierarchical clustering of the corresponding teeth, color-coded by class

RVFor[M3/P2]Merged.pdf - RV coefficients (generalization of Pearson correlation for matrices) for each pair of classes

# Text files

UnsupervisedMethods.txt - description of methods and results for the unsupervised analysis of the data

SupervisedMachineLearningMethods.txt - description of methods and results for the machine learning-based analysis of the data

SupervisedStatisticalMethods.txt - description of methods and results for the statistical method-based analysis of the data
