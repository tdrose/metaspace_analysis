######## Script for the visualization of the panels of the figure with the comparisons of the co-regulation of lipid pairs with the colocalization in the METASPACE database
# Note: the co-localization in METASPACE is here indicated with its longer but equivalent name: "Manders' co-occurrence"

#### High-level overview of work pipeline:
# Import of different required data, followed by cleaning, data orientation, and testing of the cleaned data, each time at point where their input is needed.

# Visualization of panel c and its associated statistical analysis in R.
# Visualization of panel a and b and their associated statistical analysis in R.

# Export data out of R, and do network visualization in Cytoscape.
# Extract relevant data from the network of Cytoscape, for further construction of outer circles of panel d in R, and import these again in R.

# Visualization of outer circles of panel d
# Combination of parts of panel d in Adobe Illustrator, and also combination of all figure panels and their esthetical clean-up in Adobe Illustrator.

#### Loading of packages and definition of a local convenience function
library("reshape2")

library("retistruct") # From GitHub: Install: devtools::install_github("davidcsterratt/retistruct@v0.6.4", subdir="pkg/retistruct")
library("beanplot")

library("fmsb")
library("RColorBrewer")

library("ggplot2") 
library("viridis")

library("RColorBrewer")


# Convenience function to return what is between brackets
GetStuffBetweenBrackets <- function(x){regmatches(x, gregexpr("(?<=\\().*?(?=\\))", x, perl=T))} 


#### Import major datasets needed for the analyses

# Import of the co-regulation data from the article of (Koeberlin et al., 2015, Cell)
KoeberlinCorrelations <- read.csv(file = "./InputData/Koeberlin_Snijder_Cell_2015_lipid_lipid_correlations.txt", header = TRUE, sep = "\t", as.is = TRUE)

KoeberlinCorrelations2 <- KoeberlinCorrelations[KoeberlinCorrelations[,1] != "SM C20:0", colnames(KoeberlinCorrelations) != "SM.C20.0"]
rm(KoeberlinCorrelations)

# Import the co-localization data resulting from calculation of Sergio and Theodore in METASPACE datasets: See headings for the identifiers of the datasets
Coexp_Mander_NoTreshold_Rework21022019 <- read.csv(file = "./InputData/Coexp_Mander_NoTreshold_Rework21022019.csv", header = TRUE, sep = ";", as.is = TRUE) # Needs unzipping of file first.

# Import of the analysis by Sergio of the lipid-like nature of the molecules in these data as defined by the Classifyre-based analysis
Lipid_Classes_Sergio_Subsets_210220194 <- read.csv(file = "./InputData/Lipid_Classes_Sergio_Subsets_210220194.txt", header = TRUE, sep = "\t", as.is = TRUE)

# Overview of the nature of the molecule-entries:
# 1032 lipid-like molecules; 179 other molecules; 1 (224): wrong annotation; 1 before that (223) extra column at the end

# Eliminate non-lipids from Sergio's subset
SubsetOfLipidClassesSergio <- Lipid_Classes_Sergio_Subsets_210220194[Lipid_Classes_Sergio_Subsets_210220194$SuperClass != "Other",]


# Prepare data and calculate the chemical formulas

Coexp_Mander_NoTreshold_Rework210220192 <- Coexp_Mander_NoTreshold_Rework21022019[,3:47]
rownames(Coexp_Mander_NoTreshold_Rework210220192) <- paste(Coexp_Mander_NoTreshold_Rework21022019[,1], Coexp_Mander_NoTreshold_Rework21022019[,2], sep = "_")

MList <- lapply(1:dim(Coexp_Mander_NoTreshold_Rework210220192)[2],function(y){sapply(Coexp_Mander_NoTreshold_Rework210220192[,y],function(x){if(is.na(x)){NA}else{strsplit(GetStuffBetweenBrackets(x)[[1]], ",")[[1]][1]}})})
MMatrix <- do.call("cbind", MList)

mode(MMatrix) <- "numeric"
MMeans <- cbind(Coexp_Mander_NoTreshold_Rework21022019[,1:2], rowMeans(MMatrix, na.rm = TRUE))

MMeans2 <- cbind(MMeans, cbind(sapply(MMeans[,1], function(x){x %in% SubsetOfLipidClassesSergio[,1]}),sapply(MMeans[,2], function(x){x %in% SubsetOfLipidClassesSergio[,1]})))
LipidSubsetMMeans <- MMeans2[MMeans2[,4] & MMeans2[,5],]

KoeberlinNamesSplitUp <- strsplit(KoeberlinCorrelations2[,1], " ")
LengthKoeberlinNamesSplitUp <- sapply(KoeberlinNamesSplitUp, function(x){length(x)})

KoeberlinNamesSplitUpMatrixWithMostLipids <- do.call("rbind", KoeberlinNamesSplitUp[LengthKoeberlinNamesSplitUp == 3])
KoeberlinNamesSplitUpMatrixWithSM <- do.call("rbind", KoeberlinNamesSplitUp[LengthKoeberlinNamesSplitUp == 2])

KoeberlinNamesSplitUpMatrixWithCers <- do.call("rbind", KoeberlinNamesSplitUp[LengthKoeberlinNamesSplitUp == 1])


KoeberlinNamesSplitUpMatrixWithSM2 <- cbind("SM", do.call("rbind", strsplit(gsub("^C","", KoeberlinNamesSplitUpMatrixWithSM[,2]), ":")))
KoeberlinNamesSplitUpMatrixWithSM4 <- cbind(KoeberlinNamesSplitUpMatrixWithSM, 23+as.numeric(KoeberlinNamesSplitUpMatrixWithSM2[,2]), 6, 2, 1, 48 + (as.numeric(KoeberlinNamesSplitUpMatrixWithSM2[,2])*2)-1-(2*as.numeric(KoeberlinNamesSplitUpMatrixWithSM2[,3])))

KoeberlinNamesSplitUpMatrixWithSM5 <- cbind(paste(KoeberlinNamesSplitUpMatrixWithSM4[,1], KoeberlinNamesSplitUpMatrixWithSM4[,2], sep = " "), KoeberlinNamesSplitUpMatrixWithSM4)
colnames(KoeberlinNamesSplitUpMatrixWithSM5) <- c("FullName", "Headgroup", "FattyAcid", "C", "O", "N", "P", "H")

KoeberlinNamesSplitUpMatrixWithCers2 <- do.call("rbind",strsplit(KoeberlinNamesSplitUpMatrixWithCers, "-"))
KoeberlinNamesSplitUpMatrixWithCers4 <- do.call("cbind", list(as.character(KoeberlinNamesSplitUpMatrixWithCers), do.call("rbind", strsplit(gsub("^C", "", sapply(strsplit(KoeberlinNamesSplitUpMatrixWithCers2[,2], "\\("), "[[",1)),":")), GetStuffBetweenBrackets(KoeberlinNamesSplitUpMatrixWithCers2[,2]), sapply(strsplit(KoeberlinNamesSplitUpMatrixWithCers2[,3], "\\("), "[[",1), GetStuffBetweenBrackets(KoeberlinNamesSplitUpMatrixWithCers2[,3])))


KoeberlinNamesSplitUpMatrixWithCers5 <- do.call("cbind", list(KoeberlinNamesSplitUpMatrixWithCers4,
                                                              
                                                              C = 18 + as.numeric(KoeberlinNamesSplitUpMatrixWithCers4[,2]),
                                                              O = 3+ as.numeric(KoeberlinNamesSplitUpMatrixWithCers4[,4] == "OH"),
                                                              
                                                              N = 1,
                                                              P = 0,
                                                              
                                                              H = 36 + (2*as.numeric(KoeberlinNamesSplitUpMatrixWithCers4[,6] == "2H")) + (2*as.numeric(KoeberlinNamesSplitUpMatrixWithCers4[,2]))-1-(2*as.numeric(KoeberlinNamesSplitUpMatrixWithCers4[,3]))))
mode(KoeberlinNamesSplitUpMatrixWithCers5) <- "character"

AtomVector <- c("C", "H", "N", "O", "P")
CHNOP <- KoeberlinNamesSplitUpMatrixWithCers5[1, AtomVector]

ChemicalFormulaCHNOPNo1 <- function(CHNOP){paste(sapply(1:length(CHNOP), function(x){if(as.numeric(CHNOP[x]) > 1){paste0(AtomVector[x],CHNOP[x])}else{if(as.numeric(CHNOP[x]) == 1){AtomVector[x]}else{""}}}), collapse = "" )}
KoeberlinNamesSplitUpMatrixWithCers7b <- cbind(KoeberlinNamesSplitUpMatrixWithCers5, apply(KoeberlinNamesSplitUpMatrixWithCers5[,AtomVector], 1, ChemicalFormulaCHNOPNo1))

KoeberlinNamesSplitUpMatrixWithSM7b <- cbind(KoeberlinNamesSplitUpMatrixWithSM5, apply(KoeberlinNamesSplitUpMatrixWithSM5[,AtomVector], 1, ChemicalFormulaCHNOPNo1))
KoeberlinNamesSplitUpMatrixWithMostLipids2 <- cbind(KoeberlinNamesSplitUpMatrixWithMostLipids, do.call("rbind", strsplit(gsub("^C","", KoeberlinNamesSplitUpMatrixWithMostLipids[,3]), ":")))

HeadGroupConversionMatrix <- do.call("cbind", list("HeadGroup" = c("PC", "LPC", "PE", "LPE", "PG", "LPG", "PG/BMP", "BMP", "PS", "PI", "FA", "PA", "TAG", "DAG", "VE", "VA"),
                                                   "C" = c(8, 8, 5, 5, 6, 6, 6, 6, 6, 9, 0, 3, 3, 3, 26, 20),
                                                   
                                                   "H" = c(18, 18, 12, 12, 13, 13, 13, 13, 12, 16, 1, 7, 5, 6, 50, 30),
                                                   "N" = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
                                                   
                                                   "O" = c(6, 6, 6, 6, 8, 8, 8, 8, 8, 11, 1, 6, 3, 3, 2, 1),
                                                   "P" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0)))

ConnectorPieceConversionMatrix <-cbind("Connector" = c("a", "e"),
                                       "C" = c(1,1), 
                                       
                                       "H" = c(0,2),
                                       "N" = c(0,0),
                                       
                                       "O" = c(1,0),
                                       "P" = c(0,0))

HTips <- do.call("cbind", list("HeadGroup" = c("PC", "LPC", "PE", "LPE", "PG", "LPG", "PG/BMP", "BMP", "PS", "PI", "FA", "PA", "TAG", "DAG", "VE", "VA"),
                               "HAmount" = c(2,2,2,2,2,2,2,2,2,2,1,2,3,2,0,0))) 


KoeberlinNamesSplitUpMatrixWithMostLipids4 <- cbind(KoeberlinNamesSplitUpMatrixWithMostLipids2, do.call("rbind", lapply(1:dim(KoeberlinNamesSplitUpMatrixWithMostLipids2)[1], function(z){
  
  y <- ConnectorPieceConversionMatrix[match(unlist(strsplit(KoeberlinNamesSplitUpMatrixWithMostLipids2[z,2], split = "")), ConnectorPieceConversionMatrix[,1]),-1]
  mode(y) <- "numeric"
  
  colSums(do.call("rbind", list(as.numeric(HeadGroupConversionMatrix[KoeberlinNamesSplitUpMatrixWithMostLipids2[z,1] == HeadGroupConversionMatrix[,"HeadGroup"], -1]), 
                                if(is.matrix(y)){colSums(y)}else{y},
                                
                                c(as.numeric(KoeberlinNamesSplitUpMatrixWithMostLipids2[z,4]) - nchar(KoeberlinNamesSplitUpMatrixWithMostLipids2[z,2]),
                                  2*(as.numeric(KoeberlinNamesSplitUpMatrixWithMostLipids2[z,4]) - nchar(KoeberlinNamesSplitUpMatrixWithMostLipids2[z,2])) + as.numeric(HTips[KoeberlinNamesSplitUpMatrixWithMostLipids2[z,1] == HTips[,1],2]) - 2*as.numeric(KoeberlinNamesSplitUpMatrixWithMostLipids2[z,5]), 0,0,0)
                                
  )))
})))

KoeberlinNamesSplitUpMatrixWithMostLipids5 <- cbind(KoeberlinNamesSplitUpMatrixWithMostLipids4, apply(KoeberlinNamesSplitUpMatrixWithMostLipids4[,AtomVector], 1, ChemicalFormulaCHNOPNo1))
KoeberlinNamesAndChemicalFormulas <- do.call("rbind", list(cbind(paste(KoeberlinNamesSplitUpMatrixWithMostLipids5[, 1],KoeberlinNamesSplitUpMatrixWithMostLipids5[, 2], KoeberlinNamesSplitUpMatrixWithMostLipids5[, 3], sep = " "), KoeberlinNamesSplitUpMatrixWithMostLipids5[, 11]), 
                                                           
                                                           KoeberlinNamesSplitUpMatrixWithCers7b[,c(1,12)],
                                                           KoeberlinNamesSplitUpMatrixWithSM7b[,c(1,9)]))

LipidSubsetMMeansx <- do.call("cbind", list(LipidSubsetMMeans, LipidSubsetMMeans[,1] %in% KoeberlinNamesAndChemicalFormulas[,2], LipidSubsetMMeans[,2] %in% KoeberlinNamesAndChemicalFormulas[,2]))
any(duplicated(LipidSubsetMMeansx[,1:2])) # FALSE

KoeberlinCorrelations2b07052019 <- KoeberlinCorrelations2
colnames(KoeberlinCorrelations2b07052019) <- c("lipid.lipid.corrs", KoeberlinCorrelations2b07052019$lipid.lipid.corrs)


MeltedKoeberlin <- melt(KoeberlinCorrelations2b07052019)

MeltedKoeberlin2 <- do.call("cbind", list(MeltedKoeberlin, KoeberlinNamesAndChemicalFormulas[match(MeltedKoeberlin[,1], KoeberlinNamesAndChemicalFormulas[,1]),2], KoeberlinNamesAndChemicalFormulas[match(MeltedKoeberlin[,2], KoeberlinNamesAndChemicalFormulas[,1]),2]))
colnames(MeltedKoeberlin2) <- c("FromLipid", "ToLipid", "Coregulation", "FromLipidChemicalFormula", "ToLipidChemicalFormula")


ListOfKoeberlinMatchesMETASPACE4 <- lapply(1:dim(MeltedKoeberlin2)[1], function(x){c(sapply(MeltedKoeberlin2[x,], as.character), LipidSubsetMMeansx[(LipidSubsetMMeansx[,1] == as.character(MeltedKoeberlin2[x,4])) & (LipidSubsetMMeansx[,2] == as.character(MeltedKoeberlin2[x,5])),])})

MatrixOfKoeberlinMatchesMETASPACE4 <- do.call("rbind", ListOfKoeberlinMatchesMETASPACE4)
mode(MatrixOfKoeberlinMatchesMETASPACE4) <- "character"

ComparisonCoregulationVsCooccurrence <- as.matrix(aggregate(as.numeric(MatrixOfKoeberlinMatchesMETASPACE4[,8]), list(Coregulation = as.numeric(MatrixOfKoeberlinMatchesMETASPACE4[,3])), function(x){c(mean(x, na.rm = TRUE), median(x, na.rm = TRUE), sd(x, na.rm = TRUE), sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))), sum(is.na(x)), length(x))}))
colnames(ComparisonCoregulationVsCooccurrence) <- c("co-regulation", "co-occurrence mean", "co-occurrence median", "co-occurrence sd", "co-occurrence se", "EmptyEntries", "AllEntries")

MatrixSplitsForBins <- cbind(seq(from = -1, to = 1, by = 0.05)[-41], seq(from = -1, to = 1, by = 0.05)[-1])
x <- 0.65

MatrixOfKoeberlinMatchesMETASPACE5 <- cbind(MatrixOfKoeberlinMatchesMETASPACE4, CoregulationBin = as.character(sapply(as.numeric(MatrixOfKoeberlinMatchesMETASPACE4[,3]), function(x){which((x > MatrixSplitsForBins[,1]) & (x <= MatrixSplitsForBins[,2]))})))
MatrixOfKoeberlinMatchesMETASPACE7 <- MatrixOfKoeberlinMatchesMETASPACE5[MatrixOfKoeberlinMatchesMETASPACE5[,"CoregulationBin"] != "integer(0)",]

MatrixOfKoeberlinMatchesMETASPACE8 <- MatrixOfKoeberlinMatchesMETASPACE7[order(as.numeric(MatrixOfKoeberlinMatchesMETASPACE7[,"CoregulationBin"])),]   


MMeansxx <- do.call("cbind", list(MMeans2, MMeans2[,1] %in% KoeberlinNamesAndChemicalFormulas[,2], MMeans2[,2] %in% KoeberlinNamesAndChemicalFormulas[,2]))
any(duplicated(MMeansxx[,1:2])) # FALSE


ListOfKoeberlinMatchesMETASPACE4xx <- lapply(1:dim(MeltedKoeberlin2)[1], function(x){c(sapply(MeltedKoeberlin2[x,], as.character), MMeansxx[(MMeansxx[,1] == as.character(MeltedKoeberlin2[x,4])) & (MMeansxx[,2] == as.character(MeltedKoeberlin2[x,5])),])})

MatrixOfKoeberlinMatchesMETASPACE4xx <- do.call("rbind", ListOfKoeberlinMatchesMETASPACE4xx)
mode(MatrixOfKoeberlinMatchesMETASPACE4xx) <- "character"

ComparisonCoregulationVsCooccurrencexx <- as.matrix(aggregate(as.numeric(MatrixOfKoeberlinMatchesMETASPACE4xx[,8]), list(Coregulation = as.numeric(MatrixOfKoeberlinMatchesMETASPACE4xx[,3])), function(x){c(mean(x, na.rm = TRUE), median(x, na.rm = TRUE), sd(x, na.rm = TRUE), sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))), sum(is.na(x)), length(x))}))
colnames(ComparisonCoregulationVsCooccurrencexx) <- c("co-regulation", "co-occurrence mean", "co-occurrence median", "co-occurrence sd", "co-occurrence se", "EmptyEntries", "AllEntries")

# Test for how much the lipid pairs of METASPACE cover the lipid pairs of the co-regulation dataset, and how it differs over the co-regulation: only smaller fraction not covered and each time similar fraction
pdf("./OutputFiles/ComparisonOfMETASPACECooccurrenceAndKoeberlinCoregulationDataMissingValuesFocus10052019xxboth.pdf")

plot(ComparisonCoregulationVsCooccurrencexx[, c("co-regulation", "AllEntries")], xlim = c(-1,1), type = "l", ylab = "Entries per binned co-regulation", main = "Distribution of unmapped entries (unfocused)")
lines(ComparisonCoregulationVsCooccurrence[, c("co-regulation", "EmptyEntries")], xlim = c(-1,1), col = "Blue")

lines(ComparisonCoregulationVsCooccurrencexx[, c("co-regulation", "EmptyEntries")], xlim = c(-1,1), col = "Red")
legend("bottom", c("All lipid pairs (K?berlin data)", "Lipid pairs absent from METASPACE (lipid-like subset)", "Lipid pairs absent from METASPACE (all pairs)"), fill = c("black", "blue", "red"), border = "white", cex = 0.7)

plot(ComparisonCoregulationVsCooccurrencexx[, "co-regulation"], ComparisonCoregulationVsCooccurrencexx[, "EmptyEntries"]*100/ComparisonCoregulationVsCooccurrencexx[, "AllEntries"], xlim = c(-1,1), ylim = c(0,100), type = "l", xlab = "Co-regulation", ylab = "Ratio of unmapped entries (%)", main = "Distribution of unmapped entry ratios", col = "red")
lines(ComparisonCoregulationVsCooccurrence[, "co-regulation"], ComparisonCoregulationVsCooccurrence[, "EmptyEntries"]*100/ComparisonCoregulationVsCooccurrence[, "AllEntries"], xlim = c(-1,1), col = "blue")

abline(h = sum(MatrixOfKoeberlinMatchesMETASPACE4xx[,7] == "character(0)")*100/dim(MatrixOfKoeberlinMatchesMETASPACE4xx)[1], col = "red", lty = 2)
abline(h = sum(MatrixOfKoeberlinMatchesMETASPACE4[,7] == "character(0)")*100/dim(MatrixOfKoeberlinMatchesMETASPACE4)[1], col = "blue", lty = 2)

text(x = 0.95, y = sum(MatrixOfKoeberlinMatchesMETASPACE4[,7] == "character(0)")*100/dim(MatrixOfKoeberlinMatchesMETASPACE4)[1]+3.2, labels = "total", cex = 0.85, col = "blue")
text(x = 0.95, y = sum(MatrixOfKoeberlinMatchesMETASPACE4xx[,7] == "character(0)")*100/dim(MatrixOfKoeberlinMatchesMETASPACE4xx)[1]-2.8, labels = "total", cex = 0.85, col = "red")

legend("bottom", c("Lipid-like subset of METASPACE dataset", "All pairs of METASPACE dataset"), fill = c("blue", "red"), border = "white", cex = 0.7)
dev.off()

MatrixSplitsForBins <- cbind(seq(from = -1, to = 1, by = 0.05)[-41], seq(from = -1, to = 1, by = 0.05)[-1])
x <- 0.65

MatrixOfKoeberlinMatchesMETASPACE5xx <- cbind(MatrixOfKoeberlinMatchesMETASPACE4xx, CoregulationBin = as.character(sapply(as.numeric(MatrixOfKoeberlinMatchesMETASPACE4xx[,3]), function(x){which((x > MatrixSplitsForBins[,1]) & (x <= MatrixSplitsForBins[,2]))})))
MatrixOfKoeberlinMatchesMETASPACE7xx <- MatrixOfKoeberlinMatchesMETASPACE5xx[MatrixOfKoeberlinMatchesMETASPACE5xx[,"CoregulationBin"] != "integer(0)",]

MatrixOfKoeberlinMatchesMETASPACE8xx <- MatrixOfKoeberlinMatchesMETASPACE7xx[order(as.numeric(MatrixOfKoeberlinMatchesMETASPACE7xx[,"CoregulationBin"])),] 


#### Add information about the subclass and determine the type of pair
MatrixOfKoeberlinMatchesMETASPACE8xx2 <- do.call("cbind", list(MatrixOfKoeberlinMatchesMETASPACE8xx, 
                                                               
                                                               SubclassFirst = unlist(lapply(strsplit(MatrixOfKoeberlinMatchesMETASPACE8xx[,1], " "), function(y){if(length(y) == 3){paste(y[1],y[2], sep = "_")}else{if(length(y) == 2){"SM"}else{if(length(y) == 1){"Cer"}}}})),
                                                               SubclassSecond = unlist(lapply(strsplit(MatrixOfKoeberlinMatchesMETASPACE8xx[,2], " "), function(y){if(length(y) == 3){paste(y[1],y[2], sep = "_")}else{if(length(y) == 2){"SM"}else{if(length(y) == 1){"Cer"}}}}))))

SubclassToClassConversionMatrix <- cbind(unique(lapply(strsplit(MatrixOfKoeberlinMatchesMETASPACE8xx2[,"SubclassFirst"], "_"), "[[", 1)), c("PE", "PS", "Cer", "PC", "SM", "PC","PG", "PE", "PG"))
MatrixOfKoeberlinMatchesMETASPACE8xx4 <- do.call("cbind", list(MatrixOfKoeberlinMatchesMETASPACE8xx2, 
                                                               
                                                               ClassFirst = unlist(SubclassToClassConversionMatrix[match(sapply(strsplit(MatrixOfKoeberlinMatchesMETASPACE8xx2[,"SubclassFirst"], "_"), "[[", 1), SubclassToClassConversionMatrix[,1]),2]),
                                                               ClassSecond = unlist(SubclassToClassConversionMatrix[match(sapply(strsplit(MatrixOfKoeberlinMatchesMETASPACE8xx2[,"SubclassSecond"], "_"), "[[", 1), SubclassToClassConversionMatrix[,1]),2])))

MatrixOfKoeberlinMatchesMETASPACE8xx5 <- do.call("cbind", list(MatrixOfKoeberlinMatchesMETASPACE8xx4, PairType = (MatrixOfKoeberlinMatchesMETASPACE8xx4[,"ClassFirst"] == MatrixOfKoeberlinMatchesMETASPACE8xx4[,"ClassSecond"]) + (MatrixOfKoeberlinMatchesMETASPACE8xx4[,"SubclassFirst"] == MatrixOfKoeberlinMatchesMETASPACE8xx4[,"SubclassSecond"])))



# Sideplots giving an overview at species level how the different subclasses relate to each other qua co-regulation vs. Manders' co-occurrence

MatrixOfKoeberlinMatchesMETASPACE8xx8 <- cbind(MatrixOfKoeberlinMatchesMETASPACE8xx5, SubclassLinks = paste(MatrixOfKoeberlinMatchesMETASPACE8xx5[,"SubclassFirst"], MatrixOfKoeberlinMatchesMETASPACE8xx5[,"SubclassSecond"], sep = ";"))
i <- "SM;Cer"

for(i in unique(MatrixOfKoeberlinMatchesMETASPACE8xx8[,"SubclassLinks"])){
  if(!all(is.na(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xx8[MatrixOfKoeberlinMatchesMETASPACE8xx8[,"SubclassLinks"] == i,8])))){
    
    plot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xx8[MatrixOfKoeberlinMatchesMETASPACE8xx8[,"SubclassLinks"] == i,8]), as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xx8[MatrixOfKoeberlinMatchesMETASPACE8xx8[,"SubclassLinks"] == i,3]), ylim = c(-1,1), xlim = c(0,1), main = i, ylab = "co-regulation", xlab = "Manders' co-occurrence", col = rgb(0, 0, 0, alpha = 0.2), pch = 16) 
  }}

METASPACEToKoeberlinMatchesDataframe <- data.frame(MatrixOfKoeberlinMatchesMETASPACE8xx8)
METASPACEToKoeberlinMatchesDataframe[,c(3,8)] <- as.numeric(as.character(unlist(METASPACEToKoeberlinMatchesDataframe[,c(3,8)])))

colnames(METASPACEToKoeberlinMatchesDataframe)[c(3,8)] <- c("Coregulation", "Cooccurrence")


fit <- lm(Cooccurrence ~ Coregulation, data = METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == "Cer") & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == "Cer"),])
fit2 <- lm(Coregulation ~ Cooccurrence, data = METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == "Cer") & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == "Cer"),])

TestCorPCCer <- cor.test(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == "PC_aa") & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == "Cer"), "Cooccurrence"], METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == "PC_aa") & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == "Cer"),"Coregulation"],  method = "pearson", use = "complete.obs")


QuantitativeComparisonOfSubclasses <- do.call("rbind", lapply(unique(METASPACEToKoeberlinMatchesDataframe$SubclassSecond), function(y){do.call("rbind", lapply(unique(METASPACEToKoeberlinMatchesDataframe$SubclassFirst), function(x){if(!all(is.na(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y), "Cooccurrence"]))){c(x,y,
                                                                                                                                                                                                                                                                                                                                                                                                                                       unlist(cor.test(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y), "Cooccurrence"], METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y),"Coregulation"],  method = "pearson", use = "complete.obs")[c("estimate", "p.value", "parameter","conf.int")]),
                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                       CoregulationMean = mean(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y), "Coregulation"], na.rm = TRUE),
                                                                                                                                                                                                                                                                                                                                                                                                                                       CoregulationMedian = median(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y), "Coregulation"], na.rm = TRUE),
                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                       CooccurrenceMean = mean(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y), "Cooccurrence"], na.rm = TRUE),
                                                                                                                                                                                                                                                                                                                                                                                                                                       CooccurrenceMedian = median(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y), "Cooccurrence"], na.rm = TRUE),
                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                       CoregulationRange = range(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y), "Coregulation"], na.rm = TRUE),
                                                                                                                                                                                                                                                                                                                                                                                                                                       CooccurrenceRange = range(METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y), "Cooccurrence"], na.rm = TRUE),
                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                       LinearRegressionYCooccurrenceVsXCoregulation = lm(Cooccurrence ~ Coregulation, data = METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y),])$coefficients,
                                                                                                                                                                                                                                                                                                                                                                                                                                       LinearRegressionYCoregulationVsXCooccurrence = lm(Coregulation ~ Cooccurrence, data = METASPACEToKoeberlinMatchesDataframe[(METASPACEToKoeberlinMatchesDataframe$SubclassFirst == x) & (METASPACEToKoeberlinMatchesDataframe$SubclassSecond == y),])$coefficients)
  
}}))}))
QuantitativeComparisonOfSubclassesDataframe <- data.frame(QuantitativeComparisonOfSubclasses)

QuantitativeComparisonOfSubclassesDataframe[,1:2] <- as.character(unlist(QuantitativeComparisonOfSubclassesDataframe[,1:2]))
QuantitativeComparisonOfSubclassesDataframe[,-(1:2)] <- as.numeric(as.character(unlist(QuantitativeComparisonOfSubclassesDataframe[,-(1:2)])))


QuantitativeComparisonOfSubclassesBroad <- dcast(QuantitativeComparisonOfSubclassesDataframe, V1 ~ V2, value.var = "estimate.cor")

QuantitativeComparisonOfSubclassesBroad2 <- QuantitativeComparisonOfSubclassesBroad[,-1]
rownames(QuantitativeComparisonOfSubclassesBroad2) <- QuantitativeComparisonOfSubclassesBroad[,1]

QuantitativeComparisonOfSubclassesBroad4 <- as.matrix(QuantitativeComparisonOfSubclassesBroad2)



CorelationsToMakeNetwork <- do.call("cbind", list(Corelation = QuantitativeComparisonOfSubclassesBroad4[upper.tri(QuantitativeComparisonOfSubclassesBroad4, diag = TRUE)],
                                                  
                                                  RowNames = rownames(QuantitativeComparisonOfSubclassesBroad4)[unlist(lapply(1:nrow(QuantitativeComparisonOfSubclassesBroad4),function(x){1:x}))],
                                                  ColNames = colnames(QuantitativeComparisonOfSubclassesBroad4)[unlist(lapply(1:ncol(QuantitativeComparisonOfSubclassesBroad4),function(x){rep(x,x)}))]))

# Save the correlations to make the network
write.table(CorelationsToMakeNetwork, file="./OutputFiles/CorelationsToMakeNetwork28052019.csv", sep="\t", row.names = FALSE, quote = FALSE)

QuantitativeComparisonOfSubclassesDataframe$ClassFirst <- METASPACEToKoeberlinMatchesDataframe[match(QuantitativeComparisonOfSubclassesDataframe$V1, METASPACEToKoeberlinMatchesDataframe$SubclassFirst),"ClassFirst"]
QuantitativeComparisonOfSubclassesDataframe$ClassSecond <- METASPACEToKoeberlinMatchesDataframe[match(QuantitativeComparisonOfSubclassesDataframe$V2, METASPACEToKoeberlinMatchesDataframe$SubclassSecond),"ClassSecond"]

QuantitativeComparisonOfSubclassesDataframe$TypeColor <- sapply(rowSums(cbind(QuantitativeComparisonOfSubclassesDataframe$V1 == QuantitativeComparisonOfSubclassesDataframe$V2,QuantitativeComparisonOfSubclassesDataframe$ClassFirst == QuantitativeComparisonOfSubclassesDataframe$ClassSecond)), function(x){switch(x+1, "DarkGrey", "Orange",  "Blue")})


#### Final figure for panel C, and associated statistics. Panel C was afterwards esthetically optimized in Adobe Illustrator software.
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/PairsOfLipidSubclassesDistributionCoregulationAndCooccurrenceMediansAndMeansWithLargerDots01072019WithLinearRegression2.pdf")

plot(QuantitativeComparisonOfSubclassesDataframe$CooccurrenceMean, QuantitativeComparisonOfSubclassesDataframe$CoregulationMean, xlab = "Co-occurrence mean", ylab = "Co-regulation mean", xlim = c(0,1), ylim = c(-1,1), main = "Pairs of lipid sub-classes", col = adjustcolor(sapply(QuantitativeComparisonOfSubclassesDataframe$TypeColor, function(x){switch(x, Blue = "Brown", Orange = "#E7B800", DarkGrey = "#00AFBB")}), alpha.f = 0.4), pch = sapply(QuantitativeComparisonOfSubclassesDataframe$TypeColor, function(x){switch(x, Blue = 17, Orange = 15, DarkGrey = 16)}), cex = 2.5)
legend("topleft", c("Species", "Sub-class", "Class"), pch = c(17, 15,  16), col = adjustcolor(c("brown", "#E7B800", "#00AFBB"), alpha.f = 0.64), cex = 1, bty = "n")

LmOfMeanCooccurrencesAndMeanCoregulations <- lm(formula = QuantitativeComparisonOfSubclassesDataframe$CoregulationMean ~ QuantitativeComparisonOfSubclassesDataframe$CooccurrenceMean)
points(seq(0,1,0.1), LmOfMeanCooccurrencesAndMeanCoregulations$coefficients[2]*seq(0,1,0.1) + LmOfMeanCooccurrencesAndMeanCoregulations$coefficients[1], col = "black", type = "l", lwd = 2, lty = 2)

dev.off()
cor.test(QuantitativeComparisonOfSubclassesDataframe$CooccurrenceMean, QuantitativeComparisonOfSubclassesDataframe$CoregulationMean)


#### Base of figure panels a and b, and associated statistical tests

pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/BeanPlotsForCooccurrencesForTopCoregulationsVsAllOfThem23052019.pdf")
# Make beanplots for top 60% in their blue("#50afff") per step of 5%, and general distribution in dark grey & linear regression line

z <- 0.6
MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAbove <- MatrixOfKoeberlinMatchesMETASPACE8xx[as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xx[,"Coregulation"]) > z,]

MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached <- rbind(cbind(MatrixOfKoeberlinMatchesMETASPACE8xx, NewBin = 1),
                                                                                cbind(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAbove, NewBin = as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAbove[,"CoregulationBin"]) - min(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAbove[,"CoregulationBin"])) + 2))


BData <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]), 
                  
                  main = paste0("Validation METASPACE with K?berlin data (bin means)(", z ,")"), side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                  col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                  
                  axes=F,
                  horizontal = TRUE,
                  
                  cutmin = 0,
                  cutmax = 1,
                  
                  border = FALSE,
                  overalline = "mean", what = c(FALSE,TRUE,TRUE,TRUE),
                  
                  beanlines = "mean")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))+0.5), labels = round(seq(from = z, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


BData <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]), 
                  
                  main = paste0("Validation METASPACE with K?berlin data (bin means)(", z ,")"), side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                  col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                  
                  axes=F,
                  horizontal = TRUE,
                  
                  cutmin = 0,
                  cutmax = 1,
                  
                  border = FALSE,
                  overalline = "mean", what = c(FALSE,TRUE,TRUE,TRUE),
                  
                  beanlines = "mean")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))+0.5), labels = round(seq(from = z, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


test <- data.frame(means = BData$stats[-1], bins = 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])))

f <- lm(formula = test$bins ~ test$means)
points(sapply(2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])),function(y){(y - f$coefficients[1])/f$coefficients[2]}), 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])), col = rgb(255/255,140/255,0, alpha = 0.7), type = "l", lwd = 2)

MeanCooccurenceObservedVsCalculated <- cbind(BData$stats[2:9], sapply(2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])),function(y){(y - f$coefficients[1])/f$coefficients[2]}))
colnames(MeanCooccurenceObservedVsCalculated) <- c("MeanObserved", "MeanEstimated")

sum((MeanCooccurenceObservedVsCalculated[,"MeanObserved"] - MeanCooccurenceObservedVsCalculated[,"MeanEstimated"])^2) #0.0006305962 #With X-values (also below)
sum((MeanCooccurenceObservedVsCalculated[,"MeanObserved"] - mean(MeanCooccurenceObservedVsCalculated[,"MeanObserved"]))^2) #0.02151885 (with mean of observed Y)

RSquaredOfCooccurrenceFit <- 1 - (sum((MeanCooccurenceObservedVsCalculated[,"MeanObserved"] - MeanCooccurenceObservedVsCalculated[,"MeanEstimated"])^2)/sum((MeanCooccurenceObservedVsCalculated[,"MeanObserved"] - mean(MeanCooccurenceObservedVsCalculated[,"MeanObserved"]))^2))
cor.test(BData$stats[2:9], 2:9)


# Pearson's product-moment correlation
# 
# data:  BData$stats[2:9] and 2:9
# t = 14.309, df = 6, p-value = 7.29e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9199776 0.9975012
# sample estimates:
#       cor 
# 0.9856622 


BData2 <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]), 
                   
                   main = paste0("Validation METASPACE with K?berlin data (bin medians)(", z ,")"), side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                   col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                   
                   axes=F,
                   horizontal = TRUE,
                   
                   cutmin = 0,
                   cutmax = 1,
                   
                   border = FALSE,
                   overalline = "median", what = c(FALSE,TRUE,TRUE,TRUE),
                   
                   beanlines = "median")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))+0.5), labels = round(seq(from = z, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


BData2 <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]), 
                   
                   main = paste0("Validation METASPACE with K?berlin data (bin medians)(", z ,")"), side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                   col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                   
                   axes=F,
                   horizontal = TRUE,
                   
                   cutmin = 0,
                   cutmax = 1,
                   
                   border = FALSE,
                   overalline = "median", what = c(FALSE,TRUE,TRUE,TRUE),
                   
                   beanlines = "median")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))+0.5), labels = round(seq(from = z, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


test2 <- data.frame(medians = BData2$stats[-1], bins = 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])))

f2 <- lm(formula = test2$bins ~ test2$medians)
points(sapply(2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])),function(y){(y - f2$coefficients[1])/f2$coefficients[2]}), 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])), col = rgb(255/255,140/255,0, alpha = 0.7), type = "l", lwd = 2)


beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ ifelse(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,"F","A"), 
         
         main = "Co-occurrence distribution comparison (with means)", side = "both", xlab="Manders' co-occurrence of METASPACE lipids", ylab = "Relative density", ll = 0.02, wd = 0.32,
         col = list("DarkGrey", "#50afff"), method = "overplot",
         
         axes=F,
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "mean", what = c(FALSE,TRUE,TRUE,TRUE),
         
         beanlines = "mean")


axis(1)
axis(side = 4, at = c(0.74,1.35), labels = c("All co-regulations", paste0("High co-regulations (>", z ,")")), tick = FALSE, mgp = c(3,0,0))

PreferencesOfTickPositions <- seq(-0.5, 0.5, by = 0.25)
axis(side = 2, at = PreferencesOfTickPositions+1, labels = abs(PreferencesOfTickPositions))

ks.test(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,8]),
        as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] <= 1,8]))

ks.test(na.omit(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,8])),
        na.omit(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] <= 1,8])))

ks.test(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,8]),
        as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] <= 1,8]), exact = TRUE)

# Generally: 2-side with or without NAs: D = 0.21289, p-value < 2.2*10-16, but not exact p-value possible because of ties


ks.test(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,8]),
        as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] <= 1,8]), alternative = "greater")

ks.test(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,8]),
        as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] <= 1,8]), alternative = "less")

# D^- = 0.21289, p-value < 2.2e-16
# D^+ = 0.0010896, p-value = 0.9925


beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ ifelse(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,"F","A"), 
         
         main = "Co-occurrence distribution comparison (with medians)", side = "both", xlab="Manders' co-occurrence of METASPACE lipids", ylab = "Relative density", ll = 0.02, wd = 0.32,
         col = list("DarkGrey", "#50afff"), method = "overplot",
         
         axes=F,
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "median", what = c(FALSE,TRUE,TRUE,TRUE),
         
         beanlines = "median")


axis(1)
axis(side = 4, at = c(0.74,1.35), labels = c("All co-regulations", paste0("High co-regulations (>", z ,")")), tick = FALSE, mgp = c(3,0,0))

PreferencesOfTickPositions <- seq(-0.5, 0.5, by = 0.25)
axis(side = 2, at = PreferencesOfTickPositions+1, labels = abs(PreferencesOfTickPositions))

# Make beanplots for top 65% in their blue per step of 5%, and general distribution in dark grey & linear regression line
z <- 0.65 # Alternative with nearly exactly the same code as in previous part but with a few changes: just shown here in the way that it was re-used


MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAbove <- MatrixOfKoeberlinMatchesMETASPACE8xx[as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xx[,"Coregulation"]) > z,]

MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached <- rbind(cbind(MatrixOfKoeberlinMatchesMETASPACE8xx, NewBin = 1),
                                                                                cbind(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAbove, NewBin = as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAbove[,"CoregulationBin"]) - min(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAbove[,"CoregulationBin"])) + 2))


BData <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]), 
                  
                  main = paste0("Validation METASPACE with K?berlin data (bin means)(", z ,")"), side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                  col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                  
                  axes=F,
                  horizontal = TRUE,
                  
                  cutmin = 0,
                  cutmax = 1,
                  
                  border = FALSE,
                  overalline = "mean", what = c(FALSE,TRUE,TRUE,TRUE),
                  
                  beanlines = "mean")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))+0.5), labels = round(seq(from = z, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


BData <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]), 
                  
                  main = paste0("Validation METASPACE with K?berlin data (bin means)(", z ,")"), side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                  col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                  
                  axes=F,
                  horizontal = TRUE,
                  
                  cutmin = 0,
                  cutmax = 1,
                  
                  border = FALSE,
                  overalline = "mean", what = c(FALSE,TRUE,TRUE,TRUE),
                  
                  beanlines = "mean")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))+0.5), labels = round(seq(from = z, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


test <- data.frame(means = BData$stats[-1], bins = 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])))

f <- lm(formula = test$bins ~ test$means)
points(sapply(2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])),function(y){(y - f$coefficients[1])/f$coefficients[2]}), 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])), col = rgb(255/255,140/255,0, alpha = 0.7), type = "l", lwd = 2)


BData2 <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]), 
                   
                   main = paste0("Validation METASPACE with K?berlin data (bin medians)(", z ,")"), side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                   col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                   
                   axes=F,
                   horizontal = TRUE,
                   
                   cutmin = 0,
                   cutmax = 1,
                   
                   border = FALSE,
                   overalline = "median", what = c(FALSE,TRUE,TRUE,TRUE),
                   
                   beanlines = "median")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))+0.5), labels = round(seq(from = z, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


BData2 <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]), 
                   
                   main = paste0("Validation METASPACE with K?berlin data (bin medians)(", z ,")"), side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                   col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                   
                   axes=F,
                   horizontal = TRUE,
                   
                   cutmin = 0,
                   cutmax = 1,
                   
                   border = FALSE,
                   overalline = "median", what = c(FALSE,TRUE,TRUE,TRUE),
                   
                   beanlines = "median")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"]))+0.5), labels = round(seq(from = z, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


test2 <- data.frame(medians = BData2$stats[-1], bins = 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])))

f2 <- lm(formula = test2$bins ~ test2$medians)
points(sapply(2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])),function(y){(y - f2$coefficients[1])/f2$coefficients[2]}), 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"])), col = rgb(255/255,140/255,0, alpha = 0.7), type = "l", lwd = 2)


beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ ifelse(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,"F","A"), 
         
         main = "Co-occurrence distribution comparison (with means)", side = "both", xlab="Manders' co-occurrence of METASPACE lipids", ylab = "Relative density", ll = 0.02, wd = 0.32,
         col = list("DarkGrey", "#50afff"), method = "overplot",
         
         axes=F,
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "mean", what = c(FALSE,TRUE,TRUE,TRUE),
         
         beanlines = "mean")


axis(1)
axis(side = 4, at = c(0.74,1.35), labels = c("All co-regulations", paste0("High co-regulations (>", z ,")")), tick = FALSE, mgp = c(3,0,0))

PreferencesOfTickPositions <- seq(-0.5, 0.5, by = 0.25)
axis(side = 2, at = PreferencesOfTickPositions+1, labels = abs(PreferencesOfTickPositions))


beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,8]) ~ ifelse(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached[,"NewBin"] > 1,"F","A"), 
         
         main = "Co-occurrence distribution comparison (with medians)", side = "both", xlab="Manders' co-occurrence of METASPACE lipids", ylab = "Relative density", ll = 0.02, wd = 0.32,
         col = list("DarkGrey", "#50afff"), method = "overplot",
         
         axes=F,
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "median", what = c(FALSE,TRUE,TRUE,TRUE),
         
         beanlines = "median")


axis(1)
axis(side = 4, at = c(0.74,1.35), labels = c("All co-regulations", paste0("High co-regulations (>", z ,")")), tick = FALSE, mgp = c(3,0,0))

PreferencesOfTickPositions <- seq(-0.5, 0.5, by = 0.25)
axis(side = 2, at = PreferencesOfTickPositions+1, labels = abs(PreferencesOfTickPositions))

dev.off()


# Make beanplots for top 60% in their blue per step of 5%, and general distribution in dark grey with network data & linear regression line
KoeberlinEdgesOfNetwork16052019 <- read.csv(file = "./InputData/KoeberlinEdgesOfNetwork16052019.csv", header = TRUE, sep = ",", as.is = TRUE)

KoeberlinEdgesOfNetwork160520192 <- cbind(KoeberlinEdgesOfNetwork16052019, do.call("rbind", strsplit(KoeberlinEdgesOfNetwork16052019[,"canonicalName"], " \\(interacts) ")))
colnames(KoeberlinEdgesOfNetwork160520192)[9:10] <- c("FromLipid", "ToLipid")

KoeberlinEdgesOfNetwork160520192$FromLipidWithSpaces <- gsub("_", " ", KoeberlinEdgesOfNetwork160520192[,9])
KoeberlinEdgesOfNetwork160520192$ToLipidWithSpaces <- gsub("_", " ", KoeberlinEdgesOfNetwork160520192[,10]) 


MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached2 <- MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached

ListOfNetworkMatches <- lapply(1:nrow(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached2), function(x){KoeberlinEdgesOfNetwork160520192[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached2[x,"FromLipid"] == KoeberlinEdgesOfNetwork160520192[,"FromLipidWithSpaces"] & MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached2[x,"ToLipid"] == KoeberlinEdgesOfNetwork160520192[,"ToLipidWithSpaces"], "interactionStrength"]
})

ListOfNetworkMatches[ListOfNetworkMatches == "numeric(0)"] <- NA
MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached4 <- cbind(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached2, NetworkPresence = unlist(ListOfNetworkMatches))

MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached5 <- MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached4[!is.na(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached4[,"NetworkPresence"]),]
MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached7 <- MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached5[(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached5[,"NetworkPresence"] > 0.6) & (MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached5[,"NewBin"] == 1),]

MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8 <- rbind(cbind(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached4[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached4[,"NewBin"] == 1,], NewBin2 = 1),
                                                                                 cbind(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached7, NewBin2 = as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached7[,"CoregulationBin"]) - min(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached7[,"CoregulationBin"])) + 2))

MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached7b <- MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached5[(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached5[,"NetworkPresence"] > 0.65) & (MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached5[,"NewBin"] == 1),]


MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b <- rbind(cbind(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached4[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached4[,"NewBin"] == 1,], NewBin2 = 1),
                                                                                  cbind(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached7b, NewBin2 = as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached7b[,"CoregulationBin"]) - min(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached7b[,"CoregulationBin"])) + 2))


#### Visualizations for figure panel b, and alternative visualizations included too (redone from previous entries here above)

pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/BeanPlotsForCooccurrencesForTopCoregulationsVsAllOfThemButWithNetworkSubsectionAbove06527052019.pdf")
BData4 <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]), 
                   
                   main = "Validation METASPACE with K?berlin network data (>0.65)(medians)", side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                   col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                   
                   axes=F,
                   horizontal = TRUE,
                   
                   cutmin = 0,
                   cutmax = 1,
                   
                   border = FALSE,
                   overalline = "median", what = c(FALSE,TRUE,TRUE,TRUE),
                   
                   beanlines = "median")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]))+0.5), labels = round(seq(from = 0.65, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


BData4 <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]), 
                   
                   main = "Validation METASPACE with K?berlin network data (>0.65)(medians)", side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                   col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                   
                   axes=F,
                   horizontal = TRUE,
                   
                   cutmin = 0,
                   cutmax = 1,
                   
                   border = FALSE,
                   overalline = "median", what = c(FALSE,TRUE,TRUE,TRUE),
                   
                   beanlines = "median")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]))+0.5), labels = round(seq(from = 0.65, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


test4 <- data.frame(medians = BData4$stats[-1], bins = 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"])))

f4 <- lm(formula = test4$bins ~ test4$medians)
points(sapply(2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"])),function(y){(y - f4$coefficients[1])/f4$coefficients[2]}), 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"])), col = rgb(255/255,140/255,0, alpha = 0.7), type = "l", lwd = 2)


BData4_4 <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]), 
                     
                     main = "Validation METASPACE with K?berlin network data (>0.65)(means)", side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                     col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                     
                     axes=F,
                     horizontal = TRUE,
                     
                     cutmin = 0,
                     cutmax = 1,
                     
                     border = FALSE,
                     overalline = "mean", what = c(FALSE,TRUE,TRUE,TRUE),
                     
                     beanlines = "mean")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]))+0.5), labels = round(seq(from = 0.65, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


BData4_4 <- beanplot(as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,8]) ~ as.numeric(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]), 
                     
                     main = "Validation METASPACE with K?berlin network data (>0.65)(means)", side = "second", xlab="Manders' co-occurrence of METASPACE lipids", ll = 0.02, wd = 0.53,
                     col = c("DarkGrey", as.list(rep("#50afff", length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]))-1))), method = "overplot", ylab = "Co-regulation of lipids (K?berlin data)",
                     
                     axes=F,
                     horizontal = TRUE,
                     
                     cutmin = 0,
                     cutmax = 1,
                     
                     border = FALSE,
                     overalline = "mean", what = c(FALSE,TRUE,TRUE,TRUE),
                     
                     beanlines = "mean")
axis(side = 2, at = 1, labels = "All data", las = 2, cex.lab = 0.64, cex.axis = 0.64)

axis(side = 2, at = 1.5:(length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"]))+0.5), labels = round(seq(from = 0.65, to = 1, by = 0.05), digits = 2), las = 2, cex.lab = 0.64, cex.axis = 0.64, col = "#50afff")
axis(1, cex.lab = 0.64, cex.axis = 0.64)


test4_4 <- data.frame(means = BData4_4$stats[-1], bins = 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"])))

f4_4 <- lm(formula = test4_4$bins ~ test4_4$means)
points(sapply(2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"])),function(y){(y - f4_4$coefficients[1])/f4_4$coefficients[2]}), 2:length(unique(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,"NewBin2"])), col = rgb(255/255,140/255,0, alpha = 0.7), type = "l", lwd = 2)

dev.off()


MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached10b <- do.call("cbind", list(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b, From_Lipid = gsub(" ", "_", MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,1]), To_Lipid = gsub(" ", "_", MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached8b[,2])))
colnames(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached10b) <- c("FromLipid","ToLipid","Coregulation","FromLipidChemicalFormula", "ToLipidChemicalFormula", "FromLipidChemicalFormulaMatch", "ToLipidChemicalFormulaMatch","Cooccurrence", "FromHit", "ToHit", "FromHit2","ToHit2","CoregulationBin","NewBin","NetworkPresence","NewBin2","From_Lipid","To_Lipid")

# Export data for the network visualization in Cytoscape (first set)
write.table(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached10b, file="./MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached10b03062019.csv", sep="\t", row.names = FALSE, quote = FALSE)

MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b <- MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached10b[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached10b[,"NewBin2"] == 1,]
UniqueInverseEntriesList <- lapply(1:nrow(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b), function(x){
  
  y <- which((MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b[x,1] == MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b[,2]) &
               (MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b[x,2] == MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b[,1]))
  
  z <- if(sum(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b[c(x,y),"Cooccurrence"] != "numeric(0)") != 1){
    min(c(x,y),na.rm = TRUE)
    
  }else{
    c(x,y)[MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b[c(x,y),"Cooccurrence"] != "numeric(0)"]
    
  }
  MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached11b[z,]})


# Export data for the network visualization in Cytoscape (second set)

MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached13b <- unique(do.call("rbind", UniqueInverseEntriesList)) #59292/29646==2 (exact the double here)
write.table(MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached13b, file="./MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached13bFullyUniqueEntries03062019.csv", sep="\t", row.names = FALSE, quote = FALSE)


#### Import node and edge tables

UniqueEntries03062019DefaultNodeTable <- read.csv(file = "./InputData/MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached13bFullyUniqueEntries03062019DefaultNodeTable.csv", header = TRUE, sep = ",", as.is = TRUE)
UniqueEntries03062019DefaultEdgeTable <- read.csv(file = "./InputData/MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached13bFullyUniqueEntries03062019DefaultEdgeTable.csv", header = TRUE, sep = ",", as.is = TRUE)

#### Send chemical formulas of the lipids of K?berlin dataset for use in parts outside of R
KoeberlinNetworkChemicalFormulas <- unique(c(UniqueEntries03062019DefaultNodeTable$FromLipidChemicalFormula, UniqueEntries03062019DefaultNodeTable$ToLipidChemicalFormula))

KoeberlinNetworkChemicalFormulas2 <- KoeberlinNetworkChemicalFormulas[KoeberlinNetworkChemicalFormulas != ""]
write.table(KoeberlinNetworkChemicalFormulas2, file="./KoeberlinNetworkChemicalFormulas211062019.tsv", sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


#### Import of datasets derived from visualizations in Cytoscape and further work on it for creation of outer circular projections in panel d

# XGMLL from Cytoscape network to txt and then import in R to extract x and y coordinates (Also set number of spokes/sections for final circles around panel d)
XGMLLToTxtInput <- read.csv(file = "./InputData/TabDelimitedConversionOfXGMLLOfCytoscapeNetworkInvertedColors12062019.txt", header = TRUE, sep = "\t", as.is = TRUE)

XGMLLToTxtInput2 <- XGMLLToTxtInput[XGMLLToTxtInput$x != "",]
XGMLLToTxtInput4 <- unique(XGMLLToTxtInput2[,c("label11","x","fill","y")])

XGMLLToTxtInput4$x <- as.numeric(gsub(",",".",XGMLLToTxtInput4$x))
XGMLLToTxtInput4$y <- as.numeric(gsub(",",".",XGMLLToTxtInput4$y))

MiddleX <- sum(range(XGMLLToTxtInput4$x))/2
MiddleY <- sum(range(XGMLLToTxtInput4$y))/2

EDistIn2D <- function(x1,y1,x2,y2){sqrt((y1 - y2)^2 + (x1 - x2)^2)}
XGMLLToTxtInput4$EDist <- EDistIn2D(XGMLLToTxtInput4$x, XGMLLToTxtInput4$y, MiddleX, MiddleY)

# max(XGMLLToTxtInput4$EDist) # 793.9823
# max(XGMLLToTxtInput4$EDist)*1.1

AngularTranslocation <- function(x,y,Theta){c(x*cos(Theta)-y*sin(Theta), x*sin(Theta)+y*cos(Theta))}
AngularXYPositions <- function(StartX, StartY, NumSpokes){do.call("rbind", lapply(0:(NumSpokes-1), function(z){round(AngularTranslocation(StartX,StartY,z*2*pi/NumSpokes), digits = 7)}))}

NumSpokes <- 64 
SpokesToCrossNetworkFrom0 <- AngularXYPositions(max(XGMLLToTxtInput4$EDist)*1.1,0,NumSpokes)


# Import of the edge table of the network visualization in cytoscape & preparative steps for further calculations

EdgeTable12062019 <- read.csv(file = "./InputData/MatrixOfKoeberlinMatchesMETASPACE8xxSubsetAboveWithGeneralDataAttached13bFullyUniqueEntries03062019.csv_2(1)_2(1)_6DefaultEdgesNetwork.csv", header = TRUE, sep = ",", as.is = TRUE)
EdgeTable120620192 <- cbind(EdgeTable12062019, do.call("rbind", strsplit(EdgeTable12062019$name, " \\(\\(interacts)) ")))

EdgeTable120620194 <- cbind(EdgeTable120620192, do.call("rbind", lapply(1:nrow(EdgeTable120620192), function(x){unlist(c(XGMLLToTxtInput4[XGMLLToTxtInput4[,"label11"] == EdgeTable120620192[x,15],c("x","y")], XGMLLToTxtInput4[XGMLLToTxtInput4[,"label11"] == EdgeTable120620192[x,16], c("x","y")]))})))
colnames(EdgeTable120620194)[15:20] <- c("FromLipid","ToLipid","FromX","FromY","ToX","ToY")

EdgeTable120620195 <- cbind(EdgeTable120620194, do.call("rbind", lapply(1:nrow(EdgeTable120620194), function(x){sapply(1:NumSpokes, function(z){all(is.finite(line.line.intersection(c(EdgeTable120620194[x, "FromX"], EdgeTable120620194[x,"FromY"]), c(EdgeTable120620194[x,"ToX"], EdgeTable120620194[x,"ToY"]), c(MiddleX, MiddleY), (c(MiddleX, MiddleY) + SpokesToCrossNetworkFrom0[z,]), interior.only = TRUE)))})})))
EdgeTable120620195$Cooccurrence[is.na(EdgeTable120620195$Cooccurrence)] <- 0

MedianCooccurrencesForSpokes <- sapply(1:NumSpokes, function(x){median(EdgeTable120620195[EdgeTable120620195[,as.character(x)],"Cooccurrence"])}) # up-down-up trend
CooccurrencePerSpokeLine <- do.call("rbind",lapply(1:NumSpokes, function(x){cbind(x, EdgeTable120620195[EdgeTable120620195[,as.character(x)],"Cooccurrence"])}))

QuantilesCooccurrencesForSpokes <- do.call("cbind", lapply(1:NumSpokes, function(x){quantile(EdgeTable120620195[EdgeTable120620195[,as.character(x)],"Cooccurrence"])}))
colnames(QuantilesCooccurrencesForSpokes) <- 1:NumSpokes

# Intermediate visualization alternatives: fraction stacked(binned & maybe circular), boxplot, stripchart, beanplots (maybe like flower), radarplot (maybe both medians)
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/NetworkCooccurrenceDistributionVisualizations16130620192QuantileRadarsAdded.pdf")

boxplot(CooccurrencePerSpokeLine[,2] ~ factor(CooccurrencePerSpokeLine[,1], levels = 1:NumSpokes), xlab = "Network side", ylab = "Edge co-occurrence") #, outline = TRUE, boxlty = 0, whisklty = 0, staplelty = 0)
stripchart(CooccurrencePerSpokeLine[,2] ~ factor(CooccurrencePerSpokeLine[,1], levels = 1:NumSpokes), xlab = "Edge co-occurrence", ylab = "Network side")


beanplot(CooccurrencePerSpokeLine[,2] ~ CooccurrencePerSpokeLine[,1], 
         
         main = "Edge distribution of co-occurrences over the network", side = "second", xlab="Manders' co-occurrence", ll = 0.02, wd = 0.32,
         col = "#50afff", method = "overplot", ylab = "Network side",
         
         axes=F,
         horizontal = TRUE,
         
         cutmin = 0,
         cutmax = 1,
         
         border = FALSE,
         overalline = "median", what = c(TRUE,TRUE,TRUE,TRUE),
         
         beanlines = "median")
axis(1)

axis(side = 2, at = 1:NumSpokes, labels = 1:NumSpokes)



radarchart(as.data.frame(rbind(1,0,MedianCooccurrencesForSpokes)[,c(((3*NumSpokes/4)+1):1, NumSpokes:((3*NumSpokes/4)+2))]), axistype=1, axislabcol=rgb(0,0,0,0.5), caxislabels=seq(0,100,25), vlabels = "")

radarchart(as.data.frame(rbind(1,0,QuantilesCooccurrencesForSpokes[2:4,])[,c(((3*NumSpokes/4)+1):1, NumSpokes:((3*NumSpokes/4)+2))]), axistype=1, axislabcol=rgb(0,0,0,0.5), caxislabels=seq(0,100,25), vlabels = "", pcol = colorRampPalette(c(rgb(255/255,204/255,0),rgb(0,153/255,153/255)))(3))
radarchart(as.data.frame(rbind(1,0,MedianCooccurrencesForSpokes)[,c(((3*NumSpokes/4)+1):1, NumSpokes:((3*NumSpokes/4)+2))]), axistype=1, 
           
           #custom polygon
           pcol=rgb(0.31,0.69,1,0.9) , pfcol=rgb(0.31,0.69,1,0.5) , plwd=4, 
           
           #custom the grid
           cglcol="grey", cglty=1, axislabcol=rgb(0,0,0,0.5), caxislabels=seq(0,100,25), cglwd=0.8,
           
           #custom labels
           vlcex=0.8, vlabels = "")


radarchart(as.data.frame(rbind(1,0,QuantilesCooccurrencesForSpokes[2:4,])[,c(((3*NumSpokes/4)+1):1, NumSpokes:((3*NumSpokes/4)+2))]), axistype=1, 
           
           #custom polygon
           pcol= rgb(0.31,0.69,1,0) , pfcol= rgb(0.31,0.69,1,0.32) , plwd=4, 
           
           #custom the grid
           cglcol="grey", axislabcol=rgb(0,0,0,0.5), caxislabels=seq(0,100,25), cglwd=0.1,
           
           #custom labels
           vlcex=0.8, vlabels = "")

dev.off()


CutOffRangesMatrix <- cbind(seq(0,0.9,0.1),seq(0.1,1,0.1))
x <- 1; y <- 1

CooccurrencePerSpokeLine2 <- cbind(CooccurrencePerSpokeLine, CooccurrenceBin = sapply(CooccurrencePerSpokeLine[,2], function(x){which((CutOffRangesMatrix[,1] <= x) & (x < CutOffRangesMatrix[,2]))}))
colnames(CooccurrencePerSpokeLine2) <- c("Spoke", "Cooccurrence","CooccurrenceBin")


OverviewMatrixContentCooccurrenceBins <- dcast(as.data.frame(CooccurrencePerSpokeLine2[,c(1,3)]), CooccurrenceBin ~ Spoke)

OverviewMatrixContentCooccurrenceBins2 <- as.matrix(OverviewMatrixContentCooccurrenceBins[,-1])
rownames(OverviewMatrixContentCooccurrenceBins2) <- OverviewMatrixContentCooccurrenceBins[,1]

# Intermediate figure and calculations
pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/BarplotsForBinnedDistributionsOfCooccurrencesInDirectionsOfNetwork1613062019.pdf")

barplot(OverviewMatrixContentCooccurrenceBins2, col = colorRampPalette(c(rgb(255/255,204/255,0),rgb(0,153/255,153/255)))(nrow(OverviewMatrixContentCooccurrenceBins2)), border="white", space=0.04, xlab = "Direction", ylab = "Frequency")
OverviewMatrixContentCooccurrenceBins2ColNormalized <- apply(OverviewMatrixContentCooccurrenceBins2, 2, function(y){y/sum(y)})

barplot(OverviewMatrixContentCooccurrenceBins2ColNormalized, col = colorRampPalette(c(rgb(255/255,204/255,0),rgb(0,153/255,153/255)))(nrow(OverviewMatrixContentCooccurrenceBins2)), border="white", space=0.04, xlab="Direction", ylab = "Relative frequency")
dev.off()

# OverviewMatrixContentCooccurrenceBins2m <- melt(OverviewMatrixContentCooccurrenceBins2)
# colnames(OverviewMatrixContentCooccurrenceBins2m) <- c("CooccurrenceBin","Spoke","Count")

OverviewMatrixContentCooccurrenceBins2ColNormalizedm <- melt(OverviewMatrixContentCooccurrenceBins2ColNormalized)
colnames(OverviewMatrixContentCooccurrenceBins2ColNormalizedm) <- c("CooccurrenceBin","Spoke","Count")

ColSumEdgeTable180620195 <- as.data.frame(cbind((1:NumSpokes), colSums(EdgeTable120620195[,(20+(1:NumSpokes))])))
colnames(ColSumEdgeTable180620195) <- c("Spoke", "EdgesAmount")


#### Visualization of the outer circles of panel d

pdf("C:/Users/Kevin/Documents/RData/Rest2/ACGData/FiguresByKT/LogarithmicEdgeAbundanceCircleWithGreyScale3225062019NowWithViridisOuterRingToo2.pdf")
ggplot(ColSumEdgeTable180620195) +      
  
  
  geom_bar(aes(x=as.factor(Spoke), y= 1, fill= log10(EdgesAmount)), stat="identity", alpha=0.8) +

  scale_fill_gradientn(colours = c("#FAFAFA", "black")) + 
  ylim(-32,1) +
  
  theme_minimal() +
  theme(
    
    #   legend.position = "none",
    axis.text = element_blank(),
    
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  
  coord_polar(start = ((NumSpokes/2)-1)*pi/NumSpokes, direction = 1)



ggplot(OverviewMatrixContentCooccurrenceBins2ColNormalizedm) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(Spoke), y=Count, fill= as.factor(CooccurrenceBin)), stat="identity", alpha=0.8) +
  
  scale_fill_viridis(discrete=TRUE, direction = -1) + 
  # scale_fill_manual(values = OrangeToDarkTealPallette) +
  
  
  ylim(-8,1) +
  
  theme_minimal() +
  theme(
    
    legend.position = "none",
    axis.text = element_blank(),
    
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  
  coord_polar(start = ((NumSpokes/2)-1)*pi/NumSpokes, direction = 1)
dev.off()