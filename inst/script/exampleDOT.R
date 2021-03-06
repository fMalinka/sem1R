library(sem1R)
mysem1R <- new(sem1R)
myExample <- getDatasetExample()
mysem1R$setDataset(myExample$datamatrix)
mysem1R$createCOLOntology("DOT", myExample$colOntologyPath, myExample$colOntologyDesc)
mysem1R$createROWOntology("GO", myExample$rowOntologyPath, myExample$rowOntologyDesc)
mysem1R$filterTh <- 50
mysem1R$objective <- "auc"
mysem1R$ruleDepth <- 3
mysem1R$nrules <- 2
mysem1R$featureSelectionMethod <- 0
mysem1R$minLevel <- 2
myhypothesis <- mysem1R$findDescription()

str(myhypothesis)