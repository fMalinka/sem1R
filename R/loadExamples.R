getDatasetExample <- function()
{

dotOBOpath <- system.file('extdata', 'dotOntology.obo', package="sem1R")
goOBOpath <- system.file('extdata', 'go-basic-reduced.obo', package="sem1R")
mymatrixpath <- system.file('extdata', 'dotmatrix.csv', package="sem1R")
initColpath <- system.file('extdata', 'initColDot.csv', package="sem1R")
initRowpath <- system.file('extdata', 'initRowDot_reduced.csv', package="sem1R")

mystruct <- list()
if(dotOBOpath == "" || goOBOpath == "" || mymatrixpath == "" ||
initColpath == "" || initRowpath =="")
{
mystruct <- NULL
}
else
{
mystruct[["colOntologyPath"]] <- dotOBOpath
mystruct[["rowOntologyPath"]] <- goOBOpath
mymatrix <- read.csv(mymatrixpath, stringsAsFactors = FALSE)
rownames(mymatrix) <- mymatrix[,1]
mymatrix <- mymatrix[,-1]
mystruct[["datamatrix"]] <- as.matrix(mymatrix)
initCol <- read.csv(initColpath, stringsAsFactors = FALSE, header = FALSE)
initRow <- read.csv(initRowpath, stringsAsFactors = FALSE, header = FALSE)

colList <- list()
for(icol in 1:nrow(initCol))
{
myvct <- as.character(initCol[icol,-1])
myvct <- myvct[myvct!=""]
colList <- append(colList, list(myvct))
names(colList)[icol] <- initCol[icol,1]
}
mystruct[["colOntologyDesc"]] <- colList

rowList <- list()
for(irow in 1:nrow(initRow))
{
myvct <- as.character(initRow[irow,-1])
myvct <- myvct[myvct!="" & myvct!="NA"]
rowList <- append(rowList, list(myvct))
names(rowList)[irow] <- initRow[irow,1]
}
mystruct[["rowOntologyDesc"]] <-  rowList
}

return(mystruct)
}


