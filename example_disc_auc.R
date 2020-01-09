library(sem1R)
mysem1R <- new(sem1R)

####################### FUNTIONS ###################################

PrepareOntologyDesc <- function(mydata, colCSV = "example/initsegmentFBbtWithoutComments.csv", geneASOC = "example/gene_association.fb")
{
  #http://www.geneontology.org/gene-associations/readme/fb.README
  mapFbgn2GO <- read.delim(geneASOC, header = FALSE)
  colnames(mapFbgn2GO) <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier",
                            "GO_ID", "DB:Reference", "Evidence", "From", "Aspect",
                            "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type",
                            "Taxon", "Date", "Assigned_by", "Annotation_Extension",
                            "Gene_Product_Form_ID")
  mapFbbt2GO <- read.csv(colCSV, header = FALSE, row.names = 1)
  
  #PREPARE COLUMN DESCRIPTIONS
  colDESC_df <- mapFbbt2GO[colnames(mydata),]
  colDESC <- as.list(as.data.frame(t(colDESC_df)))
  #erase empty elements
  for(i in 1:length(colDESC))
  {
    colDESC[[i]] <- levels(colDESC[[i]])[levels(colDESC[[i]]) != ""]
  }
  
  #PREPARE ROW DESCRIPTIONS
  ids <- match(mapFbgn2GO$DB_Object_ID, rownames(mydata))
  a <- mapFbgn2GO[!is.na(ids),]
  
  rowDESC <- list()
  for(i in 1:nrow(mydata))
  {
    rowDESC[[i]] <- unique(as.character(a[which(a$DB_Object_ID == rownames(mydata)[i]), 5]))
  }
  return(list(row=rowDESC, col=colDESC))
}

####################### FUNTIONS END ###################################

mydata <- read.csv("example/discMatrix.csv", header = TRUE, check.names = FALSE, row.names = 1)
mydata <- as.matrix(mydata)
#get onto Desc
ontoDesc <- PrepareOntologyDesc(mydata, geneASOC = "example/gene_association.fb", colCSV = "example/initsegmentFBbtWithoutComments.csv")

#library(lattice)
#levelplot(as.matrix(mydata))

mysem1R$setDataset(mydata)
#load Ontology
rowOntoPath <- "example/go-basic.obo"
colOntoPath <- "example/fbbt-simple.obo"
mysem1R$createCOLOntology("FBGN", colOntoPath, ontoDesc$col)
mysem1R$createROWOntology("GO", rowOntoPath, ontoDesc$row)

mysem1R$ruleFormat <- "both"

mysem1R$filterTh <- 100
mysem1R$objective <- "auc"

mysem1R$ruleDepth <- 9
mysem1R$nrules <- 10
mysem1R$featureSelectionMethod <- 0
mysem1R$minLevel <- 3
myhypothesis <- mysem1R$findDescription()

