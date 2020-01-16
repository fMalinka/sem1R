# sem1R - Concept rule learning with and ontology-based refinement operator inducing semantic rule from omics data

<img align="right" width="200" src="sem1r_logo_col.png">
<p style="text-align: justify">sem1R is a machine learning algorithm that finds interesting, hidden, and non-trivial patterns in omics data. The algorithm produces a set of semantical prediction rules that form data into clusters or biclusters, this depends on a type of ontologies. Here, we distingues between two types of ontologies: an ontology describing rows (e.g. genes) an columns (e.g. samples). Practically, for gene expression data, where rows represent genes and column represent samples, we recommend to use Gene ontology or any pathway ontologies as a row ontology. Choosing a proper column ontology is depending on a type of experiment, e.g. OBO Foundry provides almost two hundreds and many of them are domain specific so some anatomical ontologies can be used as well.

The sem1R is based on rule learning methods, where two reduction procedures were added and that make the algorithm extremely fast and efficient in comparison with traditional approach. In additional, it is relative easy to use, because all important methods are included into the package.</p>

## Getting Started
The algorithm is implemented in C++ and provided as `R` package. The following instructions will show you how to install all prerequisites and the sem1R package as well into your local machine. Afterwards, we will demonstrate the sem1R on real gene expression dataset.

### Prerequisites
We required to use R in version 3.0.2.
All prerequisites R packages that are needed for the sem1R package are the following:
`Rcpp (>= 0.12.6)`, `RcppProgress`, `RcppArmadillo (>= 0.7.800.2.0)`, and `BH (>= 1.72.0-3)`. All of these packages come from CRAN, so install them by `install.packages` function in R.

Or, for easier installation we recommend to install 'devtools' that can download and install the project instantly from gitHub using only one command.

For installing 'devtools' package run `R` and type the following:
```
install.packages("devtools")
```

### Installing
If you choosen instalation via devtools, you would go to the terminal and run `R` and then the following commands:

```
library(devtools)
install_gitHub("fmalinka/sem1R")
```
All prerequisites packages should be installed automatically.

For non-devtools users, download the sem1R package and extract it to your arbitrary folder.
Run the terminal, go to the folder that you have already choosen and install all prerequisities if you have not installed yet using `install.packages` R's function. Go to the folder and build the package.

```
cd /my/path/to/package
R CMD build .
```
The package in tar.gz format will be named as `sem1R_[version].tar.gz'. The concrete name depends on the package version.
Then, install the sem1R package.

```
R CMD INSTALL sem1R_[version].tar.gz
```

And finally, check whether the sem1R package has been installed.
Run R and load the package.
```
library(sem1R)
```


## Running the example
Running example that we present here comes from [1] and shows a gene expression over Drosophila melanogaster imaginal discs. All necessary files for this example are located at [example](example) folder.

### Data matrix
A file [discMatrix.csv](example/discMatrix.csv) contains binary information about gene expression over imaginal discs of Drosophila melanogaster. The matrix is two-dimensional where rows represent genes and columns represent samples (locations). Each dimension has own identifier, i.e. genes are described by FBgn (FlyBase) identifiers and columns by your notation. Ones in the matrix mean "expressed" and zeros mean "non-expressed" in the given positions. Obviously, process of binarization has to be done if your data are not in the binary format.

### Ontologies
Ontologies are the second type of input that has to be given to your algorithm. Ontology has to be in OBO format (https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html) and relationships of terms must be acyclic. For many interesting ontologies look at OBO Foundry (http://obofoundry.org/). In our running example, we provide two type of different ontologies. Gene ontology, located at [example/go-basic.obo](example/go-basic.obo), aims to rows of the data matrix and FBBT ontology (http://obofoundry.org/ontology/fbbt.html), located at [example/fbbt-simple.obo](example/fbbt-simple.obo), focuses on the columns.

### Connection between the data matrix and the ontologies
Now, the last step is to established a connection between our data matrix and all given ontologies. Firstly, we look at the rows which are described by the FBGN identifiers. The mapping from FBGN identifiers to Gene ontology is provided at [example/gene_association.fb](example/gene_association.fb) file. So, for each row we have to find the corresponding GO terms. The example of data matrix and the connection is shown bellow.


```
[[1]]
character(0)

[[2]]
 [1] "GO:0007442" "GO:0007440" "GO:0016348" "GO:0048617" "GO:0048619"
 [6] "GO:0046872" "GO:0003676" "GO:0000122" "GO:0045944" "GO:0045893"
[11] "GO:0045892" "GO:0048565"

[[3]]
[1] "GO:0006030" "GO:0008061" "GO:0005615" "GO:0008362" "GO:0040014"
[6] "GO:0035151"

[[4]]
 [1] "GO:0007267" "GO:0007446" "GO:0005886" "GO:0007474" "GO:0016348"
 [6] "GO:0042067" "GO:0045198" "GO:0001736" "GO:0017147" "GO:0001737"
[11] "GO:0060071" "GO:0004672" "GO:0006468" "GO:0035159" "GO:0044719"

[[5]]
 [1] "GO:0008017" "GO:0000281" "GO:0051233" "GO:0007052" "GO:0005737"
 [6] "GO:0000916" "GO:0051533" "GO:0005813" "GO:0051298" "GO:0022008"
[11] "GO:0060429" "GO:0007293"

[[6]]
 [1] "GO:0007391" "GO:0007424" "GO:0005886" "GO:0008101" "GO:0009953"
 [6] "GO:0005025" "GO:0050431" "GO:0007476" "GO:0007179" "GO:0030509"
[11] "GO:0007304" "GO:0007448" "GO:0030718" "GO:0042078" "GO:0046845"
[16] "GO:0030707" "GO:0001763" "GO:0006468" "GO:0004672" "GO:0007507"
[21] "GO:0007181" "GO:0045705" "GO:0001745" "GO:0004702" "GO:0005524"
[26] "GO:0045887" "GO:0007274" "GO:0048100" "GO:0035215" "GO:0061327"
[31] "GO:0010629" "GO:0045570" "GO:0005771" "GO:0022407" "GO:0005515"
[36] "GO:0048636" "GO:0090254" "GO:0006357" "GO:0045595" "GO:0016477"
[41] "GO:0007488" "GO:0045927" "GO:0035230" "GO:0005769" "GO:0060799"
[46] "GO:0008582" "GO:0008354" "GO:0030721" "GO:0005829"
```



```
$`1#ptn.1`
[1] "FBbt:00006033" "FBbt:00008110" "FBbt:00008111" "FBbt:00111536"

$`1#ptn.10`
[1] "FBbt:00006029" "FBbt:00008110" "FBbt:00111536" "FBbt:00111537"

$`1#ptn.11`
[1] "FBbt:00006029" "FBbt:00008110" "FBbt:00008111" "FBbt:00111536"

$`1#ptn.12`
[1] "FBbt:00006029" "FBbt:00008110" "FBbt:00111537"

$`1#ptn.13`
[1] "FBbt:00006029" "FBbt:00008111" "FBbt:00111536" "FBbt:00111537"

$`1#ptn.14`
[1] "FBbt:00008111" "FBbt:00111537"

```

Simillarly, 

### Run sem1R

```
library(sem1R)
mysem1R <- new(sem1R)
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
```

The whole source code of this example you can find at [example_disc_auc.R](example_disc_auc.R).

## Authors

* **František Malinka**

## Citation
Cite please:

Malinka, F., Zelezny, F., Klema, J. Concept rule learning with and ontology-based refinement operator inducing semantic rule from omics data. Unpublished.

Kléma, J., Malinka, F. & železný, F. Semantic biclustering for finding local, interpretable and predictive expression patterns. BMC Genomics 18, 752 (2017) doi:10.1186/s12864-017-4132-5

## Acknowledgments

* Jiri Klema for his supervision and considered comments (http://ida.felk.cvut.cz/klema/)

## License

This project is licensed under the MIT License.

