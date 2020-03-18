# sem1R - Finding semantic patterns in omics data using concept rule learning with an ontology-based refinement operator

<img align="right" width="200" src="sem1r_logo_col.png">
<p style="text-align: justify">sem1R is a machine learning algorithm that finds interesting, hidden, and non-trivial semantic patterns in omics data. The algorithm produces a set of prediction rules that form data into clusters or biclusters, this depends on a type of used ontologies (column, row, or both). Here, we make distinctions between two types of ontologies: an ontology describing rows (e.g. genes) and columns (e.g. samples). Practically, for gene expression data, where rows represent genes and column represent samples, we recommend to use Gene ontology or any pathway ontologies as a row ontology. Choosing a proper column ontology is depending on a type of experiment, e.g. OBO Foundry provides almost two hundreds ontologies and many of them are domain specific so some anatomical ontologies can be used as well. An example of gene expression dataset that adresses simultaneously column and row ontologies is DOT (Dresden Ovary Table) at http://tomancak-srv1.mpi-cbg.de/DOT/main.html.

The sem1R is based on rule learning methods, where two reduction procedures make the algorithm extremely fast and efficient in comparison with the traditional CN2 approach. In additional, it is easy to use, because all important methods are included into this R package.</p>

## Getting Started
The algorithm is implemented in C++ and provided as `R` package. The following instructions will show you how to install all prerequisites and the sem1R package as well into your local machine. Afterwards, we will demonstrate the sem1R on real gene expression dataset.

### Prerequisites
We required to use R in version 3.4.
All prerequisites R packages that are needed for the sem1R package are the following:
`Rcpp (>= 0.12.6)`, `RcppProgress`, `RcppArmadillo (>= 0.7.800.2.0)`, and `BH (>= 1.72.0-3)`. All of these packages come from CRAN, so install them by `install.packages` function in R.

Or, for easier installation we recommend to install 'devtools' that can download and install the project instantly from gitHub using only one command.

For installing 'devtools' package run `R` and type the following:
```
install.packages("devtools")
```

### Installing
If you choosen instalation via devtools, you would go to the terminal, run `R` and then the following commands:

```
>library(devtools)
>install_github("fmalinka/sem1R")
```
All prerequisites packages should be installed automatically.

Or for non-devtools users, download the sem1R package and extract it to your arbitrary folder.
Run the terminal, go to the folder and install all prerequisities using `install.packages` R's function. Finally, build the package.

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
>library(sem1R)
```


## Running the example
Running example that we present here comes from Dresden Ovary Table (DOT) located at http://tomancak-srv1.mpi-cbg.de/DOT/main.html. Since the original data matrix is to complex for a brief algorithm exhibition, we will work just with a submatrix of the original matrix in this tutorial. All necessary files for this example are stored at [inst](inst) folder.

### Data matrix
A file [dotmatrix.csv](inst/extdata/dotmatrix.csv) contains binary information about gene expression of the fruit fly adult ovary in many locations. The matrix is two-dimensional where rows represent genes and columns represent samples (locations). Each dimension has own identifier, i.e. genes are described by FBgn (FlyBase) identifiers and columns by your notation. Ones in the matrix mean "expressed" and zeros mean "non-expressed" in the given positions. Obviously, process of binarization has to be done if your data are not in the binary format. Look below how the matrix looks like.


```
            X1.8.somatic.cells X1.9.germline.cells X1.10.terminal.filament
FBgn0033019                  1                   1                       1
FBgn0263251                  1                   1                       1
FBgn0037224                  1                   1                       1
FBgn0038013                  1                   1                       1
FBgn0037358                  1                   1                       1
```

### Ontologies
Ontologies are the second type of input that has to be given to sem1R algorithm. Ontology has to be in OBO format (https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html) and relationships of terms must be acyclic (usually OBO ontologies are acyclic). For many other interesting ontologies look at OBO Foundry (http://obofoundry.org/). In our running example, we provide two type of different ontologies. Gene ontology, located at [inst/extdata/go-basic-reduced.obo](inst/extdata/go-basic-reduced.obo), aims to rows of the data matrix and DOT ontology (http://tomancak-srv1.mpi-cbg.de/cgi-bin-public/ovary_annotation_hierarchy.pl), located at [inst/extdata/dotOntology.obo](inst/extdata/dotOntology.obo), focuses on the columns.

### Connection between the data matrix and the ontologies
Now, the last step is to established an annotation, a connection between our data matrix and all given ontologies. Firstly, we look at the rows which are described by the FBGN identifiers. Result of mapping from FBGN identifiers to Gene ontology terms id is provided at [inst/extdata/initRowDot_reduced.csv](inst/extdata/initRowDot_reduced.csv) file. File showing a mapping from data matrix colums to DOT ontology can be found at [inst/extdata/initColDot.csv](inst/extdata/initColDot.csv) file.


### Run sem1R

Finally, let's run the example!

First of all, load the R library and create a new class `sem1R`. Then, we load the example data containing all necessary files described above.
```
> library(sem1R)
> mysem1R <- new(sem1R)
> myExample <- getDatasetExample()

```
Now, we load the data matrix to the sem1R class. Be sure, that the data matrix is a 'matrix' R type and has named rows and columns. It is important! Note that public methods of the class are call by $ symbol.
```
> mysem1R$setDataset(myExample$datamatrix)
```
Then, we load all ontologies. For this, use `createCOLOntology` or `createROWOntology` methods, it depends on your matrix design generally. The first argument of these methods is name of ontolgy, the second argument set up path to the corresponding obo file, and the last one is a list of vectors representing the connection between rows/columns and ontologies. For the proper format look at one of the examples (myExample$colOntologyDesc or myExample$rowOntologyDesc). When you have more than one ontology, just call the corresponding method one again. However, the name of ontology mush be unique!
```
> mysem1R$createCOLOntology("DOT", myExample$colOntologyPath, myExample$colOntologyDesc)
> mysem1R$createROWOntology("GO", myExample$rowOntologyPath, myExample$rowOntologyDesc)
```
Now, we set all algorithm parameters (see R manual).
```
> mysem1R$filterTh <- 50
> mysem1R$objective <- "auc"
> mysem1R$ruleDepth <- 3
> mysem1R$nrules <- 2
> mysem1R$featureSelectionMethod <- 0
> mysem1R$minLevel <- 2
```
If you want to check out correcness of the connection of data matrix and the ontologies, call 'mysem1R$checkRowDescription()' or 'mysem1R$checkColDescription()'.

Finally, run the algorithm and save the results!
```
> myhypothesis <- mysem1R$findDescription()
```

When it ends ...
```
[sem1R SETTINGS]
filter threshold: 50
rule depth: 3
significance threshold: 6.635
objective function: auc
number of rules: 2
featureSelectionMethod: 0
ruleFormat: both
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
```
... your final rule set will be printed on STDOUT.
```
******************************************************************
************************** FINAL RULESET *************************
===== RULE 1=====
 STATS: score 0.536225 t-score: 2182.28 POSITIVE: 23351 NEGATIVE: 11649
 RULE: GO:0044763 AND GO:0043229
 DETAILS: 
ID: GO:0044763
NAME: single-organism cellular process
DEF: "Any process that is carried out at the cellular level, occurring within a single organism." [GOC:jl]
level: 2

ID: GO:0043229
NAME: intracellular organelle
DEF: "Organized structure of distinctive morphology and function, occurring within the cell. Includes the nucleus, mitochondria, plastids, vacuoles, vesicles, ribosomes and the cytoskeleton. Excludes the plasma membrane." [GOC:go_curators]
level: 2

COVERED:
  POSITIVE:
X1.8.somatic.cells, X1.9.germline.cells, X1.11.cap.cells, X1.13.follicle.stem.cells, X1.15.interfollicular.stalk.cells, X1.18.germline.stem.cells, X1.20.presumptive.nurse.cells, X2.26.oocyte, X4.29.oocyte, X5.30.oocyte, X3.32.nurse.cells, X5.34.nurse.cells, X2.35.somatic.cells, X4.37.somatic.cells, X2.39.follicle.cells, X3.40.follicle.cells, X5.42.follicle.cells, X3.44.interfollicular.stalk.cells, X5.46.interfollicular.stalk.cells, X3.48.anterior.follicle.cells, X4.49.border.cells, X2.58.posterior.follicle.cells, X4.60.posterior.follicle.cells, X5.62.centripetally.migrating.follicle.cells, X2.64.anterior.restriction, X2.66.nurse.cells_nuclear.foci, X2.69.cytoplasmic.foci, X3.71.anterior.restriction, X5.74.anterior.restriction, X3.75.posterior.restriction, X5.77.posterior.restriction, X5.81.cortical.enrichment, X3.82.nurse.cells_nuclear.foci, X5.84.nurse.cells_nuclear.foci, X4.86.nurse.cells_perinuclear, X3.105.basal.restrictrion, X4.106.basal.restrictrion, X3.108.apical.restriction, X5.110.apical.restriction, X4.112.cytoplasmic.foci, X5.113.cytoplasmic.foci, X3.115.cytoplasmic.foci, X5.117.cytoplasmic.foci, X3.119.cytoplasmic.foci, X4.120.cytoplasmic.foci, X2.128.nuclear.foci, X4.130.nuclear.foci, X4.139.anterior.follicle.cell, X5.140.squamous.follicle.cells, X4.143.cortical.enrichment, X2.145.cortical.enrichment, X4.147.cortical.enrichment, X5.149.follicle.cells.overlaying.the.oocyte, X2.164.perinuclear, X4.166.perinuclear, X3.168.oocyte.nucleus, X5.170.oocyte.nucleus, 
  NEGATIVE:
X1.8.somatic.cells, X1.10.terminal.filament, X1.12.escort.cells, X1.14.follicle.cells, X1.17.posterior.follicle.cells, X1.19.cystoblast, X1.21.presumptive.oocyte, X4.29.oocyte, X2.31.nurse.cells, X5.34.nurse.cells, X3.36.somatic.cells, X5.38.somatic.cells, X4.41.follicle.cells, X2.43.interfollicular.stalk.cells, X4.45.interfollicular.stalk.cells, X3.48.anterior.follicle.cells, X5.50.border.cells, X3.59.posterior.follicle.cells, X5.62.centripetally.migrating.follicle.cells, X2.65.posterior.restriction, X2.67.nurse.cells_perinuclear, X2.70.apical.restriction, X5.74.anterior.restriction, X4.76.posterior.restriction, X4.80.cortical.enrichment, X4.83.nurse.cells_nuclear.foci, X3.85.nurse.cells_perinuclear, X3.105.basal.restrictrion, X5.107.basal.restrictrion, X4.109.apical.restriction, X4.112.cytoplasmic.foci, X2.114.cytoplasmic.foci, X4.116.cytoplasmic.foci, X3.119.cytoplasmic.foci, X5.121.cytoplasmic.foci, X4.130.nuclear.foci, X4.139.anterior.follicle.cell, X3.142.cortical.enrichment, X5.144.cortical.enrichment, X4.147.cortical.enrichment, X5.149.follicle.cells.overlaying.the.oocyte, X3.165.perinuclear, X3.168.oocyte.nucleus, X5.170.oocyte.nucleus, 
===== =====

```

And the structure of returned hypothesis is the following:
```
> str(myhypothesis[[1]])
List of 9
 $ ruleID     : int 1
 $ score      : num 0.536
 $ tscore     : num 2182
 $ nCoveredPOS: int 23351
 $ nCoveredNEG: int 11649
 $ rules      : chr [1:2] "GO:0044763" "GO:0043229"
 $ details    : chr [1:6] "ID: GO:0044763" "NAME: single-organism cellular process" "DEF: \"Any process that is carried out at the cellular level, occurring within a single organism.\" [GOC:jl]" "ID: GO:0043229" ...
 $ coveredPOS : chr [1:23351] "FBgn0039115,X1.8.somatic.cells" "FBgn0022238,X1.8.somatic.cells" "FBgn0262601,X1.8.somatic.cells" "FBgn0029134,X1.8.somatic.cells" ...
 $ coveredNEG : chr [1:11649] "FBgn0026737,X1.8.somatic.cells" "FBgn0003087,X1.8.somatic.cells" "FBgn0031873,X1.8.somatic.cells" "FBgn0003514,X1.8.somatic.cells" ...
```
where ruleID represents order of the induced rule, score represents quality of the rule depends on the type of evaluation function, tscore represents chi square score of the rule, positive and negative covered is a number of examples covered by the rule, rules represents a conjunction of ontological terms, details provides additional information about the terms in conjunction, and finally covered represents covered examples expressed by their position in the matrix.


## Authors

* **František Malinka**

## Citation
Cite please:

Malinka, F., Zelezny, F., Klema, J. Finding semantic patterns in omics data using concept rule learning with an ontology-based refinement operator. Unpublished.

Kléma, J., Malinka, F. & železný, F. Semantic biclustering for finding local, interpretable and predictive expression patterns. BMC Genomics 18, 752 (2017) doi:10.1186/s12864-017-4132-5

## Acknowledgments

* Jiri Klema for his supervision and considered comments (http://ida.felk.cvut.cz/klema/).

## License

This project is licensed under the MIT License.

