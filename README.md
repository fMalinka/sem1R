# sem1R - Finding semantic patterns in omics data using concept rule learning with an ontology-based refinement operator

<img align="right" width="200" src="sem1r_logo_col.png">
<p style="text-align: justify">sem1R is a machine learning algorithm that finds interesting, hidden, and non-trivial semantic patterns in omics data. The algorithm produces a set of prediction rules that form data into clusters or biclusters, this depends on a type of ontologies (column, row, or both). Here, we distingues between two types of ontologies: an ontology describing rows (e.g. genes) an columns (e.g. samples). Practically, for gene expression data, where rows represent genes and column represent samples, we recommend to use Gene ontology or any pathway ontologies as a row ontology. Choosing a proper column ontology is depending on a type of experiment, e.g. OBO Foundry provides almost two hundreds ontologies and many of them are domain specific so some anatomical ontologies can be used as well. An example of gene expression dataset that adresses simultaneously column and row ontologies is DOT (Dresden Ovary Table) at http://tomancak-srv1.mpi-cbg.de/DOT/main.html.

The sem1R is based on rule learning methods, where two reduction procedures make the algorithm extremely fast and efficient in comparison with the traditional CN2 approach. In additional, it is easy to use, because all important methods are included into this R package.</p>

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
If you choosen instalation via devtools, you would go to the terminal, run `R` and then the following commands:

```
>library(devtools)
>install_gitHub("fmalinka/sem1R")
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
Running example that we present here comes from [1] and shows a gene expression over Drosophila melanogaster imaginal discs. All necessary files for this example are located at [example](example) folder.

### Data matrix
A file [discMatrix.csv](example/discMatrix.csv) contains binary information about gene expression over imaginal discs of Drosophila melanogaster. The matrix is two-dimensional where rows represent genes and columns represent samples (locations). Each dimension has own identifier, i.e. genes are described by FBgn (FlyBase) identifiers and columns by your notation. Ones in the matrix mean "expressed" and zeros mean "non-expressed" in the given positions. Obviously, process of binarization has to be done if your data are not in the binary format.

### Ontologies
Ontologies are the second type of input that has to be given to sem1R algorithm. Ontology has to be in OBO format (https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html) and relationships of terms must be acyclic (usually OBO ontologies are acyclic). For many other interesting ontologies look at OBO Foundry (http://obofoundry.org/). In our running example, we provide two type of different ontologies. Gene ontology, located at [example/go-basic.obo](example/go-basic.obo), aims to rows of the data matrix and FBBT ontology (http://obofoundry.org/ontology/fbbt.html), located at [example/fbbt-simple.obo](example/fbbt-simple.obo), focuses on the columns.

### Connection between the data matrix and the ontologies
Now, the last step is to established an annotation, a connection between our data matrix and all given ontologies. Firstly, we look at the rows which are described by the FBGN identifiers. The mapping from FBGN identifiers to Gene ontology is provided at [example/gene_association.fb](example/gene_association.fb) file. So, for each row we have to find the corresponding GO terms. An example of subset of data matrix is shown bellow.

```
>head(mydata)
            1#ptn.1 1#ptn.10 1#ptn.11 1#ptn.12 1#ptn.13 1#ptn.14
FBgn0037992       0        0        0        1        1        0
FBgn0024244       1        1        1        1        1        0
FBgn0026077       1        1        1        1        1        1
FBgn0000658       0        1        1        1        1        0
FBgn0030241       1        1        1        1        1        1
FBgn0003716       1        1        1        1        1        1
```

The connection is represented by List data type, where each element contains a vector of ontological term ids (in jargon of OBO format, it is a term id) that are associated with the corresponding row. List can be named or unnamed, always depends on the order of elements. Bellow, we show an example of the connection on the columns level.

```
>head(ontoDesc$col)
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
.
.
.

```
And, an example of the connection for the rows level. Note that the first row, described as FBgn0037992, is not associated with any gene terms.
```
>head(ontoDesc$row)
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

For source code that generates this kind of connection for the specific dataset [discMatrix.csv](example/discMatrix.csv), see `PrepareOntologyDesc` function in file [example_disc_auc.R](example_disc_auc.R).

### Run sem1R

Finally, let's run the example!

Firstly, load the library and create a new class `sem1R`. Then, we load the data matrix and convert it to matrix data type. Be sure, that the matrix has named rows and columns. It is important!
```
>library(sem1R)
>mysem1R <- new(sem1R)
>mydata <- read.csv("example/discMatrix.csv", header = TRUE, check.names = FALSE, row.names = 1)
>mydata <- as.matrix(mydata)
```
Now, we make a connection between rows/columns and ontologies.
```
>ontoDesc <- PrepareOntologyDesc(mydata, geneASOC = "example/gene_association.fb", colCSV = "example/initsegmentFBbtWithoutComments.csv")
```
When our data are in the required format, we put them to our class. As you can see from the example above, all public methods of the package are accesible throw dollar sign. To load the dataset, use `setDataset` method. Then, load all ontologies. For this, use `createCOLOntology` or `createROWOntology` methods, it depends on your matrix design. The first argument of these methods is name, the second argument set up path to obo files, and the last one is a list of vectors representing the connection between rows/columns and ontologies. When you have more than one ontology, just call the corresponding method one again. However, the name of ontology mush be unique!
```
>mysem1R$setDataset(mydata)
>mysem1R$createCOLOntology("FBGN", "example/fbbt-simple.obo", ontoDesc$col)
>mysem1R$createROWOntology("GO", "example/go-basic.obo", ontoDesc$row)
```

Now, we set same parameters (see manual page).
```
>mysem1R$ruleFormat <- "both"
>mysem1R$filterTh <- 100
>mysem1R$objective <- "auc"
>mysem1R$ruleDepth <- 9
>mysem1R$nrules <- 10
>mysem1R$featureSelectionMethod <- 0
>mysem1R$minLevel <- 3
```

Finally, run the algorithm and save the results!
```
>myhypothesis <- mysem1R$findDescription()
```

When it ends ...
```
[sem1R SETTINGS]
filter threshold: 100
rule depth: 5
significance threshold: 6.635
objective function: auc
number of rules: 10
exhaustive test: 0
featureSelectionMethod: 0
ruleFormat: both
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
```
... your final rule set will be printed on STDOUT and results will be also returned as a List where each rule contains basic information like .
```
******************************************************************
************************** FINAL RULESET *************************
===== RULE 1=====
 STATS: score 0.62601 t-score: 4039.26 POSITIVE: 31368 NEGATIVE: 4842
 RULE: FBbt:00000169 AND FBbt:00007006
 DETAILS: 
ID: FBbt:00000169
NAME: embryonic mesothoracic segment
DEF: "Mesothoracic segment of the embryo." [FBC:DOS]
level: 7

ID: FBbt:00007006
NAME: developing material anatomical entity
DEF: 
level: 2

COVERED:
  POSITIVE:
1#ptn.1, 1#ptn.10, 1#ptn.11, 1#ptn.12, 1#ptn.13, 1#ptn.14, 1#ptn.15, 1#ptn.16, 1#ptn.18, 1#ptn.19, 1#ptn.3, 1#ptn.4, 1#ptn.5, 1#ptn.6, 1#ptn.7, 1#ptn.8, 1#ptn.9, 2#ptn.10, 2#ptn.11, 2#ptn.12, 2#ptn.13, 2#ptn.14, 2#ptn.3, 2#ptn.5, 2#ptn.6, 2#ptn.7, 2#ptn.8, 2#ptn.9, 
  NEGATIVE:
1#ptn.1, 1#ptn.17, 1#ptn.2, 2#ptn.10, 2#ptn.3, 
===== =====
```

And the structure of returned hypothesis is the following:
```
> str(myhypothesis[[1]])
List of 8
 $ ruleID         : int 1
 $ score          : num 0.626
 $ tscore         : num 4039
 $ positiveCovered: int 31368
 $ negativeCovered: int 4842
 $ rules          : chr [1:2] "FBbt:00000169" "FBbt:00007006"
 $ details        : chr [1:6] "ID: FBbt:00000169" "NAME: embryonic mesothoracic segment" "DEF: \"Mesothoracic segment of the embryo.\" [FBC:DOS]" "ID: FBbt:00007006" ...
 $ covered        : chr [1:31368] "FBgn0024244,1#ptn.1" "FBgn0026077,1#ptn.1" "FBgn0030241,1#ptn.1" "FBgn0003716,1#ptn.1" ...
```
where ruleID represents order of the induced rule, score represents quality of the rule depends on the type of evaluation function, tscore represents chi square score of the rule, positive and negative covered is a number of examples covered by the rule, rules represents a conjunction of ontological terms, details provides additional information about the terms in conjunction, and finally covered represents covered examples expressed by their position in the matrix.

The whole source code of this example you can find at [example_disc_auc.R](example_disc_auc.R).

## Authors

* **František Malinka**

## Citation
Cite please:

Malinka, F., Zelezny, F., Klema, J. Concept rule learning with and ontology-based refinement operator inducing semantic rule from omics data. Unpublished.

Kléma, J., Malinka, F. & železný, F. Semantic biclustering for finding local, interpretable and predictive expression patterns. BMC Genomics 18, 752 (2017) doi:10.1186/s12864-017-4132-5

## Acknowledgments

* Jiri Klema for his supervision and considered comments (http://ida.felk.cvut.cz/klema/).

## License

This project is licensed under the MIT License.

