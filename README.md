# sem1R - Concept single rule learning with an ontology-based refinement operator

sem1R is a machine learning algorithm that finds interesting, hidden, and non-trivial patterns in two-dimensional kind of data. The algorithm generates 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites
We required to use R in version 3.0.2.
All prerequisites R packages that are needed for the sem1R package are the following:
`Rcpp`, `RcppProgress`, `RcppArmadillo`, and `BH`. All of these packages come from CRAN, so install them by `install.packages` function in R.

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
The package in tar.gz format will be named as `sem1R_version.tar.gz'.
Then, install the sem1R package.

```
R CMD INSTALL sem1R_version.tar.gz
```

And finally, check whether the sem1R package was installed.
Run R and load the package.
```
library(sem1R)
```


## Running the tests

Explain how to run the automated tests for this system

```
Give an example
```

## Authors

* **František Malinka** - *Initial work*

## License

This project is licensed under the MIT License.

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc


