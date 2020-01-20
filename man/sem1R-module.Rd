\name{sem1R}
\alias{sem1R}
\docType{class}
\title{
Main class representing sem1R experiment.
}
\description{
This main class stores all required settings and the final induced
hypothesis.
}
\section{Objects from the Class}{
All fields and methods are accesible via dollar extractor.
}
\section{fields}{
\describe{
\item{featureSelectionMethod}{\code{numeric(1)}:
selects a feature selection method (0: FS_atLeastOne, 1: FS_onlySig,
2: FS_sigAtLeastOne)}
\item{filterTh}{\code{numeric(1)}: width of beam in Beam search}
\item{minLevel}{\code{numeric(1)}: minimum level of specificity valid
for each ontological term}
\item{nrules}{\code{numeric(1)}: maximum number of induced rules}
\item{objective}{\code{character(1)}: selects an evalution function
(acc: accuracy, f1: F1-score, auc: Area Under Curve)}
\item{ruleDepth}{\code{numeric(1)}: maximal length of rule}
\item{ruleFormat}{\code{character(1)}: description format for covered
examples (col: only column description, row: only row description,
both: column and row description)}
\item{signTH}{\code{numeric(1)}: minimum significance value of a rule}
\item{verbose}{\code{numeric(1)}: 0 is non-verbose, otherwise debug
information will be printed on STDOUT}
}
}
\section{methods}{
\describe{
\item{computeTermsEnrichment}{\code{}: compute enrichment score for
each ontological terms in provided ontologies}
\item{createCOLOntology}{\code{character, character, list}:
name of ontology, path to ontology in OBO format, list of descriptions.
Append column ontology to the Class.}
\item{createROWOntology}{\code{character, character, list}:
name of ontology, path to ontology in OBO format, list of descriptions.
Append row ontology to the Class.}
\item{findDescription}{\code{}:
runs the algorithm and finds a hypothesis}
\item{setDataset}{\code{matrix}:
load a two-dimensional binary data matrix}
}
}
\details{
The third argument of \env{createCOLOntology} and \env{createROWOntology}
method must be R list where each element is a vector of characters
associating rows/columns with the corresponding ontological terms.
}
\references{
Malinka, F., Zelezny, F., Klema, J. Concept rule learning
with and ontology-based refinement operator inducing semantics rule
from omics data. Unpublished.
}
\examples{
#new class needs to be established
myExp <- new(sem1R)
#for an example visit http://www.github.com/fmalinka/sem1r or see a vignett
}
\keyword{sem1R-class}
