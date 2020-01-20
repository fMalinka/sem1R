\name{getDatasetExample}
\alias{getDatasetExample}
\title{Get DOT example dataset}
\description{
This is a method that prepare all required inputs for running sem1R
algorithm. The input data are adjusted versions (subset) of DOT
(Dresden Ovary Table) dataset.
}
\usage{
getDatasetExample()
}
\details{
The method simply searches for a given fragment ion in an xcmsFragment
object type given a certain ppm error window
}
\value{
A list frame with the following elements:

\item{colOntologyPath}{Path to DOT column ontology in OBO format}
\item{rowOntologyPath}{Path to Gene ontology row ontology in OBO
format (simplified version)}
\item{datamatrix}{Input binary DOT matrix (subset)}
\item{colOntologyDesc}{List of column anotation}
\item{rowOntologyDesc}{List of row anotation}
}
\seealso{
Dataset and column ontology comes from
\url{http://tomancak-srv1.mpi-cbg.de/DOT/main.html}
}
\examples{
myExample <- getDatasetExample()
}
\author{Frantisek Malinka}
