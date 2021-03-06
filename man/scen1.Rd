% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{scen1}
\alias{scen1}
\title{Simulated Data - Scenario 1}
\format{A list with the following elements:
\describe{
\item{pedigrees}{A data.frame representing the pedigrees of all simulated families. It contains the
columuns \code{famid}, \code{id}, \code{dadid}, \code{momid}, and \code{sex}.}
\item{fam.nf}{A integer vector where the entry i indicates the number of individuals at family i.}
\item{phen.nf}{A list of size 100 in which each element represents a data.frame with the values of the
phenotypes X, Y, and Z for all 900 individuals.}
}}
\usage{
scen1
}
\description{
It contains 100 replicates of data simulated according to scenario 1, as described
in \insertCite{ribeiro2019family}{FamilyBasedPGMs}.
}
\details{
In this scenario, the true DAG in the total, genetic, and environmental components is
the unshielded collider \eqn{X \rightarrow Z \leftarrow Y}{X -> Z <- Y} and its moral graph is complete.

It was simulated pedigrees for 30 families, each with 30 individuals (N = 900).

The same pedigrees were used for simulating 100 replicates of family data for three Gaussian phenotypes,
namely X, Y, and Z. No covariates were used in the simulation.
}
\examples{
data(scen1)

}
\references{
{
   \insertAllCited
}
}
\keyword{datasets}
