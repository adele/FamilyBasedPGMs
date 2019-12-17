#' Simulated Data - Scenario 1
#'
#' It contains 100 replicates of data simulated according to scenario 1, as described
#' in \insertCite{ribeiro2019family}{FamilyBasedPGMs}.
#'
#' In this scenario, the true DAG in the total, genetic, and environmental components is
#' the unshielded collider \eqn{X \rightarrow Z \leftarrow Y}{X -> Z <- Y} and its moral graph is complete.
#'
#' It was simulated pedigrees for 30 families, each with 30 individuals (N = 900).
#'
#' The same pedigrees were used for simulating 100 replicates of family data for three Gaussian phenotypes,
#' namely X, Y, and Z. No covariates were used in the simulation.
#'
#' @format A list with the following elements:
#' \describe{
#' \item{pedigrees}{A data.frame representing the pedigrees of all simulated families. It contains the
#' columuns \code{famid}, \code{id}, \code{dadid}, \code{momid}, and \code{sex}.}
#' \item{fam.nf}{A integer vector where the entry i indicates the number of individuals at family i.}
#' \item{phen.nf}{A list of size 100 in which each element represents a data.frame with the values of the
#' phenotypes X, Y, and Z for all 900 individuals.}
#' }
#'
#' @examples
#' data(scen1)
#'
#' @references{
#'    \insertAllCited
#' }
"scen1"

#' Simulated Data - Scenario 2
#'
#' It contains 100 replicates of data simulated according to scenario 2, as described
#' in \insertCite{ribeiro2019family}{FamilyBasedPGMs}.
#'
#' In this scenario, the true DAG in the total, genetic, and environmental components is
#' the fork \eqn{Y \leftarrow X \rightarrow Z}{Y <- X -> Z} and its moral graph is exactly its undirected version.
#'
#' It was simulated pedigrees for 30 families, each with 30 individuals (N = 900).
#'
#' The same pedigrees were used for simulating 100 replicates of family data for three Gaussian phenotypes,
#' namely X, Y, and Z. No covariates were used in the simulation.
#'
#' @format A list with the following elements:
#' \describe{
#' \item{pedigrees}{A data.frame representing the pedigrees of all simulated families. It contains the
#' columuns \code{famid}, \code{id}, \code{dadid}, \code{momid}, and \code{sex}.}
#' \item{fam.nf}{A integer vector where the entry i indicates the number of individuals at family i.}
#' \item{phen.nf}{A list of size 100 in which each element represents a data.frame with the values of the
#' phenotypes X, Y, and Z for all 900 individuals.}
#' }
#'
#' @examples
#' data(scen2)
#'
#' @references{
#'    \insertAllCited
#' }
"scen2"

#' Simulated Data - Scenario 3
#'
#' It contains 100 replicates of data simulated according to scenario 3, as described
#' in \insertCite{ribeiro2019family}{FamilyBasedPGMs}.
#'
#' In this scenario, the true DAG in the total, genetic, and environmental components is
#' \eqn{Z \leftarrow X \rightarrow Y \rightarrow Z}{Z <- X -> Y -> Z}  and its moral graph is complete.
#' However, data for the total DAG may be unfaithful, since \eqn{\rho_{X,Y|Z}}{rho_XY|Z} is almost zero.
#'
#' It was simulated pedigrees for 30 families, each with 30 individuals (N = 900).
#'
#' The same pedigrees were used for simulating 100 replicates of family data for three Gaussian phenotypes,
#' namely X, Y, and Z. No covariates were used in the simulation.
#'
#' @format A list with the following elements:
#' \describe{
#' \item{pedigrees}{A data.frame representing the pedigrees of all simulated families. It contains the
#' columuns \code{famid}, \code{id}, \code{dadid}, \code{momid}, and \code{sex}.}
#' \item{fam.nf}{A integer vector where the entry i indicates the number of individuals at family i.}
#' \item{phen.nf}{A list of size 100 in which each element represents a data.frame with the values of the
#' phenotypes X, Y, and Z for all 900 individuals.}
#' }
#'
#' @examples
#' data(scen3)
#'
#' @references{
#'    \insertAllCited
#' }
"scen3"

#' Simulated Data - Scenario 4
#'
#' It contains 100 replicates of data simulated according to scenario 4, as described
#' in \insertCite{ribeiro2019family}{FamilyBasedPGMs}.
#'
#' In this scenario, the DAG of the total and environmental components is the unshielded collider
#' \eqn{X \rightarrow Z \leftarrow Y}{X -> Z <- Y}. However, the edge \eqn{X \rightarrow Z}{X -> Z} is absent in the genetic PGM,
#' i.e., the genetic PGM has only the edge \eqn{Z \leftarrow Y}{Z <- Y}.
#'
#' It was simulated pedigrees for 30 families, each with 30 individuals (N = 900).
#'
#' The same pedigrees were used for simulating 100 replicates of family data for three Gaussian phenotypes,
#' namely X, Y, and Z. No covariates were used in the simulation.
#'
#' @format A list with the following elements:
#' \describe{
#' \item{pedigrees}{A data.frame representing the pedigrees of all simulated families. It contains the
#' columuns \code{famid}, \code{id}, \code{dadid}, \code{momid}, and \code{sex}.}
#' \item{fam.nf}{A integer vector where the entry i indicates the number of individuals at family i.}
#' \item{phen.nf}{A list of size 100 in which each element represents a data.frame with the values of the
#' phenotypes X, Y, and Z for all 900 individuals.}
#' }
#'
#' @examples
#' data(scen4)
#'
#' @references{
#'    \insertAllCited
#' }
"scen4"

#' Simulated Data - Scenario 5
#'
#' It contains 100 replicates of data simulated according to scenario 5, as described
#' in \insertCite{ribeiro2019family}{FamilyBasedPGMs}.
#'
#' In this scenario, the true DAG in the total, genetic, and environmental components is
#' \eqn{X \rightarrow Z \leftarrow U \rightarrow W \leftarrow Y}{X -> Z <- U -> W <- Y}.
#' Since U is hidden, the true MAG is \eqn{X \rightarrow Z \leftrightarrow W \leftarrow Y}{X -> Z <-> W <- Y}.
#'
#' It was simulated pedigrees for 30 families, each with 30 individuals (N = 900).
#'
#' The same pedigrees were used for simulating 100 replicates of family data for five Gaussian phenotypes,
#' namely X, Y, Z, W, and U. No covariates were used in the simulation.
#'
#' @format A list with the following elements:
#' \describe{
#' \item{pedigrees}{A data.frame representing the pedigrees of all simulated families. It contains the
#' columuns \code{famid}, \code{id}, \code{dadid}, \code{momid}, and \code{sex}.}
#' \item{fam.nf}{A integer vector where the entry i indicates the number of individuals at family i.}
#' \item{phen.nf}{A list of size 100 in which each element represents a data.frame with the values of the
#' phenotypes X, Y, and Z for all 900 individuals.}
#' }
#'
#' @examples
#' data(scen5)
#'
#' @references{
#'    \insertAllCited
#' }
"scen5"

