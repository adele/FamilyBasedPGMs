#' @title Learning of Gaussian Undirected PGMs from Family Data
#'
#' @description A Gaussian undirected PGM and its decomposition into genetic and environmental components
#' are learned from observational family data by assigning an edge between every pair of variables such that
#' the partial correlation between the two variables in question (in the respective component) given all the other
#' variables is significantly different from zero.
#'
#' The zero partial correlation tests are derived in the work by
#' \insertCite{ribeiro2019family;textual}{FamilyBasedPGMs}. These tests are based on univariate polygenic linear mixed
#' models \insertCite{almasy1998multipoint}{FamilyBasedPGMs}, with two components of variance: the polygenic or
#' family-specific random effect, which models the phenotypic variability across the families,
#' and the environmental or subject-specific error, which models phenotypic variability after removing the
#' familial aggregation effect.
#'
#' @param phen.df A data.frame with phenotype variables of only sampled subjects.
#' Column names must be properly set with the names of the phenotypes.
#' @param covs.df A data.frame with covariates of only sampled subjects.
#' Column names must be properly set with the names of the covariates.
#' @param pedigrees A data.frame with columuns \code{famid}, \code{id}, \code{dadid},
#' \code{momid}, and \code{sex} columns for all sampled and non-sampled subjects.
#' @param sampled A logical vector in which element i indicates whether individual i was sampled or not.
#' @param fileID  A character string to be used as prefix in the filenames of
#' R objects with the partial correlation results. Note that covariates are not
#' identified in these files.
#' @param dirToSave Path to the folder you want to save the output objects.
#' @param alpha The significance level to be used in the partial correlation tests.
#' @param correction A character string indicating the correction method
#' to be used in the \code{p.adjust} function. The options are:
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none".
#' @param max_cores An integer value indicating the maximum number of CPU cores
#' to be used for parallel execution.
#' @param minK A scalar indicating the minimum dimension allowed
#' in the dimensionality reduction for confounding correction.
#' @param maxFC A scalar between 0 and 1 indicating the maximum fraction of
#' confounding allowed.
#' @param orthogonal A logical value indicating whether the transformation matrix used in the
#' confounding correction is orthogonal or not.
#' @param useGPU A logical value indicating whether GPU cores can be used for parallel execution.
#' @param debug A logical value indicating whether some debug messages can be shown.
#' @param logFile Optional file path and name to save progress and error messages.
#' If not provided and debug is True a default file is created in the dirToSave folder.
#'
#' @return Returns a list with the following elements:
#' \describe{
#' \item{pcor}{A list with the total (\code{pcor_t}), genetic (\code{pcor_g}),
#' and environmental (\code{pcor_e}) partial correlation matrices.}
#' \item{adjM}{A list with the total (\code{t}), genetic (\code{g}),
#' and environmental (\code{e}) partial adjacency matrices.}
#' \item{udg}{A list with the total (\code{t}), genetic (\code{g}),
#' and environmental (\code{e}) igraph objects representing the respective undirected graphs.}
#'}
#'
#' @examples
#' \donttest{
#' data(scen3) # available simulated datasets are scen1, scen2, scen3, and scen4
#'
#' scenario = 3 # data was simulated according to scenario 3
#'
#' fam.nf <- scen3$fam.nf
#' pedigrees <- scen3$pedigrees
#' phen.df <- scen3$phen.df[[1]] # accessing the first replicate
#' covs.df <- NULL # no covariates were used in the simulation process.
#'
#' N <- sum(fam.nf) # total number of individuals
#' sampled <- rep(1, N) # in simulated data, all individuals were sampled.
#'
#' fileID <- paste0("scen", scenario)
#' dirToSave <- paste0("./objects-UDG-", fileID, "/")
#' dir.create(dirToSave, showWarnings=FALSE)
#' alpha = 0.05
#'
#' udgs.out <- learnFamilyBasedUDGs(phen.df, covs.df, pedigrees, sampled,
#'                                  fileID, dirToSave, alpha, correction=NULL,
#'                                  max_cores=NULL, minK=10, maxFC = 0.05,
#'                                  orthogonal=TRUE, useGPU=FALSE, debug=TRUE)
#'
#' # the adjacency matrix of the learned undirected genetic PGM
#' udgs.out$adjM$g
#'
#' # the estimates, p-values, and effective sizes of the genetic partial correlations
#' udgs.out$pCor$pCor_g
#'
#' # plotting the the learned undirected genetic PGM as an `igraph` object:
#' plot(udgs.out$udg$g, vertex.size=30, vertex.color="lightblue")
#'
#' #' # the adjacency matrix of the learned undirected environmental PGM
#' udgs.out$adjM$e
#'
#' # the estimates, p-values, and effective sizes of the environmental partial correlations
#' udgs.out$pCor$pCor_e
#'
#' # plotting the the learned undirected environmental PGM as an `igraph` object:
#' plot(udgs.out$udg$e, vertex.size=30, vertex.color="lightblue")
#' }
#'
#' @references{
#'    \insertAllCited
#' }
#'
#' @importFrom Rdpack reprompt
#' @import doMC
#' @import foreach
#' @importFrom utils combn
#' @importFrom parallel detectCores
#'
#' @export learnFamilyBasedUDGs
learnFamilyBasedUDGs <- function(phen.df, covs.df, pedigrees, sampled, fileID, dirToSave,
                                 alpha=0.05, correction=NULL, max_cores=NULL,
                                 minK=10, maxFC = 0.01, orthogonal=TRUE, useGPU=FALSE,
                                 debug=TRUE, logFile=NULL)
{
  if (debug && is.null(logFile)){
    logFile = paste0(dirToSave, "learnFamilyBasedUDGs_", format(Sys.time(), "%Y%m%d-%Hh%Mm%Ss"), ".log")
    cat("Learning UDGs from Family Data...\n", file=logFile, append=FALSE)
  }

  preprocessFamilyBasedUDGs(phen.df, covs.df, pedigrees, sampled, fileID, dirToSave,
                            max_cores, minK, maxFC, orthogonal, useGPU,
                            debug, logFile)

  udgPCorM <- selectBestUDGPCorMatrices(dim(phen.df)[2], alpha, maxFC, fileID, dirToSave,
                                        savePlots=TRUE, variables=colnames(phen.df),
                                        debug=debug, logFile=logFile)

  adjMatrices <- list()
  adjMatrices$g <- getAdjMatrixFromSymmetricPCor(udgPCorM$pCor_g, alpha, correction)$adjacency.matrix
  adjMatrices$e <- getAdjMatrixFromSymmetricPCor(udgPCorM$pCor_e, alpha, correction)$adjacency.matrix
  adjMatrices$t <- getAdjMatrixFromSymmetricPCor(udgPCorM$pCor_t, alpha, correction)$adjacency.matrix

  udgList <- list()
  if (requireNamespace("igraph", quietly = TRUE)) {
    udgList$g <- igraph::as.undirected(igraph::graph.adjacency(adjMatrices$g))
    udgList$e <- igraph::as.undirected(igraph::graph.adjacency(adjMatrices$e))
    udgList$t <- igraph::as.undirected(igraph::graph.adjacency(adjMatrices$t))
  }
  return(list(pCor=udgPCorM, adjM=adjMatrices, udg=udgList))
}

#' @title Learning of Gaussian Directed Acyclic PGMs from Family Data
#'
#' @description The CPDAG representing the Markov equivalence classes to which the Gaussian directed acyclic
#' PGM belong and its decomposition into genetic and environmental components are learned
#' from observational family data by applying the IC/PC algorithm
#' \insertCite{pearl2000causality,spirtes2000causation}{FamilyBasedPGMs}
#' with the zero partial correlation tests derived in the work by
#' \insertCite{ribeiro2019family;textual}{FamilyBasedPGMs} as d-separation oracles.
#'
#' These tests are based on univariate polygenic linear mixed models
#' \insertCite{almasy1998multipoint}{FamilyBasedPGMs}, with two components of variance: the polygenic or
#' family-specific random effect, which models the phenotypic variability across the families,
#' and the environmental or subject-specific error, which models phenotypic variability after removing the
#' familial aggregation effect.
#'
#'
#' @param phen.df A data.frame with phenotype variables of only sampled subjects.
#' Column names must be properly set with the names of the phenotypes.
#' @param covs.df A data.frame with covariates of only sampled subjects.
#' Column names must be properly set with the names of the covariates.
#' @param pedigrees A data.frame with columuns \code{famid}, \code{id}, \code{dadid},
#' \code{momid}, and \code{sex} columns for all sampled and non-sampled subjects.
#' @param sampled A logical vector in which element i indicates whether individual i was sampled or not.
#' @param fileID  A character string to be used as prefix in the filenames of RData
#' objects with the partial correlation results. Note that covariates are not
#' identified in these files.
#' @param dirToSave Path to the folder you want to save the output objects.
#' @param alpha The significance level to be used in the partial correlation tests.
#' @param max_cores An integer indicating the maximum number of CPU cores
#' to be used for parallel execution.
#' @param minK A scalar indicating the minimum dimension allowed
#' in the dimensionality reduction for confounding correction.
#' @param maxFC A scalar between 0 and 1, indicating the maximum fraction of
#' confounding allowed.
#' @param orthogonal A logical value indicating whether the transformation matrix used in the
#' confounding correction is orthogonal or not.
#' @param hidden_vars A logical value indicating if the causal structure learning method should
#' account for hidden variables. The \code{rfci} algorithm is used if hidden_vars is \code{TRUE} and
#' the \code{pc} algorithm is used otherwise.
#' @param maj.rule A logical value to be used in the \code{skeleton} function, indicating whether the
#' majority rule must be applied or not.
#' @param useGPU A logical value indicating whether GPU cores can be used for parallel execution.
#' @param debug A logical value indicating whether some debug messages can be shown.
#' @param savePlots A logical value indicating whether plots for the confounding correction must be generated.
#' @param logFile Optional file path and name to save progress and error messages.
#' If not provided and debug is True a default file is created in the dirToSave folder.
#'
#' @return Returns a list with the partial correlation matrices (\code{pcor}), the
#' adjacency matrices (\code{adjM}), and with the igraph objects representing the
#' undirected PGM (\code{udg}).
#'
#' @examples
#' \donttest{
#' data(scen3) # available simulated datasets are scen1, scen2, scen3, and scen4
#'
#' scenario = 3 # data was simulated according to scenario 3
#'
#' fam.nf <- scen3$fam.nf
#' pedigrees <- scen3$pedigrees
#' phen.df <- scen3$phen.df[[1]] # accessing the first replicate
#' covs.df <- NULL # no covariates were used in the simulation process.
#'
#' N <- sum(fam.nf) # total number of individuals
#' sampled <- rep(1, N) # in simulated data, all individuals were sampled.
#'
#' fileID <- paste0("scen", scenario)
#' dirToSave <- paste0("./objects-PC-", fileID, "/")
#' dir.create(dirToSave, showWarnings=FALSE)
#'
#' alpha = 0.05
#'
#' dags <- learnFamilyBasedDAGs(phen.df, covs.df, pedigrees, sampled,
#'                              fileID, dirToSave, alpha, max_cores=NULL,
#'                              minK=10, maxFC = 0.05, orthogonal=TRUE,
#'                              hidden_vars=FALSE, maj.rule=TRUE,
#'                              useGPU=FALSE, debug=TRUE, savePlots=FALSE)
#'
#' # the adjacency matrix of the learned directed acyclic genetic PGM
#' as(dags$g, "amat")
#'
#' # plotting the the learned directed acyclic genetic PGM as an `igraph` object:
#' plot.igraph(graph.adjacency(adjM_g), vertex.size=30, vertex.color="lightblue")
#'
#' # the adjacency matrix of the learned directed acyclic environmental PGM
#' as(dags$e, "amat")
#'
#' # plotting the the learned directed acyclic environmental PGM as an `igraph` object:
#' plot.igraph(graph.adjacency(adjM_e), vertex.size=30, vertex.color="lightblue")
#' }
#'
#' @references{
#'    \insertAllCited
#' }
#'
#' @importFrom Rdpack reprompt
#' @importFrom doMC registerDoMC
#' @import foreach
#'
#'@export learnFamilyBasedDAGs
learnFamilyBasedDAGs <- function(phen.df, covs.df, pedigrees, sampled, fileID, dirToSave, alpha=0.05,
                                 max_cores=NULL, minK=10, maxFC = 0.01, orthogonal=TRUE, hidden_vars=FALSE,
                                 maj.rule=TRUE, useGPU=FALSE, debug=TRUE, savePlots=FALSE, logFile=NULL) {
  if (debug && is.null(logFile)){
    logFile = paste0(dirToSave, "learnFamilyBasedDAGs_", format(Sys.time(), "%Y%m%d-%Hh%Mm%Ss"), ".log")
    cat("Learning DAGs from Family Data...\n", file=logFile, append=FALSE)
  }

  preprocessedPvaluesCSVFile <- paste0(dirToSave, fileID, "_preprPvalues.csv")

  preprocessFamilyBasedDAGs(phen.df, covs.df,
      pedigrees, sampled, alpha,
      fileID, dirToSave, preprocessedPvaluesCSVFile,
      max_cores=max_cores, minK=minK, maxFC = maxFC,
      orthogonal=orthogonal, savePlots=savePlots, useGPU=useGPU,
      debug=debug, logFile=logFile)

  if (is.null(max_cores)) {
    max_cores <- detectCores()
  }
  registerDoMC(max_cores)

  type = "" # a hack to avoid problems when checking --as-cran
  dags <- foreach(type = c("g","e", "t")) %dopar% {
    generateFamilyDAG(phen.df, covs.df, pedigrees, sampled,
                      preprocessedPvaluesCSVFile, type, alpha, hidden_vars=hidden_vars, maj.rule=maj.rule,
                      minK=minK, maxFC = maxFC, orthogonal=orthogonal,
                      dirToSave, fileID, savePlots=savePlots, useGPU=useGPU, debug=debug, logFile=logFile)
  }
  names(dags) <- c("g","e", "t")
  save(dags, file=paste0(dirToSave, fileID, "_dags.RData"))

  return(dags)
}

#' @title Conditional Independence Test for Family Data
#'
#' @description Computes the p-value of the test of the null hypothesis of zero genetic, environmental,
#' or total partial correlation between X and Y given S, if \code{dagg}, \code{dage}, or \code{dagt}
#' is set as TRUE, respectively.
#'
#' The unconfounded estimation and significance test for the genetic and environmental partial correlation
#' coefficients are described in detail in \insertCite{ribeiro2019family;textual}{FamilyBasedPGMs}.
#' They are based on univariate polygenic linear mixed models \insertCite{almasy1998multipoint}{FamilyBasedPGMs},
#' with two components of variance: the polygenic or family-specific random effect,
#' which models the phenotypic variability across the families,
#' and the environmental or subject-specific error, which models phenotypic variability after removing the
#' familial aggregation effect.
#'
#' This function can be used as a d-separation oracle in causal structure learning algorithms such as the
#' IC/PC algorithm when the variables are normal distributed and the observations are clusterized in families.
#'
#' @param x An integer value indicating the position of variable X.
#' @param y An integer value  indicating the position of variable Y.
#' @param S An integer vector with the positions of zero or more conditioning variables in S.
#' @param suffStat A list with the elements "phen.df", "covs.df", "pedigrees",
#' "minK", "maxFC", "orthogonal", "alpha", "dirToSave", "fileID",
#' "savePlots", and "useGPU", as described in function
#' \code{\link{learnFamilyBasedDAGs}}, "dagg", "dage", "dagt", to specify whether the
#' genetic, environmental, or total conditional independence test should be conducted,
#' and also "Phi2" and "Z", both obtained by the function
#' \code{\link{calculateBlockDiagonal2PhiMatrix}} with the argument "squaredRoot" set,
#' respectively, as FALSE and TRUE.
#'
#' @return The p-value of the test.
#'
#' @references{
#'    \insertAllCited
#' }
#'
#' @importFrom Rdpack reprompt
#' @importFrom utils write.table
#'
#'@export familyBasedCITest
familyBasedCITest <- function(x, y, S, suffStat) {
  ret <- familyBasedPCorTest(x, y, S, suffStat$phen.df, suffStat$covs.df,
           suffStat$pedigrees, suffStat$Phi2, suffStat$Z,
           suffStat$minK, suffStat$maxFC, suffStat$orthogonal, suffStat$alpha,
           suffStat$dirToSave, suffStat$fileID, suffStat$savePlots, useGPU=suffStat$useGPU)

  currow <- c(x=x,y=y,S=paste(S, collapse=" "), getRetVec(ret))

  if (suffStat$dagg == TRUE) {
     pvalue = as.numeric(ret["pcorg_p.value"])
    if (!is.null(suffStat$usedPvaluesCSVFile_g)) {
      write.table(as.matrix(t(currow)), file = suffStat$usedPvaluesCSVFile_g, sep = ";", row.names=FALSE, col.names = FALSE, append=TRUE, quote=FALSE)
    }
  } else if (suffStat$dage == TRUE) {
     pvalue = as.numeric(ret["pcore_p.value"])
    if (!is.null(suffStat$usedPvaluesCSVFile_e)) {
      write.table(as.matrix(t(currow)), file = suffStat$usedPvaluesCSVFile_e, sep = ";", row.names=FALSE, col.names = FALSE, append=TRUE, quote=FALSE)
    }
  } else {
     pvalue = as.numeric(ret["pcort_p.value"])
    if (!is.null(suffStat$usedPvaluesCSVFile_t)) {
      write.table(as.matrix(t(currow)), file = suffStat$usedPvaluesCSVFile_t, sep = ";", row.names=FALSE, col.names = FALSE, append=TRUE, quote=FALSE)
    }
  }
  return(pvalue)
}

#' @title Kinship matrix of a known pedigree structure
#'
#' @description Computes the block diagonal kinship matrix with the degree of relatedness
#' between individuals all individuals of \eqn{F} families, multiplied by 2.
#' It is usually denoted as \eqn{2 \bm\Phi}, in which \eqn{\bm \Phi = diag\{\Phi^{(f)}, f=1,\ldots,F\}} and
#' each entry \eqn{\Phi^{(f)}_{ij}} is the probability that two alleles sampled at random from individuals
#' \eqn{i} and \eqn{j} of family \eqn{f} are identical by descent.
#' Thus, when the pedigree structure is known, \eqn{\Phi^{(f)}_{ii} = 1/2},
#' \eqn{\Phi^{(f)}_{ij} = 1/4} if \eqn{i} and \eqn{j} are siblings or if one of them is a parent of the other,
#' \eqn{\Phi^{(f)}_{ij} = 1/8} if one of them is a grand-parent of the other, and so on
#' \insertCite{lange2003mathematical}{FamilyBasedPGMs}.
#'
#' @param ped  A data.frame with with the pedigrees of the families, with columuns
#' \code{famid}, \code{id}, \code{dadid}, \code{momid}, and \code{sex} for all sampled
#' and non-sampled subjects.
#' @param squaredRoot a logical value indicating if the square root of the
#' kinship matrix  \eqn{2 \bm\Phi} (i.e., the Z matrix) must be computed.
#' @param sampled A logical vector in which element i indicates whether individual i was sampled or not.
#'
#' @return The kinship matrix \eqn{2 \bm\Phi} or its squared root, if squaredRoot is TRUE.
#'
#' @references{
#'    \insertAllCited
#' }
#'
#' @importFrom Rdpack reprompt
#'
#' @export calculateBlockDiagonal2PhiMatrix
calculateBlockDiagonal2PhiMatrix <- function(ped, squaredRoot=FALSE, sampled=NULL) {
  famids <- unique(ped$famid)

  n <- dim(ped)[1]
  if (!is.null(sampled)) {
    n <- length(which(sampled == 1))
  }

  Phi2 <- matrix(0, n, n) # initializing matrix corresponding to
  # 2*Phi = diag(2*Phi_i, i = 1,.., nfam)

  if (!is.null(sampled)) {
    ped <- cbind(ped, sampled=sampled)
  }

  lastid <- 0
  for (famid in famids) {
    pedf <- ped[ped$famid == famid,]
    nf <- dim(pedf)[1]
    Phi <- kinship(id=pedf$id, dadid=pedf$dadid, momid=pedf$momid, sex=pedf$sex)

    if (!is.null(sampled)) {
      sampled_ids <- which(pedf$sampled == 1)
      nf <- length(sampled_ids)
      Phi <- Phi[sampled_ids, sampled_ids]
    }

    if (nf > 0) {
      blockDiagMatrix <- 2*Phi # Phi2
      if (squaredRoot) {
        blockDiagMatrix <- matrixSqrt(blockDiagMatrix) # sqrt(Phi2)
      }

      Phi2[(lastid+1):(lastid+nf),(lastid+1):(lastid+nf)] <- blockDiagMatrix
      lastid <- lastid + nf
    }
  }

  if (!is.null(sampled)) {
    dimnames(Phi2)[[1]] <-(ped[which(sampled==1),"id"])
  } else {
    dimnames(Phi2)[[1]] <- 1:n
  }

  dimnames(Phi2)[[2]] <-  dimnames(Phi2)[[1]]

  return(Phi2)
}


