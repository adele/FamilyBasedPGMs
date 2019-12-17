#' @importFrom stats as.formula cor.test model.matrix residuals
#' @importFrom coxme lmekin ranef
#' @importFrom stats sd
#' @importFrom MASS ginv
getFamilyBasedPCorResults <- function(Xind, Yind, Sinds, phen.df, covs.df, pedigrees, Phi2, Z,
                                      minK=10, maxFC = 0.01, orthogonal=TRUE, fileID=NULL, dirToSave=NULL,
                                      useGPU=FALSE, debug=TRUE, logFile="./log.txt") {

  # pedigrees must be sampled

  #######################
  # Computing rho_X,Y|S #
  #######################

  N <- dim(phen.df)[1]

  phen.names <- colnames(phen.df)
  Xname <- phen.names[Xind]
  Yname <- phen.names[Yind]
  Snames <- phen.names[Sinds]
  if (!is.null(covs.df)) {
    Snames <- c(Snames, colnames(covs.df))
    phen.df <- cbind(phen.df, covs.df)
  }

  if (debug) {
    cat(paste0("Computing rho_", Xname, ",", Yname, "|{", paste0(Snames, collapse=","), "}"),
        file=logFile, append=TRUE, sep="\n")
  }

  ##################################
  # Checking previous computations #
  ##################################

  results_g_filepath <- ""
  if (!is.null(fileID) && !is.null(dirToSave)) {
    if (Xind < Yind) {
      results_g_filepath <- paste0(dirToSave, fileID, "_results_g_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), ".rds")
    } else {
      results_g_filepath <- paste0(dirToSave, fileID, "_results_g_", Yind, "_", Xind, "_", paste(Sinds, collapse=""), ".rds")
    }
  }

  results_g_done = FALSE
  if (results_g_filepath != "" && file.exists(results_g_filepath)) {
    results_g <- readRDS(file=results_g_filepath) # loading results_g
    ret = checkUnconfoundedEstimatesList(results_g, maxFC, N)
    if (ret == "done") { # results are complete
      results_g_done = TRUE
    } else if (ret == "incomplete") {
      if (debug) {
        cat(paste0("Updating results_g for ", fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), "\n"),
            file=logFile, append=TRUE, sep="\n")
      }
    } else if (ret == "error") {
      file.remove(results_g_filepath)
      if (debug) {
        cat(paste0("Redoing results_g for ", fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), "\n"),
            file=logFile, append=TRUE, sep="\n")
      }
    }
  }


  results_e_filepath <- ""
  if (!is.null(fileID) && !is.null(dirToSave)) {
    if (Xind < Yind) {
      results_e_filepath <- paste0(dirToSave, fileID, "_results_e_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), ".rds")
    } else {
      results_e_filepath <- paste0(dirToSave, fileID, "_results_e_", Yind, "_", Xind, "_", paste(Sinds, collapse=""), ".rds")
    }
  }


  results_e_done = FALSE
  if (results_e_filepath != "" && file.exists(results_e_filepath)) {
    results_e <- readRDS(file=results_e_filepath) # loading results_e
    ret = checkUnconfoundedEstimatesList(results_e, maxFC, N)
    if (ret == "done") { # results are complete
      results_e_done = TRUE
    } else if (ret == "incomplete") {
      if (debug) {
        cat(paste0("Updating results_e for ", fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), "\n"),
            file=logFile, append=TRUE, sep="\n")
      }
    } else if (ret == "error") {
      file.remove(results_e_filepath)
      if (debug) {
        cat(paste0("Redoing results_e for ", fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), "\n"),
            file=logFile, append=TRUE, sep="\n")
      }
    }
  }



  ###################
  # Model for X | S #
  ###################

  if (debug) {
    cat(paste0("Fitting ", Xname, " to ", "|{", paste0(Snames, collapse=","), "}"),
        file=logFile, append=TRUE, sep="\n")
  }

  formulaXS <- as.formula(paste(Xname, " ~ 1 + (1 | pedigrees$id)"))
  if (length(Sinds) > 0)
    formulaXS <-  as.formula(paste(Xname,  "~ 1 + ", paste(paste(Snames, "+ "),collapse=""), " (1 | pedigrees$id)"))

  fitFilePath <- ""
  if (!is.null(fileID) && !is.null(dirToSave)) {
    fitFilePath <- paste0(dirToSave, fileID, "_fit_", Xind, "_", paste(Sinds, collapse=""), ".rds")
  }

  if (file.exists(fitFilePath)) {
    curfit <- readRDS(fitFilePath)
  } else {
    curfit = lmekin(formulaXS, data=phen.df, varlist=Phi2, method="REML")
    if (fitFilePath != "") {
      saveRDS(curfit, file=fitFilePath) # we are always saving an object named curfit
    }
  }
  fitXS <- curfit

  sigma2XS_g <- as.numeric(fitXS$vcoef)
  if (sigma2XS_g < 1e-06) {
    sigma2XS_g = 0
  }
  sigma2XS_e <- as.numeric(fitXS$sigma)^2
  if (sigma2XS_e < 1e-06) {
    sigma2XS_e = 0
  }
  sigma2XS <- (sigma2XS_g + sigma2XS_e)

  #print(paste0("sigma2g_(", Xname, "| {", paste0(Snames, collapse=","), "}) = ", sigma2XS_g))
  #print(paste0("sigma2e_(", Xname, "| {", paste0(Snames, collapse=","), "}) = ", sigma2XS_e))

  if (!results_g_done || !results_e_done) {
    if (sigma2XS_g != 0 || sigma2XS_e != 0) {
      DXS <- diag(N) * sigma2XS_g # It is 0 if sigma2g=0
      RXS <- diag(N) * sigma2XS_e # It is 0 if sigma2e=0
      MXSinv <- (Z %*% DXS %*% t(Z) + RXS)
      MXS <- ginv(MXSinv)
      #V <- Minv * sigma2

      Sformula <- paste(paste(Snames, " + "), collapse="")
      Sformula <- as.formula(paste("~ ", Sformula, " 1"))
      S <- model.matrix(Sformula, data=phen.df)
      QXS <- MXS - MXS %*% S %*% ginv(t(S) %*% MXS %*% S) %*% t(S) %*% MXS

      temp = Z %*% DXS %*% t(Z) %*% QXS
      BgXS = temp %*% Z %*% DXS %*% t(Z) %*% t(temp) # B-A - corresponding to Zb - It is 0 if sigma2g=0
      AgXS = temp %*% RXS %*% t(temp) #- corresponding to e - It is 0 if sigma2e=0.
      BgXS = BgXS + AgXS # To be maximized in the same time that Ag is minimized


      temp = RXS %*% QXS
      BeXS = temp %*% RXS %*% t(temp) # B-A - corresponding to e - It is 0 if sigma2e=0
      AeXS = temp %*% Z %*% DXS %*% t(Z) %*% t(temp) # corresponding to Zb = It is 0 if sigma2g=0
      BeXS = BeXS + AeXS

      ranefXS <- Z %*% DXS %*% t(Z) %*% QXS %*% phen.df[,Xname]
      #head(ranef(fitXS)[[1]])
      #head(ranefXS)

      resXS <- RXS %*% QXS %*% phen.df[,Xname]
      #head(residuals(fitXS))
      #head(resXS)

      margerrXS <- MXSinv %*% QXS %*% phen.df[,Xname]
      #head(resXS + ranefXS)
      #head(margerrXS)
    }
  }

  ##################
  # Model for Y| S #
  ##################

  if (debug) {
    cat(paste0("Fitting ", Yname, " to ", "|{", paste0(Snames, collapse=","), "}"),
        file=logFile, append=TRUE, sep="\n")
  }

  formulaYS <- as.formula(paste(Yname, " ~ 1 + (1 | pedigrees$id)"))
  if (length(Sinds) > 0)
    formulaYS <-  as.formula(paste(Yname,  "~ 1 + ", paste(paste(Snames, "+ "),collapse=""), " (1 | pedigrees$id)"))

  fitFilePath <- ""
  if (!is.null(fileID) && !is.null(dirToSave)) {
    fitFilePath <- paste0(dirToSave, fileID, "_fit_", Yind, "_", paste(Sinds, collapse=""), ".rds")
  }

  if (file.exists(fitFilePath)) {
    curfit <- readRDS(fitFilePath)
  } else {
    curfit = lmekin(formulaYS, data=phen.df, varlist=Phi2, method="REML")
    if (fitFilePath != "") {
      saveRDS(curfit, file=fitFilePath) # we are always saving an object named curfit
    }
  }
  fitYS <- curfit


  sigma2YS_g <- as.numeric(fitYS$vcoef)
  if (sigma2YS_g < 1e-06) {
    sigma2YS_g = 0
  }
  sigma2YS_e <- as.numeric(fitYS$sigma)^2
  if (sigma2YS_e < 1e-06) {
    sigma2YS_e = 0
  }
  sigma2YS <- (sigma2YS_g + sigma2YS_e)

  if (!results_g_done || !results_e_done) {
    if (sigma2YS_g != 0 || sigma2YS_e != 0) {
      DYS <- diag(N) * sigma2YS_g
      RYS <- diag(N) * sigma2YS_e
      MYSinv <- (Z %*% DYS %*% t(Z) + RYS)
      MYS <- ginv(MYSinv)

      Sformula <- paste(paste(Snames, " + "), collapse="")
      Sformula <- as.formula(paste("~ ", Sformula, " 1"))
      S <- model.matrix(Sformula, data=phen.df)
      QYS <- MYS - MYS %*% S %*% ginv(t(S) %*% MYS %*% S) %*% t(S) %*% MYS

      temp = Z %*% DYS %*% t(Z) %*% QYS
      BgYS = temp %*% Z %*% DYS %*% t(Z) %*% t(temp) # B-A - corresponding to Zb
      AgYS = temp %*% RYS %*% t(temp) #- corresponding to e
      BgYS = BgYS + AgYS # To be maximized in the same time that Ag is minimized

      temp = RYS %*% QYS
      BeYS = temp %*% RYS %*% t(temp) # B-A - corresponding to e
      AeYS = temp %*% Z %*% DYS %*% t(Z) %*% t(temp) # corresponding to Zb
      BeYS = BeYS + AeYS


      ranefYS <- Z %*% DYS %*% t(Z) %*% QYS %*% phen.df[,Yname]
      #head(ranef(fitYS)[[1]])
      #head(ranefYS)

      resYS <- RYS %*% QYS %*% phen.df[,Yname]
      #head(residuals(fitYS))
      #head(resYS)

      margerrYS <- MYSinv %*% QYS %*% phen.df[,Yname]
      #head(resYS + ranefYS)
      #head(margerrYS)
    }
  }

  margerrXS = ranef(fitXS)[[1]] + residuals(fitXS)
  margerrYS = ranef(fitYS)[[1]] + residuals(fitYS)
  rho_t = cor.test(margerrXS, margerrYS)
  rho_t

  #############
  # rho_X,Y|S #
  #############

  if (!is.null(fileID) && debug) {
    cat(paste0("Processing results_g for ", fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), "\n"),
        file=logFile, append=TRUE, sep="\n")
  }

  if (!results_g_done) {
    if (sigma2XS_g == 0 || sigma2YS_g == 0) {
      # It is not possible to compute rho_g, since some standard deviantion is zero
      results_g <- list()
      results_g[[1]] <- list(k=N, p.value=1, estimate=0, sdx=0, sdy=0, FCxs=NaN, FCys=NaN)
    } else {
      # Since both sigma2XS_g and sigma2YS_g are not zero, all variables above were defined
      if (sigma2XS_e == 0 && sigma2YS_e == 0) {
        # If both sigma2XS_e and sigma2YS_e are zero, then there is no confounding in g
        rho_g = cor.test(ranefXS, ranefYS)
        results_g <- list()
        results_g[[1]] <- list(k=N, p.value=rho_g$p.value, estimate=rho_g$estimate, sdx=sd(ranefXS), sdy=sd(ranefYS), FCxs=0, FCys=0)
      } else {
        preComputedResults <- NULL
        if (results_g_filepath != "" && file.exists(results_g_filepath)) {
          results_g <- readRDS(file=results_g_filepath) # loading existing results_g
          preComputedResults <- results_g
        }
        if (debug) {
          cat(paste0("Getting unconfounded estimates for ", Xname, ",", Yname, "|{", paste0(Snames, collapse=","), "}_g"),
              file=logFile, append=TRUE, sep="\n")
        }
        results_g <- getUnconfoundedEstimates(BgXS, BgYS, AgXS, AgYS, maxFC, minK, orthogonal,
                                              ranefXS, ranefYS, preComputedResults, useGPU,
                                              debug=debug, logFile=logFile)
        if (length(results_g) == 0) {
          rho_g = cor.test(ranefXS, ranefYS)
          results_g[[1]] <- list(k=N, p.value=rho_g$p.value, estimate=rho_g$estimate, sdx=sd(ranefXS), sdy=sd(ranefYS), FCxs=0, FCys=0)
        }
      }
    }
    # saving results_g
    if (results_g_filepath != "" ) {
      if (debug) {
        cat(paste0("Saving results_g for ", fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), "\n"),
            file=logFile, append=TRUE, sep="\n")
      }
      saveRDS(results_g, file=results_g_filepath)
    }
  }

  if (!is.null(fileID) && debug) {
    cat(paste0("Processing results_e for ", fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), "\n"),
        file=logFile, append=TRUE, sep="\n")
  }

  if (!results_e_done) {
    if (sigma2XS_e == 0 || sigma2YS_e == 0) {
      # It is not possible to compute rho_e, since some standard deviantion is zero
      results_e <- list()
      results_e[[1]] <- list(k=N, p.value=1, estimate=0, sdx=0, sdy=0, FCxs=NaN, FCys=NaN)
    } else {
      # Since both sigma2XS_e and sigma2YS_e are not zero, all variables above were defined
      if (sigma2XS_g == 0 && sigma2YS_g == 0) {
        # If both sigma2XS_g and sigma2YS_g are zero, then there is no confounding in e
        rho_e = cor.test(resXS, resYS)
        results_e <- list()
        results_e[[1]] <- list(k=N, p.value=rho_e$p.value, estimate=rho_e$estimate,
                               sdx=sd(resXS), sdy=sd(resYS), FCxs=0, FCys=0)
      } else {
        preComputedResults <- NULL
        if (results_e_filepath != "" && file.exists(results_e_filepath)) {
          results_e <- readRDS(file=results_e_filepath) # loading existing results_e
          preComputedResults <- results_e
        }
        if (debug) {
          cat(paste0("Getting unconfounded estimates for ", Xname, ",", Yname, "|{", paste0(Snames, collapse=","), "}_e"),
              file=logFile, append=TRUE, sep="\n")
        }
        results_e <- getUnconfoundedEstimates(BeXS, BeYS, AeXS, AeYS, maxFC, minK,
                                              orthogonal, resXS, resYS, preComputedResults, useGPU,
                                              debug=debug, logFile=logFile)
        if (length(results_e) == 0) {
          rho_e = cor.test(resXS, resYS)
          results_e[[1]] <- list(k=N, p.value=rho_e$p.value, estimate=rho_e$estimate,
                                 sdx=sd(resXS), sdy=sd(resYS), FCxs=0, FCys=0)
        }
      }
    }
    # saving results_e
    if (results_e_filepath != "" ) {
      saveRDS(results_e, file=results_e_filepath)
      if (debug) {
        cat(paste0("Saving results_e for ", fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), "\n"),
            file=logFile, append=TRUE, sep="\n")
      }
    }
  }

  return(list(results_g=results_g, results_e=results_e,
          sigma2XS_e=sigma2XS_e, sigma2YS_e=sigma2YS_e,
          sigma2XS_g=sigma2XS_g, sigma2YS_g=sigma2YS_g,
          k_t=N, rho_t=rho_t$estimate, p.value_t=rho_t$p.value,
          sdx_t=sd(margerrXS), sdy_t=sd(margerrYS),
          Xind=Xind, Yind=Yind, Sinds=Sinds))
}



checkUnconfoundedEstimatesList <- function(results, maxFC, N) {
  redoList <- c()
  for (i in 1:length(results)) {
    curResultsNames <- names(results[[i]])
    if (is.null(curResultsNames) || length(curResultsNames) != 7 ||
      sum(curResultsNames == c("k","p.value","estimate","sdx","sdy","FCxs","FCys")) != 7) {
      redoList <- c(redoList, i)
    }
  }

  if (length(redoList) > 0) {
    return("error")
  }

  FCysVec <- as.numeric(unlist(lapply(results, '[', "FCys")))
  FCxsVec <- as.numeric(unlist(lapply(results, '[', "FCxs")))
  kVec <- as.numeric(unlist(lapply(results, '[', "k")))

  if ((length(results) == 1 && is.nan(FCysVec) && is.nan(FCxsVec)) ||
      (length(FCysVec) == length(FCxsVec) && (max(FCysVec, na.rm=T) >= maxFC || max(FCxsVec, na.rm=T) >= maxFC)) ||
      (length(FCysVec) == length(FCxsVec) && max(kVec) == N)) {
    return("done")
  } else {
    return("incomplete")
  }
}


#' @importFrom stats sd
#' @importFrom R.utils withTimeout
getUnconfoundedEstimates <- function(Bxs, Bys, Axs, Ays, maxFC, minK, orthogonal, rx, ry,
                                     preComputedResults=NULL, useGPU=FALSE, gpuCores=5,
                                     debug=TRUE, logFile="./log.txt") {
  n = dim(Bxs)[1]
  k = minK
  FCxs = FCys = 0

  if (is.null(preComputedResults)) {
    resultsList <- list()
  } else {
    resultsList <- preComputedResults
  }

  if (useGPU && requireNamespace("gpuR", quietly = TRUE)) {

    if (debug) {
      cat("\nUsing GPU in getUnconfoundedEstimates()..\n",
          file=logFile, append=TRUE, sep="\n")
    }

    curK <- minK
    kList <- c()
    while(curK <= n) {
      kList <- c(kList, curK)
      curK <- round(curK + 0.1*curK)
    }

    funGPU <- function(k) {
      gpuR::setContext(3L)
      tryCatch({
        W = mcrotate_orthogonal(Axs + Ays, Bxs + Bys, k, orthogonal, useGPU)
        FCxs <- computeFC(W,Axs,Bxs,n,useGPU)
        FCys <- computeFC(W,Ays,Bys,n,useGPU)

        rot_rx = t(W)%*% rx
        rot_ry = t(W)%*% ry
        ctest = cor.test(rot_rx, rot_ry)

        return(list(k=k, p.value=ctest$p.value, estimate=ctest$estimate, sdx=sd(rot_rx),
                    sdy=sd(rot_ry), FCxs=FCxs, FCys=FCys))
      },
      error = function(e) {
        cat(paste0("Error: ", e),
            file=logFile, append=TRUE, sep="\n")
        return(list(e))
      })
    }

    cl <- gpuR::makeCluster(gpuCores, outfile=logFile) # log can only be seen using outfile=""
    gpuR::clusterExport(cl, c("funGPU", "mcrotate_orthogonal", "mymcrotate",
                              "myginv", "computeFC", "Bxs", "Bys", "Axs", "Ays",
                              "maxFC", "minK", "orthogonal", "rx", "ry", "useGPU", "n"), envir=environment())
    gpuR::setContext(3L)
    resultsList <- gpuR::parLapply(cl,kList,funGPU)
    gpuR::stopCluster(cl)
  } else {
    nTries <- 0
    i <- length(resultsList)+1
    kVec <- as.numeric(unlist(lapply(resultsList, '[', "k")))
    while (FCxs <= maxFC && FCys <= maxFC && k <= n && nTries <= 5) {
      if (!(k %in% kVec)) {
        retValue <- tryCatch({
          W = mcrotate_orthogonal(Axs + Ays, Bxs + Bys, k, orthogonal, useGPU)
  #        W =  withTimeout({mcrotate_orthogonal(Axs + Ays, Bxs + Bys, k, orthogonal)}, timeout=120)
          if (debug) {
            cat(paste0("k:", k, "\n"),
                file=logFile, append=TRUE, sep="\n")
          }

          FCxs <- withTimeout({computeFC(W,Axs,Bxs,n,useGPU)}, timeout=60)
          if (debug) {
            cat(paste0("FCxs:", FCxs, "\n"),
                file=logFile, append=TRUE, sep="\n")
          }

          FCys <- withTimeout({computeFC(W,Ays,Bys,n,useGPU)}, timeout=60)
          if (debug) {
            cat(paste0("FCys:", FCys, "\n"),
                file=logFile, append=TRUE, sep="\n")
          }

          rot_rx = t(W)%*% rx
          rot_ry = t(W)%*% ry
          ctest = cor.test(rot_rx, rot_ry)
          if (debug) {
            cat(paste0("rho:", ctest$estimate, ". p-value:", ctest$p.value, "\n"),
                file=logFile, append=TRUE, sep="\n")
          }

          resultsList[[i]] <- list(k=k, p.value=ctest$p.value, estimate=ctest$estimate, sdx=sd(rot_rx), sdy=sd(rot_ry), FCxs=FCxs, FCys=FCys)
          i = i+1
          0
        },
        error = function(e) {
          if (debug) {
            cat(paste0("Error for k=", k, ", after ", nTries, " tries: ", e),
                file=logFile, append=TRUE, sep="\n")
          }
          return(-1)
        })
        if (retValue == -1) {
          nTries <- nTries + 1
        } else {
          nTries <- 0
        }
      }
      k = round(k + 0.1*k)
    }
  }

  return(resultsList)
}

computeFC <- function(W, A, B, n, useGPU=FALSE) {
  if (useGPU && requireNamespace("gpuR", quietly = TRUE)) {
    A <- gpuR::gpuMatrix(A)
    B <- gpuR::gpuMatrix(B)
    W <- gpuR::gpuMatrix(W)
  }
  FC <- sum(as.vector(diag(myginv(as.matrix(t(W) %*% B %*% W), useGPU) %*% (t(W) %*% A %*% W))))/n
  return(FC)
}

myginv <- function (X, useGPU = FALSE) {
  tol = sqrt(.Machine$double.eps)

  if (useGPU && requireNamespace("gpuR", quietly = TRUE)) {
    X <- gpuR::gpuMatrix(X)
  }

  Xsvd <- svd(X)
  decOrder <- order(as.vector(Xsvd$d), decreasing=TRUE)
  Xsvd.d <- as.vector(Xsvd$d)[decOrder]
  Xsvd.u <- as.matrix(Xsvd$u)[,decOrder]
  Xsvd.v <- as.matrix(Xsvd$v)[,decOrder]

  # X - Xsvd.u %*% diag(Xsvd.d) %*% t(Xsvd.v)

  Positive <- Xsvd.d > max(tol * Xsvd.d[1L], 0)

  ds <- diag(1/Xsvd.d[Positive])
  us <- as.matrix(Xsvd.u[, Positive, drop=FALSE])
  vs <- as.matrix(Xsvd.v[, Positive, drop=FALSE])

  if (useGPU && requireNamespace("gpuR", quietly = TRUE)) {
    ds <- gpuR::gpuMatrix(ds)
    us <- gpuR::gpuMatrix(us)
    vs <- gpuR::gpuMatrix(vs)
  }

  X.ginv <- as.matrix(vs %*% ds %*% t(us))
  return(X.ginv)
}


#' @importFrom Matrix rankMatrix
mymcrotate <- function(A, B, s, useGPU = FALSE) {
  r <- rankMatrix(B)
  if(is.null(s)) s <- r

  if (useGPU && requireNamespace("gpuR", quietly = TRUE)) {
    A <- gpuR::gpuMatrix(A)
    B <- gpuR::gpuMatrix(B)
  }

  B.svd <- svd(B)
  decOrder <- order(as.vector(B.svd$d), decreasing=TRUE)
  B.svd.d <- as.vector(B.svd$d)[decOrder]
  B.svd.u <- as.matrix(B.svd$u)[,decOrder]
  B.svd.v <- as.matrix(B.svd$v)[,decOrder]

  svDiag <- diag(1/sqrt(B.svd.d[1:r]))
  Tr <- B.svd.u[, 1:r]

  if (useGPU && requireNamespace("gpuR", quietly = TRUE)) {
    Tr <- gpuR::gpuMatrix(Tr)
    svDiag <- gpuR::gpuMatrix(svDiag)
  }

  A.star <- svDiag %*% t(Tr) %*% A %*% Tr %*% svDiag
 # Note that B.star == I:
 # svDiag %*% t(Tr) %*% B %*% Tr %*% svDiag )

  A.star.svd <- svd(A.star) #fast.svd(A.star)
  decOrder <- order(as.vector(A.star.svd$d), decreasing=TRUE)
  A.star.svd.d <- as.vector(A.star.svd$d)[decOrder]
  A.star.svd.u <- as.matrix(A.star.svd$u)[,decOrder]
  A.star.svd.v <- as.matrix(A.star.svd$v)[,decOrder]

  index <- seq(r, length.out = s, by = -1)
  index <- sort(index[index >= 0])

  AstarU <- A.star.svd.u[,index]
  if (useGPU && requireNamespace("gpuR", quietly = TRUE)) {
    AstarU <- gpuR::gpuMatrix(AstarU)
  }

  W <- Tr %*% svDiag %*% AstarU
  return(as.matrix(W))
}


mcrotate_orthogonal <- function (A, B, s, orthogonal, useGPU=FALSE) {
  W <- mymcrotate(A, B, s, useGPU)

  if(orthogonal) {
    # NOTE: qr.gpuMatrix(W): non-square matrix not currently supported for 'qr'
    #    if (useGPU) {
    #      require(gpuR)
    #      W <- gpuMatrix(W)
    #    }

    W.QR = qr(W)
    W.Q = qr.Q(W.QR)
    # W.R = qr.R(W.QR)
    W = W.Q
  }
  return(as.matrix(W))
}

#' @importFrom utils combn
selectBestUDGPCorMatrices <- function(phen.p, alpha, maxFC, fileID, resultsDir,
                                      savePlots=TRUE, variables=NULL, debug=TRUE, logFile="./log.txt") {

  nei_pCor_g_estimates <- matrix(NA, phen.p, phen.p)
  if (!is.null(variables) && length(variables) == phen.p) {
    colnames(nei_pCor_g_estimates) <- variables
    rownames(nei_pCor_g_estimates) <- variables
  }

  nei_pCor_g_pvalues <- nei_pCor_g_estimates
  nei_pCor_g_ks <- nei_pCor_g_estimates

  nei_pCor_e_estimates <- nei_pCor_g_estimates
  nei_pCor_e_pvalues <- nei_pCor_g_estimates
  nei_pCor_e_ks <- nei_pCor_g_estimates

  nei_pCor_t_estimates <- nei_pCor_g_estimates
  nei_pCor_t_pvalues <- nei_pCor_g_estimates

  pseq <- 1:phen.p
  xyinds <- combn(pseq,2)
  for (xyind in 1:dim(xyinds)[2]) {
    Xind <- xyinds[1,xyind]
    Yind <- xyinds[2,xyind]
    Sinds <- pseq[-c(Xind, Yind)]

    resultFile <- paste0(resultsDir, fileID, "_", Xind, "_", Yind, "_", paste(Sinds, collapse=""), ".rds")
    if (file.exists(resultFile)) {
      if (debug) {
        cat(paste0("Processing ", resultFile),
            file=logFile, append=TRUE, sep="\n")
      }
      curResult <- readRDS(resultFile)

      if (savePlots) {
        saveEdgeResultsPlot(curResult, resultsDir, fileID, Xind, Yind, Sinds, alpha, maxFC=maxFC,
                            debug=debug, logFile=logFile)
      }
      resultsList_g <- curResult$results_g
      resultsList_e <- curResult$results_e
      if (debug) {
        cat(paste0("Selecting indices for edge ", Xind, "-", Yind),
            file=logFile, append=TRUE, sep="\n")
      }
      selectedInds <- selectBestResult(curResult$results_g, curResult$results_e, curResult$rho_t,
                                       curResult$p.value_t, alpha, maxFC=maxFC, debug=debug, logFile=logFile)

      if (debug) {
        cat(paste0("done for edge ", Xind, "-", Yind),
            file=logFile, append=TRUE, sep="\n")
      }

      nei_pCor_g_estimates[Xind,Yind] <- resultsList_g[[selectedInds[1]]]$estimate
      nei_pCor_g_pvalues[Xind,Yind] <- resultsList_g[[selectedInds[1]]]$p.value
      nei_pCor_g_ks[Xind,Yind] <- resultsList_g[[selectedInds[1]]]$k

      nei_pCor_g_estimates[Yind,Xind] <- nei_pCor_g_estimates[Xind,Yind]
      nei_pCor_g_pvalues[Yind,Xind] <- nei_pCor_g_pvalues[Xind,Yind]
      nei_pCor_g_ks[Yind,Xind] <- nei_pCor_g_ks[Xind,Yind]

      nei_pCor_e_estimates[Xind,Yind] <- resultsList_e[[selectedInds[2]]]$estimate
      nei_pCor_e_pvalues[Xind,Yind] <- resultsList_e[[selectedInds[2]]]$p.value
      nei_pCor_e_ks[Xind,Yind] <- resultsList_e[[selectedInds[2]]]$k

      nei_pCor_e_estimates[Yind,Xind] <- nei_pCor_e_estimates[Xind,Yind]
      nei_pCor_e_pvalues[Yind,Xind] <- nei_pCor_e_pvalues[Xind,Yind]
      nei_pCor_e_ks[Yind,Xind] <- nei_pCor_e_ks[Xind,Yind]

      nei_pCor_t_estimates[Xind,Yind] <- curResult$rho_t
      nei_pCor_t_pvalues[Xind,Yind] <- curResult$p.value_t

      nei_pCor_t_estimates[Yind,Xind] <- nei_pCor_t_estimates[Xind,Yind]
      nei_pCor_t_pvalues[Yind,Xind] <- nei_pCor_t_pvalues[Xind,Yind]

    } else {
      if (debug) {
        cat(paste0("Could not find file: ", resultFile),
            file=logFile, append=TRUE, sep="\n")
      }
    }
  }

  return(list(pCor_g=list(estimates=nei_pCor_g_estimates, pvalues=nei_pCor_g_pvalues, k=nei_pCor_g_ks),
    pCor_e=list(estimates=nei_pCor_e_estimates, pvalues=nei_pCor_e_pvalues, k=nei_pCor_e_ks),
    pCor_t=list(estimates=nei_pCor_t_estimates, pvalues=nei_pCor_t_pvalues)))
}



#' @importFrom grDevices dev.off png
#' @importFrom graphics abline lines par plot
saveEdgeResultsPlot <- function(curResult, folderToSave, fileID, Xind, Yind, Sinds, alpha, maxFC,
                                debug=TRUE, logFile="./log.txt") {
    if (length(curResult$results_g) == 0 || length(curResult$results_e) == 0) {
      if (debug) {
        cat(paste0("Results for edge ", Xind, "-", Yind, "_", paste(Sinds, collapse=""), " are empty!"),
            file=logFile, append=TRUE, sep="\n")
      }
      return
    }

    selectedInds <- selectBestResult(curResult$results_g, curResult$results_e, curResult$rho_t,
                                     curResult$p.value_t, alpha, maxFC=maxFC, debug=debug, logFile=logFile)


    png(filename=paste0(folderToSave, fileID, "_edge_",Xind,"-",Yind, "_", paste(Sinds, collapse=""), ".png"),
        width=1280, height=768)
    par(mfcol=c(1,2), mar=c(5.1,4.1,7.1,2.1))
    for (type in c("g", "e")) {
      resultsList <- list()
      if (type == "g") {
        resultsList <- curResult$results_g

        estSigma2_XS <- curResult$sigma2XS_g
        estSigma2_YS <- curResult$sigma2YS_g
        selectedIndex <- selectedInds[1]
      } else if (type == "e") {
        resultsList <- curResult$results_e

        estSigma2_XS <- curResult$sigma2XS_e
        estSigma2_YS <- curResult$sigma2YS_e
        selectedIndex <- selectedInds[2]
      }

      if (length(resultsList) > 1) {
        pvaluesVec <- as.numeric(unlist(lapply(resultsList, '[', "p.value")))
        estimatesVec <- as.numeric(unlist(lapply(resultsList, '[', "estimate")))
        sdxVec <- as.numeric(unlist(lapply(resultsList, '[', "sdx")))
        sdyVec <- as.numeric(unlist(lapply(resultsList, '[', "sdy")))
        FCysVec <- as.numeric(unlist(lapply(resultsList, '[', "FCys")))
        FCxsVec <- as.numeric(unlist(lapply(resultsList, '[', "FCxs")))

        if (debug) {
          cat(paste0("Edge ", Xind, "-", Yind,  type, "; est: ", round(estimatesVec[selectedIndex],4),
                     " p=", round(pvaluesVec[selectedIndex],4)),
              file=logFile, append=TRUE, sep="\n")
        }

        plot(FCysVec + FCxsVec, -log10(pvaluesVec),
          col="blue", type="l", xlim=c(0,2*maxFC),ylim=c(-10,15),
          main=paste0("Edge ", Xind, "-", Yind,  type,
          "\nEstimated pCor_t: ", round(curResult$rho_t,4), " (P=", round(curResult$p.value_t,4), ")",
          "\nEstimated sigma2_X|S: ", round(estSigma2_XS,4), ", Estimated sigma2_Y|S: ", round(estSigma2_YS,4),"\nEstimated pCor_", type, " (k=", resultsList[[selectedIndex]][[1]],"): ",
          round(estimatesVec[selectedIndex],4), " (P=", round(pvaluesVec[selectedIndex],4),")"))
        lines(FCysVec + FCxsVec, estimatesVec*10, col="red")
        abline(v=(FCysVec + FCxsVec)[selectedIndex], col="green")
        abline(h=0, col="gray")
        abline(h=-log10(0.05), col="cyan")
      } else {
        if (debug) {
          cat(paste0("Edge ", Xind, "-", Yind,  type, "; est: ",
                     round(resultsList[[1]]$estimate,4), " p=", round(resultsList[[1]]$p.value,4)),
              file=logFile, append=TRUE, sep="\n")
        }
        plot(sum(resultsList[[1]]$FCxs, resultsList[[1]]$FCys, na.rm=TRUE),
           -log10(resultsList[[1]]$p.value),
          col="blue", type="l", ylim=c(-10,15), xlim=c(0,2*maxFC),
          main=paste0("Edge ", Xind, "-", Yind,  type,
          "\nEstimated pCor_t: ", round(curResult$rho_t,4), " (P=", round(curResult$p.value_t,4), ")",
          "\nEstimated sigma2_X|S: ", round(estSigma2_XS,4), ", Estimated sigma2_Y|S: ", round(estSigma2_YS,4),"\nEstimated pCor_", type, " (k=", resultsList[[1]]$k,"): ",
          round(resultsList[[1]]$estimate,4), " (P=", round(resultsList[[1]]$p.value,4),")"))
        lines(sum(resultsList[[1]]$FCxs, resultsList[[1]]$FCys, na.rm=TRUE), resultsList[[1]]$estimate*10, col="red")
        abline(h=0, col="gray")
        abline(h=-log10(0.05), col="cyan")
      }
    }
    dev.off()
}



selectBestResult <- function(resultsList_g, resultsList_e, rho_t, pvalue_t, alpha, maxFC=0.01,
                             debug=TRUE, logFile="./log.txt") {
    if (length(resultsList_g) == 0 || length(resultsList_e) == 0) {
      return(c(-1,-1))
    }

    if (debug) {
    cat(paste0("Selecting best result..."),
        file=logFile, append=TRUE, sep="\n")
    }

    FCxsVec_g <- as.numeric(unlist(lapply(resultsList_g, '[', "FCxs")))
    FCysVec_g <- as.numeric(unlist(lapply(resultsList_g, '[', "FCys")))

    validIds_g <- intersect(which(FCxsVec_g <= maxFC), which(FCysVec_g <= maxFC))
    if (length(validIds_g) == 0 && length(FCxsVec_g) >= 1 && !is.nan(FCxsVec_g[1]))
      validIds_g <- 1

    FCxsVec_g <- FCxsVec_g[validIds_g]
    FCysVec_g <- FCysVec_g[validIds_g]

    pvaluesVec_g <- as.numeric(unlist(lapply(resultsList_g, '[', "p.value")))
    if (length(validIds_g) > 0) {
      pvaluesVec_g <- pvaluesVec_g[validIds_g]
    }
    estimatesVec_g <- as.numeric(unlist(lapply(resultsList_g, '[', "estimate")))
    if (length(validIds_g) > 0) {
      estimatesVec_g <- estimatesVec_g[validIds_g]
    }

    FCxsVec_e <- as.numeric(unlist(lapply(resultsList_e, '[', "FCxs")))
    FCysVec_e <- as.numeric(unlist(lapply(resultsList_e, '[', "FCys")))

    validIds_e <- intersect(which(FCxsVec_e <= maxFC), which(FCysVec_e <= maxFC))
    if (length(validIds_e) == 0 && length(FCxsVec_e) >= 1 && !is.nan(FCxsVec_e[1]))
      validIds_e <- 1

    FCxsVec_e <- FCxsVec_e[validIds_e]
    FCysVec_e <- FCysVec_e[validIds_e]

    pvaluesVec_e <- as.numeric(unlist(lapply(resultsList_e, '[', "p.value")))
    if (length(validIds_e) > 0) {
      pvaluesVec_e <- pvaluesVec_e[validIds_e]
    }

    estimatesVec_e <- as.numeric(unlist(lapply(resultsList_e, '[', "estimate")))
    if (length(validIds_e) > 0) {
      estimatesVec_e <- estimatesVec_e[validIds_e]
    }

    n_ests_g <- length(estimatesVec_g)
    signs_g <- table(sign(estimatesVec_g))
    n_ests_e <- length(estimatesVec_e)
    signs_e <- table(sign(estimatesVec_e))

    selected_index_g = -1
    selected_index_e = -1

    if (pvalue_t > alpha) {
      # If there is a zero total partial correlation
       if ((!is.na(signs_g['-1']) && signs_g['-1'] >= sum(signs_g)/4 &&
              !is.na(signs_e['1']) && signs_e['1'] >= sum(signs_e)/4) ||
           (!is.na(signs_e['-1']) && signs_e['-1'] >= sum(signs_e)/4 &&
              !is.na(signs_g['1']) && signs_g['1'] >= sum(signs_g)/4)) {
        # If all g and e estimates have oppositve signs,
        if (min(pvaluesVec_g) > alpha && min(pvaluesVec_e) > alpha) {
          # Then or both are zero,
          if (debug)
            cat("Total is zero, with g and e zero.",
                file=logFile, append=TRUE, sep="\n")
          selected_index_g <- which(pvaluesVec_g == max(pvaluesVec_g))[1]
          selected_index_e <- which(pvaluesVec_e == max(pvaluesVec_e))[1]
        } else {
          if (debug)
            cat("Total is zero, and g and e have opposite signs.",
                file=logFile, append=TRUE, sep="\n")
          # or both effects are significant, but are cancelled out.
          selected_index_g <- getStableEstimateIdx(FCxsVec_g, FCysVec_g, pvaluesVec_g, estimatesVec_g)
          if (pvaluesVec_g[selected_index_g] > alpha) {
            selected_index_g <- which(pvaluesVec_g == min(pvaluesVec_g))[1]
          }

          selected_index_e <- getStableEstimateIdx(FCxsVec_e, FCysVec_e, pvaluesVec_e, estimatesVec_e)
          if (pvaluesVec_e[selected_index_e] > alpha) {
            selected_index_e <- which(pvaluesVec_e == min(pvaluesVec_e))[1]
          }
        }
      } else {
          # then it is expected that both g and e estimates are zero
        if ((length(pvaluesVec_g) == 1 && pvaluesVec_g[1] > alpha) ||
            (intervalPercentage(FCxsVec_g, FCysVec_g, which(pvaluesVec_g > alpha)) > 0.3)) {
          if (debug)
            cat("Total is zero, g and e have same signal, and g is zero.",
                file=logFile, append=TRUE, sep="\n")
          selected_index_g <- which(pvaluesVec_g == max(pvaluesVec_g))[1]
        } else {
          if (debug)
            cat("WARNING: Total is zero, and g is the most stable estimate.",
                file=logFile, append=TRUE, sep="\n")
          selected_index_g <- getStableEstimateIdx(FCxsVec_g, FCysVec_g, pvaluesVec_g, estimatesVec_g)
        }

        if ((length(pvaluesVec_e) == 1 && pvaluesVec_e[1] > alpha) ||
            (intervalPercentage(FCxsVec_e, FCysVec_e, which(pvaluesVec_e > alpha)) > 0.3)) {
          if (debug)
            cat("Total is zero, g and e have same signal, and e is zero.",
                file=logFile, append=TRUE, sep="\n")
          selected_index_e <- which(pvaluesVec_e == max(pvaluesVec_e))[1]
        } else {
          if (debug)
            cat("WARNING: Total is zero, and e is the most stable estimate.",
                file=logFile, append=TRUE, sep="\n")
          selected_index_e <- getStableEstimateIdx(FCxsVec_e, FCysVec_e, pvaluesVec_e, estimatesVec_e)
        }
      }
    } else {
      # If there is a non-zero total partial correlation
      # Then g or e (or both) should be non-zero
       if ((!is.na(signs_g['-1']) && signs_g['-1'] >= sum(signs_g)/5 &&
              !is.na(signs_e['1']) && signs_e['1'] >= sum(signs_e)/5) ||
           (!is.na(signs_e['-1']) && signs_e['-1'] >= sum(signs_e)/5 &&
              !is.na(signs_g['1']) && signs_g['1'] >= sum(signs_g)/5)) {
        # If all g and e estimates have oppositve signs,
        if (debug)
          cat("Total is nonzero, g and e have opposite signs.",
              file=logFile, append=TRUE, sep="\n")
        selected_index_g <- getStableEstimateIdx(FCxsVec_g, FCysVec_g, pvaluesVec_g, estimatesVec_g)
        if (pvaluesVec_g[selected_index_g] > alpha) {
          #selected_index_g <- which(pvaluesVec_g == min(pvaluesVec_g))[1]
          finalIds <- ceiling(length(estimatesVec_g)/2):length(estimatesVec_g)
          selected_index_g <- finalIds[which.min(pvaluesVec_g[finalIds])][1]
        }
        if (any(signs_g == length(estimatesVec_g)) &&
              pvaluesVec_g[selected_index_g] > alpha) {
           selected_index_g <- which.min(pvaluesVec_g)[1]
        }

        selected_index_e <- getStableEstimateIdx(FCxsVec_e, FCysVec_e, pvaluesVec_e, estimatesVec_e)
        if (pvaluesVec_e[selected_index_e] > alpha) {
          finalIds <- ceiling(length(estimatesVec_e)/2):length(estimatesVec_e)
          selected_index_e <- finalIds[which.min(pvaluesVec_e[finalIds])][1]
        }
        if (any(signs_e == length(estimatesVec_e)) &&
              pvaluesVec_e[selected_index_e] > alpha) {
           selected_index_e <- which.min(pvaluesVec_e)[1]
        }
      } else {
        if ((min(pvaluesVec_g) > alpha)) {
          if (intervalPercentage(FCxsVec_g, FCysVec_g, which(abs(estimatesVec_g) > 0.1)) > 0.20) {
            if (debug)
               cat("Total is nonzero, g zero but estimates are moderate to high",
                   file=logFile, append=TRUE, sep="\n")
            selected_index_g <- which(pvaluesVec_g == min(pvaluesVec_g))[1]
          } else {
          # g is zero: or all p-values are greater than alpha or
            # the g estimates cross the zero and most of them have p-value greater than alpha
            if (debug)
              cat("Total is nonzero, but g is zero.",
                  file=logFile, append=TRUE, sep="\n")
            selected_index_g <- which(pvaluesVec_g == max(pvaluesVec_g))[1]
          }
        }

        if ((min(pvaluesVec_e) > alpha)) {
          if (intervalPercentage(FCxsVec_e, FCysVec_e, which(abs(estimatesVec_e) > 0.1)) > 0.20) {
            if (debug)
               cat("Total is nonzero, e zero but estimates are moderate to high",
                   file=logFile, append=TRUE, sep="\n")
            selected_index_e <- which(pvaluesVec_e == min(pvaluesVec_e))[1]
          } else {
            # e is zero: or all p-values are greater than alpha or
            # the e estimates cross the zero and most of them have p-value greater than alpha
            if (debug)
              cat("Total is nonzero, but e is zero.",
                  file=logFile, append=TRUE, sep="\n")
            selected_index_e <- which(pvaluesVec_e == max(pvaluesVec_e))[1]
          }
        }

        if (selected_index_e == -1 || selected_index_g == -1) {

          if (selected_index_e != -1 && selected_index_g == -1) {
            # Since e is zero, g is probably significant.
            # We will make the selction giving preference to the most stable one.
            candidate_g <- getStableEstimateIdx(FCxsVec_g, FCysVec_g, pvaluesVec_g, estimatesVec_g)

            if (pvaluesVec_g[candidate_g] <= alpha) {
              if (debug)
                cat("Total is nonzero, g is nonzero (the most stable estimate).",
                    file=logFile, append=TRUE, sep="\n")
              selected_index_g <- candidate_g
            } else {
              if (debug)
                cat("Total is nonzero, g may be nonzero (the first significant one).",
                    file=logFile, append=TRUE, sep="\n")

              finalIds <- ceiling(length(estimatesVec_g)/2):length(estimatesVec_g)
              selected_index_g <- finalIds[which.min(pvaluesVec_g[finalIds])][1]
            }
          } else if (selected_index_g != -1 && selected_index_e == -1) {
            # Since g is zero, e is probably significant.
            # We will make the selction giving preference to the most stable one.
           candidate_e <- getStableEstimateIdx(FCxsVec_e, FCysVec_e, pvaluesVec_e, estimatesVec_e)

            if (pvaluesVec_e[candidate_e] <= alpha) {
              if (debug)
                cat("Total is nonzero, e is nonzero (the most stable estimate).",
                    file=logFile, append=TRUE, sep="\n")
              selected_index_e <- candidate_e
            } else {
              if (debug)
                cat("Total is nonzero, e may be nonzero (the first significant one).",
                    file=logFile, append=TRUE, sep="\n")
              finalIds <- ceiling(length(estimatesVec_e)/2):length(estimatesVec_e)
              selected_index_e <- finalIds[which.min(pvaluesVec_e[finalIds])][1]
            }
          } else if (selected_index_g == -1 && selected_index_e == -1) {
            # I can be more strict, because there is not a problematic coufounding,
            # but I have to select at least one componente that is significant.
            candidate_g <- getStableEstimateIdx(FCxsVec_g, FCysVec_g, pvaluesVec_g, estimatesVec_g)
            candidate_e <- getStableEstimateIdx(FCxsVec_e, FCysVec_e, pvaluesVec_e, estimatesVec_e)
            if (pvaluesVec_e[candidate_e] <= alpha ||
                pvaluesVec_g[candidate_g] <= alpha) {

              if (pvaluesVec_e[candidate_e] <= alpha) {
                if (debug)
                  cat("Total is nonzero, e is nonzero (the most stable one).",
                      file=logFile, append=TRUE, sep="\n")
                selected_index_e <- candidate_e
              } else if (length(estimatesVec_e) > 1 &&
                  max(signs_e) == sum(signs_e) &&
                  length(which(abs(estimatesVec_e) > 0.05)) == length(estimatesVec_e)) {
                selected_index_e <- which(pvaluesVec_e == min(pvaluesVec_e))[1]
              } else {
                if (debug)
                  cat("Total is nonzero, e is zero.",
                      file=logFile, append=TRUE, sep="\n")
                selected_index_e <- which(pvaluesVec_e == max(pvaluesVec_e))[1]
              }

               if (pvaluesVec_g[candidate_g] <= alpha) {
                if (debug)
                  cat("Total is nonzero, g is nonzero (the most stable one).",
                      file=logFile, append=TRUE, sep="\n")
                selected_index_g <- candidate_g
               } else if (length(estimatesVec_g) > 1 &&
                  max(signs_g) == sum(signs_g) &&
                  length(which(abs(estimatesVec_g) > 0.05)) == length(estimatesVec_g)) {
                selected_index_g <- which(pvaluesVec_g == min(pvaluesVec_g))[1]
               } else {
                if (debug)
                  cat("Total is nonzero, g is zero.",
                      file=logFile, append=TRUE, sep="\n")
                selected_index_g <- which(pvaluesVec_g == max(pvaluesVec_g))[1]
               }
            } else {
              # Then g and e have significant effect for some k (they were not set as zero before)
              # but they are not the most stable ones.
              if ((intervalPercentage(FCxsVec_e, FCysVec_e,
                  which(pvaluesVec_e <= alpha)) > 0.20) ||
                  (intervalPercentage(FCxsVec_g, FCysVec_g,
                  which(pvaluesVec_g <= alpha)) > 0.20)) {
                if (intervalPercentage(FCxsVec_e, FCysVec_e,
                  which(pvaluesVec_e <= alpha)) > 0.20) {
                  selected_index_e <- which(pvaluesVec_e == min(pvaluesVec_e))[1]
                } else {
                  selected_index_e <- which(pvaluesVec_e == max(pvaluesVec_e))[1]
                }

                if (intervalPercentage(FCxsVec_g, FCysVec_g,
                  which(pvaluesVec_g <= alpha)) > 0.20) {
                  selected_index_g <- which(pvaluesVec_g == min(pvaluesVec_g))[1]
                } else {
                  selected_index_g <- which(pvaluesVec_g == max(pvaluesVec_g))[1]
                }
              } else {
                if (debug)
                  cat("Total is nonzero, g and e are nonzero (the first significant ones).",
                      file=logFile, append=TRUE, sep="\n")
                selected_index_e <- which(pvaluesVec_e == min(pvaluesVec_e))[1]
                selected_index_g <- which(pvaluesVec_g == min(pvaluesVec_g))[1]
              }
            }
          }
        }
      }
    }

    if (selected_index_g == -1) {
      if (debug) {
        cat("WARNING: unexpected case for g!",
            file=logFile, append=TRUE, sep="\n")
      }
      selected_index_g <- getStableEstimateIdx(FCxsVec_g, FCysVec_g, pvaluesVec_g, estimatesVec_g)
    }

    if (selected_index_e == -1) {
      if (debug) {
        cat("WARNING: unexpected case for e!",
            file=logFile, append=TRUE, sep="\n")
      }
      selected_index_e <- getStableEstimateIdx(FCxsVec_e, FCysVec_e, pvaluesVec_e, estimatesVec_e)
    }

    return(c(selected_index_g, selected_index_e))
}


intervalPercentage <- function(FCxsVec, FCysVec, indices) {
  if (length(indices) == 0)
    return(0)

  sumFCs <- FCxsVec + FCysVec
  sequences <- c()
  newseq = TRUE
  for (i in 1:length(indices)) {
    if (newseq) {
      if (indices[i] > 1)
        sequences <- c(sequences, (sumFCs[indices[i]]+sumFCs[indices[i]-1])/2)
      else
        sequences <- c(sequences, sumFCs[indices[i]])

      newseq = FALSE
    }

    if ((i < length(indices) && indices[i+1] != indices[i]+1) ||
        (i == length(indices))) {
       if (indices[i] < length(FCxsVec))
        sequences <- c(sequences, (sumFCs[indices[i]]+sumFCs[indices[i]+1])/2)
      else
        sequences <- c(sequences, sumFCs[indices[i]])
      newseq = TRUE
    }
  }

  totalInterval <- 0
  i = 2
  while (i <= length(sequences)) {
    totalInterval <- totalInterval + (sequences[i] - sequences[i-1])
    i = i+2
  }
  sumFCsInterval <- sumFCs[length(sumFCs)] - sumFCs[1]
  if (sumFCsInterval == 0)
    return(0)

  return(totalInterval/sumFCsInterval)
}


# both lists must be sorted!
getNClosestPoints <- function(pointList, allPoints) {
  closestPoints <- c()
  i = 1
  for (valueId in 1:length(allPoints)) {
    if (i > length(pointList))
      break
    if (allPoints[valueId] == pointList[i]) {
      closestPoints <- c(closestPoints, allPoints[valueId])
      i = i+1
    } else if (allPoints[valueId] > pointList[i]) {
      if (abs(allPoints[valueId] - pointList[i]) < abs(allPoints[valueId-1] - pointList[i])) {
        closestPoints <- c(closestPoints, allPoints[valueId])
      } else {
        closestPoints <- c(closestPoints, allPoints[valueId-1])
      }
      i = i+1
    }
  }
  closestPoints <- c(closestPoints, allPoints[valueId])
  return(closestPoints)
}

#' @importFrom stats var
getStableEstimateIdx <- function(FCxsVec, FCysVec, pvaluesVec, estimatesVec) {
  if (length(estimatesVec) == 1) {
    return(1)
  }

  limiters <- seq((FCysVec + FCxsVec)[1], (FCysVec + FCxsVec)[length(FCysVec + FCxsVec)], length.out=6)
  limiters <- getNClosestPoints(limiters, (FCysVec + FCxsVec))
  limitersIds <- which((FCysVec + FCxsVec) %in% limiters)

  varsVec <- c()
  for (i in 1:(length(limitersIds)-1)) {
    varsVec <- c(varsVec, var(estimatesVec[limitersIds[i]:limitersIds[i+1]]))
  }

  stableIds <- which(round(varsVec,4) <= round(summary(varsVec)[4],4))
  finalStableIds <- c()
  for (id in stableIds) {
    finalStableIds <- c(finalStableIds, limitersIds[id]:limitersIds[id+1])
  }
  finalStableIds <- unique(finalStableIds)

  selected_index <- which(pvaluesVec[finalStableIds] == min(pvaluesVec[finalStableIds]))[1]
  selected_index <- finalStableIds[selected_index]

  return(selected_index)
}


#' @importFrom stats p.adjust
adjustSymmetricPValuesMatrix <- function(pvaluesMatrix, correction="fdr") {
  lowind <- which(lower.tri(pvaluesMatrix), arr.ind=TRUE)
  uppind <- which(upper.tri(pvaluesMatrix), arr.ind=TRUE)

  adjp <- p.adjust(pvaluesMatrix[lowind], method=correction)
  pvaluesMatrix[lowind] <- adjp
  pvaluesMatrix[uppind] <- adjp

  return(pvaluesMatrix)
}

getAdjMatrixFromSymmetricPCorPvalues <- function(pvaluesMatrix, alpha=0.05, correction=NULL) {
  nVertices <- dim(pvaluesMatrix)[1]

  if (!is.null(correction)) {
    pvaluesMatrix <- adjustSymmetricPValuesMatrix(pvaluesMatrix, correction)
  }

  edges <- which(pvaluesMatrix <= alpha, arr.ind=TRUE)
  adjMatrix <- matrix(0, nVertices, nVertices)
  adjMatrix[edges] <- 1
  colnames(adjMatrix) <- colnames(pvaluesMatrix)
  row.names(adjMatrix) <- colnames(pvaluesMatrix)

  return(list(adjacency.matrix=adjMatrix, adjusted.pvalues=pvaluesMatrix))
}

getAdjMatrixFromSymmetricPCor <- function(pCor, alpha, correction="fdr") {
  pvaluesMatrix <- pCor$pvalues
  diag(pvaluesMatrix) <- 1

  return(getAdjMatrixFromSymmetricPCorPvalues(pvaluesMatrix, alpha=0.05, correction))
}

#' @importFrom kinship2 kinship
getKinshipList <- function(pedigrees, sampled=NULL, debug=TRUE, logFile=NULL) {
  if (!is.null(sampled) && length(sampled) != dim(pedigrees)[1]) {
    cat("Sampled vector size and number of rows of pedigrees data frame must be equal",
        file="log_error.log", append=TRUE, sep="\n")
    return(NULL)
  }

  if (!is.null(sampled))
    pedigrees <- cbind(pedigrees, sampled=sampled)

  famids = unique(pedigrees$famid)

  Phi_list <- list()
  if (!is.null(sampled))
    sampled_Phi_list <- list()

  for (famid in famids) {
    ped <- pedigrees[pedigrees$famid == famid,]
    Phi_list[[famid]] <- kinship(id=ped$id, dadid=ped$dadid, momid=ped$momid, sex=ped$sex)

    if (!is.null(sampled)) {
      sampled_ids <- which(ped$sampled == 1)
      sampled_Phi_list[[famid]] <- Phi_list[[famid]][sampled_ids, sampled_ids]
    }
  }

  if (!is.null(sampled))
    return(list(Phi_list=Phi_list, sampled_Phi_list=sampled_Phi_list))
  else
    return(Phi_list)
}


matrixSqrt <- function(m, inverse=FALSE, tol = sqrt(.Machine$double.eps)) {
  eigDecomp <- eigen(m, symmetric=TRUE)
  Positive <- eigDecomp$values > max(tol * eigDecomp$values[1], 0)
  D = NULL
  if (inverse) {
    D = diag(1/sqrt(eigDecomp$values[Positive]))
  } else {
    D = diag(sqrt(eigDecomp$values[Positive]))
  }
  sqrtM <- eigDecomp$vectors[, Positive, drop = FALSE] %*% D %*% t(eigDecomp$vectors[, Positive, drop = FALSE])
  return(sqrtM)
}

#' @importFrom utils combn
#' @importFrom parallel detectCores
#' @import foreach
#' @import doMC
preprocessFamilyBasedUDGs <- function(phen.df, covs.df, pedigrees, sampled, fileID, dirToSave, max_cores=NULL,
                                      minK=10, maxFC = 0.01, orthogonal=TRUE, useGPU=FALSE,
                                      debug=TRUE, logFile=NULL) {

  if (debug && is.null(logFile)){
    logFile = paste0(dirToSave, "preprocessFamilyBasedUDGs_", format(Sys.time(), "%Y%m%d-%Hh%Mm%Ss"), ".log")
  }

  pseq <- 1:ncol(phen.df)
  N <- dim(phen.df)[1]
  xyinds <- combn(pseq,2)

  Phi2 <- calculateBlockDiagonal2PhiMatrix(pedigrees, FALSE, sampled=sampled)
  Z <- calculateBlockDiagonal2PhiMatrix(pedigrees, TRUE, sampled=sampled)

  sampledPed <- pedigrees
  if(!is.null(sampled)) {
    sampledPed <- sampledPed[which(sampled==1),]
  }

  joblabels <- c("x", "y", "S")
  jobqueue <- matrix(NA, 0,length(joblabels))
  colnames(jobqueue) <- joblabels

  for (xyind in 1:dim(xyinds)[2]) {
    Xind <- xyinds[1,xyind] # Xind is always less than Yind
    Yind <- xyinds[2,xyind]
    Sinds <- pseq[-c(Xind, Yind)]
    jobqueue <- rbind(jobqueue, c(Xind,Yind,paste(Sinds, collapse=" ")))
  }

  if (dim(jobqueue)[1] > 0) {
    jobqueue = as.data.frame(jobqueue, stringsAsFactors=FALSE)
    jobqueue$x <- as.numeric(jobqueue$x)
    jobqueue$y <- as.numeric(jobqueue$y)

    fitqueue = as.data.frame(rbind(cbind(rownames(jobqueue), jobqueue[,"x"], jobqueue[,"S"]), cbind(rownames(jobqueue), jobqueue[,"y"], jobqueue[,"S"])), stringsAsFactors=FALSE)
    fitqueue[,1] <- as.numeric(fitqueue[,1])
    fitqueue[,2] <- as.numeric(fitqueue[,2])

    fitqueue = fitqueue[order(fitqueue[,1]),]
    notDuplJobs <- !duplicated(fitqueue[,c(2,3)])
    xNotDuplJobs <- notDuplJobs[which(1:length(notDuplJobs) %% 2 == 1)]
    yNotDuplJobs <- notDuplJobs[which(1:length(notDuplJobs) %% 2 == 0)]

    if (is.null(max_cores)) {
      max_cores <- detectCores()
    }
    registerDoMC(max_cores)

    batches <- list()
    i = 1
    notJobs <- intersect(which(xNotDuplJobs), which(yNotDuplJobs))
    notBatches <- split(notJobs, ceiling(seq_along(notJobs)/max_cores))
    for (batch in 1:length(notBatches)) {
      batches[i] <- notBatches[batch]
      i = i + 1
    }

    xorJobs <- sort(union(intersect(which(xNotDuplJobs == FALSE), which(yNotDuplJobs)), intersect(which(yNotDuplJobs == FALSE), which(xNotDuplJobs))))
    if (length(xorJobs) > 0) {
      xorBatches <- split(xorJobs, ceiling(seq_along(xorJobs)/max_cores))
      for (batch in 1:length(xorBatches)) {
        batches[i] <- xorBatches[batch]
        i = i + 1
      }
    }

    trueJobs <- intersect(which(xNotDuplJobs == FALSE), which(yNotDuplJobs == FALSE))
    if (length(trueJobs) > 0) {
      trueBatches <- split(trueJobs, ceiling(seq_along(trueJobs)/max_cores))
      for (batch in 1:length(trueBatches)) {
        batches[i] <- trueBatches[batch]
        i = i + 1
      }
    }

    # process in parallel all jobs in jobqueue and append the output in results.
    for (batchId in 1:length(batches)) {
      jobi = 0 # a hack to avoid problems when checking --as-cran

      foreach(jobi = batches[[batchId]]) %dopar% {
        x <- as.numeric(jobqueue[jobi, "x"])
        y <- as.numeric(jobqueue[jobi, "y"])
        S <- as.numeric(unlist(strsplit(as.character(jobqueue[jobi, "S"]), " ")))
        if (debug) {
          cat(paste0(Sys.getpid(),
                     " - Checking computation for ", x, ".", y, "|{",paste(S, collapse=","), "}!"),
              file=logFile, append=TRUE, sep="\n")
        }

        processFile <- TRUE
        fileToSave <- paste0(dirToSave, fileID, "_", x, "_", y, "_", paste(S, collapse=""), ".rds")
        if (file.exists(fileToSave)) {
          if (debug) {
            cat(paste0(Sys.getpid(),
                       " - file already exists for ", x, ".", y, "|{",paste(S, collapse=","), "}..."),
                file=logFile, append=TRUE, sep="\n")
          }
          curResult <- readRDS(fileToSave)
          resultsList_g <- curResult$results_g
          resultsList_e <- curResult$results_e

          ret_g <- checkUnconfoundedEstimatesList(resultsList_g, maxFC, N)
          ret_e <- checkUnconfoundedEstimatesList(resultsList_e, maxFC, N)

          if (ret_g == "error" || ret_e == "error") {
            if (debug) {
              cat(paste0("Problem in UnconfoundedEstimatesList - removing file: ", fileToSave, "\n"),
                  file=logFile, append=TRUE, sep="\n")
              cat(paste0(Sys.getpid(),
                         " - redoing computation for ", x, ".", y, "|{",paste(S,   collapse=","), "}..."),
                  file=logFile, append=TRUE, sep="\n")
            }
            file.remove(fileToSave)
          } else if (ret_g == "done" && ret_e == "done") {
            if (debug) {
              cat(paste0(Sys.getpid(),
                         " - computation for ", x, ".", y, "|{",paste(S, collapse=","), "} is done!"),
                  file=logFile, append=TRUE, sep="\n")
            }
            processFile <- FALSE
          } else {
            if (debug) {
              cat(paste0(Sys.getpid(),
                         " - updating computation for ", x, ".", y, "|{",paste(S, collapse=","), "}..."),
                  file=logFile, append=TRUE, sep="\n")
            }
          }
        }

        if (processFile) {
          if (debug) {
            cat(paste0(Sys.getpid(),
                       " - doing computation for ", x, ".", y, "|{",paste(S, collapse=","), "}..."),
                file=logFile, append=TRUE, sep="\n")
          }
          curResult <- getFamilyBasedPCorResults(x, y, S, phen.df, covs.df, sampledPed, Phi2, Z, minK,
                                                 maxFC, orthogonal, fileID, dirToSave, useGPU, debug, logFile)
          if (!is.null(fileID) && !is.null(dirToSave)) {
            saveRDS(curResult, file=fileToSave)
          }
        }
      }
    }
  }
}

#' @importFrom utils combn write.table
#' @import foreach
preprocessFamilyBasedDAGs <- function(phen.df, covs.df, pedigrees, sampled, alpha,
    fileID, dirToSave, preprocessedPvaluesCSVFile, max_cores=NULL, minK=10, maxFC = 0.01,
    orthogonal=TRUE, savePlots=FALSE, useGPU=FALSE, debug=TRUE, logFile=NULL) {

  if (debug) {
    cat("Preprocessing pvalues for PC algorithm...",
        file=logFile, append=TRUE, sep="\n")
  }

  Phi2 <- calculateBlockDiagonal2PhiMatrix(pedigrees, FALSE, sampled=sampled)
  Z <- calculateBlockDiagonal2PhiMatrix(pedigrees, TRUE, sampled=sampled)

  sampledPed <- pedigrees
  if(!is.null(sampled)) {
    sampledPed <- sampledPed[which(sampled==1),]
  }

  csvFileName <- preprocessedPvaluesCSVFile

  # preprocessing using parallel computing
  p <- dim(phen.df)[2]
  pseq <- 1:p
  skel_t <- matrix(TRUE, p, p)
  diag(skel_t) <- FALSE
  skel_g <- skel_t
  skel_e <- skel_t

  labels <- c("x", "y", "S","pcort", "pcort_p.value",
        "pcorg", "pcorg_p.value",
        "pcore", "pcore_p.value",
        "se", "sg")

  if (!file.exists(csvFileName)) {
    write.table(as.matrix(t(labels)), file = csvFileName, sep = ";", row.names=FALSE, col.names = FALSE, append=FALSE, quote=FALSE)
  }

  orders <- 0:(p-2)

  for (iorder in orders) {
    for (pairType in c("upper", "lower")) {
      allpairs_g <- which(skel_g == TRUE, arr.ind=TRUE)
      allpairs_e <- which(skel_e == TRUE, arr.ind=TRUE)
      allpairs_t <- which(skel_t == TRUE, arr.ind=TRUE)

      allpairs <- unique(rbind(allpairs_g, allpairs_e, allpairs_t))

      remainingEdges <- dim(allpairs)[1]
      if (debug) {
        cat("Remaining", remainingEdges, "edges",
            file=logFile, append=TRUE, sep="\n")
      }

      if (remainingEdges > 0) {

        if(pairType == "upper") {
          pairs <- allpairs[which(allpairs[,1] < allpairs[,2]),,drop=FALSE]
          pairs <- pairs[order(pairs[,1]),,drop=FALSE]
        } else {
          pairs <- allpairs[which(allpairs[,1] > allpairs[,2]),,drop=FALSE]
          pairs <- pairs[order(pairs[,1]),,drop=FALSE]
        }

        joblabels <- c("x", "y", "S", "dagt", "dagg", "dage")
        jobqueue <- matrix(NA, 0,length(joblabels))
        colnames(jobqueue) <- joblabels

        results <- c()
        for (pairind in 1:dim(pairs)[1]) {
          pair <- pairs[pairind,]
          x <- as.integer(pair[1])
          y <- as.integer(pair[2])
          others <- pseq[-pair]

          dagg=sum(duplicated(rbind(allpairs_g, pair)))
          dage=sum(duplicated(rbind(allpairs_e, pair)))
          dagt=sum(duplicated(rbind(allpairs_t, pair)))

          if (iorder > 0 && iorder < (p-2)) {
            if (debug) {
              cat(paste0("Processing iorder ", iorder),
                  file=logFile, append=TRUE, sep="\n")
            }
            Ssets <- combn(others, iorder)
            for(Sind in 1:dim(Ssets)[2]) {
              S = Ssets[,Sind]
              pcorinfo <- getPCorInfo(csvFileName,x,y,S,dagg,dage,dagt, debug=debug);
              iscompl <- isPcorInfoComplete(pcorinfo, dagt, dagg, dage);
              if (length(pcorinfo) > 0 && sum(iscompl) == 3) {
                results <- rbind(results, pcorinfo)
              } else {
                jobqueue <- rbind(jobqueue, c(x,y,paste(S, collapse=" "), !iscompl))
              }
            }
          } else if(iorder == 0) {
            S = c()
            if (debug) {
              cat(paste0("Processing iorder ", iorder),
                  file=logFile, append=TRUE, sep="\n")
              cat(paste0("Checking ", x, ".", y, "|{",paste(S, collapse=","), "}..."),
                  file=logFile, append=TRUE, sep="\n")
            }
            pcorinfo <- getPCorInfo(csvFileName,x,y,S,dagg,dage,dagt, debug=debug);
            iscompl <- isPcorInfoComplete(pcorinfo, dagt, dagg, dage);
            if (length(pcorinfo) > 0 && sum(iscompl) == 3) {
              results <- rbind(results, pcorinfo)
            } else {
              jobqueue <- rbind(jobqueue, c(x,y,"", !iscompl))
            }
          } else if (iorder == (p-2)) {
            if (debug) {
              cat(paste0("Processing iorder ", iorder),
                  file=logFile, append=TRUE, sep="\n")
            }
            S = others
            pcorinfo <- getPCorInfo(csvFileName,x,y,S,dagg,dage,dagt, debug=debug);
            iscompl <- isPcorInfoComplete(pcorinfo, dagt, dagg, dage);
            if (length(pcorinfo) > 0 && sum(iscompl) == 3) {
              results <- rbind(results, pcorinfo)
            } else {
              jobqueue <- rbind(jobqueue, c(x,y,paste(S, collapse=" "), !iscompl))
            }
          }
        }

        if (dim(jobqueue)[1] > 0) {
          jobqueue = as.data.frame(jobqueue, stringsAsFactors=FALSE)
          jobqueue$x <- as.numeric(jobqueue$x)
          jobqueue$y <- as.numeric(jobqueue$y)

          fitqueue = as.data.frame(rbind(cbind(rownames(jobqueue), jobqueue[,"x"], jobqueue[,"S"]), cbind(rownames(jobqueue), jobqueue[,"y"], jobqueue[,"S"])), stringsAsFactors=FALSE)
          fitqueue[,1] <- as.numeric(fitqueue[,1])
          fitqueue[,2] <- as.numeric(fitqueue[,2])

          fitqueue = fitqueue[order(fitqueue[,1]),]
          notDuplJobs <- !duplicated(fitqueue[,c(2,3)])
          xNotDuplJobs <- notDuplJobs[which(1:length(notDuplJobs) %% 2 == 1)]
          yNotDuplJobs <- notDuplJobs[which(1:length(notDuplJobs) %% 2 == 0)]


          if (is.null(max_cores)) {
            max_cores <- detectCores()
          }
          registerDoMC(max_cores)

          batches <- list()
          i = 1
          notJobs <- intersect(which(xNotDuplJobs), which(yNotDuplJobs))
          notBatches <- split(notJobs, ceiling(seq_along(notJobs)/max_cores))
          for (batch in 1:length(notBatches)) {
            batches[i] <- notBatches[batch]
            i = i + 1
          }

          xorJobs <- sort(union(intersect(which(xNotDuplJobs == FALSE), which(yNotDuplJobs)), intersect(which(yNotDuplJobs == FALSE), which(xNotDuplJobs))))
          if (length(xorJobs) > 0) {
            xorBatches <- split(xorJobs, ceiling(seq_along(xorJobs)/max_cores))
            for (batch in 1:length(xorBatches)) {
              batches[i] <- xorBatches[batch]
              i = i + 1
            }
          }

          trueJobs <- intersect(which(xNotDuplJobs == FALSE), which(yNotDuplJobs == FALSE))
          if (length(trueJobs) > 0) {
            trueBatches <- split(trueJobs, ceiling(seq_along(trueJobs)/max_cores))
            for (batch in 1:length(trueBatches)) {
              batches[i] <- trueBatches[batch]
              i = i + 1
            }
          }

          # process in parallel all jobs in jobqueue and append the output in results.
          for (batchId in 1:length(batches)) {
            curResults <- c()
            jobi = 0 # a hack to avoid problems when checking --as-cran
            curResults <- foreach(jobi = batches[[batchId]],
                                  .combine='rbind') %dopar% {
              #for(jobi in batches[[batchId]]) {
              x <- as.numeric(jobqueue[jobi, "x"])
              y <- as.numeric(jobqueue[jobi, "y"])
              S <- as.numeric(unlist(strsplit(as.character(jobqueue[jobi, "S"]), " ")))
              dagt <- as.logical(jobqueue[jobi, "dagt"])
              dagg <- as.logical(jobqueue[jobi, "dagg"])
              dage <-as.logical(jobqueue[jobi, "dage"])

              suffStat <- list()
              suffStat$preprocessedPvaluesCSVFile <- csvFileName
              suffStat$savePlots <- savePlots
              suffStat$dagt <- as.logical(dagt)
              suffStat$dagg <- as.logical(dagg)
              suffStat$dage <- as.logical(dage)

              if (debug) {
                cat(paste0(Sys.getpid(), " - doing computation for ", x, ".", y, "|{",paste(S, collapse=","), "}..."),
                    file=logFile, append=TRUE, sep="\n")
              }

              ret <- familyBasedPCorTest(x, y, S, phen.df, covs.df, sampledPed, Phi2, Z,
                minK=minK, maxFC=maxFC, orthogonal=orthogonal, alpha=alpha,
                dirToSave, fileID, suffStat$savePlots, useGPU=useGPU, debug=debug, logFile=logFile)

              if (debug) {
                cat(paste0(Sys.getpid(), " - computation for ", x, ".", y, "|{",paste(S, collapse=","), "} is done:\n",
                           paste(ret, collapse=";")),
                    file=logFile, append=TRUE, sep="\n")
              }
              currow <- c(x=x,y=y,S=paste(S, collapse=" "), getRetVec(ret))
              currow
            }
            if (!is.null(curResults)) {
              if (length(curResults) == 11) {
                # converting the vector to a matrix
                curResultsColNames <- names(curResults)
                curResults <- matrix(curResults, ncol=11)
                colnames(curResults) <- curResultsColNames
              }
              updatePreprocessedPvaluesTable(curResults, csvFileName)

              results <- rbind(results, curResults)
            }
          }
        }

        if (is.null(results)) {
          if (length(results) > 0)
            colnames(results) <- labels

          # updating skel matrices
          rowids_g <- as.numeric(which(as.numeric(results[,"pcorg_p.value"]) > alpha))
          skel_g[cbind(as.numeric(results[rowids_g, "x"]), as.numeric(results[rowids_g, "y"]))] <- FALSE
          skel_g[lower.tri(skel_g)] <- (t(skel_g))[lower.tri(skel_g)]

          rowids_e <- as.numeric(which(as.numeric(results[,"pcore_p.value"]) > alpha))
          skel_e[cbind(as.numeric(results[rowids_e, "x"]), as.numeric(results[rowids_e, "y"]))] <- FALSE
          skel_e[lower.tri(skel_e)] <- (t(skel_e))[lower.tri(skel_e)]

          rowids_t <- as.numeric(which(as.numeric(results[,"pcort_p.value"]) > alpha))
          skel_t[cbind(as.numeric(results[rowids_t, "x"]), as.numeric(results[rowids_t, "y"]))] <- FALSE
          skel_t[lower.tri(skel_t)] <- (t(skel_t))[lower.tri(skel_t)]
        }
      }
    }
  }
}

# Function to compute the partial correlation coefficient
# It must be independent of the PC algorithm or any other approach that uses suffStat.
# pedigrees must be sampled
familyBasedPCorTest <- function(x, y, S, phen.df, covs.df, pedigrees, Phi2, Z,
        minK=10, maxFC = 0.01, orthogonal=TRUE, alpha=0.05,
        dirToSave=NULL, fileID=NULL, savePlots=FALSE, useGPU=FALSE,
        debug=TRUE, logFile="./log.txt") {

  if (debug) {
    cat(paste0("Running familyBasedPCorTest for ", fileID, "; useGPU=", useGPU),
        file=logFile, append=TRUE, sep="\n")
  }

  fileToSave <- NULL
  if (!is.null(fileID) && !is.null(dirToSave)) {
    fileToSave <- paste0(dirToSave, fileID, "_", x, "_", y, "_", paste(S, collapse=""), ".rds")
    if (x > y) {
      fileToSave <- paste0(dirToSave, fileID, "_", y, "_", x, "_", paste(S, collapse=""), ".rds")
    }
  }

  # NOTE: We have to compute this every time to check if results are complete.
  curResult <- getFamilyBasedPCorResults(x, y, S, phen.df, covs.df, pedigrees, Phi2, Z, minK, maxFC,
                                         orthogonal, fileID, dirToSave, useGPU, debug, logFile)
  results <- selectFamilyBasedPCor(curResult, alpha, maxFC, curResultFileName=fileToSave,
                                   debug=debug, logFile=logFile)

  if (savePlots && !is.null(fileID) && !is.null(dirToSave)) {
      if (debug) {
        cat("Generating result plot...",
            file=logFile, append=TRUE, sep="\n")
      }
      saveEdgeResultsPlot(curResult, dirToSave, fileID, x, y, S, alpha, maxFC, debug=debug, logFile=logFile)
  }

  return(results)
}


# curResult is the output of the getFamilyBasedPCorResults function
selectFamilyBasedPCor <- function(curResult, alpha, maxFC, curResultFileName="", debug=TRUE, logFile="./log.txt") {

  if (debug) {
    cat(paste0("Selecting from results in ", curResultFileName),
        file=logFile, append=TRUE, sep="\n")
  }
  resultsList_g <- curResult$results_g
  resultsList_e <- curResult$results_e

  selectedInds <- selectBestResult(curResult$results_g, curResult$results_e,
    curResult$rho_t, curResult$p.value_t, alpha, maxFC, debug=debug, logFile=logFile)

  pcorg <- resultsList_g[[selectedInds[1]]]$estimate
  pcorg_p.value <- resultsList_g[[selectedInds[1]]]$p.value
  sg <- resultsList_g[[selectedInds[1]]]$k

  pcore <- resultsList_e[[selectedInds[2]]]$estimate
  pcore_p.value <- resultsList_e[[selectedInds[2]]]$p.value
  se <- resultsList_e[[selectedInds[2]]]$k

  ret <- list(pcort=as.numeric(curResult$rho_t), pcort_p.value=curResult$p.value_t,
              pcorg=as.numeric(pcorg), pcorg_p.value=pcorg_p.value,
              pcore=as.numeric(pcore), pcore_p.value=pcore_p.value,
              se=se, sg=sg)

  return(unlist(ret))
}


isPcorInfoComplete <- function(pcorinfo, dagt, dagg, dage) {
  isDAGInfoComplete <- function(dag, pvalue) {
    if (!dag || (dag && !is.na(as.numeric(pvalue)))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  return(c(
    isDAGInfoComplete(dagt, pcorinfo["pcort_p.value"]),
    isDAGInfoComplete(dagg, pcorinfo["pcorg_p.value"]),
    isDAGInfoComplete(dage, pcorinfo["pcore_p.value"])))
}


#' @importFrom utils read.table write.table
getPCorInfo <- function(csvfile, x, y, S, dagg, dage, dagt,
                        debug=TRUE, logFile="./log.txt") {
  retvalue = FALSE;

  if (file.exists(csvfile)) {
    if (length(S) == 0)
      S <- ""

    preprocessedPvalues <- read.table(file=csvfile, sep=";", header=TRUE)
    narows <- which(is.na(preprocessedPvalues$S))
    if (length(narows) > 0)
      preprocessedPvalues[narows,"S"] <- ""

    Svalue <- paste(S, collapse=" ")

    xvalues <- which(preprocessedPvalues[,"x"] == x)
    yvalues <- which(preprocessedPvalues[,"y"] == y)
    Svalues <- which(preprocessedPvalues[,"S"] == Svalue)
    ind <- intersect(intersect(xvalues,yvalues), Svalues)

    if (length(ind) > 1) {
      if (debug) {
        cat(paste0("WARNING in pcorExists: repeated entries for ", x, ".", y, "|{",paste(S, collapse=","), "}"),
            file=logFile, append=TRUE, sep="\n")
      }
      maxrow <- apply(preprocessedPvalues[ind,], 2, FUN=function(x){max(x,na.rm=TRUE)})
      preprocessedPvalues[ind[1],] <- maxrow
      preprocessedPvalues <- preprocessedPvalues[-ind[-1],]
      write.table(preprocessedPvalues, file = csvfile, sep = ";",
                  row.names=FALSE, col.names = TRUE, append=FALSE, quote=FALSE)
      ind = ind[1]
    }

    return(preprocessedPvalues[ind,])
  } else {
    if (debug) {
      cat("WARNING!!! Empty return",
          file=logFile, append=TRUE, sep="\n")
    }
    return(c())
  }
}



getRetVec <- function(familyBasedPCorTestRet) {
    return(c(
    familyBasedPCorTestRet["pcort"],
    familyBasedPCorTestRet["pcort_p.value"],
    familyBasedPCorTestRet["pcorg"],
    familyBasedPCorTestRet["pcorg_p.value"],
    familyBasedPCorTestRet["pcore"],
    familyBasedPCorTestRet["pcore_p.value"],
    familyBasedPCorTestRet["se"],
    familyBasedPCorTestRet["sg"]))
}

#' @importFrom utils write.table
updatePreprocessedPvaluesTable <- function(curResults, csvFileName) {
  if (length(curResults) > 0 && nrow(curResults) > 0) {
    labels <- c("x", "y", "S","pcort", "pcort_p.value",
      "pcorg", "pcorg_p.value",
      "pcore", "pcore_p.value",
      "se", "sg")

    if (!file.exists(csvFileName)) {
      write.table(as.matrix(t(labels)), file = csvFileName, sep = ";", row.names=FALSE, col.names = FALSE, append=FALSE, quote=FALSE)
    }

    write.table(curResults, file = csvFileName, sep = ";", row.names=FALSE,
      col.names = FALSE, append=TRUE, quote=FALSE)
  }
}

#' @importFrom utils write.table
#' @importFrom pcalg pc rfci
generateFamilyDAG <- function(phen.df, covs.df, pedigrees, sampled=NULL,
    preprocessedPvaluesCSVFile, type, alpha,
    hidden_vars = FALSE, maj.rule=FALSE,
    minK=10, maxFC = 0.01, orthogonal=TRUE,
    dirToSave, fileID, savePlots, useGPU=FALSE, debug=TRUE, logFile="./log.txt") {
  labels <- c("x", "y", "S","pcort", "pcort_p.value",
        "pcorg", "pcorg_p.value",
        "pcore", "pcore_p.value",
        "se", "sg")

  # NOTE: cannot run this in parallel, since problems occur
  # when the csvfile is updated by multiple processes.

  dagg <- dage <- dagt <- FALSE
  usedPvaluesCSVFile_g <- usedPvaluesCSVFile_e <- usedPvaluesCSVFile_t <- NULL
  if (type == "g") {
    dagg <- TRUE
    usedPvaluesCSVFile_g <- paste0(dirToSave, fileID, "_usedPvalues_g.csv")
    write.table(as.matrix(t(labels)), file = usedPvaluesCSVFile_g, sep = ";", row.names=FALSE, col.names = FALSE, append=FALSE, quote=FALSE)
  } else if (type == "e") {
    dage <- TRUE
    usedPvaluesCSVFile_e <- paste0(dirToSave, fileID, "_usedPvalues_e.csv")
    write.table(as.matrix(t(labels)), file = usedPvaluesCSVFile_e, sep = ";", row.names=FALSE, col.names = FALSE, append=FALSE, quote=FALSE)
  } else {
    dagt <- TRUE
    usedPvaluesCSVFile_t <- paste0(dirToSave, fileID, "_usedPvalues_t.csv")
    write.table(as.matrix(t(labels)), file = usedPvaluesCSVFile_t, sep = ";", row.names=FALSE, col.names = FALSE, append=FALSE, quote=FALSE)
  }

  if (debug) {
    cat(paste0("Generating DAG_", type, " for ", fileID, "; useGPU=", useGPU),
        file=logFile, append=TRUE, sep="\n")
  }

  Phi2 <- calculateBlockDiagonal2PhiMatrix(pedigrees, FALSE, sampled=sampled)
  Z <- calculateBlockDiagonal2PhiMatrix(pedigrees, TRUE, sampled=sampled)

  sampledPed <- pedigrees
  if(!is.null(sampled)) {
    sampledPed <- sampledPed[which(sampled==1),]
  }

  suffStat <- list(
    phen.df=phen.df, covs.df=covs.df, pedigrees=sampledPed,
    minK=minK, maxFC=maxFC, orthogonal=orthogonal,
    alpha=alpha,
    Phi2=Phi2, Z=Z,
    dagg=dagg, dage=dage, dagt=dagt,
    maj.rule = maj.rule,
    preprocessedPvaluesCSVFile=preprocessedPvaluesCSVFile,
    usedPvaluesCSVFile_g=usedPvaluesCSVFile_g,
    usedPvaluesCSVFile_e=usedPvaluesCSVFile_e,
    usedPvaluesCSVFile_t=usedPvaluesCSVFile_t,
    savePlots=savePlots,
    dirToSave = dirToSave, fileID = fileID, useGPU=useGPU)

  if (hidden_vars == FALSE) {
    if (suffStat$maj.rule == FALSE) {
      if (debug) {
        cat("\nRunning conservative PC algorithm...",
            file=logFile, append=TRUE, sep="\n")
      }
      dag <- pc(suffStat=suffStat, indepTest = familyBasedCITest,
                alpha=alpha, labels = paste0(colnames(phen.df), "_", type), verbose = TRUE,
                skel.method="stable", solve.confl=TRUE, conservative=TRUE, u2pd="relaxed")
    } else {
      if (debug) {
        cat("\nRunning PC algorithm with majority rule ...",
            file=logFile, append=TRUE, sep="\n")
      }
      dag <- pc(suffStat=suffStat, indepTest = familyBasedCITest,
                alpha=alpha, labels = paste0(colnames(phen.df), "_", type), verbose = TRUE,
                skel.method="stable", solve.confl=TRUE, maj.rule = TRUE, u2pd="relaxed")
    }
  } else {
    if (suffStat$maj.rule == FALSE) {
      if (debug) {
        cat("\nRunning conservative RFCI algorithm...",
            file=logFile, append=TRUE, sep="\n")
      }
      dag <- rfci(suffStat=suffStat, indepTest = familyBasedCITest,
                alpha=alpha, labels = paste0(colnames(phen.df), "_", type), verbose = TRUE,
                skel.method="stable", conservative=TRUE)
    } else {
      if (debug) {
        cat("\nRunning RFCI algorithm with majority rule ...",
            file=logFile, append=TRUE, sep="\n")
      }
      dag <- rfci(suffStat=suffStat, indepTest = familyBasedCITest,
                alpha=alpha, labels = paste0(colnames(phen.df), "_", type), verbose = TRUE,
                skel.method="stable", maj.rule = TRUE)
    }
  }
  return(dag)
}


#' @importFrom methods as
printEdgeList <- function(dag) {
  nodeLabels <- dag@graph@nodes
  adjM <- as(dag, "amat")
  cat("\nAdjacency Matrix:\n")
  print(adjM)

  cat("\nType 1 Edges:\n")
  edgePairs <- which(adjM == 1, arr.ind=TRUE)
  if (dim(edgePairs)[1] > 0) {
    for (i in 1:dim(edgePairs)[1]) {
      if(adjM[edgePairs[i,2],edgePairs[i,1]] == 1) {
        # it is an undirected edge
        if (edgePairs[i,2] < edgePairs[i,1])
          cat(paste0("\n", nodeLabels[edgePairs[i,2]], " --- ",  nodeLabels[edgePairs[i,1]]))
      } else {
        # it is a directed edge (column causes row)
        cat(paste0("\n", nodeLabels[edgePairs[i,2]], " --> ",  nodeLabels[edgePairs[i,1]]))
      }
    }
  }

  cat("\n\nType 2 Edges:\n")
  edgePairs2 <- which(adjM == 2, arr.ind=TRUE)
  if (dim(edgePairs2)[1] > 0) {
    for (i in 1:dim(edgePairs2)[1]) {
      if(adjM[edgePairs2[i,2],edgePairs2[i,1]] == 2) {
        # it is an undirected edge
        if (edgePairs2[i,2] < edgePairs2[i,1])
          cat(paste0("\n", nodeLabels[edgePairs2[i,2]], " <--> ",  nodeLabels[edgePairs2[i,1]]))
      } else {
        cat(paste0("\n", nodeLabels[edgePairs2[i,2]], " --> ",  nodeLabels[edgePairs2[i,1]]))
      }
    }
  }
}














#TODO
library(igraph) # as.undirected
library(corpcor)# - pseudoinverse

##########################################
# Determining Partial Correlation Matrix #
##########################################

mypmin <- function(vec1, vec2) {
  pminvec <- c()
  for (i in 1:length(vec1)) {
    if (abs(vec1[i]) < abs(vec2[i]))
      pminvec <- c(pminvec, vec1[i])
    else
      pminvec <- c(pminvec, vec2[i])
  }
  return(pminvec)
}

mypmax <- function(vec1, vec2) {
  pmaxvec <- c()
  for (i in 1:length(vec1)) {
    if (abs(vec1[i]) > abs(vec2[i]))
      pmaxvec <- c(pmaxvec, vec1[i])
    else
      pmaxvec <- c(pmaxvec, vec2[i])
  }
  return(pmaxvec)
}

# rule can be "max" or "min"
symmetrizeMatrix <- function(amatrix, rule="max") {
  lowind <- which(lower.tri(amatrix), arr.ind=TRUE)
  uppind <- which(upper.tri(amatrix), arr.ind=TRUE)
  uppind <- uppind[order(uppind[,1]),]
  offdiagind <- rbind(lowind, uppind)
  new.values <- amatrix[offdiagind]
  if (rule == "max") {
    new.values <- mypmax((new.values[1:(length(new.values)/2)]), (new.values[(length(new.values)/2+1):length(new.values)]))
  } else {
    new.values <- mypmin((new.values[1:(length(new.values)/2)]), (new.values[(length(new.values)/2+1):length(new.values)]))
  }

  amatrix[offdiagind] <- new.values
  return(amatrix)
}

convertInvCovToPCor <- function(SigmaInv, symmetrize=FALSE) {
  pcorMatrix <- -SigmaInv
  diag(pcorMatrix) <- -diag(pcorMatrix)
  pcorMatrix <- cov2cor(pcorMatrix)
  if (symmetrize)
    pcorMatrix <- symmetrizeMatrix(pcorMatrix)
  return(pcorMatrix)
}

getPartialCorrelationMatrixFromInvCov <- function(SigmaInv, n, ncov=0, method="z", symmetrize=FALSE) {
  pcorMatrix <- convertInvCovToPCor(SigmaInv, symmetrize)

  # To test if a sample partial correlation vanishes,
  # Fisher's z-transform of the partial correlation can be used.
  # See https://en.wikipedia.org/wiki/Partial_correlation#As_conditional_independence_test


  gp <- (dim(pcorMatrix)[1]-2) + ncov # number of conditioning variables

  # TODO verificar outras estatisticas (copiado da funcao pcor da biblioteca ppcor)
  if (method == "kendall") {
    # Kendall MG, Stuart A. (1973) The Advanced Theory of Statistics,
    # Volume 2 (3rd Edition), ISBN 0-85264-215-6, Section 27.22
    statistics <- pcorMatrix/sqrt(2 * (2 * (n - gp) + 5)/(9 * (n -
                                                                 gp) * (n - 1 - gp)))
    p.values <- 2 * pnorm(-abs(statistics))
  }
  else if (method == "z") {
    # z statistic follows a Gaussian distribution
    # with zero mean and unit standard deviation.
    # The null hyphothesis of zero partial correlation is rejected if
    # stat > qnorm(1 - alpha/2)
    # The two-sided p-values of the tests are:

    # Note that log1p(x) =: log(1+x), so
    # log1p(2*r/(1-r)) = log(1 + 2*r/(1-r))
    # = log((1-r+2r)/(1-r))
    # = log(1+r) - log(1-r)
    # = log((1+r)/(1-r))
    if ((n - gp - 3) > 0) {
      statistics <- sqrt(n - gp - 3)
    } else {
      statistics <- 0
    }
    statistics <- statistics * abs(0.5 * log1p(2*pcorMatrix/(1-pcorMatrix)))

    p.values <- 2 * (1 - pnorm(abs(statistics))) # or 2 * pnorm(-abs(stat))
  }
  else { # shenskin
    # David Shenskin - 2003 - page 1013 - eq. 28.101
    # total number of variables employed in the analysis
    # when there are two predictors and one criterion variable, nu = 3
    nu <- dim(pcorMatrix)[1] + ncov # gp + 2
    statistics <- pcorMatrix * sqrt((n - nu)/(1 - pcorMatrix^2))
    p.values <- 2 * pt(-abs(statistics), (n - nu))
    # statistics <- pcorMatrix * sqrt((n - 2 - gp)/(1 - pcorMatrix^2))
    # p.values <- 2 * pt(-abs(statistics), (n - 2 - gp))
  }
  diag(p.values) <- NA

  return(list(estimates=pcorMatrix, p.values=p.values, statistics=statistics))
}

getPartialCorrelationMatrix <- function(Sigma, n, ncov=0, method="z", symmetrize=FALSE) {
  SigmaInv <- pseudoinverse(Sigma)
  return(getPartialCorrelationMatrixFromInvCov(SigmaInv, n, ncov, method, symmetrize=symmetrize))
}

############################
# Getting Undirected Graph #
############################

# It applies the correction for multiple tests using all off-diagonal values
adjustPValuesMatrix <- function(pvaluesMatrix, correction="fdr") {
  lowind <- which(lower.tri(pvaluesMatrix), arr.ind=TRUE)
  uppind <- which(upper.tri(pvaluesMatrix), arr.ind=TRUE)
  uppind <- uppind[order(uppind[,1]),]
  offdiagind <- rbind(lowind, uppind)
  adjp <- p.adjust(pvaluesMatrix[offdiagind], method=correction)
  pvaluesMatrix[offdiagind] <- adjp
  return(pvaluesMatrix)
}

getUndirectedGraphFromPCorPvalues <- function(pvaluesMatrix, alpha=0.05, correction="fdr", rule="and") {
  nVertices <- dim(pvaluesMatrix)[1]

  if (!is.null(correction)) {
    pvaluesMatrix <- adjustPValuesMatrix(pvaluesMatrix, correction)
  }

  if (rule == "and") {
    pvaluesMatrix <- symmetrizeMatrix(pvaluesMatrix, rule="max")
  } else { #rule == or
    pvaluesMatrix <- symmetrizeMatrix(pvaluesMatrix, rule="min")
  }

  edges <- which(pvaluesMatrix <= alpha, arr.ind=TRUE)
  adjMatrix <- matrix(0, nVertices, nVertices)
  adjMatrix[edges] <- 1
  colnames(adjMatrix) <- colnames(pvaluesMatrix)
  row.names(adjMatrix) <- colnames(pvaluesMatrix)

  udg <- as.undirected(graph.adjacency(adjMatrix))

  return(list(udg=udg, adj.pvalues=pvaluesMatrix))
}

getUndirectedGraphFromPCor <- function(pCor, alpha, correction="fdr") {
  pvaluesMatrix <- pCor$p.values
  diag(pvaluesMatrix) <- 1

  return(getUndirectedGraphFromPCorPvalues(pvaluesMatrix, alpha=0.05, correction))
}

getUndirectedGraphFromCov <- function(Sigma, n, ncov=0, alpha=0.05, correction="fdr", method="z", symmetrize=FALSE) {
  pCor <- getPartialCorrelationMatrix(Sigma, n, ncov, method, symmetrize)
  ret <- getUndirectedGraphFromPCor(pCor, alpha, correction)
  return(list(udg=ret$udg, adj.pvalues=ret$adj.pvalues, pCor=pCor))
}


