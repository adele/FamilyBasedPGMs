##############################################
# Auxiliary functions for simulation studies #
##############################################

# Methodology of Rebonato and Jackel (2000) to create a
# positive definite matrix out of a non-positive definite matrix.
# Reference: Rebonato and Jackel, "The most general methodology for creating a
# valid correlation matrix for risk management and option pricing purposes",
# Journal of Risk, Vol 2, No 2, 2000
fixNonPositiveDefiniteMatrixIter <- function(origMat) {
    cholStatus <- try(u <- chol(origMat), silent = TRUE)
    cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    newMat <- origMat

    iter <- 0
    while (cholError) {

        iter <- iter + 1
        cat("iteration ", iter, "\n")

        # replace -ve eigen values with small +ve number
        newEig <- eigen(newMat)
        newEig2 <- ifelse(newEig$values < 0, 1e-10, newEig$values)

        # create modified matrix eqn 5 from Brissette et al 2007,
        # inv = transp for eig vectors
        newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)

        # normalize modified matrix eqn 6 from Brissette et al 2007
        newMat <- newMat #/sqrt(diag(newMat) %*% t(diag(newMat)))

        # try chol again
        cholStatus <- try(u <- chol(newMat), silent = TRUE)
        cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    }
    return(list(newMat, iter))
}

fixNonPositiveDefiniteMatrix <- function(origMat) {
    return((fixNonPositiveDefiniteMatrixIter(origMat))[[1]])
}


##########################
# Generating the dataset #
##########################

computeSinglePhenotype <- function(Phi_list, fam.nf, phen.mu, phen.var_g, phen.var_e) {
  phen_df <- c()
  for (f in 1:length(fam.nf)) {
    if (f == 1) {
      firstIdx = 1
    } else {
      firstIdx = sum(fam.nf[1:(f-1)]) + 1
    }
    lastIdx = firstIdx + fam.nf[f] - 1

    CovYf <- kronecker(2*Phi_list[[f]], phen.var_g) + kronecker(diag(1,fam.nf[f]), phen.var_e)
    Sigma_f = fixNonPositiveDefiniteMatrix(CovYf)
    phen_f <- mvrnorm(1, mu=phen.mu[firstIdx:lastIdx], Sigma=Sigma_f)
    phen_f <- matrix(phen_f, byrow=T, ncol=1)
    phen_df <- rbind(phen_df, phen_f)
  }
  return(phen_df)
}

#' @importFrom MASS mvrnorm
generateFamilyDependentSEM <- function(betayx, betazx.y, betazy.x,
  varx_g, varx_e, vary.x_g, vary.x_e, varz.xy_g, varz.xy_e,
  fam.nf, pedigrees) {
  Phi_list <- getKinshipList(pedigrees)

  N <- sum(fam.nf)
  muX <- rep(0, N) # as.numeric(kronecker(rep(1, nf), phen.mu0))
  X <- computeSinglePhenotype(Phi_list, fam.nf, muX, varx_g, varx_e)

  muY.X <- betayx * X
  Y <- computeSinglePhenotype(Phi_list, fam.nf, muY.X, vary.x_g, vary.x_e)

  muZ.XY <- betazx.y * X + betazy.x * Y
  Z <- computeSinglePhenotype(Phi_list, fam.nf, muZ.XY, varz.xy_g, varz.xy_e)

  return(data.frame(X=X, Y=Y, Z=Z))
}

generateFamilyDependentFCI_SEM <- function(betazx.u, betazu.x, betawy.u, betawu.y,
                                           varx_g, varx_e, varu_g, varu_e, vary_g, vary_e,
                                           varz.xu_g, varz.xu_e, varw.yu_g, varw.yu_e,
                                           fam.nf, pedigrees) {
  Phi_list <- getKinshipList(pedigrees)

  N <- sum(fam.nf)

  muX <- rep(0, N) # as.numeric(kronecker(rep(1, nf), phen.mu0))
  X <- computeSinglePhenotype(Phi_list, fam.nf, muX, varx_g, varx_e)

  muU <- rep(0, N) # as.numeric(kronecker(rep(1, nf), phen.mu0))
  U <- computeSinglePhenotype(Phi_list, fam.nf, muU, varu_g, varu_e)

  muY <- rep(0, N) # as.numeric(kronecker(rep(1, nf), phen.mu0))
  Y <- computeSinglePhenotype(Phi_list, fam.nf, muY, vary_g, vary_e)

  muZ.XU <- betazx.u * X + betazu.x * U
  Z <- computeSinglePhenotype(Phi_list, fam.nf, muZ.XU, varz.xu_g, varz.xu_e)

  muW.YU <- betawy.u * Y + betawu.y * U
  W <- computeSinglePhenotype(Phi_list, fam.nf, muW.YU, varw.yu_g, varw.yu_e)

  return(data.frame(X=X, Y=Y, Z=Z, W=W, U=U))
}

getScenarioParametersFromList <- function(parametersTable, rowId) {
  betayx <- parametersTable[rowId, 1]
  betazx.y <- parametersTable[rowId, 2]
  betazy.x <- parametersTable[rowId, 3]
  varx_g <- parametersTable[rowId, 4]
  varx_e <- parametersTable[rowId, 5]
  vary.x_g <- parametersTable[rowId=6]
  vary.x_e <- parametersTable[rowId, 7]
  varz.xy_g <- parametersTable[rowId, 8]
  varz.xy_e <- parametersTable[rowId, 9]
}


getFCIScenarioParameters <- function(scenario) {
  if (scenario == 1) {

    betazx.u = 0.75
    betazu.x = 0.80
    betawy.u = 0.75
    betawu.y = 0.80
    varx_g = 0.875
    varx_e = 0.250
    varu_g = 0.425
    varu_e = 0.625
    vary_g = 0.875
    vary_e = 0.250
    varz.xu_g = 0.250
    varz.xu_e = 0.275
    varw.yu_g = 0.250
    varw.yu_e = 0.275
  }

  varx = varx_g + varx_e
  varu = varu_g + varu_e
  vary = vary_g + vary_e

  varz.xu <- varz.xu_g + varz.xu_e
  varz_g = betazx.u^2 * varx_g + betazu.x^2 * varu_g + varz.xu_g
  varz_e = betazx.u^2 * varx_e + betazu.x^2 * varu_e + varz.xu_e
  varz = varz_g + varz_e

  varw.yu <- varw.yu_g + varw.yu_e
  varw_g = betawy.u^2 * vary_g + betawu.y^2 * varu_g + varw.yu_g
  varw_e = betawy.u^2 * vary_e + betawu.y^2 * varu_e + varw.yu_e
  varw = varw_g + varw_e

  covxz = betazx.u * (varx_g + varx_e)
  covxz_g = betazx.u * (varx_g)
  covxz_e = betazx.u * (varx_e)

  rhoxz = covxz/(sqrt(varx)*sqrt(varz))
  rhoxz_g = covxz_g/(sqrt(varx_g)*sqrt(varz_g))
  rhoxz_e = covxz_e/(sqrt(varx_e)*sqrt(varz_e))

  covzw = betazu.x * betawu.y * (varu_g + varu_e)
  covzw_g = betazu.x * betawu.y * (varu_g)
  covzw_e = betazu.x * betawu.y * (varu_e)

  rhozw = covzw/(sqrt(varz)*sqrt(varw))
  rhozw_g = covzw_g/(sqrt(varz_g)*sqrt(varw_g))
  rhozw_e = covzw_e/(sqrt(varz_e)*sqrt(varw_e))

  covzu = betazu.x * (varu_g + varu_e)
  covzu_g = betazu.x * (varu_g)
  covzu_e = betazu.x * (varu_e)

  rhozu = covzu/(sqrt(varz)*sqrt(varu))
  rhozu_g = covzu_g/(sqrt(varz_g)*sqrt(varu_g))
  rhozu_e = covzu_e/(sqrt(varz_e)*sqrt(varu_e))

  covuw = betawu.y * (varu_g + varu_e)
  covuw_g = betawu.y * (varu_g)
  covuw_e = betawu.y * (varu_e)

  rhouw = covuw/(sqrt(varu)*sqrt(varw))
  rhouw_g = covuw_g/(sqrt(varu_g)*sqrt(varw_g))
  rhouw_e = covuw_e/(sqrt(varu_e)*sqrt(varw_e))

  covwy = betawy.u * (vary_g + vary_e)
  covwy_g = betawy.u * (vary_g)
  covwy_e = betawy.u * (vary_e)

  rhowy = covwy/(sqrt(varw)*sqrt(vary))
  rhowy_g = covwy_g/(sqrt(varw_g)*sqrt(vary_g))
  rhowy_e = covwy_e/(sqrt(varw_e)*sqrt(vary_e))

  betaxz = covxz/varz
  betauz = covzu/varz
  betauw = covuw/varw
  betayw = covwy/varw
  betawz = covzw/varz
  betazw = covzw/varw

  covwx.z = - betaxz * covzw - betawz * covxz + betawz * betaxz * varz
  covwx.z_g = - betaxz * covzw_g - betawz * covxz_g + betawz * betaxz * varz_g
  covwx.z_e = - betaxz * covzw_e - betawz * covxz_e + betawz * betaxz * varz_e

  covux.z = - betaxz * covzu - betauz * covxz + betauz * betaxz * varz
  covux.z_g = - betaxz * covzu_g - betauz * covxz_g + betauz * betaxz * varz_g
  covux.z_e = - betaxz * covzu_e - betauz * covxz_e + betauz * betaxz * varz_e

  varx.z = varx - 2 * betaxz * covxz + betaxz^2 * varz
  varx.z_g = varx_g - 2 * betaxz * covxz_g + betaxz^2 * varz_g
  varx.z_e = varx_e - 2 * betaxz * covxz_e + betaxz^2 * varz_e

  varu.z = varu - 2* betauz * covzu + betauz^2 * varz
  varu.z_g = varu_g - 2* betauz * covzu_g + betauz^2 * varz_g
  varu.z_e = varu_e - 2* betauz * covzu_e + betauz^2 * varz_e

  varw.z = varw - 2* betawz * covzw + betawz^2 * varz
  varw.z_g = varw_g - 2* betawz * covzw_g + betawz^2 * varz_g
  varw.z_e = varw_e - 2* betawz * covzw_e + betawz^2 * varz_e


  rhoux.z = covux.z/(sqrt(varx.z) * sqrt(varu.z))
  rhoux.z_g = covux.z_g/(sqrt(varx.z_g) * sqrt(varu.z_g))
  rhoux.z_e = covux.z_e/(sqrt(varx.z_e) * sqrt(varu.z_e))

  covuy.w = - betayw * covuw - betauw * covwy + betauw * betayw * varw
  covuy.w_g = - betayw * covuw_g - betauw * covwy_g + betauw * betayw * varw_g
  covuy.w_e = - betayw * covuw_e - betauw * covwy_e + betauw * betayw * varw_e

  covzy.w = - betayw * covzw - betazw * covwy + betazw * betayw * varw
  covzy.w_g = - betayw * covzw_g - betazw * covwy_g + betazw * betayw * varw_g
  covzy.w_e = - betayw * covzw_e - betazw * covwy_e + betazw * betayw * varw_e

  varu.w = varu - 2* betauw * covuw + betauw^2 * varw
  varu.w_g = varu_g - 2* betauw * covuw_g + betauw^2 * varw_g
  varu.w_e = varu_e - 2* betauw * covuw_e + betauw^2 * varw_e

  varz.w = varz - 2* betazw * covzw + betazw^2 * varw
  varz.w_g = varz_g - 2* betazw * covzw_g + betazw^2 * varw_g
  varz.w_e = varz_e - 2* betazw * covzw_e + betazw^2 * varw_e

  vary.w = vary - 2* betayw * covwy + betayw^2 * varw
  vary.w_g = vary_g - 2* betayw * covwy_g + betayw^2 * varw_g
  vary.w_e = vary_e - 2* betayw * covwy_e + betayw^2 * varw_e

  rhouy.w = covuy.w/(sqrt(vary.w) * sqrt(varu.w))
  rhouy.w_g = covuy.w_g/(sqrt(vary.w_g) * sqrt(varu.w_g))
  rhouy.w_e = covuy.w_e/(sqrt(vary.w_e) * sqrt(varu.w_e))

  rhozy.w = covzy.w/(sqrt(vary.w) * sqrt(varz.w))
  rhozy.w_g = covzy.w_g/(sqrt(vary.w_g) * sqrt(varz.w_g))
  rhozy.w_e = covzy.w_e/(sqrt(vary.w_e) * sqrt(varz.w_e))

  rhowx.z = covwx.z/(sqrt(varx.z) * sqrt(varw.z))
  rhowx.z_g = covwx.z_g/(sqrt(varx.z_g) * sqrt(varw.z_g))
  rhowx.z_e = covwx.z_e/(sqrt(varx.z_e) * sqrt(varw.z_e))

  return(list(betazx.u=betazx.u, betazu.x=betazu.x, betawy.u=betawy.u, betawu.y=betawu.y,
              varx_g=varx_g, varx_e=varx_e, varu_g=varu_g, varu_e=varu_e, vary_g=vary_g, vary_e=vary_e,
              varz.xu_g=varz.xu_g, varz.xu_e=varz.xu_e, varw.yu_g=varw.yu_g, varw.yu_e=varw.yu_e,
              rhoxz=rhoxz, rhoxz_g=rhoxz_g, rhoxz_e=rhoxz_e, #rho(X,Z)
              rhozu=rhozu, rhozu_g=rhozu_g, rhozu_e=rhozu_e,
              rhouw=rhouw, rhouw_g=rhouw_g, rhouw_e=rhouw_e,
              rhowy=rhowy, rhowy_g=rhowy_g, rhowy_e=rhowy_e, #rho(W,Y)
              rhowx.z=rhowx.z, rhowx.z_g=rhowx.z_g, rhowx.z_e=rhowx.z_e,
              rhoux.z=rhoux.z, rhoux.z_g=rhoux.z_g, rhoux.z_e=rhoux.z_e,
              rhouy.w=rhouy.w, rhouy.w_g=rhouy.w_g, rhouy.w_e=rhouy.w_e,
              rhozy.w=rhozy.w, rhozy.w_g=rhozy.w_g, rhozy.w_e=rhozy.w_e,
              rhozw=rhozw, rhozw_g=rhozw_g, rhozw_e=rhozw_e) #rho(Z,W)
         )
}


getScenarioParameters <- function(scenario) {
  if (scenario == 1) {

    parameters <- c(0.0, 0.75, 0.80, 0.875, 0.250, 0.425, 0.625, 0.250, 0.275)
    betayx <- parameters[1]
    betazx.y <- parameters[2]
    betazy.x <- parameters[3]
    varx_g <- parameters[4]
    varx_e <- parameters[5]
    vary.x_g <- parameters[6]
    vary.x_e <- parameters[7]
    varz.xy_g <- parameters[8]
    varz.xy_e <- parameters[9]

    pCor_t <- diag(3)
    pCor_t[1,2] = pCor_t[2,1] = -0.5539314
    pCor_t[1,3] = pCor_t[3,1] = 0.7392961
    pCor_t[2,3] = pCor_t[3,2] = 0.7492686

    pCor_g <- diag(3)
    pCor_g[1,2] = pCor_g[2,1] = -0.6112120
    pCor_g[1,3] = pCor_g[3,1] = 0.8143451
    pCor_g[2,3] = pCor_g[3,2] = 0.7218537

    pCor_e <- diag(3)
    pCor_e[1,2] = pCor_e[2,1] = -0.4949355
    pCor_e[1,3] = pCor_e[3,1] = 0.5816751
    pCor_e[2,3] = pCor_e[3,2] = 0.7698004


    #r12t, r13t, r23t, r12g, r13g, r23g, r12e, r13e, r23e
    marg_corr <- c(0.0000000,  0.5880770,  0.6060122,
                   0.0000000,  0.6966364,  0.5178755,
                   0.0000000,  0.4152274,  0.7003010)

    #r12.3t, r13.2t, r23.1t, r12.3g, r13.2g, r23.1g, r12.3e, r13.2e, r23.1e
    part_corr <- c( -0.5539314,  0.7392961,  0.7492686,
                    -0.6112120,  0.8143451,  0.7218537,
                    -0.4949355,  0.5816751,  0.7698004)

    trueAdjDAG_t <- matrix(0, 3, 3)
    trueAdjDAG_t[3,1] <- 1
    trueAdjDAG_t[3,2] <- 1

    trueAdjUDG_t <- matrix(1, 3, 3)
    diag(trueAdjUDG_t) <- 0

    trueAdjDAG_g <- trueAdjDAG_t
    trueAdjDAG_e <- trueAdjDAG_t

    trueAdjUDG_g <- trueAdjUDG_t
    trueAdjUDG_e <- trueAdjUDG_t

#    > print(c(varx, vary, varz))
#    parameters parameters parameters
#      1.125000   1.050000   1.829813
#    > print(c(varx.y, varx.z, vary.z, vary.x, varz.x, varz.y))
#    parameters parameters parameters parameters parameters parameters
#     1.1250000  0.7359361  0.6643867  1.0500000  1.1970000  1.1578125
#    >
#    > print(c(varx_g, vary_g, varz_g))
#    parameters parameters parameters
#      0.875000   0.425000   1.014188
#    > print(c(varx.y_g, varx.z_g, vary.z_g, vary.x_g, varz.x_g, varz.y_g))
#    parameters parameters parameters parameters parameters parameters
#     0.8750000  0.4854311  0.3265660  0.4250000  0.5220000  0.7421875
#    >
#    > print(c(varx_e, vary_e, varz_e))
#    parameters parameters parameters
#      0.250000   0.625000   0.815625
#    > print(c(varx.y_e, varx.z_e, vary.z_e, vary.x_e, varz.x_e, varz.y_e))
#    parameters parameters parameters parameters parameters parameters
#     0.2500000  0.2505050  0.3378208  0.6250000  0.6750000  0.4156250
#    >
#    > print(c(rhoyx, rhoxz, rhoyz))
#    parameters parameters parameters
#     0.0000000  0.5880770  0.6060122
#    > print(c(rhoyx_g, rhoxz_g, rhoyz_g))
#    parameters parameters parameters
#     0.0000000  0.6966364  0.5178755
#    > print(c(rhoyx_e, rhoxz_e, rhoyz_e))
#    parameters parameters parameters
#     0.0000000  0.4152274  0.7003010
#    >
#    > print(c(varx_g/varx, vary_g/vary, varz_g/varz))
#    parameters parameters parameters
#     0.7777778  0.4047619  0.5542576
#    > print(c(rhoyx.z, rhozx.y, rhozy.x))
#    parameters parameters parameters
#    -0.5539314  0.7392961  0.7492686
#    > print(c(rhoyx.z_g, rhozx.y_g, rhozy.x_g))
#    parameters parameters parameters
#    -0.6112120  0.8143451  0.7218537
#    > print(c(rhoyx.z_e, rhozx.y_e, rhozy.x_e))
#    parameters parameters parameters
#    -0.4949355  0.5816751  0.7698004

  } else if (scenario == 2) {

    parameters <- c(-0.8, 0.80, 0.00, 0.425, 0.600, 0.250, 0.725, 0.475, 0.275)
    betayx <- parameters[1]
    betazx.y <- parameters[2]
    betazy.x <- parameters[3]
    varx_g <- parameters[4]
    varx_e <- parameters[5]
    vary.x_g <- parameters[6]
    vary.x_e <- parameters[7]
    varz.xy_g <- parameters[8]
    varz.xy_e <- parameters[9]

    pCor_t <- diag(3)
    pCor_t[1,2] = pCor_t[2,1] = -5.139177e-01
    pCor_t[1,3] = pCor_t[3,1] = 5.859564e-01
    pCor_t[2,3] = pCor_t[3,2] = -2.596613e-16

    pCor_g <- diag(3)
    pCor_g[1,2] = pCor_g[2,1] = -0.6478292
    pCor_g[1,3] = pCor_g[3,1] = 0.4739953
    pCor_g[2,3] = pCor_g[3,2] = 0.0000000

    pCor_e <- diag(3)
    pCor_e[1,2] = pCor_e[2,1] = -0.4349207
    pCor_e[1,3] = pCor_e[3,1] = 0.6932896
    pCor_e[2,3] = pCor_e[3,2] = 0.0000000


    #r12t, r13t, r23t, r12g, r13g, r23g, r12e, r13e, r23e
    marg_corr <- c( -0.6341981,  0.6830606, -0.4331958,
                    -0.7218537,  0.6034262, -0.4355854,
                    -0.5884368,  0.7633486, -0.4491824)

    #r12.3t, r13.2t, r23.1t, r12.3g, r13.2g, r23.1g, r12.3e, r13.2e, r23.1e
    part_corr <- c( -5.139177e-01,  5.859564e-01, -2.596613e-16,
                    -0.6478292,  0.4739953,  0.0000000,
                    -0.4349207,  0.6932896,  0.0000000)

    trueAdjDAG_t <- matrix(0, 3, 3)
    trueAdjDAG_t[2,1] <- 1
    trueAdjDAG_t[3,1] <- 1

    trueAdjUDG_t <- matrix(0, 3, 3)
    trueAdjDAG_t[1,2] <- 1
    trueAdjDAG_t[1,3] <- 1
    trueAdjDAG_t[2,1] <- 1
    trueAdjDAG_t[3,1] <- 1

    trueAdjDAG_g <- trueAdjDAG_t
    trueAdjDAG_e <- trueAdjDAG_t

    trueAdjUDG_g <- trueAdjUDG_t
    trueAdjUDG_e <- trueAdjUDG_t

#    print(c(varx, vary, varz))
#    parameters parameters parameters
#         1.025      1.631      1.406
#    > print(c(varx.y, varx.z, vary.z, vary.x, varz.x, varz.y))
#    parameters parameters parameters parameters parameters parameters
#     0.6127376  0.5467639  1.3249289  0.9750000  0.7500000  1.1421521
#    >
#    > print(c(varx_g, vary_g, varz_g))
#    parameters parameters parameters
#         0.425      0.522      0.747
#    > print(c(varx.y_g, varx.z_g, vary.z_g, vary.x_g, varz.x_g, varz.y_g))
#    parameters parameters parameters parameters parameters parameters
#     0.2150680  0.2824981  0.4307988  0.2500000  0.4750000  0.6126435
#    >
#    > print(c(varx_e, vary_e, varz_e))
#    parameters parameters parameters
#         0.600      1.109      0.659
#    > print(c(varx.y_e, varx.z_e, vary.z_e, vary.x_e, varz.x_e, varz.y_e))
#    parameters parameters parameters parameters parameters parameters
#     0.3976695  0.2642657  0.8941301  0.7250000  0.2750000  0.5295085
#    >
#    > print(c(rhoyx, rhoxz, rhoyz))
#    parameters parameters parameters
#    -0.6341981  0.6830606 -0.4331958
#    > print(c(rhoyx_g, rhoxz_g, rhoyz_g))
#    parameters parameters parameters
#    -0.7218537  0.6034262 -0.4355854
#    > print(c(rhoyx_e, rhoxz_e, rhoyz_e))
#    parameters parameters parameters
#    -0.5884368  0.7633486 -0.4491824
#    >
#    > print(c(varx_g/varx, vary_g/vary, varz_g/varz))
#    parameters parameters parameters
#     0.4146341  0.3200490  0.5312945
#    > print(c(rhoyx.z, rhozx.y, rhozy.x))
#       parameters    parameters    parameters
#    -5.139177e-01  5.859564e-01 -2.596613e-16
#    > print(c(rhoyx.z_g, rhozx.y_g, rhozy.x_g))
#    parameters parameters parameters
#    -0.6478292  0.4739953  0.0000000
#    > print(c(rhoyx.z_e, rhozx.y_e, rhozy.x_e))
#    parameters parameters parameters
#    -0.4349207  0.6932896  0.0000000


  } else if (scenario == 3) {
    parameters <- c(-0.8, -0.9,  0.9,  0.6,  0.55,  0.925,  0.3,  0.25, 0.975)
    betayx <- parameters[1]
    betazx.y <- parameters[2]
    betazy.x <- parameters[3]
    varx_g <- parameters[4]
    varx_e <- parameters[5]
    vary.x_g <- parameters[6]
    vary.x_e <- parameters[7]
    varz.xy_g <- parameters[8]
    varz.xy_e <- parameters[9]

    pCor_t <- diag(3)
    pCor_t[1,2] = pCor_t[2,1] = 0.004686763
    pCor_t[1,3] = pCor_t[3,1] = -0.567485469
    pCor_t[2,3] = pCor_t[3,2] = 0.668964732

    pCor_g <- diag(3)
    pCor_g[1,2] = pCor_g[2,1] = 0.4244711
    pCor_g[1,3] = pCor_g[3,1] = -0.7658174
    pCor_g[2,3] = pCor_g[3,2] = 0.8659171

    pCor_e <- diag(3)
    pCor_e[1,2] = pCor_e[2,1] = -0.4077895
    pCor_e[1,3] = pCor_e[3,1] = -0.4348273
    pCor_e[2,3] = pCor_e[3,2] = 0.4466625

    #r12t, r13t, r23t, r12g, r13g, r23g, r12e, r13e, r23e
    marg_corr <- c(-0.6126326, -0.7592639, 0.8092371,
                   -0.5416214, -0.7821601, 0.8771786,
                   -0.7347634, -0.7364439, 0.7460787)

    #r12.3t, r13.2t, r23.1t, r12.3g, r13.2g, r23.1g, r12.3e, r13.2e, r23.1e
    part_corr <- c( 0.004686763, -0.567485469, 0.668964732,
                    0.4244711, -0.7658174, 0.8659171,
                   -0.4077895, -0.4348273, 0.4466625)

    # notation col -> row
    trueAdjDAG_t <- matrix(0, 3, 3)
    trueAdjDAG_t[3,1] <- 1
    trueAdjDAG_t[2,1] <- 1
    trueAdjDAG_t[3,2] <- 1

    trueAdjUDG_t <- matrix(1, 3, 3)
    diag(trueAdjUDG_t) <- 0

    trueAdjDAG_g <- trueAdjDAG_t
    trueAdjDAG_e <- trueAdjDAG_t

    trueAdjUDG_g <- trueAdjUDG_t
    trueAdjUDG_e <- trueAdjUDG_t

#    > print(c(varx, vary, varz))
#    [1] 1.15000 1.96100 5.23531
#    > print(c(varx.y, varx.z, vary.z, vary.x, varz.x, varz.y))
#    [1] 0.7183835 0.4870461 0.6768104 1.2250000 2.2172500 1.8068906
#    >
#    > print(c(varx_g, vary_g, varz_g))
#    [1] 0.60000 1.30900 2.57389
#    > print(c(varx.y_g, varx.z_g, vary.z_g, vary.x_g, varz.x_g, varz.y_g))
#    [1] 0.4377287 0.2341570 0.3454861 0.9250000 0.9992500 0.6045603
#    >
#    > print(c(varx_e, vary_e, varz_e))
#    [1] 0.55000 0.65200 2.66142
#    > print(c(varx.y_e, varx.z_e, vary.z_e, vary.x_e, varz.x_e, varz.y_e))
#    [1] 0.2806548 0.2528892 0.3313243 0.3000000 1.2180000 1.2023304
#    >
#    > print(c(rhoyx, rhoxz, rhoyz))
#    [1] -0.6126326 -0.7592639  0.8092371
#    > print(c(rhoyx_g, rhoxz_g, rhoyz_g))
#    [1] -0.5416214 -0.7821601  0.8771786
#    > print(c(rhoyx_e, rhoxz_e, rhoyz_e))
#    [1] -0.7347634 -0.7364439  0.7460787
#    >
#    > print(c(varx_g/varx, vary_g/vary, varz_g/varz))
#    [1] 0.5217391 0.6675166 0.4916404
#    > print(c(rhoyx.z, rhozx.y, rhozy.x))
#    [1]  0.004686763 -0.567485469  0.668964732
#    > print(c(rhoyx.z_g, rhozx.y_g, rhozy.x_g))
#    [1]  0.4244711 -0.7658174  0.8659171
#    > print(c(rhoyx.z_e, rhozx.y_e, rhozy.x_e))
#    [1] -0.4077895 -0.4348273  0.4466625
  } else if (scenario == 4) {
    ###################################################
    # Defining zero genetic variability for X.        #
    # There is an unshielded collider in the          #
    # environmental network and only the edge         #
    # Y->Z in the genetic network.                    #
    # The total network han an unshielded collider.   #
    # Note that the effect for the rhoyx.z_e is       #
    # negative even though there is no negative coef. #
    # Zero heritability for X and high heritability   #
    # for Y and Z.                                    #
    ###################################################

    betayx = 0
    betazx.y = 0.25
    betazy.x = 0.55
    varx_e =  2.5
    varx_g =  0
    vary.x_e <- 0.75
    vary.x_g <- 1.75
    varz.xy_e <- 0.25
    varz.xy_g <- 0.75

    pCor_t <- diag(3)
    pCor_t[1,2] = pCor_t[2,1] = -0.2412257
    pCor_t[1,3] = pCor_t[3,1] = 0.3676073
    pCor_t[2,3] = pCor_t[3,2] = 0.6562050

    pCor_g <- diag(3)
    pCor_g[1,2] = pCor_g[2,1] = -0.03723522
    pCor_g[1,3] = pCor_g[3,1] = 0.00000000
    pCor_g[2,3] = pCor_g[3,2] = 0.64325443

    pCor_e <- diag(3)
    pCor_e[1,2] = pCor_e[2,1] = -0.4258620
    pCor_e[1,3] = pCor_e[3,1] = 0.6201737
    pCor_e[2,3] = pCor_e[3,2] = 0.6897489

    #r12t, r13t, r23t, r12g, r13g, r23g, r12e, r13e, r23e
    marg_corr <- c( 0.0000000, 0.2858310, 0.6288281,
                    0.0000000, 0.0000000, 0.6432544,
                    0.0000000, 0.4967813, 0.5986164)


    #r12.3t, r13.2t, r23.1t, r12.3g, r13.2g, r23.1g, r12.3e, r13.2e, r23.1e
    part_corr <- c( -0.2412257, 0.3676073, 0.6562050,
                    -0.03723522, 0.00000000, 0.64325443,
                    -0.4258620, 0.6201737, 0.6897489)

    trueAdjDAG_t <- matrix(0, 3, 3)
    trueAdjDAG_t[3,1] <- 1
    trueAdjDAG_t[3,2] <- 1

    trueAdjDAG_g <- trueAdjDAG_t
    trueAdjDAG_g[3,1] <- 0

    trueAdjDAG_e <- trueAdjDAG_t

    trueAdjUDG_t <- matrix(1, 3, 3)
    diag(trueAdjUDG_t) <- 0

    trueAdjUDG_g <- trueAdjDAG_g
    trueAdjDAG_t[2,3] <- 1

    trueAdjUDG_e <- trueAdjUDG_t

    # print(c(varx, vary, varz))
    # 2.5000 2.5000 1.9125
    # print(c(varx.y, varx.z, vary.z, vary.x, varz.x, varz.y))
    # 2.500000 2.295752 1.511438 2.500000 1.756250 1.156250

    # print(c(varx_g, vary_g, varz_g))
    # 0.000000 1.750000 1.279375
    # print(c(varx.y_g, varx.z_g, vary.z_g, vary.x_g, varz.x_g, varz.y_g))
    # 0.0000000 0.1366328 1.0273159 1.7500000 1.2793750 0.7500000

    # print(c(varx_e, vary_e, varz_e))
    # 2.500000 0.750000 0.633125
    # print(c(varx.y_e, varx.z_e, vary.z_e, vary.x_e, varz.x_e, varz.y_e))
    # 2.500000 2.159119 0.484122 0.750000 0.476875 0.406250

    # print(c(rhoyx, rhoxz, rhoyz))
    # 0.0000000 0.2858310 0.6288281
    # print(c(rhoyx_g, rhoxz_g, rhoyz_g))
    # 0.0000000 0.0000000 0.6432544
    # print(c(rhoyx_e, rhoxz_e, rhoyz_e))
    # 0.0000000 0.4967813 0.5986164

    # print(c(varx_g/varx, vary_g/vary, varz_g/varz))
    # 0.0000000 0.7000000 0.6689542
    # print(c(rhoyx.z, rhozx.y, rhozy.x))
    # -0.2412257  0.3676073  0.6562050
    # print(c(rhoyx.z_g, rhozx.y_g, rhozy.x_g))
    # -0.03723522  0.00000000  0.64325443
    # print(c(rhoyx.z_e, rhozx.y_e, rhozy.x_e))
    # -0.4258620  0.6201737  0.6897489
  } else {
    return(NULL)
  }
  return(list(betayx=betayx, betazx.y=betazx.y, betazy.x=betazy.x,
              varx_g=varx_g, varx_e=varx_e,
              vary.x_g=vary.x_g, vary.x_e=vary.x_e,
              varz.xy_g=varz.xy_g, varz.xy_e=varz.xy_e,
              pCor_t=pCor_t, pCor_g=pCor_g, pCor_e=pCor_e,
              marg_corr=marg_corr, part_corr=part_corr,
              DAG_t=trueAdjDAG_t, DAG_g=trueAdjDAG_g, DAG_e=trueAdjDAG_e,
              UDG_t=trueAdjUDG_t, UDG_g=trueAdjUDG_g, UDG_e=trueAdjUDG_e))
}


