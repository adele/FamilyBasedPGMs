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

#' @importFrom MASS mvrnorm
generateFamilyDependentSEM <- function(betayx, betazx.y, betazy.x,
  varx_g, varx_e, vary.x_g, vary.x_e, varz.xy_g, varz.xy_e,
  fam.nf, pedigrees) {

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

getScenarioParametersFromList <- function(parametersTable, rowId) {
  betayx <- parametersTable[rowId, 1]
  betazx.y <- parametersTable[rowId, 2]
  betazy.x <- parametersTable[rowId, 3]
  varx_g <- parametersTable[rowId, 4]
  varx_e <- parametersTable[rowId, 5]
  vary.x_g <- parametersTable[rowId, 6]
  vary.x_e <- parametersTable[rowId, 7]
  varz.xy_g <- parametersTable[rowId, 8]
  varz.xy_e <- parametersTable[rowId, 9]
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


