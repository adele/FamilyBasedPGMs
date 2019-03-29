'%!in%' <- function(x,y)!('%in%'(x,y))

# sex: 1="male", 2="female"

#' @importFrom stats runif
simulatePedigree <- function(initialIndId, famid, nind, nInitFounderCouples, probToMarry, probChildren) {
    pedigree <- data.frame()
    getNonParents <- function(pedigree, candidates, curFounders) {
        foundersList <- pedigree$foundersid
        nonParents <- c()
        for (i in candidates) {
            if (length(which(curFounders %in% foundersList[[i]])) == 0)
                nonParents <- c(nonParents, pedigree[i,"indid"])
        }
        return(nonParents)
    }

    marry <- function(pedigree, indIDIter, groomid, brideid, probChildren) {
        pedigree[which(pedigree[,2] == groomid),"ismarried"] <- TRUE
        pedigree[which(pedigree[,2] == brideid),"ismarried"] <- TRUE
        curgeneration = pedigree[which(pedigree[,2] == brideid),"generationid"]
        curfounders =  sort(unique(c(pedigree$foundersid[[which(pedigree[,"indid"] == brideid)]], pedigree$foundersid[[which(pedigree[,"indid"] == groomid)]])))

        curfam = pedigree[which(pedigree[,2] == brideid),"famid"]

        r = runif(1, 0, 1)
        nChildren = which(probChildren == max(probChildren))[1]
        for (i in 1:length(probChildren)) {
            if (r > sum(probChildren[1:i]))
                nChildren = i+1
        }

        for (child in 1:nChildren) {
            r = runif(1,0,1)
            sex = 1
            if (r < 0.5) {
                sex = 2
            }
            pedigree <- rbind(pedigree, c(curfam, indIDIter+child, brideid, groomid, sex, FALSE, curgeneration+1, NA))
            pedigree$foundersid[[indIDIter+child]] <- curfounders
        }
        return(pedigree)
    }

    for (i in 0:(nInitFounderCouples-1)) {
        pedigree <- rbind(pedigree, c(famid, i*2+1, NA, NA, 1, FALSE, 1))
        pedigree <- rbind(pedigree, c(famid, i*2+2, NA, NA, 2, FALSE, 1))
    }
    colnames(pedigree) <- c("famid", "indid", "momid", "dadid", "sex", "ismarried", "generationid")

    listOfFounderCouples <- list()
    for (id in 1:(2*nInitFounderCouples))
        listOfFounderCouples[[id]] = id

    pedigree$foundersid <- listOfFounderCouples

    genID = 0
    while (max(pedigree$ind) < nind) {
        genID = genID + 1

        # Trying to marry every single woman with a single man
        # of the same generation.
        # A founder man is created if no single man is found.
        singleWomen <- intersect(which(pedigree[,"ismarried"] == FALSE), which(pedigree[,"sex"] == 2))
        singleWomen <- intersect(singleWomen, which(pedigree[,"generationid"] == genID))

        noSingleWomen = FALSE
        if (length(singleWomen) == 0)
            noSingleWomen = TRUE

        for (womID in singleWomen) {
            singleMen <- intersect(which(pedigree[,"ismarried"] == FALSE), which(pedigree[,"sex"] == 1))
            singleMen <- intersect(singleMen, which(pedigree[,"generationid"] == genID))
            singleMen <- getNonParents(pedigree, singleMen, pedigree[womID, "foundersid"][[1]])
            if (length(singleMen) > 0) {
                pedigree <- marry(pedigree, max(pedigree$ind), singleMen[1], womID, probChildren)
            } else {
                r = runif(1,0,1)
                if (r < probToMarry) {
                    pedigree <- rbind(pedigree, c(famid, (max(pedigree$ind)+1), NA, NA, 1, FALSE, genID, (max(pedigree$ind)+1)))
                    pedigree <- marry(pedigree, max(pedigree$ind), max(pedigree$ind), womID, probChildren)
                }
            }
        }

        # Trying to marry every single man with a single woman
        # of the same generation.
        # A founder woman is created if no single woman is found.
        singleMen <- intersect(which(pedigree[,"ismarried"] == FALSE), which(pedigree[,"sex"] == 1))
        singleMen <- intersect(singleMen, which(pedigree[,"generationid"] == genID))

        noSingleMen = FALSE
        if (length(singleMen) == 0)
            noSingleMen = TRUE

        for (manID in singleMen) {
            singleWomen <- intersect(which(pedigree[,"ismarried"] == FALSE), which(pedigree[,"sex"] == 2))
            singleWomen <- intersect(singleWomen, which(pedigree[,"generationid"] == genID))
            singleWomen <- getNonParents(pedigree, singleWomen, pedigree[manID, "foundersid"][[1]])
            if (length(singleWomen) > 0) {
                pedigree <- marry(pedigree, max(pedigree$ind), manID, singleWomen[1], probChildren)
                print("WARNING: This is unexpected!! Verify, please!")
            } else {
                r = runif(1,0,1)
                if (r < probToMarry) {
                    pedigree <- rbind(pedigree, c(famid, (max(pedigree$ind)+1), NA, NA, 2, FALSE, genID, (max(pedigree$ind)+1)))
                    pedigree <- marry(pedigree, max(pedigree$ind), manID, max(pedigree$ind), probChildren)
                }
            }
        }

        if (noSingleWomen && noSingleMen) {
            print(paste0("WARNING: No singles in generation ", genID, "!"))
            genID = genID-2
        }
    }

    pedigree = pedigree[1:nind,1:5]
    pedigree$indid = pedigree$indid + initialIndId
    pedigree$momid = pedigree$momid + initialIndId
    pedigree$dadid = pedigree$dadid + initialIndId

    colnames(pedigree) <- c("famid", "id", "momid", "dadid", "sex")

    return(pedigree)
}

#' @title Pedigree Simulation
#'
#' @param fam.nf A integer vector where the entry i indicates the number of individuals at family i.
#' @param nInitFounderCouples A integer value representing the number of founders in the first generation.
#' @param probToMarry  A scalar between 0 and 1 indicating the probability of a individual getting married.
#' Marriages only occur between individuals of the same generation.
#' @param probChildren A vector with values between 0 and 1 in which entry i indicates the probability
#' of a couple has i children.
#'
#' @examples
#' fam.nf <- rep(10, 20)  # for 20 families with 10 individuals
#' nInitFounderCouples <- 1
#' probChildren <- c(0, 0, 0.3, 0.3, 0.3, 0.1)
#' probToMarry = 0.85
#' pedigrees <- simulatePedigrees(fam.nf, nInitFounderCouples, probToMarry, probChildren)
#'
#' @export simulatePedigrees
simulatePedigrees <- function(fam.nf, nInitFounderCouples, probToMarry, probChildren) {
    pedigrees <- c()
    n <- 0
    for (famid in 1:length(fam.nf)) {
      nind = fam.nf[famid]
      curped = simulatePedigree(n, famid, nind, nInitFounderCouples, probToMarry, probChildren)
      pedigrees <- rbind(pedigrees, curped)
      n <- dim(pedigrees)[1]
    }
    return (pedigrees)
}

#' @title Plotting of Pedigree Chart
#'
#' @description Creates a chart of the pedigree of family \code{famid}.
#'
#' @param pedigrees A data.frame object with at least the columns: "famid", "id", "momid", "dadid", and "sex"
#' @param famid A integer value identifying the family for which the pedigree chart should be plotted.
#'
#' @examples
#' \donttest{
#'   data(scen1)
#'   pedigrees <- scen1$pedigrees
#'
#'   # plot the pedigree chart of the fifth simulated family
#'   plotFamilyPedigree(pedigrees, 5)
#' }
#'
#' @importFrom kinship2 pedigree
#'
#' @export plotFamilyPedigree
plotFamilyPedigree <- function(pedigrees, famid) {
    curped = pedigrees[which(pedigrees$famid == famid),]
    if (nrow(curped) > 0) {
        ped = pedigree(id=curped$id, momid=curped$momid, dadid=curped$dadid, sex=curped$sex)
        plot(ped)
    } else {
        print(paste("There is no family with id", famid))
    }
}


