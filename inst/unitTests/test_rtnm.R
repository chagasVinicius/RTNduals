# Unit tests for MBR-class methods
test_mbr <- function()
{
    data("dt4rtn", package = "RTN")
    tfs1 <- dt4rtn$tfs[c("FOXM1", "E2F2")]
    tfs2 <- dt4rtn$tfs[c("PTTG1", "RARA")]
    ##mbrPreprocess
    rmbr <- mbrPreprocess(gexp=dt4rtn$gexp, regulatoryElements1=tfs1, 
                           regulatoryElements2=tfs2, gexpIDs=dt4rtn$gexpIDs)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status["Preprocess"]=="[x]" && status[1]=="[x]" && 
                status[1]=="[x]")
    ##mbrPermutation
    rmbr <- mbrPermutation(rmbr, nPermutations=10, estimator="pearson")
    status <- mbrGet(rmbr, what="status")
    checkTrue(status["Permutation"]=="[x]" && status[2]=="[x]" && 
                status[2]=="[x]")
    ##mbrBootstrap
    rmbr <- mbrBootstrap(rmbr, estimator="pearson", nBootstrap=10)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status["Bootstrap"]=="[x]" && status[3]=="[x]" && 
                status[3]=="[x]")
    ##mbrDpiFilter
    rmbr <- mbrDpiFilter(rmbr)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status["DPI.filter"]=="[x]" && status[4]=="[x]" && 
                status[4]=="[x]")
    ##mbr.combine.TNIs
    tni1 <- mbrGet(rmbr, what="TNI1")
    tni2 <- mbrGet(rmbr, what="TNI2")
    rmbr <- tni2mbrPreprocess(tni1, tni2)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[1:4]=="[x]" && status[1:4]=="[x]" && 
                status[1:4]=="[x]")
    ##mbrAssociation
    rmbr <- mbrAssociation(rmbr, prob=0, estimator="pearson")
    status <- mbrGet(rmbr, what="status")
    motifsInformation <- mbrGet(rmbr, what="motifsInformation")
    checkTrue(status["Association"]=="[x]" && 
                  is.data.frame(motifsInformation) && 
                  (ncol(motifsInformation)==11 || ncol(motifsInformation)==10))
    ##mbr.motifs
    rmbr <- mbrDuals(rmbr)
    motifsInformation <- mbrGet(rmbr, what="motifsInformation")
    checkTrue(is.data.frame(motifsInformation))
}
