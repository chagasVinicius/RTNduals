##generic functions
##------------------------------------------------------------------------------
setGeneric("mbrPreprocess",
           function(gexp, regulatoryElements1, regulatoryElements2, 
                    verbose=TRUE,...)
             standardGeneric("mbrPreprocess"), package="RTNduals")
setGeneric("mbrPermutation",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrPermutation"), package="RTNduals")
setGeneric("mbrBootstrap",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrBootstrap"), package="RTNduals")
setGeneric("mbrDpiFilter",
           function(object, verbose = TRUE, ...)
             standardGeneric("mbrDpiFilter"), package="RTNduals")
setGeneric("mbrAssociation",
           function(object, regulatoryElements1=NULL, 
                    regulatoryElements2=NULL, minRegulonSize=50, 
                    prob=0.95, estimator='spearman', 
                    pAdjustMethod="BH", verbose=TRUE)
             standardGeneric("mbrAssociation"), package="RTNduals")
setGeneric("mbrDuals",
           function(object, supplementary.table = NULL,
                    evidenceColname=NULL, verbose = TRUE)
             standardGeneric("mbrDuals"), package="RTNduals")
setGeneric("tni2mbrPreprocess",
           function(tni1, tni2,  verbose = TRUE)
             standardGeneric("tni2mbrPreprocess"), package="RTNduals")
setGeneric("mbrGet",
           function(object, what="status")
             standardGeneric("mbrGet"), package="RTNduals")
