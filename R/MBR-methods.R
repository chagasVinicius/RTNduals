################################################################################
##########################         MBR-methods      ############################
################################################################################
##------------------------------------------------------------------------------
#########################################
###class builder
#########################################
.MBRmaker <- 
    function(class="MBR", gexp=gexp,
             regulatoryElements1=regulatoryElements1,
             regulatoryElements2=regulatoryElements2)
    {
        #---creates the 'MBR' object
        object <- newMBR(class=class, gexp=gexp,
                         regulatoryElements1=regulatoryElements1,
                         regulatoryElements2=regulatoryElements2,
                         TNI1=NULL, TNI2=NULL, 
                         testedElementsTNI1=character(),
                         testedElementsTNI2=character(),
                         dualRegulons=character(),
                         results=list(), para=list(), summary=list(),
                         status=character())
        #---status
        status <- rep('[ ]', 1, 5)
        names(status) <- c('Preprocess', 'Permutation',
                           'Bootstrap', 'DPI.filter', 
                           'Association')
        #---parameters
        sum.info.para <- list()
        sum.info.para$TNIs$perm <- NA
        sum.info.para$TNIs$boot <- NA
        sum.info.para$TNIs$dpi <- NA
        sum.info.para$MBR$association <- matrix(NA, 1, 3)
        colnames(sum.info.para$MBR$association) <- c('minRegulonSize',
                                                     'prob',
                                                     'estimator')
        rownames(sum.info.para$MBR$association) <- 'Parameter'
        #---summary motifsInformation
        sum.info.summary <- list()
        sum.info.summary$MBR$Duals <- matrix(NA, 1, 1)
        colnames(sum.info.summary$MBR$Duals) <- 'numberDuals'
        rownames(sum.info.summary$MBR$Duals) <- 'duals'
        
        #---set
        object <- .mbr.set(name="status", para=status, object=object)
        object <- .mbr.set(name="para", para=sum.info.para, object=object)
        object <- .mbr.set(name="summary", para=sum.info.summary, object=object)
        return(object)
    }




#----------------------------------------------------------
#it creates the 'MBR' class object
newMBR <- 
    function(class, gexp, regulatoryElements1, regulatoryElements2, 
             TNI1, TNI2, testedElementsTNI1, testedElementsTNI2,
             dualRegulons, results, para, summary, status)
    {
        #---checks
        if(missing(gexp)) stop("NOTE: 'gexp' is missing ", call.=FALSE)
        if(missing(regulatoryElements1)) 
            stop("NOTE: 'regulatoryElements1' is missing", call.=FALSE)
        if(missing(regulatoryElements2)) 
            stop("NOTE: 'regulatoryElements2' is missing", call.=FALSE)
        mbr.checks(name='TNI', TNI1)
        mbr.checks(name='TNI', TNI2)
        mbr.checks(name='testedElementsTNI', testedElementsTNI1)
        mbr.checks(name='testedElementsTNI', testedElementsTNI2)
        mbr.checks(name='dualRegulons', dualRegulons)
        mbr.checks(name='results', results)
        mbr.checks(name='para', para)
        mbr.checks(name='summary', summary)
        mbr.checks(name='status', status)
        mbr.checks(name='gexp', gexp)
        mbr.checks(name='regulatoryElements1', regulatoryElements1)
        mbr.checks(name='regulatoryElements2', regulatoryElements2)
        
        #---creating TNIs
        regulonsTNI1 <- new("TNI", gexp=gexp, 
                            transcriptionFactors=regulatoryElements1)
        regulonsTNI2 <- new("TNI", gexp=gexp, 
                            transcriptionFactors=regulatoryElements2)
        
        #---creating the MBR-object
        new(class, TNI1=regulonsTNI1, TNI2=regulonsTNI2,
            testedElementsTNI1=testedElementsTNI1,
            testedElementsTNI2=testedElementsTNI2,
            dualRegulons=dualRegulons, results=results,
            para=para, summary=summary, status=status)
    }

##------------------------------------------------------------------------------
#' A preprocessing function for objects of class MBR.
#'
#' @param gexp A numerical matrix, typically with mRNA and/or miRNA expression 
#' values.
#' @param regulatoryElements1 A named vector with regulatory elements listed in 
#' 'gexp' rownames.
#' @param regulatoryElements2 A named vector with regulatory elements listed in 
#' 'gexp' rownames.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed on to 
#' \code{\link[RTN:tni.preprocess]{tni.preprocess}} function.
#' @return A preprocessed 'MBR-class' object.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrPreprocess-methods
#' @aliases mbrPreprocess
#' @export

##Regulons pre-processing method
setMethod("mbrPreprocess",
          "matrix",
          function(gexp, regulatoryElements1, regulatoryElements2, 
                   verbose=TRUE,...)
          {
           ##---
           mbr.checks(name="verbose", para=verbose)  
           object <- .MBRmaker(class="MBR", gexp=gexp,
                              regulatoryElements1=regulatoryElements1,
                              regulatoryElements2=regulatoryElements2)  
           ##---get TNIs
           TNI1 <- mbrGet(object, what="TNI1")
           TNI2 <- mbrGet(object, what="TNI2")
           ##---pre-processing TNIs
           if(verbose) cat("-Preprocessing TNI objects...\n\n")
           TNI1 <- tni.preprocess(TNI1, verbose=verbose,...=...)
           TNI2 <- tni.preprocess(TNI2, verbose=verbose,...=...)
           tfs1 <- tni.get(TNI1, what="tfs")
           tfs2 <- tni.get(TNI2, what="tfs")
           mbr.checks(name="regulatoryElements", 
                      para=c(tfs1,tfs2))
           
           ##---set
           object <- .mbr.set(name="TNI1", para=TNI1, object=object)
           object <- .mbr.set(name="TNI2", para=TNI2, object=object)
           object <- .mbr.set(name="statusUpdate", 
                              para="Preprocess", object=object)

           
           return(object)
          }
)
##------------------------------------------------------------------------------
#' Inference of transcriptional networks.
#'
#' This function takes an MBR object and computes two transcriptional networks 
#' inferred 
#' by mutual information (with multiple hypothesis testing corrections).
#' 
#' @param object A preprocessed object of class \linkS4class{MBR}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed on to the 
#' \code{\link[RTN:tni.permutation]{tni.permutation}} function.
#' @return An \linkS4class{MBR} object with two mutual information matrices, 
#' one in each "TNI" slot.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrPermutation-methods
#' @aliases mbrPermutation
#' @export

## permutation
setMethod("mbrPermutation",
          "MBR",
          function(object, verbose=TRUE, ...)
          {
           ##---checks
           mbr.checks(name="object", para=object)
           mbr.checks(name="verbose", para=verbose)
           ##---get TNIs
           TNI1 <- mbrGet(object, what="TNI1")
           TNI2 <- mbrGet(object, what="TNI2")
           
           ##---permutation TNIs
           if(verbose)
               cat("-Performing permutation analysis for two TNI objects...\n\n")
           TNI1 <- tni.permutation(TNI1, verbose=verbose,...=...)
           TNI2 <- tni.permutation(TNI2, verbose=verbose,...=...)
           #---get
           tni1_summary <- tni.get(TNI1, what="summary")
           tni2_summary <- tni.get(TNI2, what="summary")
           mbr_summary <- mbrGet(object, what="summary")
           mbr_para <- mbrGet(object, what="para")
           
           #---changes
           mbr_para$TNIs$perm <- tni1_summary$para$perm
           
           mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
           mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
           
           #---set
           object <- .mbr.set(name="TNI1", para=TNI1, object=object)
           object <- .mbr.set(name="TNI2", para=TNI2, object=object)
           object <- .mbr.set(name="statusUpdate", 
                              para="Permutation", object=object)
           object <- .mbr.set(name="para", para=mbr_para, object=object)
           object <- .mbr.set(name="summary", para=mbr_summary, object=object)
           return(object)
          }
)

#' Inference of consensus transcriptional networks.
#'
#' This function takes an MBR object and computes two consensus transcriptional 
#' networks.
#'
#' @param object A processed objec of class \linkS4class{MBR} evaluated by the 
#' method \code{\link[RTNduals:mbrPermutation]{mbrPermutation}}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed to the 
#' \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} function.
#' @return An \linkS4class{MBR} object with two consensus mutual information 
#' matrices, one in each "TNI" slot.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' ##---mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrBootstrap-methods
#' @aliases mbrBootstrap
#' @export

##------------------------------------------------------------------------------
## bootstrap method
setMethod("mbrBootstrap",
          "MBR",
          function(object, verbose=TRUE, ...)
          {
           ##---checks
           mbr.checks(name="object", para=object)
           mbr.checks(name="verbose", para=verbose)
           ##---get TNIs
           TNI1 <- mbrGet(object, what="TNI1")
           TNI2 <- mbrGet(object, what="TNI2")
           
           ##---bootstrap TNIs
           if(verbose) cat("-Performing bootstrap analysis for two TNI 
                           objects...\n\n")
           TNI1 <- tni.bootstrap(TNI1, verbose=verbose,...=...)
           TNI2 <- tni.bootstrap(TNI2, verbose=verbose,...=...)
           
           #---get
           tni1_summary <- tni.get(TNI1, what="summary")
           tni2_summary <- tni.get(TNI2, what="summary")
           mbr_summary <- mbrGet(object, what="summary")
           mbr_para <- mbrGet(object, what="para")
           
           #---changes
           mbr_para$TNIs$boot <- tni1_summary$para$boot
           
           mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
           mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
           
           #---set
           object <- .mbr.set(name="TNI1", para=TNI1, object=object)
           object <- .mbr.set(name="TNI2", para=TNI2, object=object)
           object <- .mbr.set(name="statusUpdate", para="Bootstrap", object=object)
           object <- .mbr.set(name="para", para=mbr_para, object=object)
           object <- .mbr.set(name="summary", para=mbr_summary, object=object)
           return(object)
          }
)

#' A filter based on the Data Processing Inequality (DPI) algorithm.
#'
#' This function takes an MBR object and computes two transcriptional networks 
#' filtered by the data processing inequality algorithm.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by 
#' the methods
#'  \code{\link[RTNduals:mbrPermutation]{mbrPermutation}} and 
#'  \code{\link[RTNduals:mbrBootstrap]{mbrBootstrap}}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed to the 
#' \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}} function.
#' @return An \linkS4class{MBR} object with two DPI-filtered mutual information 
#' matrices, one in each "TNI" slot.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' ##---mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' ##---mbrDpiFilter
#' rmbr <- mbrDpiFilter(rmbr)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrDpiFilter-methods
#' @aliases mbrDpiFilter
#' @export

##------------------------------------------------------------------------------
## dpi filter method
setMethod("mbrDpiFilter",
          "MBR",
          function(object, verbose=TRUE, ...)
          {
           ##---checks
           mbr.checks(name="object", para=object)
           mbr.checks(name="verbose", para=verbose)
           ##---get TNIs
           TNI1 <- mbrGet(object, what="TNI1")
           TNI2 <- mbrGet(object, what="TNI2")
           
           ##---Dpi filter TNIs
           if(verbose) cat("-Applying dpi filter for two TNI objects...\n")
           TNI1 <-tni.dpi.filter(TNI1, verbose=verbose, ...=...)
           TNI2 <-tni.dpi.filter(TNI2, verbose=verbose, ...=...)
           #---get
           tni1_summary <- tni.get(TNI1, what="summary")
           tni2_summary <- tni.get(TNI2, what="summary")
           mbr_summary <- mbrGet(object, what="summary")
           mbr_para <- mbrGet(object, what="para")
           
           #---changes
           mbr_para$TNIs$dpi <- tni1_summary$para$dpi
           
           mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
           mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
           
           #---set
           object <- .mbr.set(name="TNI1", para=TNI1, object=object)
           object <- .mbr.set(name="TNI2", para=TNI2, object=object)
           object <- .mbr.set(name="statusUpdate", para="DPI.filter", object=object)
           object <- .mbr.set(name="para", para=mbr_para, object=object)
           object <- .mbr.set(name="summary", para=mbr_summary, object=object)
           return(object)
          }
)

#' Motifs analysis and inference of 'dual regulons'.
#'
#' This function takes an MBR object and compares the shared regulon 
#' targets in order to test whether regulon pairs agree on the predicted 
#' downstream effects.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the 
#' methods \code{\link[RTNduals:mbrPermutation]{mbrPermutation}}, 
#' \code{\link[RTNduals:mbrBootstrap]{mbrBootstrap}} 
#' and \code{\link[RTNduals:mbrDpiFilter]{mbrDpiFilter}}.
#' @param regulatoryElements1 An optional character vector specifying which 
#' 'TNI1' regulatory elements should be evaluated. If 'NULL' all regulatory 
#' elements will be evaluated.
#' @param regulatoryElements2 An optional character vector specifying which 
#' 'TNI2' regulatory elements should be evaluated. If 'NULL' all regulatory 
#' elements will be evaluated.
#' @param minRegulonSize A single integer or numeric value specifying the 
#' minimum number of elements in a regulon. Gene sets with fewer than this 
#' number are removed from the analysis.
#' @param prob A quantile filter applyed to the association metric used to 
#' infer 'dual regulons'.
#' @param estimator A character value specifying the estimator used in the 
#' association analysis. One of "spearman" (default), "kendall", or "pearson", 
#' can be abbreviated.
#' @param pAdjustMethod A single character value specifying the p-value 
#' adjustment method to be used (see 'p.adjust' for details).
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with two data.frames in the slot 
#' 'results' listing the inferred 'dual regulons' and a hypergeometric test 
#' for each 'dual regulon'.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' ##---mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' ##---mbrDpiFilter
#' rmbr <- mbrDpiFilter(rmbr)
#' ##---mbrAssociation
#' rmbr <- mbrAssociation(rmbr, prob=0.75)
#'
#' @import RTN 
#' @importFrom stats p.adjust phyper
#' @importFrom stats cor quantile
#'
#' @import methods
#' @docType methods
#' @rdname mbrAssociation-methods
#' @aliases mbrAssociation
#' @export

##------------------------------------------------------------------------------
##Inference of duals
setMethod("mbrAssociation",
          "MBR",
          function(object, regulatoryElements1=NULL, regulatoryElements2=NULL, 
                   minRegulonSize=30, prob=0.95, 
                   estimator='spearman', pAdjustMethod="BH", verbose=TRUE)
          {
            ##---gets
              TNI1 <- mbrGet(object, what="TNI1")
              tni1_gexp <- tni.get(TNI1, what="gexp")
              tni1_para <- tni.get(TNI1, what="para")
              TNI2 <- mbrGet(object, what="TNI2")
            ##---checks
              mbr.checks(name="minRegulonSize", para=minRegulonSize)
              mbr.checks(name="prob", para=prob)
              mbr.checks(name="estimator", para=estimator)
              mbr.checks(name="pAdjustMethod", para=pAdjustMethod)
              mbr.checks(name="verbose", para=verbose)
              if(is.null(regulatoryElements1))
                  {
                  if(verbose) 
                      cat("-Selecting regulatory elements from TNI1 object...\n")
                  regulatoryElements1 <- tni.get(TNI1, "tfs")
                  } else
                      {
                          regulatoryElements1 <- .checkRegel(TNI1, 
                                                             regulatoryElements1)
                          }
              if(is.null(regulatoryElements2))
                  {
                  if(verbose) 
                      cat("-Selecting regulatory elements from TNI2 object...\n")
                  regulatoryElements2 <- tni.get(TNI2, "tfs")
                  } else
                      {
                          regulatoryElements2 <- .checkRegel(TNI2,
                                                             regulatoryElements2)
                          }
              mbr.checks(name="regulatoryElements", 
                         para=c(regulatoryElements1, regulatoryElements2))
              mbr.checks(name="numberRegElements", 
                         para=c(regulatoryElements1, regulatoryElements2))
              ##-----get regulons
              what <- "refregulons.and.mode"
              ##if (tnet=="dpi") what <- "regulons.and.mode"
              regulons1 <- tni.get(TNI1, what=what)
              regulons2 <- tni.get(TNI2, what=what)
              
              ##-----get regulatory elements
              regulons1 <- regulons1[regulatoryElements1]
              regulons2 <- regulons2[regulatoryElements2]
              
              ##-----get regulons by min size
              size1 <- unlist(lapply(regulons1, length))
              idx <- size1 >= (minRegulonSize)
              regulons1 <- regulons1[idx]
              regulatoryElements1 <- regulatoryElements1[idx]
              size1 <- size1[idx]
              
              ##---
              size2 <- unlist(lapply(regulons2, length))
              idx <- size2 >= (minRegulonSize)
              regulons2 <- regulons2[idx]
              regulatoryElements2 <- regulatoryElements2[idx]
              size2 <- size2[idx]
              
              ##-----group regulons and regulatory element
              regulons <- c(regulons1, regulons2)
              regel <- c(regulatoryElements1, regulatoryElements2)
              
              ##-----Correlation
              tnet <- .regMatrix(regulons, regel, getNames=FALSE)
              tbmi <- tnet
              
              ##-----
              tnet <- .tni.cor(tni1_gexp, tnet, estimator=estimator, 
                               dg=0, asInteger=FALSE, 
                               mapAssignedAssociation=TRUE)
              regcor <- cor(tnet[, regulatoryElements1], 
                            tnet[, regulatoryElements2], method=estimator)
              
              ##-----
              rownames(regcor) <- names(regulatoryElements1)
              colnames(regcor) <- names(regulatoryElements2)
              
              ##-----select motifs based on 'prob' quantile
              if(verbose) 
                  cat(paste("-Getting duals in probs >", prob, "...\n", 
                            sep = ""))
              pvlist <- .motifsquantile(regcor=regcor, th=prob)
              
              if(nrow(pvlist) > 0)
                  {
                  ##-----Mutual Information
                  pvlist <- .getMIorP(pvlist, regel, tbmi, size1, size2)
                  ##----PadjustValue
                  if(!is.null(TNI1@results$adjpv))
                      {
                      tbpValue <- TNI1@results$adjpv
                      cutoff <- tni1_para$perm$pValueCutoff
                      pvlist <- .getMIorP(pvlist, regel, tbpValue, size1, size2, 
                                          mutualInformation=FALSE, cutoff=cutoff)
                      }
                  
                  ##-----Jaccard
                  pvlist <- .jcOverlap(pvlist, regel, tbmi)
                  
                  ##-----Hypergeometric
                  ##tni1 <- mbrGet(object, what="TNI1")
                  universe <- rownames(tni1_gexp)
                  hyperresults <- .mbr.hyper(pvlist=pvlist, regulons=regulons, 
                                             regel=regel, universe=universe, 
                                             pAdjustMethod=pAdjustMethod, 
                                             verbose=verbose)
                  pvlist$Hypergeometric.Adjusted.Pvalue <- hyperresults$Adjusted.Pvalue
                  pvlist$Hypergeometric.Pvalue <- hyperresults$Pvalue
                  }else
                      {
                          warning("No 'dual' has been found for the input parameters.")
                          }
              ##-----organize pvlist
              if("MI.Adjusted.Pvalue"%in%colnames(pvlist))
                  {
                  pvlist <- pvlist[,c("Regulon1","Size.Regulon1","Regulon2",
                                      "Size.Regulon2","Jaccard.coefficient", 
                                      "Hypergeometric.Pvalue",
                                      "Hypergeometric.Adjusted.Pvalue", 
                                      "MI","MI.Adjusted.Pvalue","R","Quantile")]
                  } else
                      {
                          pvlist <- pvlist[,c("Regulon1","Size.Regulon1","Regulon2",
                                              "Size.Regulon2","Jaccard.coefficient", 
                                              "Hypergeometric.Pvalue",
                                              "Hypergeometric.Adjusted.Pvalue", 
                                              "MI","R","Quantile")]
                          }
              
              ##-----para
              mbr_para <- mbrGet(object,what="para")
              sum.info.par <- c(minRegulonSize, prob, estimator)
              mbr_para$MBR$association['Parameter', ] <- sum.info.par
              
              ##---
              mbr_summary <- mbrGet(object, what="summary")
              info.summary.results <- nrow(pvlist)
              mbr_summary$MBR$Duals['duals',] <- info.summary.results
              ##-----set
              object <- .mbr.set(name="statusUpdate", 
                                 para="Association", object=object)
              object <- .mbr.set(name="para", 
                                 para=mbr_para, object=object)
              object <- .mbr.set(name="summary", 
                                 para=mbr_summary, object=object)
              object <- .mbr.set(name="testedElementsTNI1", 
                                 para=regulatoryElements1, object=object)
              object <- .mbr.set(name="testedElementsTNI2", 
                                 para=regulatoryElements2, object=object)
              object <- .mbr.set(name="dualRegulons", 
                                 para=rownames(pvlist), object=object)
              object <- .mbr.set(name="motifsInformation", 
                                 para=pvlist, object=object)
              object <- .mbr.set(name="hypergeometricResults", 
                                 para=hyperresults, object=object)    
              return(object)
              }
          )


#' A summary for results from the MBR methods.
#'
#' This function lists the inferred 'dual regulons' and, if available, 
#' adds external evidences.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the 
#' method \code{\link[RTNduals:mbrAssociation]{mbrAssociation}}.
#' @param supplementary.table An optional 'data.frame' with three columns 
#' representing 
#' (1) regulatory elements of 'TNI1', (2) regulatory elements of 'TNI2', and 
#' (3) external evidences between the regulatory elements.
#' @param evidenceColname A single character value specifying a column in 
#' the 'supplementary.table'.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with a data.frame in the slot 'results' 
#' listing the input additional evidences.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' ##---mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' ##---mbrDpiFilter
#' rmbr <- mbrDpiFilter(rmbr)
#' ##---mbrAssociation
#' rmbr <- mbrAssociation(rmbr, prob=0.75)
#' ##---a 'toy' table with supplementary evidences
#' motifsInformation <- mbrGet(rmbr, what="motifsInformation")
#' n <- nrow(motifsInformation)
#' supplementaryTable <- motifsInformation[1:n,c("Regulon1","Regulon2")]
#' supplementaryTable$ToyEvidence <- rnorm(n)
#' ##---mbrDuals
#' rmbr <- mbrDuals(rmbr, supplementary.table = supplementaryTable, 
#' evidenceColname = "ToyEvidence")
#' ##---motifsInformation with 'Evidence'
#' motifsInformation <- mbrGet(rmbr, what="motifsInformation")
#' head(motifsInformation)
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import methods
#' @docType methods
#' @rdname mbrDuals-methods
#' @aliases mbrDuals
#' @export

##------------------------------------------------------------------------------
##organize duals
setMethod( "mbrDuals",
           "MBR",
           function(object, supplementary.table=NULL, evidenceColname, 
                    verbose=TRUE)
           {
            ##---checks
            mbr.checks(name="object", para=object)
            ##---
            motifsInformation <- mbrGet(object, what="motifsInformation")
            if(is.null(dim(motifsInformation)))
             stop("'motifsInformation' seems null!")
            if(verbose) cat("-Sorting by the R value...\n")
            idx <- sort(abs(motifsInformation[,"R"]), decreasing=TRUE, 
                        index.return=TRUE)
            motifsInformation <- motifsInformation[idx$ix, ]
            
            #---set
            object <- .mbr.set(name="motifsInformation", 
                               para=motifsInformation, object=object)
            object <- .mbr.set(name="dualRegulons", 
                               para=rownames(motifsInformation), object=object)
            if(!is.null(supplementary.table))
            {
             ##---checks
             if(missing(evidenceColname)) 
              stop("'evidenceColname' should be a character value present in 
                   colnames of supplementary.table!")
             mbr.checks(name="supplementary.table", para=supplementary.table)
             mbr.checks(name="uniqueInput", para=supplementary.table)
             mbr.checks(name="evidenceColname", para=evidenceColname)
             
             ##---consistency
             if(verbose) 
              cat("-Checking the 'supplementary.table' consistency...\n")
             supplementary.table <- .consisSuppTable(object, 
                                                     supplementary.table, 
                                                     evidenceColname, 
                                                     verbose=verbose)
             ##---find duals
             object <- .checkLoops(object, supplementary.table, 
                                   evidenceColname, verbose=verbose)
            }
            return(object)
           }
)

#' A preprocessing function for objects of class MBR.
#'
#' This function merges two TNI class objects and creates one MBR class object.
#'
#' @param tni1 A 'TNI' class object.
#' @param tni2 Another 'TNI' class object
#' @param verbose A single logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object.
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' ##---mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' rmbr <- mbrDpiFilter(rmbr)
#' ##---tni2mbrPreprocess
#' tni1 <- mbrGet(rmbr, what="TNI1")
#' tni2 <- mbrGet(rmbr, what="TNI2")
#' rmbr <- tni2mbrPreprocess(tni1, tni2)
#'
#' @import methods
#' @docType methods
#' @rdname tni2mbrPreprocess-methods
#' @aliases tni2mbrPreprocess
#' @export

##------------------------------------------------------------------------------
##Combine two TNIs produced separately
setMethod("tni2mbrPreprocess",
          "TNI",
          function (tni1,  tni2,  verbose=TRUE)
          {
           if(missing(tni1)) stop("NOTE: 'tni1' is missing ", call.=FALSE)
           if(missing(tni2)) stop("NOTE: 'tni2' is missing ", call.=FALSE)        
           mbr.checks (name='tni', para=tni1)
           mbr.checks (name='tni', para=tni2)
           .combineTNIs (tni1=tni1, tni2=tni2, verbose=verbose)
           #---get
           gexp <- tni.get(tni1, what="gexp")
           regulatoryElements1 <- tni.get(tni1, what="tfs")
           regulatoryElements2 <- tni.get(tni2, what="tfs")
           ##---- creates MBR object
           object <- .MBRmaker(class="MBR", gexp=gexp,
                               regulatoryElements1=regulatoryElements1,
                               regulatoryElements2=regulatoryElements2)
           #---TNIs Update
           object <- .mbr.set(name="TNI1", para=tni1, object=object)
           object <- .mbr.set(name="TNI2", para=tni2, object=object)
           #---statu update
           tni_status <- tni.get(tni1, what="status")
           status <- names(tni_status[tni_status=="[x]"])
           object <- .mbr.set(name="statusUpdate", para=status, object=object)
           #---get Updates
           tni1_summary <- tni.get(tni1, what="summary")
           tni2_summary <- tni.get(tni2, what="summary")
           mbr_summary <- mbrGet(object, what="summary")
           mbr_para <- mbrGet(object, what="para")
           
           ##---permutation
           
           mbr_para$TNIs$perm <- tni1_summary$para$perm
           mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
           mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
           ##---bootstrap
           mbr_para$TNIs$boot <- tni1_summary$para$boot
           mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
           mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
           ##---summary dpi.filter
           mbr_para$TNIs$dpi <- tni1_summary$para$dpi
           mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
           mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
           colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
           #---set
           object <- .mbr.set(name="para", para=mbr_para, object=object)
           object <- .mbr.set(name="summary", para=mbr_summary, object=object)
           return (object)
          }
)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod( "show",
           "MBR",
           function(object)
           {
            cat("an MBR (Motifs Between Regulons) object:\n")
            message("--status:")
            print(object@status, quote=FALSE)
           }
)

#' Get information from individual slots in MBR object.
#' 
#' Get information from individual slots in an MBR object and any available 
#' results from previous analysis.
#' 
#' @param object A preprocessed object of class \linkS4class{MBR}
#' @param what a single character value specifying which information should be 
#' retrieved from the slots. Options: "TNI1", "TNI2", "testedElementsTNI1", 
#' "testedElementsTNI2", "dualRegulons", "results", "para", "summary", 
#' "status", "motifsInformation" and "hyperResults"
#' @return A slot content from a object of class 'MBR' \linkS4class{MBR} object
#' @examples
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' ##---mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' ##---get the 'TNI1' slot using 'mbrGet'
#' tni1 <- mbrGet(rmbr, what="TNI1")
#' 
#' @import methods
#' @docType methods
#' @rdname mbrGet-methods
#' @aliases mbrGet
#' @export
##------------------------------------------------------------------------------
##get slots from MBR object
setMethod( "mbrGet",
           "MBR", 
           function(object, what="status")
           {
            ##---check input arguments
            mbr.checks(name="object", para=object)
            mbr.checks(name="mbrGet", para=what)
            ##---Association options any change needs update!
            optsAssoci <- c("testedElementsTNI1", "testedElementsTNI2", 
                            "dualRegulons", "motifsInformation", "results", 
                            "hyperResults")
            ##---get query
            if(what=="TNI1")
            {
             query <- object@TNI1
            }
            else if(what=="TNI2")
            {
             query <- object@TNI2
            }
            else if(what=="para")
            {
             query <- object@para
            }
            else if(what=="summary")
            {
             query <- object@summary
            }
            else if(what=="status")
            {
             query <- object@status
            }
            else if(what%in%optsAssoci)
            {
             if(object@status["Association"] != "[x]")
             {
              warning("NOTE: input 'object' needs 'mbrAssociation' 
                      evaluation!")
              query <- NULL
             } else {
              if(what=="testedElementsTNI1")
              {
               query <- object@testedElementsTNI1
              }
              else if(what=="testedElementsTNI2")
              {
               query <- object@testedElementsTNI2
              }
              else if(what=="dualRegulons")
              {
               query <- object@dualRegulons
              }
              else if(what=="results")
              {
               query <- object@results
              }
              else if(what=="motifsInformation")
              {
               query <- object@results$motifsInformation
              }
              else if(what=="hyperResults")
              {
               query <- object@results$hypergeometricResults
              }
             }
            }
            return(query)
           }
)
