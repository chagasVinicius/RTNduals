
################################################################################
##########################    plot dual regulons    ############################
################################################################################

#' Plot shared target clouds between dual regulons.
#'
#' This function plots the shared target clouds between a regulon pair.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by 
#' the method \code{\link[RTNduals:mbrAssociation]{mbrAssociation}}.
#' @param names.motifs A vector with 'dual regulon' indentifiers from the 
#' 'motifsInformation' table.
#' @param filepath A character string indicating the file path where the plot 
#' should be saved.
#' @param alpha  The alpha transparency, a number in [0,1].
#' @param lncols A vector of length 2 indicating the colors of the negative 
#' and positive target clouds, respectively.
#' @param lwd  Line width, a decimal value (between 0 and 1).
#' @param estimator A character string indicating the association metric. 
#' One of "spearman" (default), "kendall", or "pearson", can be abbreviated.
#' @return A plot with the shared target clouds between dual regulons.
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
#' ##---mbrDuals
#' rmbr <- mbrDuals(rmbr)
#' ##---
#' dual <- mbrGet(rmbr, what="dualRegulons")[1]
#' mbrPlotDuals(rmbr, names.motifs=dual)
#'
#' @import graphics
#' @importFrom grDevices adjustcolor dev.off pdf colorRampPalette
#' @importFrom graphics abline axis par plot.new plot.window points title 
#' legend
#' @export

##------------------------------------------------------------------------------
mbrPlotDuals <- function(object, names.motifs = NULL, filepath=NULL, 
                           alpha=0.80,lncols=c("darkgreen","darkorange3"), 
                           lwd=0.70, estimator="spearman")
{
  ##----check object class
  mbr.checks(name="object", para=object)
  mbr.checks(name="estimator", para=estimator)
  
  rtni <- .merge.tnis(object)
  rtni_para <- tni.get(rtni, what="para")
  estimator <- rtni_para$perm$estimator
  motifstb <- mbrGet(object, what="motifsInformation")
  motifstb <- .namesMotifs.check(motifstb, names.motifs)
  
  res <- apply(motifstb, 1, function (mtfs)
  {
    mtfs <- as.character(mtfs)
    reg1 <- mtfs[1]
    reg2 <- mtfs[2]
    rval <- as.numeric(mtfs[3])
    labelMotif <- paste(reg1, reg2, sep=".vs.")
    if(!is.null(filepath))
    {
      file <- paste(filepath, labelMotif, sep="")
    }else
    {
      file <- NULL
    }
    .tni.plot.greement(rtni=rtni, estimator=estimator, 
                       duals=c(reg1, reg2), corVal=rval,file=file, 
                       alpha=alpha, lwd=lwd, lncols=lncols)
  })
}

##------------------------------------------------------------------------------
##subfunction for 'mbrPlotDuals'
.tni.plot.greement<-function(rtni,duals,corVal,file=NULL,
                             lncols=c("blue","red"), 
                             bgcols=lncols, lwd=0.70, alpha=0.80, 
                             estimator='spearman', sharedTargets=TRUE, 
                             mapAssignedAssociation=TRUE)
{
  ##---
  tfs <- tni.get(rtni, "tfs")
  idx1 <- match(duals, names(tfs))
  idx2 <- match(duals, tfs)
  idxcheck<-which(is.na(idx1))
  idx1[idxcheck]<-idx2[idxcheck]
  duals<-tfs[idx1]
  ##---
  reftnet <- tni.get(rtni, "refnet")
  gexp <- tni.get(rtni, "gexp")
  tnet<-reftnet[,duals]
  xy<-.tni.cor(gexp,tnet,asInteger=FALSE,estimator=estimator, 
               mapAssignedAssociation=mapAssignedAssociation)
  if(sharedTargets)
  {
    idx<-rowSums(tnet!=0)==2
    tnet<-tnet[idx,]
    xy<-xy[idx,]
  } else
  {
    idx<-rowSums(xy!=0)>=1
    tnet<-tnet[idx,]
    xy<-xy[idx,]
  }
  ##---
  xlab=paste(names(duals)[1],"targets (R)")
  ylab=paste(names(duals)[2],"targets (R)")
  xlim=c(-1.0,1.0)
  ylim=c(-1.0,1.0)
  bgcols[1]<-colorRampPalette(c(lncols[1],"white"))(30)[15]
  #bgcols[2]<-"white"
  bgcols[2]<-colorRampPalette(c(lncols[2],"white"))(30)[15]
  #bgcols[4]<-"white"
  bgcols<-adjustcolor(bgcols,alpha.f=alpha)
  ##---plot
  if(!is.null(file))
  {
    pdf(file=paste(file,".pdf",sep=""), height=3, width=3)
  }
  par(mgp=c(2.2, 0.5, 0),mar=c(3.5, 3.5, 1, 1) + 0.1)
  plot.new()
  plot.window(ylim=xlim,xlim=ylim)
  axis(2,cex.axis=1,las=1,tcl=-0.15,lwd=2)
  axis(1,cex.axis=1,las=1,tcl=-0.15,lwd=2)
  title(xlab=xlab,ylab=ylab,cex.lab=1)
  
  if(corVal<0)
  {
    ##---negative Dual
    tpp<-xy[(sign(tnet[, 1])==1 & sign(tnet[, 2])==-1),]
    points(tpp, col=lncols[1], pch=21, cex=0.7, bg=bgcols[1], lwd=lwd)
    
    tpp<-xy[sign(tnet[, 1])==-1 & sign(tnet[, 2])==1,]
    points(tpp,col=lncols[1],pch=21,cex=0.7,bg="white", lwd=lwd)
  }else
  {
    ##---positive Dual
    tpp<-xy[rowSums(sign(tnet))==2, ]
    points(tpp,col=lncols[2],pch=21,cex=0.7,bg="white", lwd=lwd)
    
    tpp<-xy[rowSums(sign(tnet))==-2, ]
    points(tpp,col=lncols[2],pch=21,cex=0.7,bg=bgcols[2], lwd=lwd)
  }
  
  ##---legend
  legend("topright", legend=paste("R=", corVal, sep=" "), bty="n")
  if(!is.null(file))
  {
    dev.off()
    cat(paste("File '", paste(file,".pdf",sep=""),"' generated!\n\n", 
              sep=""))
  }
  ##---report
  colnames(xy)<-paste(names(duals),"(R)",sep="")
  nms<-rownames(xy)
  annot<-rtni@annotation[nms,]
  report<-cbind(annot,format(round(xy,3)))
  invisible(report)
}


##subfunction for 'mbrPlotDuals'
.merge.tnis <- function (object)
{
  TNI1 <- mbrGet(object, "TNI1")
  TNI2 <- mbrGet(object, "TNI2")
  elreg1 <- tni.get(TNI1, "tfs")
  elreg2 <- tni.get(TNI2, "tfs")
  elregs <- c (elreg1, elreg2)
  rtni_merge <-
    new ("TNI",
         gexp = tni.get(TNI1, "gexp"),
         transcriptionFactors = elregs)
  rtni_merge@annotation <- object@TNI1@annotation
  rtni_merge@para <- tni.get(TNI1, "para")
  #---
  mirmt <- tni.get(TNI2, "refnet") [, elreg2]
  rtni_merge@results$tn.ref <- cbind (tni.get(TNI1, "refnet") [, elreg1], 
                                      mirmt)
  mirmt <- tni.get(TNI2, "tnet") [, elreg2]
  rtni_merge@results$tn.dpi <- cbind (tni.get(TNI1, "tnet") [, elreg1], 
                                      mirmt)
  rtni_merge@status [1:4] <- "[x]"
  return (rtni_merge)
}

##subfunction for 'mbrPlotDuals'
.namesMotifs.check <- function(motifstb, names.motifs)
{
  if (!is.null (names.motifs))
  {
    ##----checks names.motifs
    if(sum(names.motifs%in%rownames(motifstb)) == 0) 
      stop("-NOTE: 'names.motifs' should be in 
           '@results$motifsInformation!' \n")
    if(sum(names.motifs%in%rownames(motifstb)) != 
       length(names.motifs)) 
      stop ("Not all motifs names are available! \n")
    ##----
    motifstb <- motifstb[names.motifs, 
                         c("Regulon1","Regulon2", "R")]
  } else
  {
    motifstb <- motifstb[, c("Regulon1", "Regulon2", "R")]
  }
  motifstb[, 3] <- round(motifstb[, 3], 2)
  return(motifstb)
}
