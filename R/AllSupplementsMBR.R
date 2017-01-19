################################################################################
#################  Internal functions for RTNduals-methods  ###################
################################################################################
##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##it checks if the input regulatory elements are in the 'TNI' annotation
.checkRegel <- function(tni, regulatoryElements)
{
  #---check regulatoryElements
  tp <- sapply(colnames(tni@annotation), function(i)
  {
    sum(regulatoryElements%in%tni@annotation[, i])
  })
  colid <- names(tp[which.max (tp)])
  idx <- which(tni@annotation[, colid]%in%regulatoryElements)
  if(length(idx) < length(regulatoryElements))
  {
    warning("Not all 'regulatory elements' are available in the 'TNI' 
            annotation!" )
  }
  regulatoryElements <- tni@annotation[idx,]
  idx <- match(rownames(regulatoryElements), tni@transcriptionFactors)
  regulatoryElements <- tni@transcriptionFactors[idx]
  return (regulatoryElements)
  }

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##creates a matrix of regulons (tnet) for analysis (can be numeric or factors)
.regMatrix <- function(regulons, regel, getNames = TRUE, factors = FALSE)
{
  targets <- unique(unlist(lapply(regulons, names)))
  targets <- unique(c(regel, targets))
  xmat <- matrix(0, nrow=length(targets), ncol=length(regel))
  rownames(xmat) <- targets
  colnames(xmat) <- regel
  for(i in regel)
  {
    regs <- regulons[[i]]
    if(factors == TRUE) xmat[names(regs), i] <- as.character(regs)
    else {xmat[names(regs), i] <- regs}
  }
  if(getNames) colnames(xmat) <- names(regel)
  if (factors == TRUE) xmat <- as.data.frame(xmat, stringAsFactors=factors)
  return(xmat)
}

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##gets significant associations inferred for motifs using the 'prob' parameter
.motifsquantile <- function(regcor, th=0.99)
{
  th <- th*100
  cormat <- .cutoffquantile(regcor)
  pmat <- cormat>th
  coord <- which(pmat,arr.ind=TRUE)
  rnames <- rownames(pmat)[coord[, 1]]
  cnames <- colnames(pmat)[coord[, 2]]
  corvalues <- regcor[coord]
  qvalues <- cormat[coord]
  qvalues <- data.frame(Regulon1=rnames,Regulon2=cnames,R=corvalues, 
                        Quantile = qvalues, stringsAsFactors=FALSE)
  rownames(qvalues) <- paste(qvalues$Regulon1, qvalues$Regulon2, sep="~")
  return(qvalues)
}

#internal function '.motifsquantile'
.cutoffquantile <- function (regcor)
{
  cormat <- abs(regcor)
  cormat <- as.numeric(cut(cormat, breaks=quantile(cormat, (0:100)/100)))
  cormat <- matrix(cormat, ncol=ncol(regcor), nrow=nrow(regcor), 
                   dimnames=list(rownames(regcor), colnames(regcor)))
  return(cormat)
}

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##it gets the jaccard information between the 'duals'
.jcOverlap <- function(pvlist, regel, tnet)
{
  #-----
  regEle1 <- unique(pvlist[, "Regulon1"])
  regEle1 <- regel[regEle1]
  
  regEle2 <- unique(pvlist[, "Regulon2"])
  regEle2 <- regel[regEle2]
  
  #----(jaccard)
  ##jcagree <- .jc.overlap(regEle1, regEle2, tnet, overlap="agreement")
  ##jcdisagree <- .jc.overlap(regEle1, regEle2, tnet, overlap="disagreement")
  jcall <- .jc.overlap(regEle1, regEle2, tnet, overlap="all")
  if(!is.null(dim(jcall)))
  {
    tb <- as.matrix(pvlist[, c("Regulon1", "Regulon2")]) 
    jcinf <- apply(tb, 1, function(x)
    {
      reg <- x[1]
      reg <- regel[reg]
      t <- x[2]
      t <- regel[t]
      jcAll <- jcall[t, reg]
    })
    
  }
  else
  {
    jcinf <- jcall
  }
  pvlist <- cbind(pvlist, Jaccard.coefficient=jcinf)
  return(pvlist)
}

#internal function '.jcOverlap'
.jc.overlap <- function(regEle1, regEle2, tnet, overlap="all"){
  #overlap %in% c("all","agreement","disagreement")
  if(overlap == "all"){
    tnet[tnet != 0] <- 1
    jc <- function(x, xmat){
      c <- x+xmat
      a <- colSums(c == 2)
      b <- colSums(c > 0)
      b[b == 0] <- 1
      a/b
    }
  } else {
    if(overlap == "agreement"){
      ov <- (+1)
    } else {
      ov <- (-1)
    }
    tnet [tnet > 0] <- 1; tnet[tnet < 0] <- -1
    jc <- function(x, xmat){
      c <- x*xmat
      a <- colSums(c == ov)
      c <- abs(x) + abs(xmat)
      b <- colSums(c != 0)
      b[b == 0] <- 1
      a/b
    }
  }
  amap <- apply(tnet[, regEle1, drop=FALSE], 2, jc, 
                xmat=tnet[, regEle2, drop=FALSE])
  return(amap)
}


##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##gets mutual information between 'duals' reg. elements 
##and the p-value when available
.getMIorP <- function(pvlist, regel, tbmi, size1, size2, 
                      mutualInformation=TRUE, cutoff=0)
{
  tb <- as.matrix(pvlist[, c("Regulon1", "Regulon2")])
  mutinf <- apply(tb, 1, function(x)
  {
    reg <- x[1]
    reg <- regel[reg]
    t <- x[2]
    t <- regel[t]
    s1 <- size1[reg]
    s2 <- size2[t]
    mi <- tbmi[t, reg]
    mi <- abs(mi)
    cbind(mi, round(s1), round(s2))
  })
  if(mutualInformation)
  {
    pvlist <- cbind(pvlist, MI=mutinf[1, ], Size.Regulon1=mutinf[2, ],
                    Size.Regulon2=mutinf[3, ])
    pvlist <- pvlist[, c("Regulon1", "Size.Regulon1", "Regulon2", 
                         "Size.Regulon2", "MI", "R", "Quantile")]
    ##----MI
    pvlist <- pvlist[which(pvlist$MI != 0), ]
  }
  else
  {
    pvlist <- cbind(pvlist, MI.Adjusted.Pvalue=mutinf[1, ])
    pvlist$MI.Adjusted.Pvalue[
      pvlist$MI.Adjusted.Pvalue<cutoff] <- paste("<", cutoff, sep="")
    pvlist <- pvlist[, c("Regulon1", "Size.Regulon1", "Regulon2", 
                         "Size.Regulon2", "MI", "MI.Adjusted.Pvalue", 
                         "R", "Quantile")]
  }
  return (pvlist)
}


##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##it takes the correlation between regulons
.tni.cor<-function(x, tnet, estimator="pearson",dg=0, asInteger=TRUE, 
                   mapAssignedAssociation=TRUE){
  tfs<-colnames(tnet)
  tar<-rownames(tnet)
  ids<-unique(c(tfs,setdiff(tar,tfs)))
  x=x[ids,]
  x=t(x)
  #--
  pcorm=cor(x[,tfs],x[,tar], method=estimator,use="complete.obs")
  if(asInteger){
    pcorm[pcorm<0]=-1
    pcorm[pcorm>0]=1
  }
  if(length(tfs)>1)diag(pcorm[,tfs])=dg
  #--
  pcorm<-t(pcorm)
  colnames(pcorm)<-tfs
  if(mapAssignedAssociation)pcorm[tnet==0]=0
  pcorm
}

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##This function takes two gene sets (i.e. regulons), a vector 
##containing the size of the gene universe, and compute the number of 
##genes expected to occur in both regulons, the actual observed overlap, 
##and the pvalue from a hypergeometric test.
.mbr.hyper <- function(pvlist, regulons, regel, universe, pAdjustMethod, 
                       verbose=TRUE)
{
  if(verbose) cat("-Hypergeometric analysis...\n")
  regpairstb <- pvlist[, c("Regulon1", "Regulon2")]
  if(verbose) pb<-txtProgressBar(style=3)
  res <- NULL
  for(i in 1:nrow(regpairstb))
  {
    if(verbose) setTxtProgressBar(pb, i/nrow(regpairstb))
    vecpairs <- as.character(regpairstb[i, ])
    ##---
    reg1 <- vecpairs[1]
    reg1 <- regel[reg1]
    reg2 <- vecpairs[2]
    reg2 <- regel[reg2]
    ##---
    regulon1 <- names(regulons[[reg1]])
    regulon2 <- names(regulons[[reg2]])
    tmp <- .regulon.hyper(regulon1=regulon1, universe=universe, 
                          regulon2=regulon2)
    res <- rbind(res, tmp)
  }
  if(verbose) close(pb)
  results <- cbind(regpairstb, res)
  adjPvals <- p.adjust(results[, "Pvalue"], method = pAdjustMethod)
  results <- cbind(results, adjPvals)
  colnames(results)[ncol(results)] <- "Adjusted.Pvalue"
  results <- results[order(results[, "Pvalue"]), , drop=FALSE]
  
}

#internal function for '.mbr.hyper'
##it makes the hypergeometric test.
.regulon.hyper <- function(regulon1, universe, regulon2) {
  ##number of genes in universe
  N <- length(universe)			
  ##remove genes from gene set that are not in universe			
  regulon1 <- intersect(regulon1, universe) 
  regulon2 <- intersect(regulon2, universe)
  ##size of gene set	
  m <- length(regulon1) 							
  Nm <- N-m	
  ##regulon2 in gene set
  overlap <- intersect(regulon1, regulon2) 	
  ##number of hits between regulons		
  k <- length(overlap) 							
  n <- length(regulon2)	
  HGTresults <- phyper(k-1, m, Nm, n, lower.tail = FALSE)
  ex <- (n/N)*m
  if(m == 0 | n == 0) HGTresults <- 1
  hyp.vec <- c(N, m, n, ex, k, HGTresults)
  names(hyp.vec) <- c("Universe.Size", "R1.Size", "R2.Size", 
                      "Expected.Overlap", "Observed.Overlap", "Pvalue")
  return(hyp.vec)
}
##------------------------------------------------------------------------------
#internal function for 'tni2mbrPreprocess'
##it takes two 'TNIs' objects produced separetely and creates a 'MBR' object
.combineTNIs <- function (tni1, tni2, verbose = TRUE)
{
  ## checks gexp consistency
  v1 <- as.numeric(tni1@gexp)
  v2 <- as.numeric(tni2@gexp)
  if (verbose)
    cat("-Checking expression matrix consistency...\n")
  if (sum(v1 - v2) != 0)
    stop("The TNIs should use the same expression matrix.")
  
  ## checks parameter consistency
  if (verbose)
    cat("-Checking parameter consistency...\n")
  if(!all(tni1@para %in% tni2@para))
  {
    ## gives feedback on which parameters are wrong
    if (verbose)
    {
      idx <- which(!(tni1@para %in% tni2@para))
      
      for (i in idx)
      {
        p1 <- unlist(tni1@para[[i]])
        p2 <- unlist(tni2@para[[i]])
        
        ind <- which(!(p1 %in% p2))
        
        if (i == idx[1])
          cat ("- Parameter differences:\n")
        
        for (j in ind)
        {
          line <- c(p1[j], p2[j])
          print(line)
        }
      }
    }
    stop("TNIs were not computed using the same parameters.")
  }
  if (verbose)
    cat("-Checking whether regulatory elements are unique to each TNI...\n")
  if (any(tni1@transcriptionFactors %in% tni2@transcriptionFactors))
    stop ("There are regulatory elements listed in both TNIs.")
  
  ## checks whether both TNIs have undergone all methods in RTN from
  ## Permutation to DPI filter
  if (verbose)
    cat("-Checking if all TNI methods are completed...\n")
  
  if(any(tni1@status[1:4] != "[x]") || any(tni2@status[1:4] != "[x]"))
  {
    ## gives feedback on which methods were not run
    if (verbose)
    {
      cat("TNI1: ")
      print(tni1@status)
      cat("TNI2: ")
      print(tni2@status)
    }
    stop("Both TNIs must be evaluated by the RTN pipeline up to 
         the DPI filter.")
  }
  }

##------------------------------------------------------------------------------
#internal function for 'mbrDuals'
##it checks the consistency of supplementary.table
.consisSuppTable <- function(object, supplementary.table, evidenceColname, 
                             verbose)
{
  ##-----checks
  if(!evidenceColname%in%colnames(supplementary.table)) 
    stop("'evidenceColname' should be a valid colname in 'supplementary.table", 
         call.=FALSE)
  ##---
  idx <- which(colnames(supplementary.table)%in%evidenceColname)
  if(idx!=3) 
    stop("'evidenceColname' should be the third column in 'supplementary.table", 
         call.=FALSE)
  ##---
  motifsInformation <- mbrGet(object, what="motifsInformation")
  colnms <- colnames(motifsInformation)
  if(evidenceColname%in%colnms)
  {
    cat("-This evidence was already computed, overwrite the information...\n")
  }
  ##-----
  ##-----Calcules the consistency
  idx <- which(colnames(supplementary.table)%in%evidenceColname)
  tni1 <- mbrGet(object, what="TNI1"); annot <- tni1@annotation 
  annot <- as.matrix(annot)
  tmp <- as.matrix(supplementary.table); ttmp <- tmp[,-idx]
  colnames(tmp) <- colnames(supplementary.table)
  ##---consistency between supplementary.table and annotation
  consc <- (sum(ttmp%in%annot)/prod(dim(ttmp)))*100
  if(consc<90) 
    warning(paste("Only",paste(round(consc,2),"%",sep=""), 
                  "of 'supplementary.table' is listed in the annotation!\n"), 
            call.=FALSE)
  if(consc>90 & verbose) 
    cat(paste("-",paste(round(consc,2),"%",sep=""), 
              "of 'supplementary.table' is listed in the annotation!\n"))
  ##-----
  return(tmp)
}
##------------------------------------------------------------------------------
#internal function for 'mbrDuals'
##it checks the evidences for 'duals' in 'supplementary.table'
.checkLoops <- function (object, supplementary.table, evidenceColname, 
                         verbose=TRUE)
{
  motifsInformation <- mbrGet(object, what="motifsInformation")
  motifsInformation[, evidenceColname] <- NA
  ##---
  if(verbose) 
    cat("-Checking whether evidences in 'supplementary.table' support 
        the existence of the inferred duals...\n")
  if(verbose)pb<-txtProgressBar(style=3)
  x <- 0
  ##---
  regs1 <- mbrGet(object, what="testedElementsTNI1")
  regs2 <- mbrGet(object, what="testedElementsTNI2")
  regs <- c(regs1,regs2)
  for (i in 1:length(regs1))
  {
    rg1 <- regs1[i]
    idx1 <- which((supplementary.table[,1] %in% rg1) | 
                    (supplementary.table[,1] %in% names(rg1)))
    idx2 <- which((supplementary.table[,2] %in% rg1) | 
                    (supplementary.table[,2] %in% names(rg1)))
    tpev <- supplementary.table[c(idx1,idx2),];
    for(j in 1:length(regs))
    {
      rg <- regs[j]
      idx <- which((tpev%in%rg) | (tpev%in%names(rg)))
      tpev[idx] <- names(rg)
    }
    if(!is.null(dim(tpev)))
    {
      tpev <- rbind(tpev, tpev[, c(2,1,3)])
    }else
    {
      tpev <- rbind(tpev, tpev[c(2,1,3)]) 
    }
    nmsDuals <- paste(tpev[,1], tpev[,2],sep="~")
    rownames(tpev) <- nmsDuals
    nms <- rownames(motifsInformation)
    ids <- nms[nms%in%nmsDuals] 
    motifsInformation[ids,evidenceColname] <- tpev[ids,evidenceColname]
    x <- x+1
    if(verbose)setTxtProgressBar(pb, x/length(regs1))
  }
  if(verbose)close(pb)
  object <- .mbr.set(name="motifsInformation", para=motifsInformation, object=object)
  return (object)
}
##------------------------------------------------------------------------------
#".mbr.set" internal function
##it setts the slots of a MBR object
.mbr.set <- 
    function(name, para, object)
    {
        if(name=="para")
        {
            object@para <- para
        }
        else if(name=="summary")
        {
            object@summary <- para
        }
        else if(name=="status")
        {
            object@status <- para
        }
        else if(name=="TNI1")
        {
            object@TNI1 <- para
        }
        else if(name=="TNI2")
        {
            object@TNI2 <- para
        }
        else if(name=="statusUpdate")
        {
            object@status[para] <- "[x]"
        }
        else if(name=="testedElementsTNI1")
        {
            object@testedElementsTNI1 <- para
        }
        else if(name=="testedElementsTNI2")
        {
            object@testedElementsTNI2 <- para
        }
        else if(name=="dualRegulons")
        {
            object@dualRegulons <- para
        }
        else if(name=="motifsInformation")
        {
            object@results$motifsInformation <- para
        }
        else if(name=="hypergeometricResults")
        {
            object@results$hypergeometricResults <- para
        }
        
        return(object)
    }
