
##------------------------------------------------------------------------------
##This function is used for argument checking
mbr.checks <- function(name, para)
{
  if(name == "TNI")
  {
    if(!is.null(para))
        stop("'TNIs' should be NULL")
  }
  if(name=="testedElementsTNI")
  {
    if((!is(para, "character") && !length(para)==0))
        stop("'testedElementsTNIs' should be an empty 'character' vector")
  }
  if(name=="dualRegulons")
  {
    if((!is(para, "character") && !length(para)==0))
        stop("'dualRegulons' should be an empty 'character' vector")
  }
  if(name=="results")
  {
    if((!is(para, "list") && !length(para)==0))
        stop("'results' should be an empty 'list'")
  }
  if(name=="para")
  {
    if((!is(para, "list") && !length(para)==0))
        stop("'para' should be an empty 'list'")
  }
  if(name=="summary")
  {
    if((!is(para, "list") && !length(para)==0))
        stop("'summary' should be an empty 'list'")
  }
  if(name=="status")
  {
    if((!is(para, "character") && !length(para)==0))
        stop("'status' should be an empty 'character' vector")
  }
  if(name == "gexp")
  {
    if(!is.matrix(para) || !is.numeric(para[1, ]))
      stop("'gexp' should be a numeric matrix with genes on rows and 
           samples on cols!", call.=FALSE)
    if(is.null(rownames(para)) || is.null(colnames(para)) || 
       length(unique(rownames(para))) < length(rownames(para)) || 
       length(unique(colnames(para))) < length(colnames(para)))
      stop("the 'gexp' matrix should be named on rows and cols (unique names)", 
           call. = FALSE)
  }
  
  ##---
  else if(name == "regulatoryElements1" || name == "regulatoryElements2")
  {
    if(!(is.character(para) || is.numeric(para)) || 
       any(is.na(para)) || any(para == "" || is.null(para)))
      stop("'regulatoryElements1 and 2' should be a character vector, 
           without 'NA' or empty names!", call.=FALSE)
    if(length(unique(para)) < length(para))
      stop("'regulatoryElements1 and 2' should have unique identifiers!",
           call.=FALSE)
  }
  
  ##---
  else if(name=="object")
  {
    if(class(para) != 'MBR')
      stop("'object' should be a 'MBR' class object", call.=FALSE)
  }
  
  ##---
  else if(name=="verbose")
  {
    if(!is.logical(para))
      stop("'verbose' should be a logical value!", call.=FALSE)
  }
  
  ##---
  else if(name == "minRegulonSize")
  {
    if(!is.singleNumber(para) || !para>0)
      stop("'minRegulonSize' should be numeric and >0", call.=FALSE)
  }
  
  ##---
  else if(name == "prob")
  {
    if(!is.singleNumber(para) || (!para>=0) && (!para<1))
      stop("'prob' should be a numeric value >=0 and <1!", call.=FALSE)
  }
  
  ##---
  else if(name == "estimator")
  {
    if(!para %in% c("spearman", "kendall", "pearson"))
      stop("'estimator' should be one of 'spearman', 'kendall', 'pearson' !", 
           call.=FALSE)
  }
  
  ##---
  else if(name == "regulatoryElements")
  {
    if(sum(duplicated(para)) > 0)
    {
      dupl <- para[which(duplicated(para))]
      dupl <- paste(dupl, " ", sep = "")
      stop("All regulatory elements should be unique!
           The following regulatory elements are duplicated:\n", 
           dupl, call. = FALSE)
    }
  }
  
  ##---
  else if(name=="numberRegElements")
  {
    if(length(para) < 3)
      stop("At least 3 regulatory elements (regulatoryElement1 + regulatoryElement2) 
           need to be inputted!")
  }
  
  ##---
  else if(name=="pAdjustMethod") 
  {
    if(!is.character(para) || length(para)!=1 || 
       !(para %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", 
                     "fdr", "none")))
      stop("'pAdjustMethod' should be any one of 'holm','hochberg','hommel',
           'bonferroni', 'BH','BY','fdr' and 'none'!",call.=FALSE)
  }
  
  ##---
  else if(name == "supplementary.table")
  {
    if(!class(para) == "data.frame" || !ncol(para)==3 || !nrow(para)>=1 || 
       is.null(dim(para)))
      stop("'supplementary.table' should be a 'data.frame' class object with 
           3 columns (Regulon1, Regulon2, Evidence)!", call.=FALSE)
  }
  
  ##---
  else if(name == "evidenceColname")
  {
    if(!is.character(para) || length(para)>1 || is.null(para) || 
       !is.singleString(para))
      stop("'evidenceColname' should be a character value present in colnames 
           of supplementary.table!", call.=FALSE)
  }
  
  ##---
  else if(name=="uniqueInput")
  {
    nms <- paste(para[,1], para[,2],sep="~")
    dupliNms <- c(paste(para[,1],para[,2],sep="~"),
                  paste(para[,2], para[,1],sep="~"))
    if(sum(duplicated(dupliNms))>0)
    {
      dupliNms <- dupliNms[duplicated(dupliNms)]
      ids <- unique(nms[nms%in%dupliNms])
      stop(c(paste("The possible pairs in 'supplementary.table' should be unique!", 
                   "The following pairs are duplicated\n"), paste(ids, " ")), 
           call.=FALSE)
    }
  }
  
  ##---
  else if(name == "tni")
  {
    if((class(para) != 'TNI') || is.null(para))
    {
      stop("'tni1' and 'tni2' should be TNI-class objects!", call.=FALSE)
    }
  }
  
  ##---
  else if(name=="mbr_get")
  {
    opts <- c("TNI1", "TNI2", "testedElementsTNI1", "testedElementsTNI2", 
              "dualRegulons", "results", "para", "summary", "status", 
              "motifsInformation", "hyperResults")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", 
                 paste(opts,collapse = ", ") ) ,call.=FALSE)
  }
  
  }

##------------------------------------------------------------------------------
.txtcollapse<-function(vec){
  paste("'",paste(vec[-length(vec)], collapse = "', '"),
        "'"," and '",vec[length(vec)],"'!", sep=""
  )
}

##------------------------------------------------------------------------------
is.singleNumber<-function(para){
  (is.integer(para) || is.numeric(para)) && length(para)==1L && !is.na(para)
}
is.singleInteger<-function(para){
  lg <- (is.integer(para) || is.numeric(para)) && length(para)==1L && 
    !is.na(para)
  if(lg) lg <- (para / ceiling(para)) == 1
  return(lg)
}
is.singleString<-function(para){
  is.character(para) && length(para) == 1L && !is.na(para)
}
is.singleLogical<-function(para){
  is.logical(para) && length(para) == 1L && !is.na(para)
}
all.binaryValues<-function(para){
  all( para %in% c(0, 1, NA) )
}
all.integerValues<-function(para){
  lg <- ( all(is.integer(para)) || all(is.numeric(para)) ) && 
    !any(is.na(para))
  if(lg) lg <- all ( (para / ceiling(para)) == 1 )
  return(lg)
}
all.characterValues<-function(para){
  all(is.character(para)) && !any(is.na(para))
}


