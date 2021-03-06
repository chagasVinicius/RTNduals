% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{MBR-class}
\alias{MBR}
\alias{MBR-class}
\title{MBR objects}
\description{
MBR: an S4 class for co-regulation analysis and inference of 'dual 
regulons'.
}
\details{
The MBR class is a container for results from the MBR methods. 
The class slots are used to store information of different transcriptional networks, 
regulator annotation, infered 'dual regulons' and parameters used in the analysis. 
All the information is stored in nine slots.
}
\section{Slots}{

\describe{
\item{\code{TNI1}}{a 'TNI' object created by the RTN package.}

\item{\code{TNI2}}{another 'TNI' object created by the RTN package.}

\item{\code{testedElementsTNI1}}{regulatory elements listed in the TNI1.}

\item{\code{testedElementsTNI2}}{regulatory elements listed in the TNI2.}

\item{\code{dualRegulons}}{all possible 'duals regulons' computed by 
\code{\link[RTNduals:mbrAssociation]{mbrAssociation}}}

\item{\code{results}}{a list, results from the MBR methods.}

\item{\code{para}}{a list, parameters used in the MBR methods.}

\item{\code{summary}}{a list, summary for 'para' and 'results'.}

\item{\code{status}}{a character vector specifying the status of the MBR object 
based on the available methods.}
}}
\section{Constructor}{


There are two constructors to create an MBR object, users can opt for one of the following: 
(1) \code{\link[RTNduals:mbrPreprocess]{mbrPreprocess}}; 
(2) \code{\link[RTNduals:tni2mbrPreprocess]{tni2mbrPreprocess}}.
\itemize{
\item {(1)}: {It is used to create an MBR object without any pre-computed transcriptional network.}
\item {(2)}: {It is used to create an MBR object using available transcriptional networks.} 
}
}

