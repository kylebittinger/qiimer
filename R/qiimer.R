#' qiimer: Read QIIME output files and create plots
#' 
#' @name qiimer
#' @docType package
#' @import pheatmap
NULL


#' Sample dataset from murine gut microbiome
#' 
#' The \code{relmbeta} dataset is taken from a mouse study where wild-type 
#' and RELMbeta knockout mice were fed either a normal or high-fat diet.  
#' The diet was observed to have a pronounced effect on the gut microbiome 
#' composition.  The genotype also had an effect, but less so.
#' 
#' The bacterial 16S rDNA gene was sequenced using 454 FLX technology, 
#' producing about 26k reads.  The reads were processed via the standard QIIME
#' workflow, \code{pick_de_novo_otus.py}, which clustered the reads into 776 
#' operational taxonomic units (OTUs).  QIIME version 1.7.0 was used for the 
#' analysis.
#' 
#' The OTU table was filtered to remove OTUs appearing in only one sample.
#' Following this, a single rarefaction was performed at a level of 500 reads
#' per sample.  The unweighted UniFrac distance was then computed for each pair 
#' of samples.
#' 
#' The \code{relmbeta} dataset contains a series of related data objects.
#' \itemize{
#'   \item The \code{relmbeta} data frame lists the sample IDs, genotypes, 
#'   and dietary assignments for the mice.
#'   
#'   \item The matrix \code{relmbeta_counts} contains the number of reads 
#'   observed in each OTU, after rarefaction to 500 reads per sample.  
#'   OTUs are listed in the rows and samples are listed in the columns.
#'   
#'   \item \code{relmbeta_dist} is an object of class \code{"dist"}, containing
#'   the unweighted UniFrac distances between samples.
#'   
#'   \item \code{relmbeta_alpha} is a data frame containing the shannon 
#'   diversity for each sample at a level of 10-500 sequences per sample.
#'   
#'   \item The character vector \code{relmbeta_assignments} contains taxonomic
#'   assignments for each OTU in the study.
#'
#'   \item The \code{relmbeta_biom} object is a representation of the BIOM file
#'   produced by QIIME.  The \code{biom} package must be installed to use this
#'   object.
#' }
#' 
#' @name relmbeta
#' @docType data
#' @keywords datasets
#' @format \code{relmbeta} is a data frame with 20 rows and 3 variables.
#'   \code{relmbeta_counts} is a matrix with 337 rows and 20 columns.
#'   \code{relmbeta_dist} is a \code{"dist"} object with 20 rows.
#'   \code{relmbeta_alpha} is a data frame with 2200 rows and 4 columns.
#'   \code{relmbeta_assignments} is a character vector of lengh 337.
#'   \code{relmbeta_biom} is a \code{"biom"} object.
#' @source Hildebrandt et al. High-fat diet determines the composition of the 
#'   murine gut microbiome independently of obesity. \emph{Gastroenterology} 
#'   \strong{137}, 1716 (2009).
NULL

