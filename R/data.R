#' Example community data
#'
#' A data frame with site names as row names, species names as column names,
#' cells are the abundance of each species at each site.
#'
#' @docType data
#' @keywords datasets
#' @name comm_a
#' @format A data frame with 15 sites and 15 species.
"comm_a"

#' Example community data
#'
#' A data frame with site names as row names, species names as column names,
#' cells are the abundance of each species at each site.
#'
#' @docType data
#' @keywords datasets
#' @name comm_b
#' @format A data frame with 15 sites and 9 species.
"comm_b"

#' Example phylogeny
#'
#' A phylogeny with more species than the community data.
#'
#' @docType data
#' @keywords datasets
#' @name phylotree
#' @format Newick format.
"phylotree"

#' Example environmental data
#'
#' A data frame of site environmental variables.
#'
#' @docType data
#' @keywords datasets
#' @name envi
#' @format A data frame with 15 sites and 4 variables: sand proportion,
#' canopy shade proportion, precipitation, and minimum temperature.
"envi"

#' Example species traits data
#'
#' A data frame of species functional traits.
#'
#' @docType data
#' @keywords datasets
#' @name traits
#' @format A data frame with 18 species and 3 variables: sla,
#' height, and seed dispersal mode.
"traits"

#' Phylogeny and community data from an Oldfield ecosystem in Southern Ontario, Canada
#'
#' A list containing a phylogeny for XX species of Oldfield forbs, as well as a
#' presence / absence dataset for their occurrence across several locations in
#' Southern Ontario see Dinnage (2009) for details. Sites each had two plots which
#' experienced a different treatment each; either they has been disturbed (ploughed 
#' 1 or 2 years previously), or they were a control plot (undisturbed in recent records).
#'
#' @format A list with two elements:
#' \describe{
#'   \item{\code{phy}}{A phylogeny in \code{ape}'s \code{phy} format}
#'   \item{\code{data}}{A data.frame containing data on the occurrence of the species in \code{phy}}
#' }
#' oldfield$data is a data.frame with 1786 rows, and the following 7 columns:
#' \describe{
#'   \item{\code{site_orig}}{integer. Site ID number.}
#'   \item{\code{habitat_type}}{character. Plot treatment: disturbed or undisturbed.}
#'   \item{\code{sp}}{character. Species name using underscore to separate binomial names (to match phylogeny).}
#'   \item{\code{abundance}}{integer. Recorded abundance of species in plot.}
#'   \item{\code{disturbance}}{integer. Whether the plot was disturbed or not. 0 or 1. 0 for undisturbed, 1 for disturbed}
#'   \item{\code{site_orig}}{character. A unique site descriptor concatenating the site number with the disturbance treatment.}
#'   \item{\code{pres}}{integer. Species presence or absence in plot. 0 or 1. 0 for absent, 1 for present}
#' }
"oldfield"
