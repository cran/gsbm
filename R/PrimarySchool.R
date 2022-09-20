#' Network of interactions within a primary school in the course of a day
#' 
#' This network, collected and analyzed by J. Stehle et al. in "High-resolution measurements of face-to-face contact patterns in a primary school", records physical interactions between 226 children and 10 teachers within a primary school over the course of a day. The network data was collected using a system of sensors worn by the participants. This system records the duration of interactions between two individuals facing each other at a maximum distance of one and a half meters.
# 'The duration of these interactions varies between 20 seconds and two and a half hours. We consider that a physical interaction has been observed if the corresponding interaction duration is greater than one minute. If an interaction of less than one minute is observed, we consider that this observation may be erroneous, and treat the corresponding data as missing.
#'
#' @docType data
#'
#' @usage data(PrimarySchool)
#'
#' @format A list with 2 attributes:
#' \describe{
#'   \item{A}{adjacency matrix of the graph. A binary matrix encoding
#'   2490 undirected connections between 236 nodes, with 7054 missing entries}
#'   \item{class}{vector indicating the class of the node if the corresponding individual is a child, and otherwise that it belongs to the group of teachers.}
#' }
#' @source \url{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0023176}
"PrimarySchool"