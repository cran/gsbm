#' Character network from "Les miserables" novel
#'
#' A dataset containing Les Misérables characters network, encoding interactions
#' between characters of Victor Hugo's novel. Two characters are connected
#' whenever they appear in the same chapter. This dataset was first created
#' by Donald Knuth as part of the Stanford Graph Base.
#' (https://people.sc.fsu.edu/~jburkardt/datasets/sgb/sgb.html). It contains 77 nodes corresponding to characters of the novel, and 254
#'
#' @docType data
#'
#' @usage data(les_miserables)
#'
#' @format A list with 2 attributes:
#' \describe{
#'   \item{A}{ adjacency matrix of the graph. A binary matrix encoding
#'   254 connections between 77 nodes}
#'   \item{names}{ a vector giving the names of the characters corresponding to the
#'   rows and columns of the adjacency matrix}
#' }
#' @source \url{https://people.sc.fsu.edu/~jburkardt/datasets/sgb/sgb.html}
"les_miserables"
