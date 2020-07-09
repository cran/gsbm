## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load data----------------------------------------------------------------
library(Matrix)
library(gsbm)
library(igraph)
library(RColorBrewer)

data(PrimarySchool)
summary(PrimarySchool)

