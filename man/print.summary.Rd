\name{print.object}
\alias{print.EM}
\alias{summary.EM}
\alias{print.search}
\alias{summary.search}
\title{Functions for Printing or Summarizing Objects}
\description{
  \code{EM} and \code{search} classes  for printing and summarizing objects.
}
\usage{
\method{print}{EM}(x, ...)
\method{summary}{EM}(object, ...)
\method{print}{search}(x, ...)
\method{summary}{search}(object, ...)
}
\arguments{
  \item{x}{an object with the 'EM' (or 'search') class attributes.}
  \item{object}{an object with the 'EM' (or 'search') class attributes.}
  \item{\dots}{other possible options.}
}
\details{
  Some useful functions for printing and summarizing results.
}
\references{
Melnykov, V. (2016) Model-Based Biclustering of Clickstream Data, Computational Statistics and Data Analysis, 93, 31-45.

Melnykov, V. (2016) ClickClust: An R Package for Model-Based Clustering of Categorical Sequences, Journal of Statistical Software, 74, 1-34.
}
\author{ Melnykov, V.}

\seealso{
  \code{click.EM}.
}

