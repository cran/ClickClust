\name{ClickClust-package}
\alias{ClickClust-package}
\alias{ClickClust-package}
\docType{package}
\title{
Model-based clustering of categorical sequences 
}
\description{
The package runs finite mixture modeling and model-based clustering for categorical sequences
}
\details{
\tabular{ll}{
Package: \tab ClickClust\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-04\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab no\cr
}
Function 'click.EM' runs the EM algorithm for finite mixture models with Markov model components.

}
\author{ Volodymyr Melnykov

Maintainer: Volodymyr Melnykov <vmelnykov@cba.ua.edu>
}
\references{
Melnykov, V. (2014) Model-based biclustering of clickstream data, under revision.
}
\examples{


set.seed(123)

n.seq <- 50

p <- 5
K <- 2
mix.prop <- c(0.3, 0.7)


TP1 <- matrix(c(0.20, 0.10, 0.15, 0.15, 0.40,
                0.20, 0.20, 0.20, 0.20, 0.20,
                0.15, 0.10, 0.20, 0.20, 0.35,
                0.15, 0.10, 0.20, 0.20, 0.35,
                0.30, 0.30, 0.10, 0.10, 0.20), byrow = TRUE, ncol = p)

TP2 <- matrix(c(0.15, 0.15, 0.20, 0.20, 0.30,
                0.20, 0.10, 0.30, 0.30, 0.10,
                0.25, 0.20, 0.15, 0.15, 0.25,
                0.25, 0.20, 0.15, 0.15, 0.25,
                0.10, 0.30, 0.20, 0.20, 0.20), byrow = TRUE, ncol = p)


TP <- array(rep(NA, p * p * K), c(p, p, K))
TP[,,1] <- TP1
TP[,,2] <- TP2


# DATA SIMULATION

A <- click.sim(n = n.seq, int = c(10, 50), alpha = mix.prop, gamma = TP)
C <- click.read(p, A$S)


# EM ALGORITHM

click.EM(X = C$X, K = 2)


}
