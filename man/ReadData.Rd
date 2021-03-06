\name{click.read}
\alias{click.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reading sequences of visited states }
\description{
Prepares sequences of visited states for running the EM algorithm.
}
\usage{
click.read(S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S }{ list of numeric sequences }
}
\details{
Prepares sequences of visited states for running the EM algorithm by means of the click.EM() function.
}
\value{
  \item{X }{ dataset array (p x p x n) (p - # of states, n - # of sequences) }
  \item{y }{ vector of initial states (length n) }
}
\references{
Melnykov, V. (2016) Model-Based Biclustering of Clickstream Data, Computational Statistics and Data Analysis, 93, 31-45.

Melnykov, V. (2016) ClickClust: An R Package for Model-Based Clustering of Categorical Sequences, Journal of Statistical Software, 74, 1-34.
}
\author{ Melnykov, V.}

\seealso{click.sim, click.EM}
\examples{

set.seed(123)

n.seq <- 20

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
C <- click.read(A$S)
C$X
C$y

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ Markov model }% __ONLY ONE__ keyword per line
