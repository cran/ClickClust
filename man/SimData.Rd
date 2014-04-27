\name{click.sim}
\alias{click.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulating sequences of visited states }
\description{
Simulates sequences of visited states.
}
\usage{
click.sim(n, int = c(5, 100), alpha, beta = NULL, gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n }{ number of sequences }
  \item{int }{ interval defining the lower and upper bounds for the length of sequences }
  \item{alpha }{ vector of mixing proportions (length K)}
  \item{beta }{ matrix of initial state probabilities (K x p)}
  \item{gamma }{ array of K p x p transition probability matrices (p x p x K)}
}
\details{
Simulates 'n' sequences of visited states according to the following mixture model parameters: 'alpha' - mixing proportions, 'beta' - initial state probabilities, 'gamma' - transition probability matrices. If the matrix 'beta' is not provided, all initial states are assumed to be equal to 1 / p.
}
\value{
  \item{S }{ list of simulated sequences}
  \item{id }{ true classification of simulated sequences}
}
\references{
Melnykov, V. (2014) Model-based biclustering of clickstream data, under revision.
}
\author{ Melnykov, V.}

\seealso{click.read, click.EM}
\examples{


# SPECIFY MODEL PARAMETERS

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
A


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ Markov model }% __ONLY ONE__ keyword per line
