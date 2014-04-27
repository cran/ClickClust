\name{click.forward}
\alias{click.forward}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Forward search for equivalent states }
\description{
Runs forward search to detect blocks of equivalent states.
}
\usage{
click.forward(X, K, eps = 10^-10, r = 100, iter = 5, min.gamma = 10^-3,
  scale.const = 1.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset array (p x p x n) }
  \item{K }{ number of mixture components }
  \item{eps }{ tolerance level }
  \item{r }{ number of restarts for initialization}
  \item{iter }{ number of iterations for each short EM run}
  \item{min.gamma }{ lower bound for transition probabilities }
  \item{scale.const }{ scaling constant for avoiding numerical issues }
}
\details{
Runs forward search to detect blocks of equivalent states. States i and j are called equivalent if their behavior expressed in terms of transition probabilities is identical, i.e., the probabilities of leaving i and j to visit another state h are the same as well as the probabilities of coming to i and j from another state h are the same; this condition should hold for all mixture components.
Notation: p - number of states, n - sample size, K - number of mixture components, d - number of equivalence blocks.
}
\value{
  \item{z }{ matrix of posterior probabilities (n x K)}
  \item{alpha }{ vector of mixing proportions (length K)}
  \item{gamma }{ array of transition probabilities (d x d x K)}
  \item{states }{ detected equivalence blocks (length p)}
  \item{logl }{ log likelihood value}
  \item{BIC }{ Bayesian Information Criterion}
  \item{id }{ classification vector (length n)}
}
\references{
Melnykov, V. (2014) Model-based biclustering of clickstream data, under revision.
}
\author{ Melnykov, V.}

\seealso{backward.search, click.EM}
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


# FORWARD SEARCH

click.forward(X = C$X, K = 2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ forward search }
\keyword{ EM algorithm }
\keyword{ Markov model }% __ONLY ONE__ keyword per line
