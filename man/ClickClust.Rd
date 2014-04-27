\name{click.EM}
\alias{click.EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ EM algorithm for mixtures of Markov models }
\description{
Runs the EM algorithm for finite mixture models with Markov model components.
}
\usage{
click.EM(X, y = NULL, K, eps = 10^-10, r = 100, iter = 5, min.beta = 10^-3,
  min.gamma = 10^-3, scale.const = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset array (p x p x n) }
  \item{y }{ vector of initial states (length n) }
  \item{K }{ number of mixture components }
  \item{eps }{ tolerance level }
  \item{r }{ number of restarts for initialization}
  \item{iter }{ number of iterations for each short EM run}
  \item{min.beta }{ lower bound for initial state probabilities }
  \item{min.gamma }{ lower bound for transition probabilities }
  \item{scale.const }{ scaling constant for avoiding numerical issues }
}
\details{
Runs the EM algorithm for finite mixture models with first order Markov model components. The function returns estimated mixing proportions 'alpha' and transition probabilty matrices 'gamma'. If initial states 'y' are not provided, initial state probabilities 'beta' are not estimated and assumed to be equal to 1 / p. In this case, the total number of estimated parameters is given by M = K - 1 + K * p * (p - 1). Otherwise, initial state probabilities 'beta' are also estimated and the total number of parameters is M = K - 1 + K * (p - 1) + K * p * (p - 1).
Notation: p - number of states, n - sample size, K - number of mixture components, d - number of equivalence blocks.
}
\value{
  \item{z }{ matrix of posterior probabilities (n x K)}
  \item{id }{ classification vector (length n)}
  \item{alpha }{ vector of mixing proportions (length K)}
  \item{beta }{ matrix of initial state probabilities (K x p)}
  \item{gamma }{ array of transition probabilities (p x p x K)}
  \item{logl }{ log likelihood value}
  \item{BIC }{ Bayesian Information Criterion}
}
\references{
Melnykov, V. (2014) Model-based biclustering of clickstream data, under revision.
}
\author{ Melnykov, V.}

\seealso{click.plot, click.forward, click.backward}
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


# EM ALGORITHM (without initial state probabilities)

N2 <- click.EM(X = C$X, K = 2)
N2$BIC


# EM ALGORITHM (with initial state probabilities)

M2 <- click.EM(X = C$X, y = C$y, K = 2)
M2$BIC

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ Markov model }
\keyword{ click-plot}% __ONLY ONE__ keyword per line
