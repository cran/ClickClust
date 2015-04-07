\name{click.predict}
\alias{click.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prediction of future state visits}
\description{
Calculates the transition probability matrix associated with the M-step transition.
}
\usage{
click.predict(M = 1, gamma, pr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M }{ number of transition steps (M = 1 by default) }
  \item{gamma }{ array of transition probabilities (p x p x K)}
  \item{pr }{ vector of probabilities associated with components (length K)}
}
\details{
Returns a transition probability matrix associated with the M-step transition. If the vector pr is not specified, all components are assumed equally likely.
}
\references{
Melnykov, V. (2014) Model-based biclustering of clickstream data, accepted by Computational Statistics & Data Analysis.
}
\author{ Melnykov, V.}

\seealso{click.EM}
\examples{

set.seed(123)

n.seq <- 200

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


# EM ALGORITHM

M2 <- click.EM(X = C$X, y = C$y, K = 2)


# Assuming component probabilities given by mixing proportions, predict the next state 

click.predict(M = 1, gamma = M2$gamma, pr = M2$alpha)

# For the last location in the first sequence, predict the three-step transition
# location, given corresponding posterior probabilities

click.predict(M = 3, gamma = M2$gamma, pr = M2$z[1,])[A$S[[1]][length(A$S[[1]])],]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prediction }
\keyword{ EM algorithm }
\keyword{ Markov model }% __ONLY ONE__ keyword per line
