\name{click.var}
\alias{click.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Variance-covariance matrix estimation}
\description{
Estimates the variance-covariance matrix for model parameter estimates.
}
\usage{
click.var(X, y = NULL, alpha, beta = NULL, gamma, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset array (p x p x n) }
  \item{y }{ vector of initial states (length n) }
  \item{alpha }{ vector of mixing proportions (length K)}
  \item{beta }{ matrix of initial state probabilities (K x p)}
  \item{gamma }{ array of transition probabilities (p x p x K)}
  \item{z }{ matrix of posterior probabilities (n x K)}
}
\details{
Returns an estimated variance-covariance matrix for model parameter estimates.
}
\references{
Melnykov, V. (2014) Model-based biclustering of clickstream data, under revision.
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
C <- click.read(p, A$S)


# EM ALGORITHM

M2 <- click.EM(X = C$X, y = C$y, K = 2)


# VARIANCE ESTIMATION

V <- click.var(X = C$X, y = C$y, alpha = M2$alpha, beta = M2$beta,
               gamma = M2$gamma, z = M2$z)

# 95% confidence intervals for all model parameters

Estimate <- c(M2$alpha[-K], as.vector(t(M2$beta[,-p])),
              as.vector(apply(M2$gamma[,-p,], 3, t)))

Lower <- Estimate - qnorm(0.975) * sqrt(diag(V))
Upper <- Estimate + qnorm(0.975) * sqrt(diag(V))

cbind(Estimate, Lower, Upper)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ variance estimation }
\keyword{ EM algorithm }
\keyword{ Markov model }% __ONLY ONE__ keyword per line
