\name{click.plot}
\alias{click.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot of the obtained clustering solution }
\description{
Constructs a click-plot for the clustering solution.
}
\usage{
click.plot(X, y = NULL, file = NULL, id, margs = c(0.1, 0.1, 0.1, 0.1),
  font.cex = 2, font.col = "black", cell.cex = 1, cell.lwd = 1.3,
  cell.col = "red", sep.lwd = 1.3, sep.col = "magenta", obs.lwd = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset array (p x p x n) }
  \item{y }{ vector of initial states (length n) }
  \item{file }{ name of the output pdf-file }
  \item{id }{ estimated classification vector (length n) }
  \item{margs }{ plot margins }
  \item{font.cex }{ magnification of labels }
  \item{font.col }{ color of labels }
  \item{cell.cex }{ magnification of cells }
  \item{cell.lwd }{ width of cell frames }
  \item{cell.col }{ color of cell frames }
  \item{sep.lwd }{ width of separator lines }
  \item{sep.col }{ color of separator lines }
  \item{obs.lwd }{ width of observation lines }
}
\details{
Constructs a click-plot for the provided clustering solution. If the parameter 'file' is specified, the constructed plot will be saved in the pdf-file with the name 'file'.
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


# CONSTRUCT CLICK-PLOT

click.plot(X = C$X, y = C$y, file = NULL, id = M2$id)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ click-plot }
\keyword{ EM algorithm }
\keyword{ Markov model }% __ONLY ONE__ keyword per line
