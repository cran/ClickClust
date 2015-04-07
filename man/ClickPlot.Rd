\name{click.plot}
\alias{click.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot of the obtained clustering solution }
\description{
Constructs a click-plot for the clustering solution.
}
\usage{
click.plot(X, y = NULL, file = NULL, id, states = NULL, marg = 1,
  font.cex = 2, font.col = "black", cell.cex = 1, cell.lwd = 1.3,
  cell.col = "black", sep.lwd = 1.3, sep.col = "black",
  obs.lwd = NULL, colors = c("lightcyan", "pink", "darkred"),
  col.levels = 8, legend = TRUE, leg.cex = 1.3, top.srt = 0,
  frame = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset array (p x p x n) }
  \item{y }{ vector of initial states (length n) }
  \item{file }{ name of the output pdf-file }
  \item{id }{ classification vector (length n) }
  \item{states }{ vector of state labels (length p) }
  \item{marg }{ plot margin value (for the left and top) }
  \item{font.cex }{ magnification of labels }
  \item{font.col }{ color of labels }
  \item{cell.cex }{ magnification of cells }
  \item{cell.lwd }{ width of cell frames }
  \item{cell.col }{ color of cell frames }
  \item{sep.lwd }{ width of separator lines }
  \item{sep.col }{ color of separator lines }
  \item{obs.lwd }{ width of observation lines }
  \item{colors }{ edge colors for interpolation }
  \item{col.levels }{ number of colors obtained by interpolation }
  \item{legend }{ legend of color hues }
  \item{leg.cex }{ magnification of legend labels }
  \item{top.srt }{ rotation of state names in the top }
  \item{frame }{ frame around the plot }
}
\details{
Constructs a click-plot for the provided clustering solution. Click-plot is a graphical display representing relative transition frequencies for the partitioning specified via the parameter 'id'. If the parameter 'file' is specified, the constructed plot will be saved in the pdf-file with the name 'file'. If the width of observation lines 'obs.lwd' is not specified, median colors will be used for all cell segments.
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


# CONSTRUCT CLICK-PLOT

click.plot(X = C$X, y = C$y, file = NULL, id = M2$id)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ click-plot }
\keyword{ EM algorithm }
\keyword{ Markov model }% __ONLY ONE__ keyword per line
