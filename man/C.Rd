\name{C}
\docType{data}
\alias{C}
\title{Dataset: simulated dataset}

\description{ This dataset is used to run the backward state selection procedure (results in "B3").
}

\usage{data(utilityB3)}
\details{ Original dataset used to illustrate the utility of backward selection. }
\references{
Melnykov, V. (2014) Model-based biclustering of clickstream data, accepted by Computational Statistics & Data Analysis.

Melnykov, V. (2015) Model-based biclustering of clickstream data, under review at Journal of Statistical Software.
}
\seealso{help(B3)}
\examples{

data(utilityB3)

dev.new(width = 11, height = 11)
click.plot(X = C$X, id = B3$id, colors = c("lightyellow", "red", "darkred"), col.levels = 10)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dataset}% __ONLY ONE__ keyword per line
