\name{B3}
\docType{data}
\alias{B3}
\title{Dataset: result of backward state selection}

\description{ These data demonstrate the result of the backward state selection procedure obtained for the dataset "C".
}

\usage{data(utilityB3)}
\details{ Results of the backward state selection procedure assuming three components are provided for the dataset "C".}
\references{
Melnykov, V. (2016) Model-Based Biclustering of Clickstream Data, Computational Statistics and Data Analysis, 93, 31-45.
Melnykov, V. (2016) ClickClust: An R Package for Model-Based Clustering of Categorical Sequences, Journal of Statistical Software, 74, 1-34.
}
\seealso{help(C, package = "ClickClust")}
\examples{

data(utilityB3)

dev.new(width = 11, height = 11)
click.plot(X = C$X, id = B3$id, colors = c("lightyellow", "red", "darkred"), col.levels = 10)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dataset}% __ONLY ONE__ keyword per line
