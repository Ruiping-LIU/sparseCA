\name{ColorSound}
\alias{ColorSound}
\docType{data}
\title{
Toy example
}
\description{
Toy example in sparse CA
}
\usage{data("ColorSound")}
\format{
  A contingency table on the following 2 variables.
  \describe{
    \item{\code{Color}}{a numeric vector}
    \item{\code{Sound}}{a numeric vector}
  }
}
\details{
Data of Table where J= 9 pieces of music were presented to 22 participants who were asked to associate one of I =10 colors to each piece of music.
}
\source{
Abdi, H. and Béra, M. (2014). Correspondence Analysis. In Encyclopedia of Social Network Analysis and Mining, pages 275–284. Springer New York, New York, NY.
}
\references{
Ruiping Liu, Ndeye Niang Keita, Gilbert Saporta, Huiwen Wang. (2022). Sparse Correspondence Analysis for Contingency Tables. Advances in Data Analysis and Classification. https://doi.org/10.1007/s11634-022-00531-5
}
\examples{
data(ColorSound)
## maybe str(ColorSound) ; plot(ColorSound) ...
}
\keyword{datasets}
