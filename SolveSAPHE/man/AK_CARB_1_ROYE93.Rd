\encoding{latin1}
\name{AK_CARB_1_ROYE93}
\alias{AK_CARB_1_ROYE93}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First dissociation constant of carbonic acid in mol/kg-SW, from Roy et al.}
\description{Calculates first dissociation constant of carbonic acid in mol/kg-SW on the total pH-scale, following Roy et al. (1993)}
\usage{AK_CARB_1_ROYE93(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   Total}

\value{
  \item{K1}{First dissociation constant of carbonic acid in mol/kg-SW}
}

\references{Roy et al. (1993) -- also Handbook (1994), Millero (1979) pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_CARB_1_ROYE93(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
