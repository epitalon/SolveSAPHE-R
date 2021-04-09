\encoding{latin1}
\name{AK_AMMO_1_YAMI95}
\alias{AK_AMMO_1_YAMI95}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dissociation constant of ammonium in sea-water [mol/kg-SW]}
\description{Function returns the dissociation constant of ammonium in sea-water [mol/kg-SW]}
\usage{AK_AMMO_1_YAMI95(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   SWS}

\value{
  \item{KNH}{Dissociation constant of ammonium in mol/kg-SW}
}

\references{Yao and Millero (1995), Millero (1995) for pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_AMMO_1_YAMI95(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
