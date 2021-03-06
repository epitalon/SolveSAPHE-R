\encoding{latin1}
\name{AK_CARB_2_MILL95}
\alias{AK_CARB_2_MILL95}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Second dissociation constant of carbonic acid in mol/kg-SW from Millero et al. (1995)}
\description{Calculates second dissociation constant of carbonic acid in mol/kg-SW on the SWS pH-scale.}
\usage{AK_CARB_2_MILL95(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   SWS}

\value{
  \item{K2}{Second dissociation constant of carbonic acid in mol/kg-SW}
}

\references{Millero (1995, eq 51 -- ln K2(COM)), Millero (1979) pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_CARB_2_MILL95(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
