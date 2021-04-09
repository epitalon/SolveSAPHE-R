\encoding{latin1}
\name{AK_W_MILL95}
\alias{AK_W_MILL95}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Water dissociation constant Kw in (mol/kg-SW)^2}
\description{Calculates water dissociation constant Kw in (mol/kg-SW)^2 on the SWS pH-scale}
\usage{AK_W_MILL95(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   SWS}

\value{
  \item{Kw}{Dissociation constant of water in (mol/kg-SW^2}
}

\references{Millero (1995) for value at p_bar = 0, Millero (pers. comm. 1996) for pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_W_MILL95(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
