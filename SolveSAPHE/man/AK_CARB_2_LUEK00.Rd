\encoding{latin1}
\name{AK_CARB_2_LUEK00}
\alias{AK_CARB_2_LUEK00}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Second dissociation constant of carbonic acid in mol/kg-SW from Luecker et al.}
\description{Calculates second dissociation constant of carbonic acid in mol/kg-SW on the total pH-scale, following Luecker et al. (2000)}
\usage{AK_CARB_2_LUEK00(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   Total}

\value{
  \item{K2}{Second dissociation constant of carbonic acid in mol/kg-SW}
}

\references{Luecker et al. (2000) -- also Handbook (2007), Millero (1979) pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_CARB_2_LUEK00(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
