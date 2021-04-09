\encoding{latin1}
\name{AK_H2S_1_MILL95}
\alias{AK_H2S_1_MILL95}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dissociation constant of hydrogen sulfide in sea-water}
\description{Function returns the dissociation constant of hydrogen sulfide in sea-water}
\usage{AK_H2S_1_MILL95(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   SWS (according to Yao and Millero, 1995, p. 82: "refitted if necessary")

Total (according to Lewis and Wallace, 1998)

We stick to SWS here for the time being

The fits from Millero (1995) and Yao and Millero (1995)
derive from Millero et al. (1988), with all the coefficients
multiplied by -ln(10)
}

\value{
  \item{KHS}{Dissociation constant of hydrogen sulfide in mol/kg-SW}
}

\references{Millero et al. (1988) (cited by Millero (1995), Millero (1995) for pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_H2S_1_MILL95(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
