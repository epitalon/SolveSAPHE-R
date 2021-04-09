\encoding{latin1}
\name{AK_HF_PEFR87}
\alias{AK_HF_PEFR87}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dissociation constant of hydrogen fluoride in sea-water [mol/kg-SW]}
\description{Function returns the dissociation constant of hydrogen fluoride [mol/kg-SW]}
\usage{AK_HF_PEFR87(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   Total (according to Handbook, 2007)}

\value{
  \item{KHF}{Dissociation constant of hydrogen fluoride in mol/kg-SW}
}

\references{Perez and Fraga (1987), Millero (1995) for pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_HF_PEFR87(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
