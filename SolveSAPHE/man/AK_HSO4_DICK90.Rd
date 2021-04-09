\encoding{latin1}
\name{AK_HSO4_DICK90}
\alias{AK_HSO4_DICK90}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dissociation constant of hydrogen sulfate (bisulfate) in sea-water }
\description{Function returns the dissociation constant of hydrogen sulfate [mol/kg-SW]}
\usage{AK_HSO4_DICK90(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   Total (according to Handbook, 2007)}

\value{
  \item{KSO}{Dissociation constant of hydrogen sulfate in mol/kg-SW}
}

\references{Dickson (1990) -- also Handbook (2007), Millero (1995) for pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_HSO4_DICK90(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
