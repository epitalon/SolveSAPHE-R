\encoding{latin1}
\name{AK_BORA_DICK90}
\alias{AK_BORA_DICK90}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boric acid dissociation constant KB in mol/kg-SW}
\description{Calculates boric acid dissociation constant in mol/kg-SW on the total pH-scale}
\usage{AK_BORA_DICK90(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   Total}

\value{
  \item{KB}{Dissociation constant of boric acid in mol/kg-SW}
}

\references{Dickson (1990, eq. 23) -- also Handbook (2007, eq. 37), Millero (1979) pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_BORA_DICK90(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
