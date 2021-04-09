\encoding{latin1}
\name{AK_PHOS_1_MILL95}
\alias{AK_PHOS_1_MILL95}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First dissociation constant of phosphoric acid (H3PO4) in seawater}
\description{Calculates first dissociation constant of phosphoric acid on the SWS pH-scale}
\usage{AK_PHOS_1_MILL95(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\details{pH scale:   SWS}

\value{
  \item{KP1}{First dissociation constant of phosphoric acid in mol/kg-SW}
}

\references{Yao and Millero (1995), Millero (1995) for pressure correction}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_PHOS_1_MILL95(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
