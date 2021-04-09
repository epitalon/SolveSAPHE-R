\encoding{latin1}
\name{AK_SILI_1_MILL95}
\alias{AK_SILI_1_MILL95}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{First dissociation constant of sillicic acid (H4SiO4) in seawater}
\description{Calculates first dissociation constant of sillicic acid on the SWS pH-scale}
\usage{AK_SILI_1_MILL95(t_k, s)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
}

\details{pH scale:   SWS (according to Dickson et al, 2007)
    
No pressure correction available
}

\value{
  \item{KSi}{First dissociation constant of sillicic acid in mol/kg-SW}
}

\references{Yao and Millero (1995) cited by Millero (1995)}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_SILI_1_MILL95(t_k=298,s=35)
}
\keyword{utilities}
