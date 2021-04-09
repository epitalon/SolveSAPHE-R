\encoding{latin1}
\name{ACVT_HSWS_O_HFREE}
\alias{ACVT_HSWS_O_HFREE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conversion ratio H_SWS/H_free from free pH-scale to SWS scale}
\description{Function returns the ratio H_SWS/H_free as a function of temperature, salinity and pressure}
\usage{ACVT_HSWS_O_HFREE(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\value{
  \item{ratio}{Ratio H_SWS/H_free}
}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  ACVT_HSWS_O_HFREE(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
