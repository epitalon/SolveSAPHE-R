\encoding{latin1}
\name{ACVT_HTOT_O_HFREE}
\alias{ACVT_HTOT_O_HFREE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conversion ratio H_Tot/H_free from free pH-scale to Total scale}
\description{Function returns the ratio H_Tot/H_free as a function of temperature, salinity and pressure}
\usage{ACVT_HTOT_O_HFREE(t_k, s, p_bar)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
	\item{p_bar}{pressure in bar}
}

\value{
  \item{ratio}{Ratio H_Tot/H_free}
}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  ACVT_HTOT_O_HFREE(t_k=298,s=35, p_bar=0)
}
\keyword{utilities}
