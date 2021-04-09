\encoding{latin1}
\name{AK_CARB_0_WEIS74}
\alias{AK_CARB_0_WEIS74}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Henry's constant K0 in (mol/kg-SW)/atmosphere}
\description{Calculates Henry's constant (K0) based on Weiss (1979) formulation}
\usage{AK_CARB_0_WEIS74(t_k, s)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{t_k}{temperature in Kelvin}
	\item{s}{Salinity in psu}
}

\details{currently no pressure correction}

\value{
  \item{K0}{Henry's constant mol/(kg/atm)}
}

\references{Weiss R. F., 1974 Carbon dioxide in water and seawater: the solubility of a non-ideal gas. \emph{Marine Chemistry} \bold{2}, 203-215.}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{
  AK_CARB_0_WEIS74(t_k=298,s=35)
}
\keyword{utilities}
