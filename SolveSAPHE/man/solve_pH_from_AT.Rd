\encoding{latin1}
\name{solve_pH_from_AT}
\alias{solve_pH_from_AT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Solver for the total alkalinity-pH equations}
\description{Determines [H+] from Total alkalinity and dissolved total elements in sea water. Universal and robust algorithm from Munhoven (2013) with Newton-
Raphson iterations}
\usage{solve_pH_from_AT(p_alktot, p_dicvar, p_bortot, p_po4tot, p_siltot, 
                 p_nh4tot, p_h2stot, p_so4tot, p_flutot, p_pHscale, p_dicsel,
                 p_askVal=FALSE, p_dissoc, p_temp=18, p_sal=35, p_pres=0, p_hini)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{p_alktot}{Total alkalinity (mol/kg)}
	\item{p_dicvar}{Value of a carbonate system related variable : DIC, [CO2*], [HCO3-] or [CO3--] (mol/kg)
	
See below parameter p_dicsel}
	\item{p_bortot}{Total boron concentration (mol/kg)}
	\item{p_po4tot}{Total phosphate concentration (mol/kg)}
	\item{p_siltot}{Total silicate concentration (mol/kg)}
	\item{p_nh4tot}{Total ammonia concentration (mol/kg)}
	\item{p_h2stot}{Total sulfide concentration (mol/kg)}
	\item{p_so4tot}{Total sulphate concentration (mol/kg)}
	\item{p_flutot}{Total fluor concentration (mol/kg)}
	\item{p_pHscale}{Chosen pH scale: "T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale}
	\item{p_dicsel}{Carbonate variable selector (default = DIC). See parameter p_dicvar above. Values are:
	
"DIC"   :   p_dicvar = DIC (Dissolved Inorganic Carbon)

"CO2"   :   p_dicvar = [CO2*]

"HCO3"  :   p_dicvar = [HCO3-]

"CO3"   :   p_dicvar = [CO3--]}

	\item{p_askVal}{Optional boolean - set to TRUE if you want this function to return error on alkalinity, along with [H+] concentration value. Default is FALSE}
	\item{p_dissoc}{Named list of all dissociation constants. 
The list is optional but, if given, it should contain all members listed below excepted K2_Sil, which is itself optional. 
Member names are :

K1_DIC  :  First dissociation constant of carbonic acid (mol/kg) on chosen scale

K2_DIC  :  Second dissociation constant of carbonic acid (mol/kg) on chosen scale

K_BT    :  Dissociation constant of boric acid (mol/kg) on chosen scale

K1_PO4  :  First dissociation constant of phosphoric acid (mol/kg) on chosen scale

K2_PO4  :  Second dissociation constant of phosphoric acid (mol/kg) on chosen scale

K3_PO4  :  third dissociation constant of phosphoric acid (mol/kg) on chosen scale

K_Sil   :  First dissociation constant of sillicic acid (mol/kg) on chosen scale

K2_Sil  :  Second dissociation constant of sillicic acid (mol/kg) on chosen scale. It is optional.
If K2_Sil is absent from the list, then SiO2(OH)2 ion is not considered in the alkalinity equation.
Only SiO(OH)3 ion is.

K_NH4   :  Dissociation constant of ammonium (mol/kg) on chosen scale

K_H2S   :  Dissociation constant of hydrogen sulfide (mol/kg) on chosen scale

K_HSO4  :  Dissociation constant of hydrogen sulfate (mol/kg) on free scale

K_HF    :  Dissociation constant of hydrogen fluoride (mol/kg) on free scale

K_H2O   :  Dissociation constant of water (mol/kg) on chosen scale

Note that all dissociation constants shall be expressed in chosen pH scale except K_HF and K_HSO4 which shall be in free scale.

If the list is not given, these constants, excepted K2_Sil, will be calculated.}

	\item{p_temp}{Temperature in degree Celsius, to compute dissociation constants when p_dissoc is not given}
	\item{p_sal}{Salinity, in pratical salinty unit (psu), to compute dissociation constants when p_dissoc is not given}
	\item{p_pres}{Pressure, in bars, to compute dissociation constants when p_dissoc is not given}
	\item{p_hini}{Optional initial value of [H+] concentration
	
If p_dicsel = "CO3", a vector of two initial values are expected since there may be two solutions for pH.
Else, one initial value is expected}
}

\details{Formulations used when calculating dissociation constants:
  \itemize{
    \item{Carbonate if Total pH scale: }{Luecker et al. (2000) -- also Handbook (2007)}
    \item{Carbonate if SWS or Free pH scale: }{Millero et al. (1995) -- also Handbook (2007)}
    \item{Boric acid: }{Dickson (1990, eq. 2  3) -- also Handbook (2007, eq. 37)}
    \item{Phosphoric acid: }{Yao and Millero (1995)}
    \item{Silicic acid: }{Yao and Millero (1995) cited by Millero (1995)}
    \item{Ammonium: }{Yao and Millero (1995)}
    \item{Hydrogen sulfide: }{Millero et al. (1988) (cited by Millero (1995)}
    \item{Hydrogen sulfate: }{Dickson (1990) -- also Handbook (2007)}
    \item{Fluoric acid if Total pH scale: }{Perez and Fraga (1987)}
    \item{Fluoric acid if SWS or Free pH scale: }{Dickson and Riley (1979)}
    \item{Water: }{Millero (1995)}
  }

This function does not support vectors as arguments, only scalar values.
}

\value{If p_dicsel is "CO3", there may be one or two solutions for [H+], else there is only one. In case there are two solutions, both are returned.
  \itemize{
    \item If you set p_askVal to TRUE, the function returns a data frame containing the following columns:

    zh    :    [H+] concentration value(s) in the chosen pH scale

    val   :    Error(s) on total alkalinity, that is the deviation between total alkalinity calculated from [H+] and given total alkalinity

    \item If you set p_askVal to FALSE (default), the function returns only [H+] concentration value(s) in the chosen pH scale.
    
  }
}

\references{
Munhoven G. Munhoven G. (2013) Mathematics of the total alkalinity-pH equation - pathway
to robust and universal solution algorithms: the SolveSAPHE package v1.0.1. \emph{Geoscientif. Model Dev.}, 1367-1388
}

\author{
Guy Munhoven and Jean-Marie Epitalon
}

\examples{

## Compute [H+] from Alkalinity total and DIC, on total pH scale
p_dissoc <- list()
p_dissoc$K1_DIC <- 1.421828e-06
p_dissoc$K2_DIC <- 1.081555e-09
p_dissoc$K_BT   <- 2.526573e-09
p_dissoc$K1_PO4 <- 0.02408434
p_dissoc$K2_PO4 <- 1.076024e-06
p_dissoc$K3_PO4 <- 1.600484e-09
p_dissoc$K_Sil  <- 4.071935e-10
p_dissoc$K_NH4  <- 5.380823e-10
p_dissoc$K_H2S  <- 3.087264e-07
p_dissoc$K_HSO4 <- 0.1003021    # on free scale
p_dissoc$K_HF   <- 0.00176441   # on free scale
p_dissoc$K_H2O  <- 5.97496e-14
solve_pH_from_AT(p_alktot=2.5e-3, p_dicvar=2e-3, p_bortot=0.0004157, p_po4tot=0, p_siltot=0, 
                 p_nh4tot=0, p_h2stot=0, p_so4tot=0.0282, p_flutot=6.832e-05, p_pHscale="T",
                 p_dicsel="DIC", p_dissoc=p_dissoc)

## Giving inital [H+] value and asking for final error on alkalinity, on seawater pH scale
result <- solve_pH_from_AT(p_alktot=0.00234, p_dicvar=0.001936461, p_bortot=0.0004157, p_po4tot=0,
                p_siltot=0, p_nh4tot=0, p_h2stot=0, p_so4tot=0.0282, p_flutot=6.832e-05, 
                p_pHscale="SWS",p_dicsel="DIC", p_dissoc=p_dissoc, p_askVal=TRUE, p_hini=1.e-8)
H <- result$zh
error <- result$val

## Compute [H+] from Alkalinity total and CO3, on total free scale
H <- solve_pH_from_AT(p_alktot=0.00240, p_dicvar=2.775481e-04, p_bortot=0.0004157, p_po4tot=0,
                p_siltot=0, p_nh4tot=0, p_h2stot=0, p_so4tot=0.0282, p_flutot=6.832e-05, 
                p_pHscale="F", p_dicsel="CO3", FALSE, p_dissoc=p_dissoc)
H1 <- H[1]
if (H[2] != 1) print ( c("Second solution : ", H[2]))
}
\keyword{utilities}
