#
#    Copyright  2013, 2014, 2020, 2021 Guy Munhoven
#
#    This file is part of SolveSAPHE v. 2

#    SolveSAPHE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SolveSAPHE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with SolveSAPHE.  If not, see <http://www.gnu.org/licenses/>.
#


#===============================================================================
equation_at <- function (p_alktot, p_h,       p_dicvar, p_bortot,
                         p_po4tot, p_siltot,  p_nh4tot, p_h2stot,
                         p_so4tot, p_flutot,  
                         p_dicsel="DIC", api, p_deriv=FALSE)
#===============================================================================
{
    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alktot       :  Alkalinity total
    #  p_h            :  pH
    #  p_dicvar       :  Carbonate value
    #  p_bortot       :  Dissolved bore total
    #  p_po4tot       :  Dissolved Phosphate total
    #  p_siltot       :  Dissolved Silicium total
    #  p_nh4tot       :  Dissolved Ammonium total
    #  p_h2stot       :  Dissolved Hydrogene sulfide total
    #  p_so4tot       :  Dissolved Sulfate total
    #  p_flutot       :  Dissolved Fluor total
    #  p_dicsel       :  Carbonate variable selector (default = DIC)
    #                      p_dicsel = "DIC":   p_dicvar = DIC
    #                      p_dicsel = "CO2":   p_dicvar = [CO2*]
    #                      p_dicsel = "HCO3":  p_dicvar = [HCO3-]
    #                      p_dicsel = "CO3":   p_dicvar = [CO3--]
    #  api            :  Pi factors of dissociation constants
    #  p_deriv        :  Request to compute derivative dAlk/dH

    #-------------------------------------------------------------------------------

    # H2O - OH
    zalk_wat   = api$api1_wat/p_h - p_h/api$aphscale

    result =   ANW(p_h, p_dicvar, p_bortot, p_po4tot, p_siltot,
                    p_nh4tot, p_h2stot, p_so4tot, p_flutot,
                    p_dicsel, api, p_deriv)
                            
    z_anw = result[1]
    zdiffalk = z_anw + zalk_wat - p_alktot
    
    # Do we need to compute derivative dAlk/dH ?
    if (p_deriv) 
    {
        z_deriv_anw = result[2]
        z_deriveqn =  z_deriv_anw - api$api1_wat/p_h^2 - 1./api$aphscale

        # return alkalinity difference and its derivative
        return (c(zdiffalk, z_deriveqn))

    }
    else
    {
        # return alkalinity difference
        return (zdiffalk)
    }
}

