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
ANW <- function (p_h, p_dicvar, p_bortot,
                      p_po4tot, p_siltot,
                      p_nh4tot, p_h2stot,
                      p_so4tot, p_flutot,
                      p_dicsel, api, p_deriv=FALSE)
#===============================================================================
{


    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_h            :  pH
    #  p_dicvar       :  Carbonate value
    #  p_bortot       :  Dissolved bore total
    #  p_po4tot       :  Dissolved Phosphate total
    #  p_siltot       :  Dissolved Silicium total
    #  p_nh4tot       :  Dissolved Ammonium total
    #  p_h2stot       :  Dissolved Hydrogene sulfide total
    #  p_so4tot       :  Dissolved Sulfate total
    #  p_flutot       :  Dissolved Fluor total
    #  p_dicsel       :  Carbonate variable selector
    #                      p_dicsel = "DIC":   p_dicvar = DIC
    #                      p_dicsel = "CO2":   p_dicvar = [CO2*]
    #                      p_dicsel = "HCO3":  p_dicvar = [HCO3-]
    #                      p_dicsel = "CO3":   p_dicvar = [CO3--]
    #  api            :  Pi factors of dissociation constants
    #  p_deriv        :  Request to compute derivative dAlk/dH


    #-------------------------------------------------------------------------------


    # H2CO3 - HCO3 - CO3 : n=2, m=0
    znumer_dic = 2.*api$api2_dic + p_h*       api$api1_dic
    
    if (p_dicsel == "DIC")
    {
        # DIC is the control variable of the carbonate system
        zdenom_dic =     api$api2_dic + p_h*(      api$api1_dic + p_h)
    }
    else if (p_dicsel == "CO2")
    {
        # [CO2] is the control variable of the carbonate system
        zdenom_dic =                p_h*                  p_h
    }
    else if (p_dicsel == "HCO3")
    {
        # [HCO3] is the control variable of the carbonate system
        zdenom_dic =                p_h*       api$api1_dic
    }
    else if (p_dicsel == "CO3")
    {
        # [CO3] is the control variable of the carbonate system
        zdenom_dic =     api$api2_dic
    }
    zalk_dic   = p_dicvar * (znumer_dic/zdenom_dic)

    # B(OH)3 - B(OH)4 : n=1, m=0
    znumer_bor =       api$api1_bor
    zdenom_bor =       api$api1_bor + p_h
    zalk_bor   = p_bortot * (znumer_bor/zdenom_bor)

    # H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
    znumer_po4 = 3.*api$api3_po4 + p_h*(2.*api$api2_po4 + p_h* api$api1_po4)
    zdenom_po4 =       api$api3_po4 + p_h*(      api$api2_po4 + p_h*(api$api1_po4 + p_h))
    zalk_po4   = p_po4tot * (znumer_po4/zdenom_po4 - 1.) # Zero level of H3PO4 = 1

    # H4SiO4 - H3SiO4 - H2SiO4 : n=2, m=0
    # if api$api2_sil is zero, then it excludes H2SiO4, downgrading to :  H4SiO4 - H3SiO4 : n=1, m=0
    znumer_sil = 2.*api$api2_sil + p_h*       api$api1_sil
    zdenom_sil =    api$api2_sil + p_h*(      api$api1_sil + p_h)
    zalk_sil   = p_siltot * (znumer_sil/zdenom_sil)

    # NH4 - NH3 : n=1, m=0
    znumer_nh4 =       api$api1_nh4
    zdenom_nh4 =       api$api1_nh4 + p_h
    zalk_nh4   = p_nh4tot * (znumer_nh4/zdenom_nh4)

    # H2S - HS : n=1, m=0
    znumer_h2s =       api$api1_h2s
    zdenom_h2s =       api$api1_h2s + p_h
    zalk_h2s   = p_h2stot * (znumer_h2s/zdenom_h2s)

    # HSO4 - SO4 : n=1, m=1
    znumer_so4 =       api$api1_so4
    zdenom_so4 =       api$api1_so4 + p_h
    zalk_so4   = p_so4tot * (znumer_so4/zdenom_so4 - 1.)

    # HF - F : n=1, m=1
    znumer_flu =       api$api1_flu
    zdenom_flu =       api$api1_flu + p_h
    zalk_flu   = p_flutot * (znumer_flu/zdenom_flu - 1.)


    z_anw =   zalk_dic + zalk_bor + zalk_po4 + zalk_sil +
              zalk_nh4 + zalk_h2s + zalk_so4 + zalk_flu


    # Do we need to compute derivative dAlk/dH ?
    if (p_deriv) 
    {

        # H2CO3 - HCO3 - CO3 : n=2
        if (p_dicsel == "DIC")
        {
            # DIC is the control variable of the carbonate system
            zdnumer_dic = api$api1_dic*api$api2_dic 
                + p_h*(4.*api$api2_dic + p_h*       api$api1_dic)
            zdalk_dic   = -p_dicvar*(zdnumer_dic/zdenom_dic^2)
        }
        else if (p_dicsel == "CO2")
        {
            # [CO2] is the control variable of the carbonate system
            zdnumer_dic = 4.*api$api2_dic   
                + p_h*(2.*api$api1_dic - p_h)
            zdalk_dic   = -p_dicvar*(zdnumer_dic/p_h^3)
        }
        else if (p_dicsel == "HCO3")
        {
            # [HCO3] is the control variable of the carbonate system
            zdalk_dic   = -p_dicvar*(2.*api$api2_dic/p_h^2)
        }
        else if (p_dicsel == "CO3")
        {
            # [CO3] is the control variable of the carbonate system
            zdalk_dic   =  p_dicvar/api$api2_dic
        }

        # B(OH)3 - B(OH)4 : n=1
        zdnumer_bor = api$api1_bor
        zdalk_bor   = -p_bortot*(zdnumer_bor/zdenom_bor^2)

        # H3PO4 - H2PO4 - HPO4 - PO4 : n=3
        zdnumer_po4 = api$api2_po4 * api$api3_po4 + p_h * (4.*api$api1_po4 * api$api3_po4     +
                                        p_h * (9.*api$api3_po4 + api$api1_po4 * api$api2_po4  +
                                        p_h * (4.*api$api2_po4                                +
                                        p_h *     api$api1_po4)))
        zdalk_po4   = -p_po4tot * (zdnumer_po4/zdenom_po4^2)

        # H4SiO4 - H3SiO4 - H2SiO4 : n=2, m=0
        # if api$api2_sil is zero, then H2SiO4 does not account, downgrading to :  H4SiO4 - H3SiO4 : n=1, m=0
        zdnumer_sil = api$api1_sil*api$api2_sil 
            + p_h*(4.*api$api2_sil + p_h*       api$api1_sil)
        zdalk_sil   = -p_siltot*(zdnumer_sil/zdenom_sil^2)

        # NH4 - NH3 : n=1
        zdnumer_nh4 = api$api1_nh4
        zdalk_nh4   = -p_nh4tot * (zdnumer_nh4/zdenom_nh4^2)

        # H2S - HS : n=1
        zdnumer_h2s = api$api1_h2s
        zdalk_h2s   = -p_h2stot * (zdnumer_h2s/zdenom_h2s^2)

        # HSO4 - SO4 : n=1
        zdnumer_so4 = api$api1_so4
        zdalk_so4   = -p_so4tot * (zdnumer_so4/zdenom_so4^2)

        # HF - F : n=1
        zdnumer_flu = api$api1_flu
        zdalk_flu   = -p_flutot * (zdnumer_flu/zdenom_flu^2)

        z_derivanw =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil +
                       zdalk_nh4 + zdalk_h2s + zdalk_so4 + zdalk_flu 

        return (c(z_anw, z_derivanw ))
    }
    else
        return (z_anw)
}
