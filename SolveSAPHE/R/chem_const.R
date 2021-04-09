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




# --------------------------------------------------------
# List of subroutines for the chemical constants (PRIVATE)
# --------------------------------------------------------

# AK_CARB_0_WEIS74
# AK_CARB_1_MILL95, AK_CARB_2_MILL95
# AK_CARB_1_LUEK00, AK_CARB_2_LUEK00
# AK_CARB_1_ROYE93, AK_CARB_2_ROYE93
# AK_BORA_DICK90
# AK_W_MILL95
# AK_PHOS_1_MILL95, AK_PHOS_2_MILL95, AK_PHOS_3_MILL95
# AK_SILI_1_MILL95
# AK_H2S_1_MILL95
# AK_AMMO_1_YAMI95
# ABETA_HF_DIRI79
# AK_HF_PEFR87
# AK_HSO4_DICK90
# A_IONSTRENGTH_SALIN
# ASP_CALC_MUCC83, ASP_ARAG_MUCC83
# A_BTOT_SALIN, A_CATOT_SALIN, A_FTOT_SALIN, A_SO4TOT_SALIN
# ACVT_HSWS_O_HTOT, ACVT_HTOT_O_HFREE, ACVT_HSWS_O_HFREE
# A_RHOSW1_MUNH97, A_RHOSW2_MUNH97



# --------------------------------------
# Parameters for usage within the module
# --------------------------------------

# Gas constant
# ------------

#gasconst_bar_cm3_o_mol_k = 83.14510 # libthdyct
gasconst_bar_cm3_o_mol_k = 83.14472 # Handbook (2007)

# 0 degrees centigrade in Kelvin
# ------------------------------

t_k_zerodegc = 273.15 # Handbook (2007)



#=======================================================================
compute_dissoc_constants <- function (p_temp, p_sal, p_pres, p_pHscale="T")
#=======================================================================
#
# This function computes and returns dissociation constants
# of main chemical elements dissolved in sea water
#
# Input parameters :
#   p_temp         :  temperature in degree Celsius
#   p_sal          :  salinity in psu (Practical Salinity Unit)
#   p_pres         :  pressure in bar
#   p_pHscale      :  pH scale: "T" for Total scale (default), "SWS" for SeaWater, "F" for Free
#
# Output : this function  returns a named list of all dissociation constants. Member names are :
#
# K1_DIC  :  First dissociation constant of carbonic acid (mol/kg) on chosen scale
# K2_DIC  :  Second dissociation constant of carbonic acid (mol/kg) on chosen scale
# K_BT    :  Dissociation constant of boric acid (mol/kg) on chosen scale
# K1_PO4  :  First dissociation constant of phosphoric acid (mol/kg) on chosen scale
# K2_PO4  :  Second dissociation constant of phosphoric acid (mol/kg) on chosen scale
# K3_PO4  :  third dissociation constant of phosphoric acid (mol/kg) on chosen scale
# K_Sil   :  Dissociation constant of sillicic acid (mol/kg) on chosen scale
# K_NH4   :  Dissociation constant of ammonium (mol/kg) on chosen scale
# K_H2S   :  Dissociation constant of hydrogen sulfide (mol/kg) on chosen scale
# K_HSO4  :  Dissociation constant of hydrogen sulfate (mol/kg) on free scale
# K_HF    :  Dissociation constant of hydrogen fluoride (mol/kg) on free scale
# K_H2O   :  Dissociation constant of water (mol/kg) on chosen scale
#
# Note :  K2_Sil, second dissociation constant of sillicic acid is not computed.
{
    t_k = p_temp + t_k_zerodegc
    
    dissoc <- list()
    if (p_pHscale == "SWS")
    {
        zcvt_hsws_o_htot  = ACVT_HSWS_O_HTOT(t_k, p_sal, p_pres)

        dissoc$K1_DIC <- AK_CARB_1_MILL95(t_k, p_sal, p_pres)
        dissoc$K2_DIC <- AK_CARB_2_MILL95(t_k, p_sal, p_pres)
        dissoc$K_BT   <- AK_BORA_DICK90(t_k, p_sal, p_pres) * zcvt_hsws_o_htot
        dissoc$K1_PO4 <- AK_PHOS_1_MILL95(t_k, p_sal, p_pres)
        dissoc$K2_PO4 <- AK_PHOS_2_MILL95(t_k, p_sal, p_pres)
        dissoc$K3_PO4 <- AK_PHOS_3_MILL95(t_k, p_sal, p_pres)
        dissoc$K_Sil  <- AK_SILI_1_MILL95(t_k, p_sal)
        dissoc$K_NH4  <- AK_AMMO_1_YAMI95(t_k, p_sal, p_pres)
        dissoc$K_H2S  <- AK_H2S_1_MILL95(t_k, p_sal, p_pres)
        dissoc$K_HSO4 <- AK_HSO4_DICK90(t_k, p_sal, p_pres)
        dissoc$K_HF   <- 1. / ABETA_HF_DIRI79(t_k, p_sal, p_pres)
        dissoc$K_H2O  <- AK_W_MILL95(t_k, p_sal, p_pres)
    }
    else if (p_pHscale == "F")
    {
        zcvt_hfree_o_htot = 1. / ACVT_HTOT_O_HFREE(t_k, p_sal, p_pres)
        zcvt_hfree_o_hsws = 1. / ACVT_HSWS_O_HFREE(t_k, p_sal, p_pres)

        dissoc$K1_DIC <- AK_CARB_1_MILL95(t_k, p_sal, p_pres) * zcvt_hfree_o_hsws
        dissoc$K2_DIC <- AK_CARB_2_MILL95(t_k, p_sal, p_pres) * zcvt_hfree_o_hsws
        dissoc$K_BT   <- AK_BORA_DICK90(t_k, p_sal, p_pres)   * zcvt_hfree_o_htot 
        dissoc$K1_PO4 <- AK_PHOS_1_MILL95(t_k, p_sal, p_pres) * zcvt_hfree_o_hsws
        dissoc$K2_PO4 <- AK_PHOS_2_MILL95(t_k, p_sal, p_pres) * zcvt_hfree_o_hsws
        dissoc$K3_PO4 <- AK_PHOS_3_MILL95(t_k, p_sal, p_pres) * zcvt_hfree_o_hsws
        dissoc$K_Sil  <- AK_SILI_1_MILL95(t_k, p_sal)         * zcvt_hfree_o_hsws
        dissoc$K_NH4  <- AK_AMMO_1_YAMI95(t_k, p_sal, p_pres) * zcvt_hfree_o_hsws
        dissoc$K_H2S  <- AK_H2S_1_MILL95(t_k, p_sal, p_pres)  * zcvt_hfree_o_hsws
        dissoc$K_HSO4 <- AK_HSO4_DICK90(t_k, p_sal, p_pres)
        dissoc$K_HF   <- 1. / ABETA_HF_DIRI79(t_k, p_sal, p_pres)
        dissoc$K_H2O  <- AK_W_MILL95(t_k, p_sal, p_pres)      * zcvt_hfree_o_hsws
    }
    else
    {
        zcvt_htot_o_hsws  = 1. / ACVT_HSWS_O_HTOT(t_k, p_sal, p_pres)
        zcvt_htot_o_hfree = ACVT_HTOT_O_HFREE(t_k, p_sal, p_pres)

        dissoc$K1_DIC <- AK_CARB_1_LUEK00(t_k, p_sal, p_pres)
        dissoc$K2_DIC <- AK_CARB_2_LUEK00(t_k, p_sal, p_pres)
        dissoc$K_BT   <- AK_BORA_DICK90(t_k, p_sal, p_pres)
        dissoc$K1_PO4 <- AK_PHOS_1_MILL95(t_k, p_sal, p_pres) * zcvt_htot_o_hsws
        dissoc$K2_PO4 <- AK_PHOS_2_MILL95(t_k, p_sal, p_pres) * zcvt_htot_o_hsws
        dissoc$K3_PO4 <- AK_PHOS_3_MILL95(t_k, p_sal, p_pres) * zcvt_htot_o_hsws
        dissoc$K_Sil  <- AK_SILI_1_MILL95(t_k, p_sal)         * zcvt_htot_o_hsws
        dissoc$K_NH4  <- AK_AMMO_1_YAMI95(t_k, p_sal, p_pres) * zcvt_htot_o_hsws
        dissoc$K_H2S  <- AK_H2S_1_MILL95(t_k, p_sal, p_pres)  * zcvt_htot_o_hsws
        dissoc$K_HSO4 <- AK_HSO4_DICK90(t_k, p_sal, p_pres)
        dissoc$K_HF   <- AK_HF_PEFR87(t_k, p_sal, p_pres) / zcvt_htot_o_hfree
        dissoc$K_H2O  <- AK_W_MILL95(t_k, p_sal, p_pres)      * zcvt_htot_o_hsws
    }
    
    return (dissoc)
}


#=======================================================================
 AK_CARB_0_WEIS74 <- function (t_k, s)
#=======================================================================
{
    # Function calculates K0 in (mol/kg-SW)/atmosphere

    # References: Weiss (1979) [(mol/kg-SW)/atm]
    # pH scale  : N/A
    # Note      : currently no pressure correction


    # ------------------
    # Argument variables
    # ------------------

    #     s      : salinity
    #     t_k    : temperature in K

    zt_k_o_100 = t_k/100.

    AK_CARB_0_WEIS74 = exp( -60.2409 + 93.4517/zt_k_o_100
                    + 23.3585*log(zt_k_o_100)      
                    + (   0.023517 - 0.023656*zt_k_o_100
                      + 0.0047036*zt_k_o_100*zt_k_o_100)*s )


    return (AK_CARB_0_WEIS74)
}




#=======================================================================
AK_CARB_1_MILL95 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates first dissociation constant of carbonic acid
    # in mol/kg-SW on the SWS pH-scale.

    # References: Millero (1995, eq 50 -- ln K1(COM))
    #             Millero (1982) pressure correction
    # pH scale:   SWS


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in Kelvin
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zsqrts         : square root of salinity
    #     zds            : salinity-34.8
    #     zln_kc1_p0     : ln(K_C1) at p_bar = 0
    #     zln_kc1_pp     : pressure correction for p_bar /= 0


    # ln(K_C1) value at p_bar = 0

    zsqrts     = sqrt(s)

    zln_kc1_p0 = (    2.18867 - 2275.0360/t_k - 1.468591*log(t_k)
                + (-0.138681 -   9.33291/t_k)*zsqrts
                +  0.0726483*s
                - 0.00574938*s*zsqrts)


    # Pressure correction

    zt_degc    = t_k - t_k_zerodegc
    zds        = s - 34.8
    zrt        = gasconst_bar_cm3_o_mol_k * t_k

    zdvi       =  -25.50 - 0.151*zds + 0.1271*zt_degc
    zdki       = ( -3.08 - 0.578*zds + 0.0877*zt_degc)*1.0E-03

    zln_kc1_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_C1 value

    AK_CARB_1_MILL95 = exp( zln_kc1_p0 + zln_kc1_pp )

    return (AK_CARB_1_MILL95)
}




#=======================================================================
AK_CARB_2_MILL95 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates second dissociation constant K1
    # in mol/kg-SW on the SWS pH-scale.

    # References: Millero (1995, eq 51 -- ln K2(COM))
    #             Millero (1979) pressure correction
    # pH scale:   SWS


    # Argument variables
    # ------------------

    #     t_k    : temperature in Kelvin
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zsqrts         : square root of salinity
    #     zds            : salinity-34.8
    #     zln_kc2_p0     : ln(K_C2) at p_bar = 0
    #     zln_kc2_pp     : pressure correction for p_bar /= 0

    # ln(K_C2) value at p_bar = 0

    zsqrts     = sqrt(s)

    zln_kc2_p0 = (    -0.84226 - 3741.1288/t_k - 1.437139*log(t_k)
                + (-0.128417 -  24.41239/t_k)*zsqrts
                +  0.1195308*s
                - 0.00912840*s*zsqrts)


    # Pressure correction

    zt_degc    = t_k - t_k_zerodegc
    zds        = s - 34.8
    zrt        = gasconst_bar_cm3_o_mol_k * t_k

    zdvi       =  -15.82 + 0.321*zds - 0.0219*zt_degc
    zdki       =  ( 1.13 - 0.314*zds - 0.1475*zt_degc)*1.0E-03

    zln_kc2_pp =  (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_C2 value

    AK_CARB_2_MILL95  = exp( zln_kc2_p0 + zln_kc2_pp )

    return (AK_CARB_2_MILL95)
}




#=======================================================================
AK_CARB_1_LUEK00 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates first dissociation constant of carbonic acid
    # in mol/kg-SW on the Total pH-scale.

    # References: Luecker et al. (2000) -- also Handbook (2007)
    #             Millero (1979) pressure correction
    # pH scale:   Total


    # Argument variables
    # ------------------

    #     t_k    : temperature in Kelvin
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zds            : salinity-34.8
    #     zlog10_kc1_p0  : log_10(k_C1) at p_bar = 0
    #     zln_kc1_pp     : pressure correction for p_bar /= 0

    # log_10(K_C1) value at p_bar = 0

    zlog10_kc1_p0 = (  61.2172 - 3633.86/t_k - 9.67770*log(t_k)
                    + s*(0.011555 - s*0.0001152))


    # Pressure correction

    zt_degc    = t_k - t_k_zerodegc
    zds        = s - 34.8
    zrt        = gasconst_bar_cm3_o_mol_k * t_k

    zdvi       =  -25.50 - 0.151*zds + 0.1271*zt_degc
    zdki       = ( -3.08 - 0.578*zds + 0.0877*zt_degc)*1.0E-03

    zln_kc1_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_C1 value

    AK_CARB_1_LUEK00 = 10.**zlog10_kc1_p0 * exp(zln_kc1_pp)

    return (AK_CARB_1_LUEK00)
}




#=======================================================================
AK_CARB_2_LUEK00 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates second dissociation constant K1
    # in mol/kg-SW on the Total pH-scale.

    # References: Luecker et al. (2000) -- also Handbook (2007)
    #             Millero (1979) pressure correction
    # pH scale:   Total


    # Argument variables
    # ------------------

    #     t_k    : temperature in Kelvin
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zsqrts         : square root of salinity
    #     zds            : salinity-34.8
    #     zlog10_kc2_p0  : log_10(K_C2) at p_bar = 0
    #     zln_kc2_pp     : pressure correction for p_bar /= 0


    # log_10(K_C2) value at p_bar = 0

    zlog10_kc2_p0 = (-25.9290 - 471.78/t_k + 3.16967*log(t_k)
                    + s*(0.01781 - s*0.0001122))


    # Pressure correction

    zt_degc    = t_k - t_k_zerodegc
    zds        = s - 34.8
    zrt        = gasconst_bar_cm3_o_mol_k * t_k

    zdvi       =  -15.82 + 0.321*zds - 0.0219*zt_degc
    zdki       =  ( 1.13 - 0.314*zds - 0.1475*zt_degc)*1.0E-03

    zln_kc2_pp =  (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_C2 value

    AK_CARB_2_LUEK00  = 10.**zlog10_kc2_p0 *exp(zln_kc2_pp)

    return (AK_CARB_2_LUEK00)
}







#=======================================================================
AK_CARB_1_ROYE93 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates first dissociation constant of carbonic acid
    # in mol/kg-SW on the Total pH-scale.

    # References: Roy et al. (1993) -- also Handbook (1994)
    #             Millero (1979) pressure correction
    # pH scale  : Total
    # Note      : converted here from mol/kg-H2O to mol/kg-SW


    # Argument variables
    # ------------------

    #     t_k    : temperature in Kelvin
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zds            : salinity-34.8
    #     zln_kc1_p0     : ln(k_C1) at p_bar = 0
    #     zln_kc1_pp     : pressure correction for p_bar /= 0

    # ln(K_C1) value at p_bar = 0

    zsqrts     = sqrt(s)
    zcvt_to_kgsw    = ACVT_KGH2O_O_KGSW(s)

    zln_kc1_p0 = (-2307.1266/t_k +    2.83655 - 1.5529413*log(t_k)
                + (-4.0484/t_k - 0.20760841)*zsqrts
                + 0.08468345*s
                - 0.00654208*zsqrts*s)


    # Pressure correction

    zt_degc    = t_k - t_k_zerodegc
    zds        = s - 34.8
    zrt        = gasconst_bar_cm3_o_mol_k * t_k

    zdvi       =  -25.50 - 0.151*zds + 0.1271*zt_degc
    zdki       = ( -3.08 - 0.578*zds + 0.0877*zt_degc)*1.0E-03

    zln_kc1_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_C1 value

    AK_CARB_1_ROYE93 = exp(zln_kc1_p0 + zln_kc1_pp) * zcvt_to_kgsw

    return (AK_CARB_1_ROYE93)
}



#=======================================================================
AK_CARB_2_ROYE93 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates second dissociation constant K1
    # in mol/kg-SW on the Total pH-scale.

    # References: Roy et al. (1993) -- also Handbook (1994)
    #             Millero (1979) pressure correction
    # pH scale  : Total
    # Note      : converted here from mol/kg-H2O to mol/kg-SW


    # Argument variables
    # ------------------

    #     t_k    : temperature in Kelvin
    #     s      : salinity
    #     p_bar  : applied pressure in bar


    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zsqrts         : square root of salinity
    #     zds            : salinity-34.8
    #     zln_kc2_p0     : ln(K_C2) at p_bar = 0
    #     zln_kc2_pp     : pressure correction for p_bar /= 0

    # ln(K_C2) value at p_bar = 0

    zsqrts     = sqrt(s)
    zcvt_to_kgsw    = ACVT_KGH2O_O_KGSW(s)

    zln_kc2_p0 =   (-3351.6106/t_k - 9.226508 - 0.2005743*log(t_k)
                + ( -23.9722/t_k - 0.106901773)*zsqrts
                +  0.1130822*s
                - 0.00846934*zsqrts*s)


    # Pressure correction

    zt_degc    = t_k - t_k_zerodegc
    zds        = s - 34.8
    zrt        = gasconst_bar_cm3_o_mol_k * t_k

    zdvi       =  -15.82 + 0.321*zds - 0.0219*zt_degc
    zdki       =  ( 1.13 - 0.314*zds - 0.1475*zt_degc)*1.0E-03

    zln_kc2_pp =  (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_C2 value

    AK_CARB_2_ROYE93  = exp(zln_kc2_p0 + zln_kc2_pp) * zcvt_to_kgsw

    return (AK_CARB_2_ROYE93)
}



#=======================================================================
AK_BORA_DICK90 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates boric acid dissociation constant KB
    # in mol/kg-SW on the total pH-scale.

    # References: Dickson (1990, eq. 23) -- also Handbook (2007, eq. 37)
    #             Millero (1979) pressure correction
    # pH scale  : total


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in Kelvin
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt          : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc      : temperature in degrees Celsius
    #     zdvi         : volume change for ionization
    #     zdki         : compressibility change for ionization
    #     zsqrts       : square root of salinity
    #     zds          : salinity-34.8
    #     zln_kb_p0    : K_b at p_bar = 0
    #     zln_kb_pp    : pressure correction for p_bar /= 0


    # ln(K_B) value at p_bar = 0

    zsqrts     = sqrt(s)

    zln_kb_p0  = (( -8966.90
                            + zsqrts*( -2890.53
                            + zsqrts*(  -77.942
                            + zsqrts*(    1.728 - 0.0996*zsqrts)))) / t_k
                +  148.0248 + zsqrts*(137.1942 + zsqrts*1.62142)
                + (-24.4344 + zsqrts*(-25.085 - zsqrts*0.2474)) * log(t_k)
                + 0.053105*zsqrts*t_k)


    # Pressure correction

    zt_degc   = t_k - t_k_zerodegc
    zds       = s - 34.8
    zrt       = gasconst_bar_cm3_o_mol_k * t_k

    zdvi      = -29.48 + 0.295*zds + 0.1622*zt_degc - 0.002608*zt_degc*zt_degc
    zdki      = (-2.84 + 0.354*zds)*1.0E-03

    zln_kb_pp =  (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_B value

    AK_BORA_DICK90   = exp( zln_kb_p0 + zln_kb_pp )

    return (AK_BORA_DICK90)
}



#=======================================================================
AK_W_MILL95 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates water dissociation constant Kw in (mol/kg-SW)^2

    # References: Millero (1995) for value at p_bar = 0
    #             Millero (pers. comm. 1996) for pressure correction
    # pH scale  : SWS


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt        : R*t_k
    #     zt_degc    : temperature in degrees Celsius
    #     zdvi       : volume change for ionization
    #     zdki       : compressibility change for ionization
    #     zln_kw_p0  : ln(K_w) at p_bar = 0
    #     zln_kw_pp  : pressure correction for p_bar /= 0

    # ln(K_w) value at p_bar = 0

    zln_kw_p0 =  (  148.9802
                - 13847.26/t_k
                -  23.6521*log(t_k)
                + ( -5.977 + 118.67/t_k + 1.0495*log(t_k))*sqrt(s)
                - 0.01615*s)


    # Pressure correction

    zt_degc = t_k - t_k_zerodegc
    zrt     = gasconst_bar_cm3_o_mol_k * t_k

    zdvi    =  -20.02 + 0.1119*zt_degc - 0.1409E-02*zt_degc*zt_degc
    zdki    = ( -5.13 + 0.0794*zt_degc)*1.0E-03

    zln_kw_pp =  (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_w value

    AK_W_MILL95 = exp( zln_kw_p0 + zln_kw_pp )


    return (AK_W_MILL95)
}



#=======================================================================
AK_PHOS_1_MILL95 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the first dissociation constant
    # of phosphoric acid (H3PO4) in seawater

    # References: Yao and Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : SWS


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt          : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc      : temperature in degrees Celsius
    #     zdvi         : volume change for ionization
    #     zdki         : compressibility change for ionization
    #     zln_kp1_p0   : ln(K_p1) at p_bar = 0
    #     zln_kp1_pp   : pressure correction for p_bar /= 0

    # ln(K_P1) for p_bar = 0

    zln_kp1_p0 = (     115.54 - 4576.752/t_k - 18.453*log(t_k)
                + ( 0.69171 -  106.736/t_k)* sqrt(s)
                + (-0.01844 -  0.65643/t_k)*s )


    # Pressure correction

    zt_degc   = t_k - t_k_zerodegc
    zrt       = gasconst_bar_cm3_o_mol_k * t_k

    zdvi      =  -14.51 + 0.1211*zt_degc - 0.321E-03*zt_degc*zt_degc
    zdki      = ( -2.67 + 0.0427*zt_degc)*1.0E-03

    zln_kp1_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final value of K_P1

    AK_PHOS_1_MILL95 = exp(zln_kp1_p0 + zln_kp1_pp)

    return (AK_PHOS_1_MILL95)
}


#=======================================================================
AK_PHOS_2_MILL95 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the second dissociation constant
    # of phosphoric acid (H3PO4) in seawater

    # References: Yao and Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : SWS


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt          : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc      : temperature in degrees Celsius
    #     zdvi         : volume change for ionization
    #     zdki         : compressibility change for ionization
    #     zln_kp2_p0   : ln(K_P2) at p_bar = 0
    #     zln_kp2_pp   : pressure correction for p_bar /= 0

    # ln(K_P2) for p_bar = 0

    zln_kp2_p0 = (  172.1033
                - 8814.715/t_k
                -   27.927*log(t_k)
                + (  1.3566 -  160.340/t_k)*sqrt(s)
                + (-0.05778 +  0.37335/t_k)*s )


    # Pressure correction

    zt_degc    = t_k - t_k_zerodegc
    zrt        = gasconst_bar_cm3_o_mol_k * t_k

    zdvi       =  -23.12 + 0.1758*zt_degc -2.647E-03*zt_degc*zt_degc
    zdki       = ( -5.15 +   0.09*zt_degc)*1.0E-03

    zln_kp2_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_P2 value

    AK_PHOS_2_MILL95  = exp( zln_kp2_p0 + zln_kp2_pp )

    return (AK_PHOS_2_MILL95)
}


#=======================================================================
AK_PHOS_3_MILL95 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the third dissociation constant
    # of phosphoric acid (H3PO4) in seawater

    # References: Yao and Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : SWS

    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt          : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc      : temperature in degrees Celsius
    #     zdvi         : volume change for ionization
    #     zdki         : compressibility change for ionization
    #     zln_kp3_p0   : ln(K_P3) at p_bar = 0
    #     zln_kp3_pp   : pressure correction for p_bar /= 0


    # ln(K_P3) for p_bar = 0

    zln_kp3_p0 =   (    -18.126 -  3070.75/t_k
                + ( 2.81197 + 17.27039/t_k)*sqrt(s)
                + (-0.09984 - 44.99486/t_k)*s )


    # Pressure correction

    zt_degc   = t_k - t_k_zerodegc
    zrt       = gasconst_bar_cm3_o_mol_k * t_k

    zdvi      =  -26.57 + 0.2020*zt_degc -3.042E-03*zt_degc*zt_degc
    zdki      = ( -4.08 + 0.0714*zt_degc)*1.0E-03

    zln_kp3_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_P3 value

    AK_PHOS_3_MILL95 = exp( zln_kp3_p0 + zln_kp3_pp )

    return (AK_PHOS_3_MILL95)
}


#=======================================================================
AK_SILI_1_MILL95 <- function (t_k, s)
#=======================================================================
{
    # Function returns the first dissociation constant
    # of silicic acid (H4SiO4) in seawater

    # References: Yao and Millero (1995) cited by Millero (1995)
    # pH scale  : SWS (according to Dickson et al, 2007)
    # Note      : No pressure correction available
    # Note      : converted here from mol/kg-H2O to mol/kg-sw

    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity

    # ---------------
    # Local variables
    # ---------------

    #     zcvt_to_kgsw: fraction of pure water in 1 kg seawater at salinity s
    #     zionst      : ionic strength [mol/kg-H2O]
    #     zln_ksi1_p0 : ln(K_Si1) at p_bar = 0
    #     zln_ksi1_pp : pressure correciotn for p_bar /= 0

    # K_Si1 value at p_bar = 0

    zcvt_to_kgsw = ACVT_KGH2O_O_KGSW(s)
    zionst       = A_IONSTRENGTH_SALIN(s)/zcvt_to_kgsw  # mol/kg-H2O ##

    zln_ksi1_p0  = (    117.40 -  8904.2/t_k - 19.334 * log(t_k)
                + ( 3.5913 -  458.79/t_k) * sqrt(zionst)
                + (-1.5998 +  188.74/t_k) * zionst
                + (0.07871 - 12.1652/t_k) * zionst*zionst )


    # Pressure correction : currently none

    zln_ksi1_pp = 0.


    # Final value

    AK_SILI_1_MILL95 = exp( zln_ksi1_p0 + zln_ksi1_pp ) * zcvt_to_kgsw

    return (AK_SILI_1_MILL95)
}


#=======================================================================
AK_H2S_1_MILL95 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the dissociation constant of hydrogen sulfide in sea-water


    # References: Millero et al. (1988) (cited by Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : - SWS (according to Yao and Millero, 1995, p. 82: "refitted if necessary")
    #             - Total (according to Lewis and Wallace, 1998)
    # Note      : we stick to SWS here for the time being
    # Note      : the fits from Millero (1995) and Yao and Millero (1995)
    #             derive from Millero et al. (1988), with all the coefficients
    #             multiplied by -ln(10)

    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zt_degc      : temperature in degrees Celsius
    #     zrt          : R*t_k, R in bar*cm3/(mol*K)
    #     zdvi         : volume change for ionization
    #     zdki         : compressibility change for ionization
    #     zln_kh2s_p0  : ln(K_H2S) at p_bar = 0
    #     zln_kh2s_pp  : pressure correction for p_bar /= 0


    # K_H2S value at p_bar = 0
    # ------------------------

    zln_kh2s_p0  =  (  225.838
                    - 13275.3/t_k
                    - 34.6435 * log(t_k)
                    +  0.3449*sqrt(s)
                    -  0.0274*s )


    # Pressure correction
    # -------------------

    zt_degc      = t_k - t_k_zerodegc
    zrt          = gasconst_bar_cm3_o_mol_k * t_k

    zdvi         =  -14.80 + zt_degc*(0.0020 - zt_degc*0.400E-03)
    zdki         = (  2.89 + zt_degc*0.054)*1.0E-03

    zln_kh2s_pp  = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_H2S value
    # -----------------

    AK_H2S_1_MILL95 = exp( zln_kh2s_p0 + zln_kh2s_pp )

    return (AK_H2S_1_MILL95)
}


#=======================================================================
AK_AMMO_1_YAMI95 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the dissociation constant
    # of ammonium in sea-water [mol/kg-SW]


    # References: Yao and Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : SWS

    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zt_degc      : temperature in degrees Celsius
    #     zrt          : R*t_k, R in bar*cm3/(mol*K)
    #     zdvi         : volume change for ionization
    #     zdki         : compressibility change for ionization
    #     zln_knh4_p0  : ln(K_NH4) at p_bar = 0
    #     zln_knh4_pp  : pressure correction for p_bar /= 0


    # K_NH4 value at p_bar = 0
    # ------------------------

    zln_knh4_p0  = (   -0.25444 -  6285.33/t_k + 0.0001635*t_k
                + ( 0.46532 - 123.7184/t_k) * sqrt(s)
                + (-0.01992 +  3.17556/t_k) * s )  


    # Pressure correction
    # -------------------

    zt_degc      = t_k - t_k_zerodegc
    zrt          = gasconst_bar_cm3_o_mol_k * t_k

    zdvi         =  -26.43 + zt_degc*(0.0889 - zt_degc*0.905E-03)
    zdki         = ( -5.03 + zt_degc*0.0814)*1.0E-03

    zln_knh4_pp  = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final K_NH4 value
    # -----------------

    AK_AMMO_1_YAMI95  = exp( zln_knh4_p0 + zln_knh4_pp )

    return (AK_AMMO_1_YAMI95)
}


#=======================================================================
ACVT_KGH2O_O_KGSW <- function (s)
#=======================================================================
{
    # Function returns the mass of pure water in one kg of seawater
    # of salinity s

    # References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
    #             "Handbook (2007)" -- Handbook (2007)
    # pH scale:   N/A


    #ACVT_KGH2O_O_KGSW = 1. - 0.0010049*s # libthdyct
    ACVT_KGH2O_O_KGSW = 1. - 0.001005*s # Handbook (2007)

    return (ACVT_KGH2O_O_KGSW)
}


#=======================================================================
A_IONSTRENGTH_SALIN <- function (s)
#=======================================================================
{
    # Function calculates ionic strength in mol/kg-SW, for given salinity.

    # References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
    #             "Handbook (2007)" -- Handbook (2007)
    # pH scale:   N/A


    #A_IONSTRENGTH_SALIN = (0.019920*s) # libthdyct
    A_IONSTRENGTH_SALIN = (0.019924*s) # Handbook (2007)

    return (A_IONSTRENGTH_SALIN)
}



#=======================================================================
ABETA_HF_DIRI79 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates association constant \beta_{HF} [(mol/kg-SW)^{-1}]
    # in (mol/kg-SW)^{-1}, where
    #   \beta_{HF} = \frac{ [HF] }{ [H^{+}] [F^{-}] }

    # References: Dickson and Riley (1979)
    #             Millero (1995) for pressure correction
    # pH scale  : free
    # Note      : converted here from mol/kg-H2O to mol/kg-SW

    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zionst         : ionic strength [mol/kg-H2O]
    #     zcvt_to_kgsw   : mass of pure water in 1kg of seawater as a fct. of salinity
    #     zln_bhf_p0     : \beta_HF at p_bar = 0
    #     zln_khf_pp     : pressure correction for k_HF = 1/\beta_HF at p_bar /= 0


    # \beta_HF at p_bar = 0
    # ---------------------

    zcvt_to_kgsw    = ACVT_KGH2O_O_KGSW(s)
    zionst          = A_IONSTRENGTH_SALIN(s)/zcvt_to_kgsw 

    zln_bhf_p0      = -1590.2/t_k + 12.641 - 1.525*sqrt(zionst)


    # Pressure correction
    # -------------------

    zt_degc      = t_k - t_k_zerodegc
    zrt          = gasconst_bar_cm3_o_mol_k * t_k

    zdvi         =   -9.78 + zt_degc*(-0.0090 - zt_degc*0.942E-03)
    zdki         = ( -3.91 + zt_degc*0.054)*1.0E-03

    zln_khf_pp   = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final \beta_HF value
    # --------------------
    #  notice that  ln(k_HF(P)) = ln(k_HF(0)) + zln_khf_pp
    #         <=>  -ln(\beta_HF(P)) = -ln(\beta_HF(0)) + zln_khf_pp
    #         <=>   ln(\beta_HF(P)) =  ln(\beta_HF(0)) - zln_khf_pp

    ABETA_HF_DIRI79 = exp(zln_bhf_p0 - zln_khf_pp ) / zcvt_to_kgsw

    return (ABETA_HF_DIRI79)
}



#=======================================================================
AK_HF_PEFR87 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function calculates dissociation constant for hydrogen fluoride
    # in mol/kg-SW

    # References: Perez and Fraga (1987)
    #             Millero (1995) for pressure correction
    # pH scale  : Total (according to Handbook, 2007)

    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zln_khf_p0     : ln(K_HF) at p_bar = 0
    #     zln_khf_pp     : pressure correction for p_bar /= 0


    # ln(K_HF) at p_bar = 0

    zln_khf_p0   = 874./t_k - 9.68 + 0.111*sqrt(s)


    # Pressure correction

    zt_degc      = t_k - t_k_zerodegc
    zrt          = gasconst_bar_cm3_o_mol_k * t_k

    zdvi         =   -9.78 + zt_degc*(-0.0090 - zt_degc*0.942E-03)
    zdki         = ( -3.91 + zt_degc*0.054)*1.0E-03

    zln_khf_pp   = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final value of K_HF

    AK_HF_PEFR87 = exp( zln_khf_p0 + zln_khf_pp )

    return (AK_HF_PEFR87)
}



#=======================================================================
AK_HSO4_DICK90 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the dissociation constant of hydrogen sulfate (bisulfate)

    # References: Dickson (1990) -- also Handbook (2007)
    #             Millero (1995) for pressure correction
    # pH scale  : free
    # Note      : converted here from mol/kg-H2O to mol/kg-SW


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt            : R*t_k, R in bar*cm3/(mol*K)
    #     zt_degc        : temperature in degrees Celsius
    #     zdvi           : volume change for ionization
    #     zdki           : compressibility change for ionization
    #     zionst         : ionic strength in mol/-kg-H2O
    #     zsqrti         : square root og ion strength
    #     zcvt_to_kgsw   : mass of pure water in 1kg of seawater as a fct. of salinity
    #     zln_khso4_p0   : K_HSO4 at p_bar = 0
    #     zln_khso4_pp   : pressure correction for p_bar /= 0


    # ln(K_HSO4) at p_bar = 0

    zcvt_to_kgsw = ACVT_KGH2O_O_KGSW(s)
    zionst       = A_IONSTRENGTH_SALIN(s)/zcvt_to_kgsw
    zsqrti       = sqrt(zionst)

    zln_khso4_p0 = (   -4276.1/t_k + 141.328 -  23.093*log(t_k)
                + (-13856./t_k +  324.57 -  47.986*log(t_k)) * zsqrti
                + ( 35474./t_k -  771.54 + 114.723*log(t_k)) * zionst
                - (  2698./t_k)*zsqrti * zionst
                + (  1776./t_k)*zionst*zionst )


    # Pressure correction

    zt_degc      = t_k - t_k_zerodegc
    zrt          = gasconst_bar_cm3_o_mol_k * t_k

    zdvi         =  -18.03 + zt_degc*(0.0466 + zt_degc*0.316E-03)
    zdki         = ( -4.53 + zt_degc*0.0900)*1.0E-03

    zln_khso4_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # ln(K_HSO4) at p_bar = 0

    AK_HSO4_DICK90 = zcvt_to_kgsw * exp( zln_khso4_p0 + zln_khso4_pp )

    return (AK_HSO4_DICK90)
}



#=======================================================================
ASP_CALC_MUCC83 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns stoechiometric solubility product
    # of calcite in seawater

    # References: Mucci (1983)
    #             Millero (1995) for pressure correction
    # pH scale  : N/A
    # Units     : (mol/kg-SW)^2

    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt          : R*t_k, R in bar*cm3/(mol*K)
    #     zsqrts       : square root of salinity
    #     zt_degc      : temperature in degrees Celsius
    #     zdvi         : volume change for ionization
    #     zdki         : compressibility change for ionization
    #     zln_kp1_p0   : ln(K_p1) at p_bar = 0
    #     zln_kp1_pp   : pressure correction for p_bar /= 0


    zsqrts    = sqrt(s)

    # log10(Ksp_Calc) for p_bar = 0
    zlog10_kspcalc_p0 = ( -171.9065 -   0.077993*t_k
                +  2839.319/t_k + 71.595*log10(t_k)
                + (-0.77712 +  0.0028426*t_k + 178.34/t_k)*zsqrts
                -   0.07711*s
                + 0.0041249*s*zsqrts )


    # Pressure correction
    zt_degc   = t_k - t_k_zerodegc
    zrt       = gasconst_bar_cm3_o_mol_k * t_k

    zdvi      =  -48.76 + 0.5304*zt_degc
    zdki      = (-11.76 + 0.3692*zt_degc)*1.0E-03

    zln_kspcalc_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final value of Ksp_Calc

    ASP_CALC_MUCC83 = 10.**(zlog10_kspcalc_p0) * exp(zln_kspcalc_pp)

    return (ASP_CALC_MUCC83)
}



#=======================================================================
ASP_ARAG_MUCC83 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns stoechiometric solubility product
    # of aragonite in seawater

    # References: Mucci (1983)
    #             Millero (1979) for pressure correction
    # pH scale  : N/A
    # Units     : (mol/kg-SW)^2


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar

    # ---------------
    # Local variables
    # ---------------

    #     zrt          : R*t_k, R in bar*cm3/(mol*K)
    #     zsqrts       : square root of salinity
    #     zt_degc      : temperature in degrees Celsius
    #     zdvi         : volume change for ionization
    #     zdki         : compressibility change for ionization
    #     zln_kp1_p0   : ln(K_p1) at p_bar = 0
    #     zln_kp1_pp   : pressure correction for p_bar /= 0


    zsqrts    = sqrt(s)

    # log10(Ksp_Arag) for p_bar = 0
    zlog10_ksparag_p0 = ( -171.945   - 0.077993*t_k
                +   2903.293/t_k + 71.595*log10(t_k)
                + (-0.068393 +  0.0017276*t_k + 88.135/t_k)*zsqrts
                -    0.10018*s
                +  0.0059415*s*zsqrts )


    # Pressure correction
    zt_degc   = t_k - t_k_zerodegc
    zrt       = gasconst_bar_cm3_o_mol_k * t_k

    zdvi      =  -48.76 + 0.5304*zt_degc  + 2.8
    zdki      = (-11.76 + 0.3692*zt_degc)*1.0E-03

    zln_ksparag_pp = (-zdvi + zdki*p_bar/2.)*p_bar/zrt


    # Final value of Ksp_Arag

    ASP_ARAG_MUCC83 = 10.**(zlog10_ksparag_p0) * exp(zln_ksparag_pp)

    return (ASP_ARAG_MUCC83)
}






#=======================================================================
A_BTOT_SALIN <- function (s)
#=======================================================================
{
    # Function returns total borate concentration in mol/kg-SW
    # given the salinity of a sample

    # References: Uppström (1974), cited by  Dickson et al. (2007, chapter 5, p 10)
    #             Millero (1982) cited in Millero (1995)
    # pH scale  : N/A


    A_BTOT_SALIN = 0.000416*(s/35.)

    return (A_BTOT_SALIN)
}


#=======================================================================
A_CATOT_SALIN <- function (s)
#=======================================================================
{
    # Function returns total calcium concentration in mol/kg-SW
    # given the salinity of a sample

    # References: Culkin and Cox (1966)
    # pH scale  : N/A

    # (g Ca/kg)/Cl_permil values
    # Culkin (1967):                         0.0213
    # Culkin and Cox (DSR 13 1966):          0.02126 +/- 0.00004 (stdev)
    # Riley and Tongudai (Chem Geol 2 1967): 0.02128 +/- 0.00006 (stdev)
    # Handbook (2007):                       0.02127
    #    with reference to Riley and Tongudai (1967) (???)

    # Here:
    # (g Ca/kg)/Cl_permil = 0.02127
    # (g Ca)/(mol Ca)     = 40.078
    #  Cl_permil          = S/1.80655
    # mol Ca/kg = (0.02127/40.078) * (35/1.80655)
    #A_CATOT_SALIN = 0.010282*(s/35.)
    A_CATOT_SALIN = (0.02127/40.078) * (s/1.80655)


    return (A_CATOT_SALIN)
}


#=======================================================================
A_FTOT_SALIN <- function (s)
#=======================================================================
{
    # Function returns total calcium concentration in mol/kg-SW
    # given the salinity of a sample

    # References: Culkin (1965) (???)
    # pH scale  : N/A
    A_FTOT_SALIN = 0.000068*(s/35.)


    return (A_FTOT_SALIN)
}


#=======================================================================
A_SO4TOT_SALIN <- function(s)
#=======================================================================
{
    # Function returns total sulfate concentration in mol/kg-SW
    # given the salinity of a sample

    # References: Morris, A.W. and Riley, J.P. (1966) quoted in Handbook (2007)
    # pH scale  : N/A


    #A_SO4TOT_SALIN = 0.028234*(s/35.) # in libthdyct and Thesis
    #A_SO4TOT_SALIN = 0.02824*(s/35.)                # Handbook (2007, chap 6, p 10, tab 2, col 3)
    A_SO4TOT_SALIN = (0.1400/96.062)*(s/1.80655)  # Handbook (2007, chap 6, p 10)

    return (A_SO4TOT_SALIN)
}



#=======================================================================
ACVT_HSWS_O_HTOT <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the ratio H_SWS/H_Tot as a function of salinity s

    # Reference:  Munhoven
    # pH scale:   all

    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar


    # ---------------
    # Local variables
    # ---------------

    #     zso4_tot: total sulfate concentration in mol/kg-SW
    #     zf_tot  : total fluoride concentration in mol/kg-SW

    #-----------------------------------------------------------------------


    zso4_tot = A_SO4TOT_SALIN(s)
    zf_tot   = A_FTOT_SALIN(s)


    ACVT_HSWS_O_HTOT = (1. +  (zf_tot*ABETA_HF_DIRI79(t_k, s, p_bar))
                            /(1. + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar)))

    return (ACVT_HSWS_O_HTOT)
}



#=======================================================================
ACVT_HTOT_O_HFREE <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the ratio H_Tot/H_free as a function of salinity s

    # Reference:  Munhoven
    # pH scale:   N/A


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar


    # ---------------
    # Local variables
    # ---------------

    #     zso4_tot: total sulfate concentration in mol/kg-SW

    #-----------------------------------------------------------------------


    zso4_tot = A_SO4TOT_SALIN(s)


    ACVT_HTOT_O_HFREE = 1. + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar)

    return (ACVT_HTOT_O_HFREE)
}



#=======================================================================
ACVT_HSWS_O_HFREE <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns the ratio H_SWS/H_free as a function
    # of salinity s

    # Reference:  Munhoven
    # pH scale:   N/A


    # ------------------
    # Argument variables
    # ------------------

    #     t_k    : temperature in K
    #     s      : salinity
    #     p_bar  : applied pressure in bar


    # ---------------
    # Local variables
    # ---------------

    #     zso4_tot: total sulfate concentration in mol/kg-SW
    #     zf_tot  : total fluoride concentration in mol/kg-SW

    #-----------------------------------------------------------------------


    zso4_tot = A_SO4TOT_SALIN(s)
    zf_tot   = A_FTOT_SALIN(s)


    ACVT_HSWS_O_HFREE = ( 1. + zf_tot*ABETA_HF_DIRI79(t_k, s, p_bar)
                            + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar) )

    return (ACVT_HSWS_O_HFREE)
}


#=======================================================================
A_RHOSW1_MUNH97 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns first order approximation of \rho in (kg-SW)/(m^3-SW)

    # References: Munhoven (1997)
    #             after EOS80 (UNESCO, 1981, 1983)

    # ------------------
    # Argument variables
    # ------------------

    #     s      : salinity
    #     tk     : temperature in K
    #     p_bar  : depth in m


    # ---------------
    # Local variables
    # ---------------

    s0 =  35.5
    t_k0 = 285.16
    p_bar0 = 300.0



    A_RHOSW1_MUNH97 = ( 1039.9044 + 0.77629393*(s-s0)
                                - 0.19692738*(t_k-t_k0)      )

    return (A_RHOSW1_MUNH97)
}


#=======================================================================
A_RHOSW2_MUNH97 <- function (t_k, s, p_bar)
#=======================================================================
{
    # Function returns second order approximation of \rho in (kg-SW)/(m^3-SW)

    # References: Munhoven (1997)
    #             after EOS80 (UNESCO, 1981, 1983)

    # ------------------
    # Argument variables
    # ------------------

    #     s      : salinity
    #     tk     : temperature in K
    #     p_bar  : depth in m


    # ---------------
    # Local variables
    # ---------------

    s0 =  35.5
    t_k0 = 285.16
    p_bar0 = 300.0


    A_RHOSW2_MUNH97 =( 1040.0145
                    + 0.77629393*(s-s0)
                    - 0.25013591*(t_k-t_k0)
                    + 4.2026266E-02*(p_bar-p_bar0)
                    - 4.7473116E-03*(t_k-t_k0)*(t_k-t_k0)
                    - 4.7974224E-06*(p_bar-p_bar0)*(p_bar-p_bar0)
                    - 2.1404592E-04*(t_k-t_k0)*(p_bar-p_bar0) )

    return (A_RHOSW2_MUNH97)
}




#=======================================================================
CHECKCONSTANTS <- function()
#=======================================================================
{

    # ---------------
    # Local variables
    # ---------------

    #     s      : salinity
    #     tk     : temperature in K
    #     p_bar  : applied pressure in bar



    print ('Checking constant values generated from MOD_CHEMCONST')
    print ('')
    print (' % indicates checking against the Handbook (1994);')
    print (' $ indicates checking against the Lewis and Wallace (1998);')
    print (' * indicates checking against the Handbook (2007);')
    print ('   target values are quoted in brackets')
    print ('')
    print ('')
    print (' For S = 35, P = 0 and T/K = 298.15:')
    print ('')


    s     = 35.
    p_bar = 0.
    t_k   = 298.15


    zkc0 = AK_CARB_0_WEIS74(t_k, s)
    print ('') 
    print ('K_0 -- Weiss (1974)')
    print ('===================')
    print ('')
    print (c('   K_0            :', zkc0))
    print (c('   ln(K_0)        :', log(zkc0)))
    print (c('   pK_0           :', -log10(zkc0)))
    print ('') 


    zkhso4 = AK_HSO4_DICK90(t_k, s, p_bar)

    print ('') 
    print ('K_HSO4 -- Dickson (1990) -- pH_free')
    print ('===================================')
    print ('')
    print (c('   K_HSO4         :', zkhso4))
    print (c('   ln(K_HSO4)     :', log(zkhso4)))
    print (c('   pK_HSO4        :', -log10(zkhso4)))
    print ('') 


    zkb = AK_BORA_DICK90(t_k, s, p_bar)

    print ('') 
    print ('K_b -- Dickson (1990) -- pH_tot')
    print ('===============================')
    print ('')
    print (c('   K_b            :', zkb))
    print (c('   ln(K_b)        :', log(zkb)))
    print (c('   pK_b           :', -log10(zkb)))
    print ('') 




    zkc1 = AK_CARB_1_LUEK00(t_k, s, p_bar)

    print ('') 
    print ('K_1 -- Luecker et al (2000) -- pH_tot')
    print ('=====================================')
    print ('')
    print (c('   K_1            :', zkc1))
    print (c('   ln(K_1)        :', log(zkc1)))
    print (c('   pK_1           :', -log10(zkc1)))
    print ('') 


    zkc2 = AK_CARB_2_LUEK00(t_k, s, p_bar)

    print ('') 
    print ('K_2 -- Luecker et al (2000) -- pH_tot')
    print ('=====================================')
    print ('')
    print (c('   K_2            :', zkc2))
    print (c('   ln(K_2)        :', log(zkc2)))
    print (c('   pK_2           :', -log10(zkc2)))
    print ('') 




    zkc1 = AK_CARB_1_ROYE93(t_k, s, p_bar)

    print ('') 
    print ('K_1 -- Roy et al (1993) -- pH_tot')
    print ('=================================')
    print ('')
    print (c('   K_1            :', zkc1))
    print (c('   ln(K_1)        :', log(zkc1)))
    print (c('   pK_1           :', -log10(zkc1)))
    print ('') 


    zkc2 = AK_CARB_2_ROYE93(t_k, s, p_bar)

    print ('') 
    print ('K_2 -- Roy et al (1993) -- pH_tot')
    print ('=================================')
    print ('')
    print (c('   K_2            :', zkc2))
    print (c('   ln(K_2)        :', log(zkc2)))
    print (c('   pK_2           :', -log10(zkc2)))
    print ('') 




    zkhf = AK_HF_PEFR87(t_k, s, p_bar)

    print ('') 
    print ('K_HF -- Perez and Fraga (1987) -- pH_tot')
    print ('========================================')
    print ('')
    print (c('   K_HF           :', zkhf))
    print (c('   ln(K_HF)       :', log(zkhf)))
    print (c('   pK_HF          :', -log10(zkhf)))
    print ('') 


    zkp1 = AK_PHOS_1_MILL95(t_k, s, p_bar)

    print ('') 
    print ('K_P1 -- Millero (1995) -- pH_SWS')
    print ('================================')
    print ('')
    print (c('   K_P1           :', zkp1))
    print (c('   ln(K_P1)       :', log(zkp1)))
    print (c('   pK_1           :', -log10(zkp1)))
    print ('') 


    zkp2 = AK_PHOS_2_MILL95(t_k, s, p_bar)

    print ('') 
    print ('K_P2 -- Millero (1995) -- pH_SWS')
    print ('================================')
    print ('')
    print (c('   K_2            :', zkp2))
    print (c('   ln(K_P2)       :', log(zkp2)))
    print (c('   pK_2           :', -log10(zkp2)))
    print ('') 


    zkp3 = AK_PHOS_3_MILL95(t_k, s, p_bar)

    print ('') 
    print ('K_P3 -- Millero (1995) -- pH_SWS')
    print ('================================')
    print ('')
    print (c('   K_P3           :', zkp3))
    print (c('   ln(K_P3)       :', log(zkp3)))
    print (c('   pK_P3          :', -log10(zkp3)))
    print ('') 


    zksi1 = AK_SILI_1_MILL95(t_k, s)

    print ('') 
    print ('K_Si1 -- Millero (1995) -- pH_SWS')
    print ('=================================')
    print ('')
    print (c('   K_Si1          :', zksi1))
    print (c('   ln(K_Si1)      :', log(zksi1)))
    print (c('   pK_Si1         :', -log10(zksi1)))
    print ('') 


    zkw = AK_W_MILL95(t_k, s, p_bar)

    print ('') 
    print ('K_w -- Millero (1995) -- pH_SWS')
    print ('===============================')
    print ('')
    print (c('   K_w            :', zkw))
    print (c('   ln(K_w)        :', log(zkw)))
    print (c('   pK_w           :', -log10(zkw)))
    print ('') 


    zkh2s = AK_H2S_1_MILL95(t_k, s, p_bar)

    print ('') 
    print ('K_H2S -- Millero (1995) -- pH_SWS')
    print ('=================================')
    print ('')
    print (c('   K_H2S          :', zkh2s))
    print (c('   ln(K_H2S)      :', log(zkh2s)))
    print (c('   pK_H2S         :', -log10(zkh2s)))
    print ('')


    zknh4 = AK_AMMO_1_YAMI95(t_k, s, p_bar)

    print ('') 
    print ('K_NH4 -- Yao and Millero (1995) -- pH_SWS')
    print ('=========================================')
    print ('')
    print (c('   K_NH4          :', zknh4))
    print (c('   ln(K_NH4)      :', log(zknh4)))
    print (c('   pK_NH4         :', -log10(zknh4)))
    print ('')

    return()
}
