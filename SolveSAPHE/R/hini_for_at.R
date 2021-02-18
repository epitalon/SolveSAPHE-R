
# NaN for [H^+] results
pp_hnan = -1.

#===============================================================================
hini_ACB_DIC <- function (p_alkcb, p_dictot, p_bortot, p_hinf, p_hsup, api)
#===============================================================================

# Function returns the root for the 2nd order approximation of the
# DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
# around the local minimum, if it exists.

# Returns * p_hsup if p_alkcb <= 0
#         * p_hinf if p_alkcb >= 2*p_dictot + p_bortot
#         * the geometric mean of p_hinf and p_hsup if
#           0 < p_alkcb < 2*p_dictot + p_bortot
#           but the 2nd order approximation does not have a solution

{
    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alkcb,           :   Alkalinity from Carbonates + Bore
    #  p_dictot           :   DIC total
    #  p_bortot           :   Bore total
    #  p_hinf, p_hsup     :   infimum and supremum of pH
    #  api                :   Pi factors of dissociation constants

    #-------------------------------------------------------------------------------


    if (p_alkcb <= 0.) 
    {
        z_hini = p_hsup
    } 
    else if (p_alkcb >= (2.*p_dictot + p_bortot)) 
    {
        z_hini = p_hinf
    }
    else
    {
        zca = p_dictot / p_alkcb
        zba = p_bortot / p_alkcb

        # Coefficients of the cubic polynomial
        za2 = api$api1_bor*(1. - zba) + api$api1_dic*(1.-zca)
        za1 = api$api1_dic*api$api1_bor*(1. - zba - zca) + api$api2_dic*(1. - (zca+zca))
        za0 = api$api2_dic*api$api1_bor*(1. - zba - (zca+zca))


        # Taylor expansion around the minimum

        zd = za2*za2 - 3.*za1          # Discriminant of the quadratic equation
                                        # for the minimum close to the root

        if (zd > 0.)
        {              # If the discriminant is positive, i.e.,
                       # if the cubic has two distinct extrema
            zsqrtd = sqrt(zd)

            if (za2 < 0) {
            zhmin =        (-za2 + zsqrtd) / 3.
            } 
            else {
            zhmin = -za1 / ( za2 + zsqrtd)
            }

            z_hini = zhmin + sqrt(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)

            # make h_ini is within [p_hinf,p_hsup]
            z_hini = max(min(p_hsup, z_hini), p_hinf)

        } 
        else
        {
            z_hini = sqrt(p_hinf*p_hsup)
        }
    }

    return (z_hini)
}

#===============================================================================
hini_ACB_CO2 <- function (p_alkcb, p_co2, p_bortot, p_hinf, p_hsup, api)
#===============================================================================

# Function returns the root for the 2nd order approximation of the
# DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
# around the local minimum, if it exists.

# Returns * p_hsup if p_alkcb <= 0
#         * the geometric mean of p_hinf and p_hsup if 0 < p_alkcb
#           but the 2nd order approximation does not have a solution

{

    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alkcb,           :   Alkalinity from Carbonates + Bore
    #  p_pco2             :   [CO2*] concentration
    #  p_bortot           :   Bore total
    #  p_hinf, p_hsup     :   infimum and supremum of pH
    #  api                :   Pi factors of dissociation constants

    #-------------------------------------------------------------------------------


    if (p_alkcb <= 0.) 
    {
        z_hini = p_hsup
    }
    else
    {
        zac = p_alkcb  / p_co2
        zbc = p_bortot / p_co2

        # Coefficients of the cubic polynomial
        za3 = zac
        za2 = -api$api1_dic + api$api1_bor * (zac - zbc)
        za1 = -(api$api1_dic*api$api1_bor + 2. * api$api2_dic)
        za0 = -2.*api$api2_dic*api$api1_bor


        # Taylor expansion around the minimum

        zd = za2*za2 - 3.*za3*za1   # Discriminant of the quadratic equation
                                    # for the minimum close to the root

        if (zd > 0.)
        {              # If the discriminant is positive, i.e.,
                       # if the cubic has two distinct extrema,
                       
            zsqrtd = sqrt(zd)               # locate the minimum.

            if (za2 < 0)
            {
                zhmin = (-za2 + zsqrtd)/3.
            }
            else
            {
                zhmin = -za1/(za2 + zsqrtd)
            }

            zpcb_hmin = za0 + zhmin*(za1 + zhmin*(za2 + zhmin*za3))

            if (zpcb_hmin < 0.)
            {     # If the minimum has a negative ordinate,
                  # it can be used to derive a H_0

                z_hini = zhmin + sqrt(-zpcb_hmin/zsqrtd)

                # Check if this H_ini is compatible
                # with the Alk_CB bounds. We have to
                # check the sign of the auxiliary polynomial 
                za2 = 2. - zac + zbc
                za1 = 2. * api$api1_dic
                za0 = 2. * api$api2_dic

                zpcb_hini = za0 + z_hini*(za1 + z_hini*za2)

                if (zpcb_hini > 0.)
                {
                    # H_ini *is* compatible with Alk_CB:
                    # only make sure it is within [p_hinf,p_hsup].
                    z_hini = max(min(p_hsup, z_hini), p_hinf)

                }
                else
                {
                    # H_ini *is not* compatble with Alk_CB:
                    # use the geometric mean of p_hinf and p_hsup instead.
                    z_hini = sqrt(p_hinf*p_hsup)
                }

            }
            else
            {
                z_hini = sqrt(p_hinf*p_hsup)
            }

        }
        else
        {
            # If the cubic does not have distinct
            # extrema, use the geometric mean of p_hinf and p_hsup

            z_hini = sqrt(p_hinf*p_hsup)
        }

    }

    return (z_hini)
}


#===============================================================================
hini_ACBW_HCO3 <- function (p_alkcbw, p_hco3, p_bortot, p_hinf, p_hsup, api)
#===============================================================================

# Function returns the root for the 2nd order approximation of the
# DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
# around the local minimum, if it exists.

# Returns * p_hsup if p_alkcbw <= 0
#         * the geometric mean of p_hinf and p_hsup if 0 < p_alkcb 
#           but the 2nd order approximation does not have a solution

{

    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alkcb,           :   Alkalinity from Carbonates + Bore
    #  p_phco3            :   [HCO3-] concentration
    #  p_bortot           :   Bore total
    #  p_hinf, p_hsup     :   infimum and supremum of pH
    #  api                :   Pi factors of dissociation constants

    #-------------------------------------------------------------------------------


    if (p_alkcbw <= 0.)
    {
        z_hini = p_hsup
    } 
    else 
    {
        zac = p_alkcbw / p_hco3
        zbc = p_bortot / p_hco3

        # Coefficients of the cubic polynomial
        za3 = 1. / (api$aphscale * p_hco3)
        za2 = zac + api$api1_bor * za3 - 1.
        za0 = -(api$api1_wat / p_hco3 + 2.*api$api2_dic/api$api1_dic)
        za1 = api$api1_bor * (zac - zbc - 1.) + za0
        za0 = api$api1_bor * za0

        # Taylor expansion around the minimum

        zd = za2*za2 - 3.*za3*za1   # Discriminant of the quadratic equation
                                    # for the extrema

        if (zd > 0.) 
        {              # If the discriminant is positive, i.e.,
                       # if the cubic has two distinct extrema
            zsqrtd = sqrt(zd)

            if (za2 < 0)
            {
                zhmin =        (-za2 + zsqrtd ) / (3. * za3)
            }
            else 
            {
                zhmin = -za1 / ( za2 + zsqrtd )
            }

            # Here z_pcbw_hmin = P_CBW(H_min)
            zpcbw_hmin = za0 + zhmin * (za1 + zhmin * (za2 + zhmin*za3))

            if (zpcbw_hmin < 0.)
            {
                z_hini = zhmin + sqrt(-zpcbw_hmin/zsqrtd)
            } 
            else 
            {

                if (za2 < 0.)
                {
                    zhmax = -za1 / ( za2 - zsqrtd )
                }
                else
                {
                    zhmax =        (-za2 - zsqrtd ) / (3. * za3)
                }

                # Here z_pcbw_hmax = P_CBW(H_max) - a_0
                zpcbw_hmax = zhmax * (za1 + zhmax * (za2 + zhmax*za3))

                z_hini = zhmax * (1. - sqrt((zpcbw_hmax + za0)/zpcbw_hmax))

            }
            
            # Still have to make sure that z_hini is within [p_hinf,p_hsup]
            z_hini = max(min(p_hsup, z_hini), p_hinf)

        }
        else
        {   # If the cubic does not have distinct
            # extrema, use the geometric mean of
            # H_min and p_max for H_ini

            z_hini = sqrt(p_hinf*p_hsup)

        }

    }

    return (z_hini)
}


#===============================================================================
hini_ACBW_CO3 <- function (p_alkcbw, p_co3, p_bortot, p_hinf, p_hsup, k_nroots, api)
#===============================================================================

# Function returns the root for the 2nd order approximation of the
# DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
# around the local minimum, if it exists.

# Returns * p_hsup if p_alkcbw <= 0
#         * the geometric mean of p_hinf and p_hsup if 0 < p_alkcb 
#           but the 2nd order approximation does not have a solution

{
    #--------------------#
    # Argument variables #
    #--------------------#

    #  p_alkcbw,          :   Alkalinity from Carbonates + Bore + Water
    #  p_pco3             :   [CO3--] concentration
    #  p_bortot           :   Bore total
    #  p_hinf, p_hsup     :   infimum and supremum of pH (both are 2 element arrays)
    #  k_nroots           :   number of roots to calculate
    #  api                :   Pi factors of dissociation constants

    #-------------------------------------------------------------------------------


    if (p_alkcbw <= 0.)
    {

        # [XXX] Not sure about this ##

        z_hini1 = p_hsup[1]
        z_hini2 = p_hsup[2]

        #~ } else if (p_alkcb >= (2.*p_dictot + p_bortot)) {

        #~   z_hini = p_hinf

    } else {

        zgamma = p_co3*(api$api1_dic/api$api2_dic) - 1. / api$aphscale

        zac = p_alkcbw / p_co3
        zbc = p_bortot / p_co3

                                            # Coefficients of the cubic polynomial
        za3 = zgamma / p_co3
        za2 = -(api$api1_bor * za3 + zac - 2.)
        za0 = api$api1_wat / p_co3            # ... provisionally
        za1 = -api$api1_bor * (zac - zbc - 2.) + za0
        za0 = api$api1_bor * za0              # ... finally


        if (zgamma < 0.) {          # k_nroots = 1 always in this case

            zd = za2*za2 - 3.*za1*za3

            if (zd <= 0.) {

                z_hini1 = sqrt(p_hinf[1]*p_hsup[1])
                z_hini2 = pp_hnan

            } else {

                zsqrtd = sqrt(zd)

                zhmax = - ( za2 + zsqrtd) / (3. * za3)    # Since za3 < 0, we
                zhmin =   (-za2 + zsqrtd) / (3. * za3)    # have H_min < H_max 

                zpcbw_hmax = za0 + zhmax*(za1 + zhmax*(za2 + zhmax * za3))

                if (zpcbw_hmax >= 0)
                {
                    z_hini1 = zhmax + sqrt(zpcbw_hmax/zsqrtd)
                } 
                else
                {
                    # Here: zpcbw_hmin = P_CBW(H_min) - 
                    zpcbw_hmin = zhmin * (za1 + zhmin*(za2 + zhmin * za3))
                    z_hini1    = zhmin * (1. - sqrt(zpcbw_hmin/(zpcbw_hmin + za0)))
                }
                
                # and make sure that z_hini1 is within [p_hinf[1],p_hsup[1]]
                z_hini1 = max(min(p_hsup[1], z_hini1), p_hinf[1])

                z_hini2 = pp_hnan

            }


        } else if (zgamma > 0.) {

            if (k_nroots == 0) 
            {
                z_hini1 = pp_hnan
                z_hini2 = pp_hnan
            }
            else if (k_nroots == 1)                         # H = H_tan
            {
                z_hini1 = p_hinf[1]           # p_hinf[1] = p_hsup[1] = H_tan
                z_hini2 = pp_hnan             # in this case
            }
            else if (k_nroots == 2)
            {
                zd = za2*za2 - 3.*za1*za3

                if (zd > 0.) { 

                    zsqrtd = sqrt(zd)

                    zhmax = -( za2 + zsqrtd) / (3. * za3)    # Since za3 > 0, we
                    zhmin =  (-za2 + zsqrtd) / (3. * za3)    # have H_max < H_min

                    zpcbw_hmin = za0 + zhmin*(za1 + zhmin*(za2 + zhmin * za3))

                    if ((zhmin <= 0.) || (zpcbw_hmin >= 0.)) {

                        z_hini1 = sqrt(p_hinf[1] * p_hsup[1])
                        z_hini2 = sqrt(p_hinf[2] * p_hsup[2])

                    } 
                    else
                    {

                        zhifl = -za2 / (3. * za3)
                        zpcbw_hifl = za0 + zhifl*(za1 + zhifl*(za2 + zhifl * za3))

                        if (zhifl > 0.)
                        {

                            if (zpcbw_hifl > 0.) {

                                z_hini1 = zhmin - (zhmin - zhifl) *
                                                sqrt(zpcbw_hmin/(zpcbw_hmin - zpcbw_hifl))

                            }
                            else if (zpcbw_hifl < 0.) {

                                zpcbw_hmax = za0 + zhmax*(za1 + zhmax*(za2 + zhmax * za3))

                                z_hini1 = zhmax + (zhifl - zhmax) *
                                                sqrt(zpcbw_hmax/(zpcbw_hmax - zpcbw_hifl))

                            }
                            else {

                                z_hini1 = zhifl

                            }

                        } 
                        else {

                                z_hini1 = zhmin * (1. - sqrt(zpcbw_hmin/(zpcbw_hmin - za0)))

                        }

                        z_hini2 = zhmin + sqrt(-zpcbw_hmin/zsqrtd)

                        # and make sure that z_hini1 is within [p_hinf[1],p_hsup[1]]
                        z_hini1 = max(min(p_hsup[1], z_hini1), p_hinf[1])
                        # and z_hini2 within [p_hinf[2],p_hsup[2]]
                        z_hini2 = max(min(p_hsup[2], z_hini2), p_hinf[2])

                    }


                }
                else 
                {   # If the cubic does not have distinct
                    # extrema, use the respective geometric
                    # mean of H_inf and H_sup for H_ini
                    z_hini1 = sqrt(p_hinf[1]*p_hsup[1])
                    z_hini2 = sqrt(p_hinf[2]*p_hsup[2])

                }

            }

        } else {  # zgamma == 0.

            # The cubic equation degenerates to a quadratic

            if (k_nroots == 0) {

                z_hini1 = pp_hnan
                z_hini2 = pp_hnan

            } else {  # k_nroots = 1

                zd = za1*za1 - 4.*za2*za0

                if (za1 > 0.)
                {
                    z_hini1 =                -( za1 + sqrt(zd) ) / (za2 + za2)
                }
                else
                {
                    z_hini1 = - (za0 + za0) / ( za1 - sqrt(zd) )
                }
                
                # Still have to make sure that z_hini1 is within [p_hinf[1],p_hsup[1]]
                z_hini1 = max(min(p_hsup[1], z_hini1), p_hinf[1])

                z_hini2 = pp_hnan

            }

        }

    }

    return (c(z_hini1, z_hini2))
}
